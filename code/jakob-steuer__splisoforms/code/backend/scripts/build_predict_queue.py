#!/usr/bin/env python3
"""
build_predict_queue.py
======================
Phase 5 of the GENCODE v47 reference build.

Builds the HPC prediction queue for AlphaFold-3 on Leonardo:
  * one FASTA per unique protein_md5
  * a tab-separated manifest with target IDs, priorities, and source tags
  * SHA-256 checksums so the tarball can be verified at the cluster end
  * matching ``prediction_jobs`` rows in the DB so we can track lifecycle
    (queued → submitted → predicting → predicted → ingested)

Scope (per the AFDB Option-B decision)
--------------------------------------
**Everything gets a fresh AF3 prediction**, including the rows currently
backed by AFDB monomers, so the final structure pool is uniform under
the same model.  De-duplication is by protein_md5: identical proteins
get a single prediction job regardless of how many transcripts encode
them.

The full target set is the union of
  1. every GENCODE v47 protein-coding transcript whose protein passes the
     sequence filters (30 ≤ length ≤ 2700, no X/B/Z/U/* residues), and
  2. every novel DB isoform (no ``ensembl_transcript_id``) that passes
     the same filters.

Priority tiers
--------------
  1  MANE Select transcripts                          (HPC first wave)
  2  Other GENCODE protein-coding transcripts         (alt GENCODE)
  3  Novel DB isoforms not matched to GENCODE         (PB.* / etc.)

The ``queued_for`` column on prediction_jobs records the upstream tag:
``mane_select`` / ``alt_gencode`` / ``novel_db``.

Output layout
-------------
  <out-dir>/predict_queue_v47/
      manifest.tsv          target_id, protein_md5, protein_length,
                            priority_tier, queued_for, source_count
      fastas/<md5>.fasta    one FASTA per unique md5 (header = target_id)
      SHA256SUMS            sha256 over every fasta file
      stats.json            roll-up counts, byte sizes, tier breakdown

Usage
-----
    cd backend
    python scripts/build_predict_queue.py --dry-run
    python scripts/build_predict_queue.py
    python scripts/build_predict_queue.py --tar
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import logging
import re
import sys
import tarfile
from collections import defaultdict
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Set

sys.path.insert(0, str(Path(__file__).parent.parent))

from _backfill_io import Progress, print_summary, setup_logging  # noqa: E402
from sqlalchemy.dialects.postgresql import insert as pg_insert   # noqa: E402

from app.database import SessionLocal                            # noqa: E402
from app.models.bio import Isoform, PredictionJob                # noqa: E402

log = logging.getLogger("build_predict_queue")


# ── Constants ────────────────────────────────────────────────────────────────

GENCODE_VERSION = "v47"
REPO_ROOT       = Path(__file__).resolve().parent.parent.parent
DEFAULT_GENCODE_TSV = REPO_ROOT / "data" / "gencode" / GENCODE_VERSION / "gencode_v47.tsv"
DEFAULT_OUT_DIR     = REPO_ROOT / "data" / "gencode" / GENCODE_VERSION
QUEUE_DIRNAME       = "predict_queue_v47"

MIN_LEN = 30
MAX_LEN = 2700
NON_CANONICAL_AA = set("XBZU*")
AA_LINE_WIDTH = 60  # FASTA line length

PRIORITY_MANE  = 1
PRIORITY_ALT   = 2
PRIORITY_NOVEL = 3

QUEUED_FOR_MANE  = "mane_select"
QUEUED_FOR_ALT   = "alt_gencode"
QUEUED_FOR_NOVEL = "novel_db"


# ── Candidate record ─────────────────────────────────────────────────────────

@dataclass
class Candidate:
    target_id:     str
    protein_md5:   str
    protein_seq:   str
    protein_len:   int
    priority_tier: int
    queued_for:    str
    sources:       Set[str] = field(default_factory=set)  # debug breadcrumbs

    @property
    def fasta(self) -> str:
        header = f">{self.target_id} md5={self.protein_md5} len={self.protein_len} tier={self.priority_tier}"
        body = "\n".join(self.protein_seq[i:i + AA_LINE_WIDTH]
                         for i in range(0, len(self.protein_seq), AA_LINE_WIDTH))
        return f"{header}\n{body}\n"


# ── Filters ──────────────────────────────────────────────────────────────────

def _seq_ok(seq: str) -> bool:
    if not seq:
        return False
    if not (MIN_LEN <= len(seq) <= MAX_LEN):
        return False
    if set(seq) & NON_CANONICAL_AA:
        return False
    return True


# ── GENCODE side ─────────────────────────────────────────────────────────────

def _promote(existing: Candidate, new: Candidate) -> Candidate:
    """
    When the same md5 is reached from multiple sources, keep the *lowest*
    priority tier (MANE > alt > novel) and record the union of source tags.
    """
    if new.priority_tier < existing.priority_tier:
        new.sources |= existing.sources
        return new
    existing.sources |= new.sources
    return existing


def collect_gencode_candidates(tsv_path: Path) -> Dict[str, Candidate]:
    """Load Phase-0 TSV and emit MANE-first / alt-second candidates by md5."""
    log.info("loading GENCODE TSV: %s", tsv_path)
    by_md5: Dict[str, Candidate] = {}
    n_seen = n_filtered = 0
    with open(tsv_path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            n_seen += 1
            seq = (row.get("protein_sequence") or "").strip().upper()
            if not _seq_ok(seq):
                n_filtered += 1
                continue
            md5 = row["protein_md5"]
            tx_name = row["transcript_name"] or row["ensembl_transcript_id"]
            # The Phase-0 fetcher emits the boolean as "1"/"0".
            is_mane = row.get("is_mane_select", "") in ("1", "true", "True", "TRUE")
            cand = Candidate(
                target_id=tx_name,
                protein_md5=md5,
                protein_seq=seq,
                protein_len=len(seq),
                priority_tier=PRIORITY_MANE if is_mane else PRIORITY_ALT,
                queued_for=QUEUED_FOR_MANE if is_mane else QUEUED_FOR_ALT,
                sources={f"gencode:{row['ensembl_transcript_id']}"},
            )
            if md5 in by_md5:
                by_md5[md5] = _promote(by_md5[md5], cand)
            else:
                by_md5[md5] = cand
    log.info("  → %d GENCODE rows seen, %d filtered by length / AA filter, %d unique md5",
             n_seen, n_filtered, len(by_md5))
    return by_md5


# ── DB-novel side ────────────────────────────────────────────────────────────

def add_novel_candidates(by_md5: Dict[str, Candidate], db,
                          skip_existing_af3: bool) -> Dict[str, int]:
    """
    For every DB isoform without an Ensembl mapping that hasn't already been
    captured by GENCODE (md5-collision-safe), add a Tier-3 candidate.

    When `skip_existing_af3` is True (the sensible default), md5s that
    already have at least one isoform with ``structure_provenance =
    'af3_local_v1'`` are *not* queued: those structures came from the same
    AF3 workflow as the new GENCODE batch, so re-predicting is gratuitous.
    AFDB-legacy-only novels (no AF3 yet) are still queued — they need a
    fresh AF3 model for uniformity with the GENCODE pool.
    """
    log.info("scanning DB for novel-only isoforms …")
    rows = (
        db.query(Isoform.transcript_id, Isoform.sequence,
                 Isoform.protein_md5, Isoform.structure_provenance)
        .filter(Isoform.ensembl_transcript_id.is_(None))
        .all()
    )
    # md5 → set of structure_provenance values across all isos sharing this md5
    md5_provenance: Dict[str, set] = defaultdict(set)
    for _, _, md5, prov in rows:
        if md5 and prov:
            md5_provenance[md5].add(prov)

    n_added = n_skipped = n_collided = n_pruned_af3 = 0
    for transcript_id, sequence, protein_md5, _prov in rows:
        if not protein_md5 or not _seq_ok(sequence or ""):
            n_skipped += 1
            continue
        if protein_md5 in by_md5:
            # Already covered by a GENCODE candidate — just tag the breadcrumb.
            by_md5[protein_md5].sources.add(f"db:{transcript_id}")
            n_collided += 1
            continue
        if skip_existing_af3 and "af3_local_v1" in md5_provenance.get(protein_md5, set()):
            n_pruned_af3 += 1
            continue
        cand = Candidate(
            target_id=transcript_id,
            protein_md5=protein_md5,
            protein_seq=sequence.strip().upper(),
            protein_len=len(sequence.strip()),
            priority_tier=PRIORITY_NOVEL,
            queued_for=QUEUED_FOR_NOVEL,
            sources={f"db:{transcript_id}"},
        )
        by_md5[protein_md5] = cand
        n_added += 1
    log.info("  → +%d novel md5 candidates (+%d md5-collisions with GENCODE; "
             "%d pruned as already AF3-predicted; %d skipped by length/AA filter)",
             n_added, n_collided, n_pruned_af3, n_skipped)
    return {"added": n_added, "collided": n_collided,
            "pruned_af3": n_pruned_af3, "skipped": n_skipped}


# ── Disk output ──────────────────────────────────────────────────────────────

def _safe_target_for_filename(target_id: str, md5: str) -> str:
    """target_id can contain '-' '.' etc.  We name FASTAs by md5 to avoid any
    sanitisation drift — the manifest is the source-of-truth mapping."""
    return f"{md5}.fasta"


def write_queue(out_dir: Path, candidates: List[Candidate]) -> Dict[str, object]:
    """
    Writes:
        <out_dir>/manifest.tsv
        <out_dir>/fastas/<md5>.fasta  (one per md5)
        <out_dir>/SHA256SUMS
        <out_dir>/stats.json
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    fasta_dir = out_dir / "fastas"
    # Drop any pre-existing FASTAs so a rebuild with smaller scope (e.g. after
    # toggling --include-existing-af3) doesn't leave orphans behind.
    if fasta_dir.exists():
        n_existing = sum(1 for _ in fasta_dir.iterdir())
        if n_existing:
            log.info("clearing %d existing FASTAs from %s", n_existing, fasta_dir)
            for p in fasta_dir.iterdir():
                p.unlink()
    fasta_dir.mkdir(exist_ok=True)

    manifest_path = out_dir / "manifest.tsv"
    sums_path     = out_dir / "SHA256SUMS"
    stats_path    = out_dir / "stats.json"

    tier_counts = defaultdict(int)
    total_bytes = 0
    sha_lines: List[str] = []

    with open(manifest_path, "w", newline="") as mf:
        writer = csv.writer(mf, delimiter="\t")
        writer.writerow([
            "target_id", "protein_md5", "protein_length",
            "priority_tier", "queued_for", "source_count", "fasta_filename",
        ])
        with Progress(len(candidates), desc="write fastas") as pbar:
            for c in candidates:
                fname = _safe_target_for_filename(c.target_id, c.protein_md5)
                fpath = fasta_dir / fname
                data = c.fasta.encode()
                fpath.write_bytes(data)
                sha = hashlib.sha256(data).hexdigest()
                sha_lines.append(f"{sha}  fastas/{fname}")
                writer.writerow([
                    c.target_id, c.protein_md5, c.protein_len,
                    c.priority_tier, c.queued_for, len(c.sources), fname,
                ])
                tier_counts[c.priority_tier] += 1
                total_bytes += len(data)
                pbar.update(1)

    sums_path.write_text("\n".join(sha_lines) + "\n")

    stats = {
        "gencode_version": GENCODE_VERSION,
        "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds").replace("+00:00", "Z"),
        "n_candidates": len(candidates),
        "tier_counts": {
            "mane_select":  tier_counts[PRIORITY_MANE],
            "alt_gencode":  tier_counts[PRIORITY_ALT],
            "novel_db":     tier_counts[PRIORITY_NOVEL],
        },
        "total_fasta_bytes": total_bytes,
        "min_len": MIN_LEN,
        "max_len": MAX_LEN,
    }
    stats_path.write_text(json.dumps(stats, indent=2) + "\n")
    log.info("wrote manifest, %d FASTAs (%d bytes), checksums, stats",
             len(candidates), total_bytes)
    return stats


def make_tar(out_dir: Path, tar_path: Path) -> None:
    log.info("creating tar: %s", tar_path)
    with tarfile.open(tar_path, "w:gz") as tar:
        tar.add(out_dir, arcname=out_dir.name)


# ── DB-side: prediction_jobs ────────────────────────────────────────────────

def upsert_prediction_jobs(db, candidates: List[Candidate]) -> int:
    """
    Insert (or no-op upsert) a prediction_jobs row per candidate.  Existing
    rows with a non-failed status are left untouched.  Conflicts on target_id
    are resolved by keeping the existing row — re-runs are idempotent.
    """
    n_inserted = 0
    BATCH = 500
    rows = []
    for c in candidates:
        rows.append(dict(
            target_id=c.target_id,
            protein_md5=c.protein_md5,
            protein_sequence=c.protein_seq,
            protein_length=c.protein_len,
            priority_tier=c.priority_tier,
            queued_for=c.queued_for,
            status="queued",
            gencode_version=GENCODE_VERSION,
        ))
        if len(rows) >= BATCH:
            stmt = pg_insert(PredictionJob).values(rows)
            stmt = stmt.on_conflict_do_nothing(index_elements=["target_id"])
            result = db.execute(stmt)
            n_inserted += result.rowcount or 0
            rows.clear()
    if rows:
        stmt = pg_insert(PredictionJob).values(rows)
        stmt = stmt.on_conflict_do_nothing(index_elements=["target_id"])
        result = db.execute(stmt)
        n_inserted += result.rowcount or 0
    db.commit()
    return n_inserted


# ── Orchestrator ─────────────────────────────────────────────────────────────

def run(dry_run: bool, gencode_tsv: Path, out_root: Path,
        make_tarball: bool, skip_existing_af3: bool,
        reset_prediction_jobs: bool) -> None:
    logfile = setup_logging("build_predict_queue")

    # 1. Build candidate pool from GENCODE
    by_md5 = collect_gencode_candidates(gencode_tsv)

    # 2. Add novel DB-only candidates
    db = SessionLocal()
    try:
        if reset_prediction_jobs and not dry_run:
            n_dropped = db.query(PredictionJob).delete()
            db.commit()
            log.info("reset prediction_jobs: deleted %d existing rows", n_dropped)

        novel_stats = add_novel_candidates(by_md5, db, skip_existing_af3)
        n_novel = novel_stats["added"]

        # 3. Sort: priority_tier ASC, then protein_length ASC (small first
        #    finishes earlier on the HPC and keeps tier-1 throughput high).
        candidates = sorted(
            by_md5.values(),
            key=lambda c: (c.priority_tier, c.protein_len, c.target_id),
        )

        # 4. Write queue layout
        queue_dir = out_root / QUEUE_DIRNAME
        if not dry_run:
            stats = write_queue(queue_dir, candidates)
            n_inserted = upsert_prediction_jobs(db, candidates)
            log.info("inserted %d new prediction_jobs rows", n_inserted)
            if make_tarball:
                tar_path = out_root / f"{QUEUE_DIRNAME}.tar.gz"
                make_tar(queue_dir, tar_path)
        else:
            stats = {
                "n_candidates": len(candidates),
                "tier_counts": _tier_counts(candidates),
                "would_write_to": str(queue_dir),
            }
            n_inserted = 0
            log.info("DRY RUN — no files written, no DB rows inserted")

        print_summary(
            "PREDICT QUEUE BUILD COMPLETE" + (" (DRY RUN)" if dry_run else ""),
            counts=[
                ("MANE Select (tier 1)",   stats["tier_counts"]["mane_select"]   if not dry_run else _tier_counts(candidates)[PRIORITY_MANE]),
                ("Alt GENCODE (tier 2)",   stats["tier_counts"]["alt_gencode"]   if not dry_run else _tier_counts(candidates)[PRIORITY_ALT]),
                ("Novel DB-only (tier 3)", stats["tier_counts"]["novel_db"]      if not dry_run else _tier_counts(candidates)[PRIORITY_NOVEL]),
                ("Total queue size",       len(candidates)),
                ("Novel md5s added",       n_novel),
                ("Novel md5s pruned (already AF3)", novel_stats["pruned_af3"]),
                ("DB rows inserted",       n_inserted),
            ],
            total=len(candidates) or 1,
            logfile=logfile,
            notes=[
                f"output: {queue_dir}",
                "filters: 30≤len≤2700, no X/B/Z/U/*",
                f"skip_existing_af3: {skip_existing_af3}",
                f"reset_prediction_jobs: {reset_prediction_jobs}",
            ],
        )
    finally:
        db.close()


def _tier_counts(cands: List[Candidate]) -> Dict[int, int]:
    out = defaultdict(int)
    for c in cands:
        out[c.priority_tier] += 1
    return out


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build the AF3 prediction queue for GENCODE v47 + novel DB isoforms.",
        formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__,
    )
    parser.add_argument("--gencode-tsv", type=Path, default=DEFAULT_GENCODE_TSV)
    parser.add_argument("--out-dir",     type=Path, default=DEFAULT_OUT_DIR,
                        help="Parent dir where predict_queue_v47/ is written.")
    parser.add_argument("--tar",         action="store_true",
                        help="Also produce predict_queue_v47.tar.gz")
    parser.add_argument("--dry-run",     action="store_true")
    parser.add_argument("--include-existing-af3", action="store_true",
                        help="Include novel md5s that already have a structure_provenance='af3_local_v1' "
                             "prediction.  Default: skip them (their structures came from the same AF3 "
                             "workflow as this batch, so re-prediction is gratuitous).")
    parser.add_argument("--reset-prediction-jobs", action="store_true",
                        help="TRUNCATE the prediction_jobs table before inserting fresh rows.  Use this "
                             "when rebuilding the queue with different scope (e.g. after toggling "
                             "--include-existing-af3) to avoid leaving stale rows behind.")
    args = parser.parse_args()

    if not args.gencode_tsv.exists():
        log.error("GENCODE TSV not found: %s", args.gencode_tsv)
        log.error("Run scripts/fetch_gencode_reference.py first.")
        sys.exit(1)

    run(
        dry_run=args.dry_run, gencode_tsv=args.gencode_tsv,
        out_root=args.out_dir, make_tarball=args.tar,
        skip_existing_af3=not args.include_existing_af3,
        reset_prediction_jobs=args.reset_prediction_jobs,
    )


if __name__ == "__main__":
    main()
