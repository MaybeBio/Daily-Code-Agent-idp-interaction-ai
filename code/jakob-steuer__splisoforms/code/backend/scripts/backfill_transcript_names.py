#!/usr/bin/env python3
"""
backfill_transcript_names.py — populate `isoforms.transcript_name` (display name).

`transcript_id` is the SOURCE identifier and never changes (ENST00000397752,
PB.1969727307.p1). The PacBio IDs carry no meaning, so this derives a readable name
alongside it:

  GENCODE isoform (transcript_id LIKE 'ENST%')
      → the Ensembl `transcript_name` from the GENCODE GTF, e.g.  MET-201

  Novel isoform WITH a canonical_match_id (its closest GENCODE isoform by identity)
      → <matched_transcript_name>-NOV-<first 6 of protein_md5>, e.g.  MET-201-NOV-a1b2c3
        This revives the scheme from the retired rename_to_gencode.py, which named novel
        rows after the closest known transcript. Deterministic: protein_md5 is populated on
        100% of isoforms, so re-running reproduces identical names.

  Novel isoform with NO canonical_match_id
      → left NULL; the UI falls back to transcript_id.

Naming is order-dependent: GENCODE names are resolved first, because novel names are built
FROM them. Run AFTER backfill_canonical_match.py (which supplies canonical_match_id) — a
novel isoform whose match is unset simply gets no name.

Dependencies: GENCODE GTF (--gtf), DB tables `isoforms`.
Idempotent; --force overwrites existing names.

Usage:
    python backend/scripts/backfill_transcript_names.py [--dry-run] [--force]
"""
import argparse
import gzip
import logging
import os
import re
import sys
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from _backfill_io import Progress, print_summary, setup_logging  # noqa: E402
from app.database import SessionLocal  # noqa: E402
from app.models.bio import Isoform, Gene  # noqa: E402

log = logging.getLogger("backfill_transcript_names")

# backend/scripts/ -> backend -> repo -> db -> project root (holds references/)
DEFAULT_GTF = (
    Path(__file__).resolve().parents[4]
    / "references" / "gencode_v47" / "gencode.v47.annotation.gtf.gz"
)

NOV_TAG = "NOV"
MD5_CHARS = 6

_TX_ID = re.compile(r'transcript_id "([^"]+)"')
_TX_NAME = re.compile(r'transcript_name "([^"]+)"')


def load_gtf_names(gtf_path: Path) -> dict[str, str]:
    """Map versionless ENST -> Ensembl transcript_name (ENST00000397752 -> MET-201)."""
    if not gtf_path.exists():
        raise SystemExit(f"GENCODE GTF not found: {gtf_path}\nPass --gtf <path>.")

    names: dict[str, str] = {}
    with gzip.open(gtf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            # Only transcript rows carry transcript_name; skip exon/CDS/UTR lines.
            fields = line.split("\t", 3)
            if len(fields) < 3 or fields[2] != "transcript":
                continue
            m_id, m_name = _TX_ID.search(line), _TX_NAME.search(line)
            if not (m_id and m_name):
                continue
            # GTF carries versions (ENST00000397752.8); the DB stores them versionless.
            names[m_id.group(1).split(".")[0]] = m_name.group(1)
    return names


def run(gtf_path: Path, dry_run: bool, force: bool) -> None:
    log.info("Parsing GENCODE GTF …")
    gtf_names = load_gtf_names(gtf_path)
    log.info(f"  {len(gtf_names):,} ENST -> transcript_name")

    db = SessionLocal()
    try:
        isoforms = db.query(Isoform).all()
        by_id = {i.id: i for i in isoforms}
        gene_sym = dict(db.query(Gene.id, Gene.symbol).all())  # gene_id -> symbol (fallback base)

        stats = {"gencode": 0, "novel": 0, "gene_fallback": 0, "no_match": 0, "skipped": 0, "unmapped_enst": 0}

        # ── Pass 1: GENCODE isoforms adopt the Ensembl transcript_name ──────────────
        pbar = Progress(len(isoforms), "names")
        for iso in isoforms:
            if iso.transcript_name and not force:
                stats["skipped"] += 1
                pbar.update(1)
                continue
            if iso.transcript_id.startswith("ENST"):
                enst = iso.transcript_id.split(".")[0]
                name = gtf_names.get(enst)
                if name:
                    iso.transcript_name = name
                    stats["gencode"] += 1
                else:
                    stats["unmapped_enst"] += 1
            pbar.update(1, status=f"gencode={stats['gencode']}")
        pbar.close()

        # ── Pass 2: novel isoforms are named after their closest GENCODE match ──────
        # Runs second on purpose: it reads the names Pass 1 just assigned.
        for iso in isoforms:
            if iso.transcript_id.startswith("ENST"):
                continue
            if iso.transcript_name and not force:
                continue
            match = by_id.get(iso.canonical_match_id) if iso.canonical_match_id else None
            base = match.transcript_name if match else None
            # Fallback: novel isoforms too divergent to match a GENCODE transcript still
            # get a readable, gene-anchored name (<GENE>-NOV-<md5>) instead of showing the
            # raw PB./assembly id. Requires only the gene symbol + protein_md5.
            via_fallback = False
            if not base:
                base = gene_sym.get(iso.gene_id)
                via_fallback = True
            if not base or not iso.protein_md5:
                stats["no_match"] += 1
                continue
            iso.transcript_name = f"{base}-{NOV_TAG}-{iso.protein_md5[:MD5_CHARS]}"
            stats["gene_fallback" if via_fallback else "novel"] += 1

        if dry_run:
            db.rollback()
            log.info("--dry-run: rolled back, nothing written.")
        else:
            db.commit()

        print_summary("TRANSCRIPT NAME BACKFILL", [
            ("GENCODE named (Ensembl name)", stats["gencode"]),
            ("Novel named (-NOV-, matched)", stats["novel"]),
            ("Novel named (-NOV-, gene fallback)", stats["gene_fallback"]),
            ("Novel with no name (no gene symbol)", stats["no_match"]),
            ("ENST not in GTF", stats["unmapped_enst"]),
            ("Already named (skipped)", stats["skipped"]),
        ])
    finally:
        db.close()


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gtf", type=Path, default=DEFAULT_GTF, help="GENCODE annotation GTF (.gz)")
    ap.add_argument("--dry-run", action="store_true", help="Report without writing.")
    ap.add_argument("--force", action="store_true", help="Overwrite existing names.")
    args = ap.parse_args()

    setup_logging("backfill_transcript_names")
    run(args.gtf, args.dry_run, args.force)


if __name__ == "__main__":
    main()
