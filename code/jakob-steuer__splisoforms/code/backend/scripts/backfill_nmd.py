#!/usr/bin/env python3
"""
Gold-standard NMD backfill (50-nt rule + authoritative biotype).

This replaces the heuristic-only ``backfill_nmd_parallel.py`` in the backfill
chain.  It supplies ``analyze_isoform`` with the inputs it needs to run the real
50-nucleotide rule, instead of falling back to the protein-length heuristic
(which is why the DB was previously ~all ``unknown`` / ``predicted_translatable``
with zero decay calls).

Per-source CDS provenance (transcript-space CDS coords, cds_end = last nt of the
stop codon):
  * ccRCC   TransDecoder ``*.transdecoder.pep`` headers — the ORF whose ``len:``
            matches the isoform's protein length wins (multi-ORF transcripts).
  * GENCODE ``gencode_v47.tsv`` genomic ``cds_start_min``/``cds_end_max`` projected
            onto the isoform's ``exons_json``; validated against protein length.
  * BC      TransDecoder pep, when present (``--bc-transdecoder-dir``); otherwise
            heuristic fallback.

Authoritative biotype (optional, ``--gtf``): GENCODE ``transcript_biotype ==
nonsense_mediated_decay`` is Ensembl's own curated NMD call → ``nmd_target``
directly.  Only matters once NMD-biotype transcripts are ingested (they are not
in the protein-coding-only reference set yet — see build_source_metadata).

The heavy lifting (map building) runs once in the parent, is pickled, and loaded
by each worker — same parallel harness as the other backfills.

Usage:
    python backend/scripts/backfill_nmd.py [--workers 8] [--force] \
        [--gtf data/sources/gencode_v47/gencode.v47.annotation.gtf.gz] \
        [--bc-transdecoder-dir data/sources/breast_cancer/transdecoder]
"""
import argparse
import glob
import gzip
import hashlib
import json
import logging
import os
import pickle
import re
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Dict, List, Optional, Tuple

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from _backfill_io import Progress, print_summary, setup_logging, setup_worker_logging
from _sources_registry import SOURCES_DIR
from sqlalchemy.orm import sessionmaker

from app.analysis.nmd_prediction import analyze_isoform
from app.database import SessionLocal, create_worker_engine
from app.models.bio import Isoform, NmdPrediction

log = logging.getLogger("backfill_nmd")

# ─────────────────────────────────────────────────────────────────────────────
# TransDecoder pep parsing (matches the real header format for this dataset)
#   >PB.2149.37.p1 GENE...  ORF type:complete (+),score=44.40 len:165 PB.2149.37:822-1319(+)
# The last token carries the transcript id and CDS coords (1-based, stop-inclusive).
# ─────────────────────────────────────────────────────────────────────────────
_LAST_TOKEN = re.compile(r"([\w.]+):(\d+)-(\d+)\([+-]\)\s*$")
_ORF_TYPE   = re.compile(r"ORF type:(\S+)")
_ORF_LEN    = re.compile(r"\blen:(\d+)\b")


def parse_transdecoder_pep(paths: List[str]) -> Dict[str, List[dict]]:
    """{transcript_id: [{cds_start, cds_end, orf_len, orf_type}, ...]}."""
    out: Dict[str, List[dict]] = {}
    for path in paths:
        try:
            with open(path) as fh:
                for line in fh:
                    if not line.startswith(">"):
                        continue
                    line = line.rstrip()
                    m = _LAST_TOKEN.search(line)
                    if not m:
                        continue
                    tid = m.group(1)
                    ot  = _ORF_TYPE.search(line)
                    ln  = _ORF_LEN.search(line)
                    out.setdefault(tid, []).append({
                        "cds_start": int(m.group(2)),
                        "cds_end":   int(m.group(3)),
                        "orf_len":   int(ln.group(1)) if ln else None,
                        "orf_type":  ot.group(1) if ot else None,
                    })
        except OSError as e:
            log.warning("could not read pep file %s: %s", path, e)
    return out


_PBID_HEAD = re.compile(r"^(PB\.[\d.]+)\|")


def parse_author_pep_by_md5(paths: List[str]) -> Dict[str, List[Tuple[str, int, int, Optional[str]]]]:
    """{md5(protein): [(pbid, cds_start, cds_end, orf_type), ...]} from the ccRCC
    author pep.

    These are the author-curated ORFs that were actually folded (data/sources/
    ccrcc/author_pep/), so keying by the protein md5 lets us look up the CDS
    coordinates of the *folded* protein directly by the isoform's protein_md5
    (== structure md5). md5 excludes any trailing '*' to match the ingested seq.
    Every ORF is kept (identical proteins from different transcripts can have
    different 5' UTRs and hence different CDS coordinates — review NMD-6); the
    caller picks the entry whose pbid matches the isoform's own transcript."""
    out: Dict[str, List[Tuple[str, int, int, Optional[str]]]] = {}

    def _add(hdr: Optional[str], seq_parts: List[str]) -> None:
        if hdr is None:
            return
        protein = "".join(seq_parts).rstrip("*")
        m = _LAST_TOKEN.search(hdr)
        if not m or not protein:
            return
        md5 = hashlib.md5(protein.encode()).hexdigest()
        ot = _ORF_TYPE.search(hdr)
        pb = _PBID_HEAD.match(hdr)
        out.setdefault(md5, []).append(
            (pb.group(1) if pb else "", int(m.group(2)), int(m.group(3)),
             ot.group(1) if ot else None))

    for path in paths:
        try:
            with open(path) as fh:
                hdr: Optional[str] = None
                seq: List[str] = []
                for line in fh:
                    line = line.rstrip()
                    if line.startswith(">"):
                        _add(hdr, seq)
                        hdr, seq = line[1:], []
                    else:
                        seq.append(line)
                _add(hdr, seq)
        except OSError as e:
            log.warning("could not read author pep %s: %s", path, e)
    return out


# ─────────────────────────────────────────────────────────────────────────────
# Genomic → spliced-transcript coordinate projection (GENCODE CDS)
# ─────────────────────────────────────────────────────────────────────────────
def project_genomic_to_transcript(exons: List[List[int]], strand: str,
                                  gpos: int) -> Optional[int]:
    """1-based transcript offset of a genomic position, or None if not in an exon."""
    ordered = sorted(exons, key=lambda e: e[0], reverse=(strand == "-"))
    offset = 0
    for s, e in ordered:
        if s <= gpos <= e:
            return offset + ((e - gpos) if strand == "-" else (gpos - s)) + 1
        offset += e - s + 1
    return None


# ─────────────────────────────────────────────────────────────────────────────
# Map builders (run once in the parent)
# ─────────────────────────────────────────────────────────────────────────────
def build_ccrcc_maps(db) -> Tuple[Dict[str, Tuple[int, int]], Dict[str, str], dict]:
    """cds_map + orf_type_map for ccRCC, keyed on the isoform's protein md5.

    The CDS coordinates come from the author-curated pep (author_pep/), the ORFs
    that were actually folded — matched by the ingested protein md5 (== structure
    md5), NOT by transcript id + length against the later transdecoder re-run
    (which carries different ORF calls; see build_source_metadata / task ccRCC
    md5 mismatch). This stays correct once the sequences are the folded proteins."""
    pep_dir = SOURCES_DIR / "ccrcc" / "author_pep"
    pep_files = sorted(glob.glob(str(pep_dir / "*.fa"))) + sorted(glob.glob(str(pep_dir / "*.pep")))
    by_md5 = parse_author_pep_by_md5(pep_files)
    log.info("ccRCC: parsed %d author ORFs (by md5) from %d file(s)", len(by_md5), len(pep_files))

    cds_map: Dict[str, Tuple[int, int]] = {}
    orf_type_map: Dict[str, str] = {}
    stats = {"matched": 0, "matched_by_pbid": 0, "no_md5_match": 0}

    from sqlalchemy import text
    rows = db.execute(text(
        "SELECT transcript_id, legacy_transcript_id, md5(rtrim(sequence,'*')) AS seq_md5 "
        "FROM isoforms WHERE source_id=2 AND sequence IS NOT NULL AND sequence<>''"
    )).fetchall()

    for tid, legacy_tid, seq_md5 in rows:
        entries = by_md5.get(seq_md5)
        if not entries:
            stats["no_md5_match"] += 1
            continue
        # transcript_id may carry an "@<sample>" suffix (sample-local PB ids);
        # the bare pbid is what the pep header carries.
        pbid = (legacy_tid or tid).split("@")[0]
        entry = next((e for e in entries if e[0] == pbid), None)
        if entry is not None:
            stats["matched_by_pbid"] += 1
        else:
            entry = entries[0]
        cds_map[tid] = (entry[1], entry[2])
        if entry[3]:
            orf_type_map[tid] = entry[3]
        stats["matched"] += 1

    log.info("ccRCC: cds coords resolved for %(matched)d isoforms "
             "(by own pbid=%(matched_by_pbid)d, no_md5_match=%(no_md5_match)d)", stats)
    return cds_map, orf_type_map, stats


def _load_gencode_tsv() -> Dict[str, dict]:
    """ensembl_transcript_id → {cds_min, cds_max, strand, plen} from gencode_v47.tsv."""
    tsv = SOURCES_DIR / "gencode_v47" / "gencode_v47.tsv"
    out: Dict[str, dict] = {}
    with open(tsv) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            try:
                cmin = f[idx["cds_start_min"]]; cmax = f[idx["cds_end_max"]]
                if not cmin or not cmax:
                    continue
                out[f[idx["ensembl_transcript_id"]]] = {
                    "cds_min": int(cmin), "cds_max": int(cmax),
                    "strand":  f[idx["strand"]],
                    "plen":    int(f[idx["protein_length"]]) if f[idx["protein_length"]] else None,
                }
            except (KeyError, ValueError, IndexError):
                continue
    return out


def build_gencode_cds_map(db) -> Tuple[Dict[str, Tuple[int, int]], dict]:
    """Project GENCODE genomic CDS onto exons_json → transcript-space cds coords."""
    tsv = _load_gencode_tsv()
    cds_map: Dict[str, Tuple[int, int]] = {}
    stats = {"matched": 0, "no_tsv": 0, "no_exons": 0, "off_frame": 0, "unprojectable": 0}

    from sqlalchemy import text
    rows = db.execute(text(
        "SELECT transcript_id, exons_json, length(sequence) AS plen "
        "FROM isoforms WHERE source_id=1 AND exons_json IS NOT NULL AND exons_json<>''"
    )).fetchall()

    for tid, exons_json, plen in rows:
        meta = tsv.get(tid)
        if not meta:
            stats["no_tsv"] += 1
            continue
        try:
            exons = json.loads(exons_json)
        except (json.JSONDecodeError, TypeError):
            stats["no_exons"] += 1
            continue
        strand = meta["strand"]
        # + strand: translation start = cds_min, coding end = cds_max
        # - strand: translation start = cds_max, coding end = cds_min
        g_start = meta["cds_min"] if strand == "+" else meta["cds_max"]
        g_end   = meta["cds_max"] if strand == "+" else meta["cds_min"]
        tx_start = project_genomic_to_transcript(exons, strand, g_start)
        tx_end_coding = project_genomic_to_transcript(exons, strand, g_end)
        if tx_start is None or tx_end_coding is None:
            stats["unprojectable"] += 1
            continue
        coding_len = tx_end_coding - tx_start + 1
        # Validate against protein length (CDS excl stop = plen*3). Guard against
        # frame/coordinate drift; the TSV CDS span excludes the stop codon.
        expected = (plen or 0) * 3
        if expected and abs(coding_len - expected) > 3:
            stats["off_frame"] += 1
            continue
        cds_end = tx_end_coding + 3  # include the stop codon (50-nt rule expects it)
        cds_map[tid] = (tx_start, cds_end)
        stats["matched"] += 1

    log.info("GENCODE: cds coords projected for %(matched)d isoforms "
             "(no_tsv=%(no_tsv)d, off_frame=%(off_frame)d, unprojectable=%(unprojectable)d)", stats)
    return cds_map, stats


def build_bc_maps(db) -> Tuple[Dict[str, Tuple[int, int]], Dict[str, str], dict]:
    """BC (source_id=3) via TransDecoder pep if present, else empty (heuristic)."""
    bc_dir = SOURCES_DIR / "breast_cancer" / "transdecoder"
    pep_files = sorted(glob.glob(str(bc_dir / "*.pep"))) if bc_dir.exists() else []
    if not pep_files:
        log.info("BC: no TransDecoder pep files found locally → heuristic fallback")
        return {}, {}, {"matched": 0, "no_pep_dir": 1}
    orfs = parse_transdecoder_pep(pep_files)
    cds_map: Dict[str, Tuple[int, int]] = {}
    orf_type_map: Dict[str, str] = {}
    stats = {"matched": 0, "no_pep": 0, "no_len_match": 0}
    from sqlalchemy import text
    rows = db.execute(text(
        "SELECT transcript_id, length(rtrim(sequence,'*')) AS plen "
        "FROM isoforms WHERE source_id=3 AND sequence IS NOT NULL AND sequence<>''"
    )).fetchall()
    for tid, plen in rows:
        cands = orfs.get(tid)
        if not cands:
            stats["no_pep"] += 1
            continue
        match = next((o for o in cands if o["orf_len"] == plen), None)
        if match is None:
            stats["no_len_match"] += 1
            continue
        cds_map[tid] = (match["cds_start"], match["cds_end"])
        if match["orf_type"]:
            orf_type_map[tid] = match["orf_type"]
        stats["matched"] += 1
    log.info("BC: cds coords resolved for %(matched)d isoforms", stats)
    return cds_map, orf_type_map, stats


def build_biotype_map_from_db(db) -> Dict[str, str]:
    """{transcript_id: 'nonsense_mediated_decay'} for isoforms whose ingested
    transcript_metadata.biotype marks them as NMD.

    This is dataset-agnostic: any source whose metadata.tsv carries
    ``biotype = nonsense_mediated_decay`` (e.g. a GTF import that keeps both
    protein_coding and nonsense_mediated_decay transcripts) is picked up
    automatically — ingest already stores metadata.biotype verbatim. Empty on the
    current protein-coding-only reference set; populated once such transcripts are
    ingested."""
    from sqlalchemy import text
    rows = db.execute(text(
        "SELECT i.transcript_id FROM isoforms i "
        "JOIN transcript_metadata tm ON tm.isoform_id = i.id "
        "WHERE tm.biotype = 'nonsense_mediated_decay'"
    )).fetchall()
    out = {r[0]: "nonsense_mediated_decay" for r in rows}
    log.info("biotype: %d ingested nonsense_mediated_decay transcripts", len(out))
    return out


def build_biotype_map_from_gtf(gtf_path: Optional[str]) -> Dict[str, str]:
    """Optional supplement: NMD-biotype ids straight from a GENCODE GTF, for the
    case where biotype was not stored at ingest. Prefer the DB source above."""
    if not gtf_path or not os.path.exists(gtf_path):
        return {}
    tid_re = re.compile(r'transcript_id "([^".]+)')
    out: Dict[str, str] = {}
    opener = gzip.open if gtf_path.endswith(".gz") else open
    with opener(gtf_path, "rt") as fh:
        for line in fh:
            if '\ttranscript\t' not in line or \
               'transcript_type "nonsense_mediated_decay"' not in line:
                continue
            m = tid_re.search(line)
            if m:
                out[m.group(1)] = "nonsense_mediated_decay"
    log.info("biotype: %d nonsense_mediated_decay transcripts from GTF", len(out))
    return out


# ─────────────────────────────────────────────────────────────────────────────
# Parallel worker plumbing
# ─────────────────────────────────────────────────────────────────────────────
_worker_session_factory = None
_CDS_MAP: Dict[str, Tuple[int, int]] = {}
_ORF_TYPE_MAP: Dict[str, str] = {}
_BIOTYPE_MAP: Dict[str, str] = {}


def _worker_init(logfile: str, maps_pickle: str):
    setup_worker_logging(logfile)
    global _CDS_MAP, _ORF_TYPE_MAP, _BIOTYPE_MAP
    with open(maps_pickle, "rb") as fh:
        m = pickle.load(fh)
    _CDS_MAP, _ORF_TYPE_MAP, _BIOTYPE_MAP = m["cds"], m["orf_type"], m["biotype"]


def _get_worker_session():
    global _worker_session_factory
    if _worker_session_factory is None:
        engine = create_worker_engine(pool_size=2, max_overflow=1)
        _worker_session_factory = sessionmaker(autocommit=False, autoflush=False, bind=engine)
    return _worker_session_factory()


def process_batch(isoform_ids: List[int], force: bool):
    db = _get_worker_session()
    success = failed = skipped = 0
    errors: list = []
    try:
        for iso_id in isoform_ids:
            try:
                iso = db.query(Isoform).filter(Isoform.id == iso_id).first()
                if not iso:
                    continue
                res = analyze_isoform(
                    db, iso, force=force,
                    cds_map=_CDS_MAP, nmd_biotype_map=_BIOTYPE_MAP, orf_type_map=_ORF_TYPE_MAP,
                )
                if res:
                    success += 1
                else:
                    skipped += 1
            except Exception as e:
                db.rollback()
                failed += 1
                errors.append((iso_id, f"{type(e).__name__}: {e}"))
                log.error("isoform_id=%d %s", iso_id, e)
        return success, skipped, failed, errors
    finally:
        db.close()


def main():
    p = argparse.ArgumentParser(description="NMD backfill — 50-nt rule + biotype")
    p.add_argument("--workers", type=int, default=min(os.cpu_count() or 4, 8))
    p.add_argument("--batch-size", type=int, default=200)
    p.add_argument("--force", action="store_true",
                   help="Re-analyze even if a prediction exists (default: only missing)")
    p.add_argument("--limit", type=int, help="Cap isoforms processed (debug)")
    p.add_argument("--gtf", type=str, default=None,
                   help="Optional GENCODE GTF to supplement the DB nonsense_mediated_decay "
                        "biotype map (normally unnecessary — biotype is read from the DB).")
    args = p.parse_args()

    logfile = setup_logging("backfill_nmd")
    print("[backfill_nmd] gold-standard NMD (50-nt rule + biotype)")

    # ── Build maps once ──────────────────────────────────────────────────────
    db = SessionLocal()
    cds_map: Dict[str, Tuple[int, int]] = {}
    orf_type_map: Dict[str, str] = {}

    cc_cds, cc_orf, _ = build_ccrcc_maps(db)
    cds_map.update(cc_cds); orf_type_map.update(cc_orf)

    gc_cds, _ = build_gencode_cds_map(db)
    cds_map.update(gc_cds)
    # GENCODE protein_coding reference transcripts are complete ORFs by definition;
    # surface that instead of a null "Unknown" in the NMD panel.
    for tid in gc_cds:
        orf_type_map.setdefault(tid, "complete")

    bc_cds, bc_orf, _ = build_bc_maps(db)
    cds_map.update(bc_cds); orf_type_map.update(bc_orf)

    # Biotype: authoritative NMD calls for any ingested nonsense_mediated_decay
    # transcript (dataset-agnostic; reads transcript_metadata.biotype). A GTF may
    # be passed to supplement it if biotype was not stored at ingest.
    biotype_map = build_biotype_map_from_db(db)
    if args.gtf:
        biotype_map.update(build_biotype_map_from_gtf(args.gtf))

    # Scope: with --force, every isoform; otherwise only those lacking a prediction.
    q = db.query(Isoform.id).filter(Isoform.sequence.isnot(None), Isoform.sequence != "")
    if not args.force:
        from sqlalchemy import select
        q = q.filter(Isoform.id.not_in(select(NmdPrediction.isoform_id)))
    if args.limit:
        q = q.limit(args.limit)
    isoform_ids = [r[0] for r in q.all()]
    db.close()

    total = len(isoform_ids)
    print(f"  cds coords: {len(cds_map):,} | NMD-biotype: {len(biotype_map):,} | "
          f"orf_type: {len(orf_type_map):,}")
    print(f"  scope: {total:,} isoforms | {args.workers} workers | force={args.force}")
    if total == 0:
        print("Nothing to process.")
        return

    # Pickle the maps for the workers.
    maps_pickle = os.path.join(tempfile.gettempdir(), f"nmd_maps_{os.getpid()}.pkl")
    with open(maps_pickle, "wb") as fh:
        pickle.dump({"cds": cds_map, "orf_type": orf_type_map, "biotype": biotype_map}, fh)

    batches = [isoform_ids[i:i + args.batch_size] for i in range(0, total, args.batch_size)]
    success = skipped = failed = 0
    errors: list = []
    try:
        with Progress(total, desc="NMD") as pbar, \
             ProcessPoolExecutor(max_workers=args.workers, initializer=_worker_init,
                                 initargs=(str(logfile), maps_pickle)) as ex:
            futs = {ex.submit(process_batch, b, args.force): len(b) for b in batches}
            for fut in as_completed(futs):
                n = futs[fut]
                try:
                    s, k, f, errs = fut.result()
                except Exception as e:
                    s, k, f, errs = 0, 0, n, [(-1, f"batch failure: {e}")]
                success += s; skipped += k; failed += f; errors.extend(errs)
                pbar.update(n, status=f"ok={success} skip={skipped} err={failed}")
    finally:
        try:
            os.remove(maps_pickle)
        except OSError:
            pass

    print_summary(
        "NMD BACKFILL COMPLETE",
        counts=[("Isoforms", total), ("✓ Predicted", success),
                ("– Skipped", skipped), ("✗ Errors", failed)],
        total=total, logfile=logfile, errors=errors,
    )


if __name__ == "__main__":
    main()
