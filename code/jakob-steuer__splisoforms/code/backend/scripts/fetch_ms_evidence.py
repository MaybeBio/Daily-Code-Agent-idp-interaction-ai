#!/usr/bin/env python3
"""
fetch_ms_evidence.py
====================
Phase 4 of the GENCODE v47 reference build.

Overlays mass-spectrometry evidence from PeptideAtlas (canonical proteome)
and TransCODE Consortium (dark-proteome / peptideins) onto every isoform
in the DB, setting ``isoforms.ms_evidence = TRUE`` when at least one
*distinguishing* tryptic peptide of the isoform is observed in the
combined evidence set.

A peptide is *distinguishing* when it is present in this isoform's tryptic
digest but absent from the gene's MANE Select reference protein.  The
flag therefore answers "is this isoform-specific protein form observed
by MS?", not the weaker "is this protein observed at all".  Use
``--mode any`` to fall back to "any tryptic peptide of this isoform is
in the evidence set" if you want the broader flag.

This script makes **no** structural assumptions; it operates purely on
protein sequence.

Data sources (download manually before running)
-----------------------------------------------
  PeptideAtlas Human All-Distinct-Peptide list
      https://db.systemsbiology.net/sbeams/cgi/PeptideAtlas/buildDetails
      → "All Distinct Peptide Mappings" → save as TSV.

  TransCODE Consortium peptidein peptide table
      Supplementary table accompanying the 2026 Nature paper.
      → save as TSV with one peptide sequence per row.

Both inputs are flat text — one peptide per line — and the script can
accept multiple files (``--ms-evidence path1.tsv path2.tsv``); their
union is the evidence set.

Inputs / outputs
----------------
  --ms-evidence            list of peptide-set files (one peptide per line
                           or a TSV with a "peptide" column)
  --mode                   distinguishing (default) | any
  --min-len / --max-len    tryptic peptide length cutoffs (defaults 7..30)
  --miscleavages           tryptic miscleavage count (default 0)
  data/gencode/v47/ms_evidence_report.tsv (output report)
  backend/logs/fetch_ms_evidence.log

Algorithm
---------
1. Build the evidence set (uppercase, strip non-AA).
2. For each gene, find the MANE Select isoform; compute its tryptic peptide
   set.  Skip the gene if it has no MANE Select isoform yet (those are
   the alt-only genes that came from BC ingestion — handled by Phase 5
   once predictions land).
3. For each isoform in the gene, compute its peptide set, subtract the
   MANE peptide set to get *distinguishing* peptides, intersect with the
   evidence set.  If non-empty → ms_evidence = TRUE.

Usage
-----
    cd backend
    python scripts/fetch_ms_evidence.py \
        --ms-evidence data/peptideatlas/distinct_peptides.tsv \
                      data/transcode/peptideins.tsv \
        --mode distinguishing

    python scripts/fetch_ms_evidence.py --ms-evidence ... --dry-run
"""

from __future__ import annotations

import argparse
import logging
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Set, Tuple

sys.path.insert(0, str(Path(__file__).parent.parent))

from _backfill_io import Progress, print_summary, setup_logging  # noqa: E402

from app.database import SessionLocal                            # noqa: E402
from app.models.bio import Gene, Isoform                         # noqa: E402

log = logging.getLogger("fetch_ms_evidence")


GENCODE_VERSION = "v47"
REPO_ROOT = Path(__file__).resolve().parent.parent.parent
DEFAULT_REPORT = REPO_ROOT / "data" / "gencode" / GENCODE_VERSION / "ms_evidence_report.tsv"

# Only canonical amino acids; anything else interrupts a tryptic fragment.
_CANONICAL_AA = set("ACDEFGHIKLMNPQRSTVWY")
_AA_STRIP_RE  = re.compile(r"[^A-Z]")


# ── Tryptic digest ────────────────────────────────────────────────────────────

def tryptic_digest(
    seq: str, *, min_len: int = 7, max_len: int = 30, miscleavages: int = 0,
) -> Set[str]:
    """
    Standard in-silico trypsin: cleaves C-terminal to K/R unless the next
    residue is P.  Returns the set (not multiset) of peptides whose length
    is within [min_len, max_len], optionally with up to `miscleavages`
    skipped cleavage sites.

    Leading initiator-methionine is *not* removed here.  PeptideAtlas
    distributes peptides as observed, with or without the leading Met
    cleaved; checking both forms is left as a future refinement.
    """
    s = _AA_STRIP_RE.sub("", seq.upper())
    if not s:
        return set()

    # Find cleavage sites: positions right after K/R not followed by P.
    cuts: List[int] = [0]
    for i in range(len(s) - 1):
        if s[i] in "KR" and s[i + 1] != "P":
            cuts.append(i + 1)
    cuts.append(len(s))

    out: Set[str] = set()
    n = len(cuts) - 1
    # Generate fragments with 0..miscleavages skipped sites.
    for skip in range(0, miscleavages + 1):
        for i in range(0, n - skip):
            frag = s[cuts[i]: cuts[i + 1 + skip]]
            if min_len <= len(frag) <= max_len and set(frag).issubset(_CANONICAL_AA):
                out.add(frag)
    return out


# ── Evidence-file ingestion ───────────────────────────────────────────────────

def load_evidence(paths: List[Path]) -> Set[str]:
    """
    Load peptide sequences from one or more flat-text files.  Accepts a TSV
    with a header column named "peptide" / "sequence" / "Peptide", or a
    plain newline-delimited list.

    Strips quotes and any non-amino-acid characters; uppercases.  Empty
    after stripping → discarded.
    """
    out: Set[str] = set()
    for path in paths:
        log.info("loading MS evidence: %s", path)
        with open(path) as fh:
            header_line = fh.readline().rstrip("\n")
            if "\t" in header_line:
                cols = header_line.split("\t")
                # Heuristic column pick
                pep_col = None
                for cand in ("peptide", "sequence", "Peptide", "Sequence", "peptide_sequence"):
                    if cand in cols:
                        pep_col = cols.index(cand)
                        break
                if pep_col is None:
                    # No recognisable column; treat first line as data, first column = peptide
                    pep_col = 0
                    parts = header_line.split("\t")
                    p = _clean_peptide(parts[pep_col])
                    if p:
                        out.add(p)
                else:
                    pass  # header consumed
                for line in fh:
                    parts = line.rstrip("\n").split("\t")
                    if pep_col >= len(parts):
                        continue
                    p = _clean_peptide(parts[pep_col])
                    if p:
                        out.add(p)
            else:
                # plain peptide-per-line (header_line was actually first peptide)
                p = _clean_peptide(header_line)
                if p:
                    out.add(p)
                for line in fh:
                    p = _clean_peptide(line)
                    if p:
                        out.add(p)
    log.info("  → %d unique peptides total", len(out))
    return out


def _clean_peptide(s: str) -> str:
    s = s.strip().strip('"').upper()
    s = _AA_STRIP_RE.sub("", s)
    if len(s) < 5 or not set(s).issubset(_CANONICAL_AA):
        return ""
    return s


# ── Main pass ─────────────────────────────────────────────────────────────────

def run(
    dry_run: bool, ms_paths: List[Path], report_path: Path, mode: str,
    min_len: int, max_len: int, miscleavages: int,
) -> None:
    logfile = setup_logging("fetch_ms_evidence")
    log.info("mode=%s  min_len=%d  max_len=%d  miscleavages=%d",
             mode, min_len, max_len, miscleavages)

    evidence = load_evidence(ms_paths)
    if not evidence:
        log.error("evidence set is empty — aborting.")
        sys.exit(1)

    db = SessionLocal()
    try:
        log.info("loading isoforms + genes …")
        isoforms: List[Isoform] = (
            db.query(Isoform)
            .filter(Isoform.sequence.isnot(None))
            .filter(Isoform.sequence != "")
            .filter(Isoform.sequence != "UNKNOWN")
            .order_by(Isoform.gene_id, Isoform.id)
            .all()
        )
        log.info("  → %d isoforms with non-empty sequence", len(isoforms))

        # Group by gene_id for MANE-relative digestion.
        by_gene: Dict[int, List[Isoform]] = defaultdict(list)
        for iso in isoforms:
            by_gene[iso.gene_id].append(iso)

        # For each gene: digest the MANE Select if available, then per-isoform.
        report: List[Tuple[int, str, int, int, int, str]] = []
        # (isoform_id, transcript_id, n_pep_total, n_pep_distinguishing,
        #  n_pep_hit, decision)
        n_with_ms       = 0
        n_without_ms    = 0
        n_no_mane       = 0
        n_short_seq     = 0

        with Progress(len(by_gene), desc="ms-overlay") as pbar:
            for gene_id, group in by_gene.items():
                mane = next((i for i in group if i.is_mane_select), None)
                mane_peps: Set[str] = set()
                if mane is not None:
                    mane_peps = tryptic_digest(
                        mane.sequence, min_len=min_len, max_len=max_len,
                        miscleavages=miscleavages,
                    )

                for iso in group:
                    peps = tryptic_digest(
                        iso.sequence, min_len=min_len, max_len=max_len,
                        miscleavages=miscleavages,
                    )
                    if not peps:
                        n_short_seq += 1
                        report.append((iso.id, iso.transcript_id, 0, 0, 0, "no_tryptic_peps"))
                        continue

                    if mode == "distinguishing":
                        if mane is None:
                            # Without a MANE Select baseline we cannot decide
                            # which peptides are isoform-specific.  Flag and skip.
                            n_no_mane += 1
                            report.append((iso.id, iso.transcript_id, len(peps), 0, 0, "no_mane_baseline"))
                            continue
                        # MANE Select itself: skip — by definition not isoform-distinguishing.
                        if iso.id == mane.id:
                            report.append((iso.id, iso.transcript_id, len(peps), 0, 0, "is_mane_select"))
                            continue
                        distinguishing = peps - mane_peps
                    else:
                        distinguishing = peps  # "any" mode

                    hit = distinguishing & evidence
                    has_evidence = bool(hit)
                    if not dry_run:
                        iso.ms_evidence = has_evidence
                    if has_evidence:
                        n_with_ms += 1
                        report.append((iso.id, iso.transcript_id, len(peps),
                                       len(distinguishing), len(hit), "evidence_yes"))
                    else:
                        n_without_ms += 1
                        report.append((iso.id, iso.transcript_id, len(peps),
                                       len(distinguishing), 0, "evidence_no"))

                pbar.update(1, status=f"yes={n_with_ms} no={n_without_ms}")

        if not dry_run:
            db.commit()
            log.info("committed ms_evidence flags: %d TRUE, %d FALSE",
                     n_with_ms, n_without_ms)

        report_path.parent.mkdir(parents=True, exist_ok=True)
        with open(report_path, "w") as fh:
            fh.write("isoform_id\ttranscript_id\tn_pep_total\tn_pep_distinguishing\tn_pep_hit\tdecision\n")
            for r in report:
                fh.write("\t".join(str(x) for x in r) + "\n")
        log.info("report written: %s", report_path)

        print_summary(
            "MS EVIDENCE OVERLAY COMPLETE" + (" (DRY RUN)" if dry_run else ""),
            counts=[
                ("ms_evidence = TRUE",      n_with_ms),
                ("ms_evidence = FALSE",     n_without_ms),
                ("skipped (no MANE baseline)", n_no_mane),
                ("skipped (no tryptic peps)",  n_short_seq),
            ],
            total=(n_with_ms + n_without_ms + n_no_mane + n_short_seq) or 1,
            logfile=logfile,
            notes=[
                f"evidence peptides: {len(evidence)}",
                f"report TSV: {report_path}",
            ],
        )

    finally:
        db.close()


def main() -> None:
    parser = argparse.ArgumentParser(
        description="MS-evidence overlay: flag isoforms whose distinguishing "
                    "tryptic peptides are observed by PeptideAtlas / TransCODE.",
        formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__,
    )
    parser.add_argument("--ms-evidence", type=Path, nargs="+", required=True,
                        help="One or more peptide-set TSV / text files.")
    parser.add_argument("--mode", choices=("distinguishing", "any"),
                        default="distinguishing")
    parser.add_argument("--min-len", type=int, default=7,
                        help="Minimum tryptic peptide length (default 7).")
    parser.add_argument("--max-len", type=int, default=30,
                        help="Maximum tryptic peptide length (default 30).")
    parser.add_argument("--miscleavages", type=int, default=0,
                        help="Allowed tryptic miscleavages (default 0).")
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    for p in args.ms_evidence:
        if not p.exists():
            log.error("MS evidence file not found: %s", p)
            sys.exit(1)

    run(
        dry_run=args.dry_run, ms_paths=args.ms_evidence, report_path=args.report,
        mode=args.mode, min_len=args.min_len, max_len=args.max_len,
        miscleavages=args.miscleavages,
    )


if __name__ == "__main__":
    main()
