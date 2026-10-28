"""
Backfill PTM site accessibility analysis for isoforms using programmatic retrieval.

Usage:
    python backend/scripts/backfill_ptm.py [--force] [--gene SYMBOL] [--limit N]
                                            [--retry-errors] [--errors-file PATH]

Strategy:
    - Uses IsoProDB mapping (TranscriptMetadata) to identify target UniProt accessions.
    - Programmatically retrieves PTM features from the UniProtKB REST API.
    - No longer requires local PhosphoSitePlus files.

Serial by design: UniProt REST has per-IP rate limits, so parallel hammering is
counterproductive.  Use --retry-errors after a run to re-attempt failures.
"""
import argparse
import logging
import os
import sys
import traceback
from typing import List, Tuple

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from _backfill_io import Progress, print_summary, setup_logging
from app.analysis.ptm_accessibility import analyze_isoform
from app.database import SessionLocal
from app.models import bio as bio_models

log = logging.getLogger("backfill_ptm")

# Default path for the errors file so --retry-errors can find it automatically.
_DEFAULT_ERRORS_FILE = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "ptm_backfill_errors.txt"
)


def main():
    parser = argparse.ArgumentParser(description="Backfill PTM accessibility for all isoforms")
    parser.add_argument("--force",        action="store_true", help="Re-analyze even if results exist")
    parser.add_argument("--gene",         type=str,  help="Only process isoforms for this gene symbol")
    parser.add_argument("--limit",        type=int,  default=None, help="Limit number of isoforms to process")
    parser.add_argument("--retry-errors", action="store_true",
                        help=f"Only process transcript IDs listed in --errors-file (default: {_DEFAULT_ERRORS_FILE})")
    parser.add_argument("--errors-file",  type=str,  default=_DEFAULT_ERRORS_FILE,
                        help="Path to save (or read) the list of failed transcript IDs")
    args = parser.parse_args()

    logfile = setup_logging("backfill_ptm")

    db = SessionLocal()
    query = db.query(bio_models.Isoform)
    query = query.filter(bio_models.Isoform.sequence.isnot(None), bio_models.Isoform.sequence != "UNKNOWN")

    scope: list[str] = []
    if args.retry_errors:
        if not os.path.exists(args.errors_file):
            print(f"[backfill_ptm] No errors file at {args.errors_file}. Nothing to retry.")
            return
        with open(args.errors_file) as fh:
            tids = [line.split("\t")[0].strip() for line in fh
                    if line.strip() and not line.startswith("#")]
        if not tids:
            print("[backfill_ptm] Errors file is empty. Nothing to retry.")
            return
        query = query.filter(bio_models.Isoform.transcript_id.in_(tids))
        args.force = True
        scope.append(f"retry-errors=True ({len(tids)} ids)")
    else:
        if args.gene:
            query = query.join(bio_models.Gene, bio_models.Isoform.gene_id == bio_models.Gene.id) \
                         .filter(bio_models.Gene.symbol == args.gene)
            scope.append(f"gene={args.gene}")
        if args.limit:
            query = query.limit(args.limit)
            scope.append(f"limit={args.limit}")
    if args.force:
        scope.append("force=True")

    # Gene-ordered so consecutive isoforms share a canonical structure (SASA memo)
    # and a UniProt accession (site cache) — see ptm_accessibility / ptm_data.
    isoforms = query.order_by(bio_models.Isoform.gene_id, bio_models.Isoform.id).all()
    total = len(isoforms)

    print("[backfill_ptm] PTM accessibility backfill (UniProt REST)")
    print(f"  scope: {total:,} isoforms" + (f"  |  {' | '.join(scope)}" if scope else ""))

    if total == 0:
        print("Nothing to process.")
        db.close()
        return

    success = 0
    skipped = 0
    errors: List[Tuple[str, str]] = []

    with Progress(total, desc="PTM") as pbar:
        for iso in isoforms:
            tid = iso.transcript_id
            try:
                results = analyze_isoform(db, iso, force=args.force)
                if results:
                    success += 1
                else:
                    skipped += 1
            except Exception as exc:
                short = f"{type(exc).__name__}: {exc}"
                errors.append((tid, short))
                # Full traceback only to the logfile (DEBUG)
                log.error("isoform=%s  %s\n%s", tid, short, traceback.format_exc())
                try:
                    db.rollback()
                except Exception:
                    pass
            pbar.update(1, status=f"ok={success} skip={skipped} err={len(errors)}")

    # Persist / clean up the errors file so --retry-errors works
    if errors:
        with open(args.errors_file, "w") as fh:
            fh.write("# transcript_id\terror\n")
            for tid, msg in errors:
                fh.write(f"{tid}\t{msg}\n")
        retry_hint = [(f"Errors file: {args.errors_file}", "")]
    else:
        if os.path.exists(args.errors_file):
            os.remove(args.errors_file)
        retry_hint = None

    db.close()
    try:
        from app.analysis.ptm_data import flush_sites_cache
        flush_sites_cache()
    except Exception as exc:                 # cache is an optimisation, never fatal
        log.warning("could not flush PTM site cache: %s", exc)

    print_summary(
        "PTM BACKFILL COMPLETE",
        counts=[
            ("Isoforms processed", total),
            ("✓ Analyzed",        success),
            ("– Skipped",         skipped),
            ("✗ Errors",          len(errors)),
        ],
        total=total,
        logfile=logfile,
        errors=errors,
    )
    if errors:
        print(f"  Error list: {args.errors_file}  (re-run with --retry-errors)")
    from _backfill_io import exit_status
    sys.exit(exit_status(total, len(errors)))


if __name__ == "__main__":
    main()
