#!/usr/bin/env python3
"""
Backfill per-residue structural impact (DSSP) for existing isoforms.

Usage:
    python backend/scripts/backfill_structural_impact.py [--force] [--gene SYMBOL]

Runs DSSP secondary structure assignment on all isoforms with structures.
With --force, re-analyzes even if results already exist.
With --gene, limits to isoforms of the specified gene symbol.
"""
import argparse
import logging
import os
import sys

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from _backfill_io import Progress, print_summary, setup_logging
from app.database import SessionLocal
from app.models.bio import ExonStructuralImpact, Gene, Isoform

log = logging.getLogger("backfill_structural_impact")


def main():
    parser = argparse.ArgumentParser(description="Backfill structural impact analysis")
    parser.add_argument("--force", action="store_true", help="Re-analyze even if results exist")
    parser.add_argument("--gene", type=str, help="Limit to a specific gene symbol")
    parser.add_argument("--limit", type=int, help="Limit total isoforms to process")
    args = parser.parse_args()

    logfile = setup_logging("backfill_structural_impact")

    db = SessionLocal()
    query = db.query(Isoform)
    if args.gene:
        query = query.join(Gene, Isoform.gene_id == Gene.id).filter(Gene.symbol == args.gene)
    if args.limit:
        query = query.limit(args.limit)
    isoforms = query.all()

    total = len(isoforms)
    scope = []
    if args.gene:   scope.append(f"gene={args.gene}")
    if args.force:  scope.append("force=True")
    if args.limit:  scope.append(f"limit={args.limit}")
    print(f"[backfill_structural_impact] DSSP structural impact backfill")
    print(f"  scope: {total:,} isoforms" + (f"  |  {' | '.join(scope)}" if scope else ""))

    if total == 0:
        print("Nothing to process.")
        db.close()
        return

    from app.analysis.structural_impact import analyze_isoform

    success = 0
    skipped = 0
    failed = 0
    errors: list[tuple[str, str]] = []

    with Progress(total, desc="StructImpact") as pbar:
        for isoform in isoforms:
            try:
                if not args.force:
                    existing = db.query(ExonStructuralImpact).filter(
                        ExonStructuralImpact.isoform_id == isoform.id
                    ).first()
                    if existing:
                        skipped += 1
                        pbar.update(1, status=f"ok={success} skip={skipped} err={failed}")
                        continue

                results = analyze_isoform(db, isoform, force=args.force)
                if results:
                    success += 1
                else:
                    skipped += 1
            except Exception as e:
                failed += 1
                msg = f"{type(e).__name__}: {e}"
                errors.append((isoform.transcript_id, msg))
                log.error("isoform=%s  %s", isoform.transcript_id, msg)
            pbar.update(1, status=f"ok={success} skip={skipped} err={failed}")

    db.close()

    print_summary(
        "STRUCTURAL IMPACT BACKFILL COMPLETE",
        counts=[
            ("Isoforms processed", total),
            ("✓ Analyzed",        success),
            ("– Skipped",         skipped),
            ("✗ Errors",          failed),
        ],
        total=total,
        logfile=logfile,
        errors=errors,
    )


if __name__ == "__main__":
    main()
