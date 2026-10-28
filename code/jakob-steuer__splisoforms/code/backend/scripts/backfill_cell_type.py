#!/usr/bin/env python3
"""
Backfill cell-type expression data from Human Protein Atlas for existing genes.

Usage:
    python backend/scripts/backfill_cell_type.py [--force] [--gene SYMBOL]

Fetches single-cell type expression data from HPA for each gene.
With --force, re-fetches even if data already exists.
With --gene, limits to a specific gene symbol.
"""
import os
import sys
import argparse
import time

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from sqlalchemy import not_

from app.database import SessionLocal
from app.models.bio import Gene, CellTypeExpression


def main():
    parser = argparse.ArgumentParser(description="Backfill HPA cell-type expression")
    parser.add_argument("--force", action="store_true", help="Re-fetch even if results exist")
    parser.add_argument("--gene", type=str, help="Limit to a specific gene symbol")
    parser.add_argument("--limit", type=int, help="Max genes to process")
    args = parser.parse_args()

    db = SessionLocal()

    query = db.query(Gene)
    if args.gene:
        query = query.filter(Gene.symbol == args.gene)
    else:
        # Exclude PB.* (PacBio IDs) and 'unknown' gene symbols — these will never
        # resolve to an Ensembl gene ID and would waste API calls.
        query = query.filter(
            not_(Gene.symbol.like("PB.%")),
            Gene.symbol != "unknown",
            Gene.symbol.isnot(None),
        )

    if args.limit:
        query = query.limit(args.limit)

    genes = query.all()
    print(f"Found {len(genes)} genes to process (PB.* and unknown symbols excluded)")

    from app.analysis.cell_type_specificity import analyze_gene

    success = 0
    skipped = 0
    failed = 0

    for i, gene in enumerate(genes, 1):
        try:
            if not args.force:
                existing = db.query(CellTypeExpression).filter(
                    CellTypeExpression.gene_id == gene.id
                ).first()
                if existing:
                    skipped += 1
                    continue

            results = analyze_gene(db, gene, force=args.force)
            if results:
                success += 1
                if success % 50 == 0:
                    print(f"  [{i}/{len(genes)}] Fetched {success} genes...")
            else:
                skipped += 1

            # Rate-limit HPA API calls
            time.sleep(0.2)

        except Exception as e:
            failed += 1
            print(f"  Error on {gene.symbol}: {e}")

    db.close()
    print(f"\nDone: {success} fetched, {skipped} skipped, {failed} failed")


if __name__ == "__main__":
    main()
