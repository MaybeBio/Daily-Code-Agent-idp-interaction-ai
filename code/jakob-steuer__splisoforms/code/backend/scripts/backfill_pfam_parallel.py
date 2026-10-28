#!/usr/bin/env python3
"""
High-throughput parallel backfill for Pfam domains.
Utilizes multiple CPU cores to process the ~80,000 isoform dataset.

Usage:
    python backend/scripts/backfill_pfam_parallel.py --workers 12 --batch-size 50
"""
import argparse
import logging
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import List

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from _backfill_io import Progress, print_summary, setup_logging, setup_worker_logging
from sqlalchemy import select
from sqlalchemy.orm import sessionmaker

from app.analysis.pfam_domains import analyze_isoform
from app.database import SessionLocal, create_worker_engine
from app.models.bio import Isoform, PfamDomain

log = logging.getLogger("backfill_pfam_parallel")

# Module-level worker session factory (initialized per-process)
_worker_session_factory = None


def _get_worker_session():
    """Create a process-local engine and session factory on first call."""
    global _worker_session_factory
    if _worker_session_factory is None:
        worker_engine = create_worker_engine(pool_size=2, max_overflow=1)
        _worker_session_factory = sessionmaker(
            autocommit=False, autoflush=False, bind=worker_engine,
        )
    return _worker_session_factory()


def process_batch(isoform_ids: List[int], force: bool = False):
    """Worker: returns (success, skipped, failed, errors)."""
    db = _get_worker_session()
    success = 0
    skipped = 0
    failed = 0
    errors: list[tuple[int, str]] = []
    try:
        for iso_id in isoform_ids:
            try:
                isoform = db.query(Isoform).filter(Isoform.id == iso_id).first()
                if not isoform:
                    continue
                if not force:
                    existing = db.query(PfamDomain).filter(PfamDomain.isoform_id == iso_id).first()
                    if existing:
                        skipped += 1
                        continue
                results = analyze_isoform(db, isoform, force=force)
                if results:
                    db.commit()
                    success += 1
                else:
                    skipped += 1
            except Exception as e:
                db.rollback()
                failed += 1
                msg = f"{type(e).__name__}: {e}"
                errors.append((iso_id, msg))
                log.error("isoform_id=%d  %s", iso_id, msg)
        return success, skipped, failed, errors
    finally:
        db.close()


def main():
    parser = argparse.ArgumentParser(description="Parallel Pfam Backfill")
    parser.add_argument("--workers",    type=int, default=min(os.cpu_count() or 4, 8),
                        help="Number of parallel workers (default: min(cpu_count, 8))")
    parser.add_argument("--batch-size", type=int, default=50, help="Isoforms per batch")
    parser.add_argument("--force",      action="store_true", help="Re-analyze even if results exist")
    parser.add_argument("--limit",      type=int, help="Limit total isoforms to process")
    args = parser.parse_args()

    logfile = setup_logging("backfill_pfam_parallel")

    db = SessionLocal()
    # Gene-ordered so a batch holds siblings of the same gene: the canonical is
    # then scanned once per worker (scan memo in pfam_domains) instead of once
    # per isoform, which halves the work.
    query = db.query(Isoform.id).order_by(Isoform.gene_id, Isoform.id)
    if not args.force:
        processed = select(PfamDomain.isoform_id)
        query = query.filter(Isoform.id.not_in(processed))
    if args.limit:
        query = query.limit(args.limit)
    isoform_ids = [r[0] for r in query.all()]
    db.close()

    total = len(isoform_ids)
    scope = [f"{args.workers} workers", f"batch {args.batch_size}"]
    if args.force: scope.append("force=True")
    if args.limit: scope.append(f"limit={args.limit}")
    print("[backfill_pfam_parallel] Pfam domain backfill")
    print(f"  scope: {total:,} isoforms  |  {' | '.join(scope)}")

    if total == 0:
        print("Nothing to process.")
        return

    batches = [isoform_ids[i:i + args.batch_size] for i in range(0, total, args.batch_size)]

    success = skipped = failed = 0
    errors: list[tuple[int, str]] = []

    with Progress(total, desc="Pfam") as pbar, \
         ProcessPoolExecutor(max_workers=args.workers,
                             initializer=setup_worker_logging,
                             initargs=(str(logfile),)) as executor:
        futures = {executor.submit(process_batch, b, args.force): len(b) for b in batches}
        for future in as_completed(futures):
            batch_n = futures[future]
            try:
                s, k, f, errs = future.result()
            except Exception as e:
                s, k, f, errs = 0, 0, batch_n, [(-1, f"batch failure: {e}")]
            success += s
            skipped += k
            failed  += f
            errors.extend(errs)
            pbar.update(batch_n, status=f"ok={success} skip={skipped} err={failed}")

    print_summary(
        "PFAM BACKFILL COMPLETE",
        counts=[
            ("Isoforms processed", total),
            ("✓ With domains",     success),
            ("– Skipped",          skipped),
            ("✗ Errors",           failed),
        ],
        total=total,
        logfile=logfile,
        errors=errors,
    )


if __name__ == "__main__":
    main()
