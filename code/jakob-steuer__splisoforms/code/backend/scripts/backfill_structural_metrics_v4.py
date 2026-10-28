#!/usr/bin/env python3
"""
backfill_structural_metrics_v4.py — per-isoform structural metrics.

For every isoform with a structure: per-residue pLDDT series + secondary
structure (PyDSSP, via structural_impact.analyze_isoform, which also writes
exon_structural_impacts), PAE mean + rigid PAE blocks (PAEService).

Parallel: the per-isoform work is independent, so `--workers N` fans the target
ids out over a ProcessPoolExecutor (each worker owns its own engine/session —
same harness as backfill_nmd.py). Serial throughput is ~4 isoforms/s (~9 h for
the whole DB); 10 workers bring that to ~1 h.

Skip-existing: without --force only isoforms lacking a StructuralMetrics row or
lacking plddt_series are processed, so a stopped run resumes where it left off.

Usage:
    python backend/scripts/backfill_structural_metrics_v4.py [--workers 10] [--force]
                                                             [--genes A,B] [--limit N]
"""
import argparse
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import List, Tuple

# One BLAS/Accelerate thread per process: with N worker processes each spinning
# up a full thread pool for numpy (PyDSSP einsums, PAE matrices) the workers
# thrash and the parallel run is no faster than serial (observed: 10 workers at
# 830 % CPU, 3.5 iso/s — the serial rate). Must be set BEFORE numpy is imported;
# spawned workers re-import this module, so they inherit it too.
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np

# Add backend directory to Python path
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
sys.path.append(os.path.dirname(os.path.realpath(__file__)))

from app.database import SessionLocal, create_worker_engine
from app.models import bio as bio_models
from app.services.pae_service import PAEService
try:
    from app.analysis import structural_impact
except ImportError:
    structural_impact = None
from app.analysis import portable_structural_impact


def standardize_pae_path(pae_url: str) -> str:
    if not pae_url: return pae_url
    if "backend/static/" in pae_url:
        pae_url = "/static/" + pae_url.split("backend/static/")[-1].lstrip('/')
    elif "http" in pae_url:
        if "/static/" in pae_url:
            pae_url = "/static/" + pae_url.split("/static/")[-1].lstrip('/')
    elif not pae_url.startswith("/static/"):
        pae_url = "/static/" + pae_url.lstrip('/')
    return pae_url


def _process_one(db, pae_service, iso_id: int, force: bool) -> None:
    isoform = db.query(bio_models.Isoform).get(iso_id)
    if not isoform:
        return

    metrics = isoform.structural_metrics
    if not metrics:
        metrics = bio_models.StructuralMetrics(isoform_id=isoform.id)
        # Attach via the relationship, NOT only db.add(): the one-to-one
        # `isoform.structural_metrics` was already loaded as None above, and
        # structural_impact.analyze_isoform() reads that attribute — a bare
        # db.add() leaves it None, analyze_isoform creates a SECOND row for
        # the same isoform_id, and the UNIQUE constraint rolls the whole
        # isoform back (review STR-1: 124,687/124,687 errors, exit 0).
        isoform.structural_metrics = metrics
        db.add(metrics)
        db.flush()

    # 2. Deep hydration: pLDDT series + secondary structure (+ exon_structural_impacts)
    if not metrics.plddt_series or force:
        if structural_impact:
            structural_impact.analyze_isoform(db, isoform, **({"force": True} if force else {}))
        else:
            portable_structural_impact.analyze_isoform_portable(db, isoform)
        db.flush()

    # 3. PAE: path standardisation, mean, rigid blocks
    struct = db.query(bio_models.Structure).filter(bio_models.Structure.isoform_id == isoform.id).first()
    pae_url = metrics.pae_json_path or (struct.af3_pae_url if struct else None)
    pae_url = standardize_pae_path(pae_url)
    if pae_url:
        metrics.pae_json_path = pae_url
        if not metrics.pae_blocks or force:
            matrix = pae_service.load_pae_matrix(pae_url)
            if matrix is not None:
                blocks = pae_service.find_rigid_blocks(matrix)
                metrics.pae_blocks = json.dumps(blocks)
                metrics.pae_mean = float(np.mean(matrix))


# ── worker plumbing ──────────────────────────────────────────────────────────
_worker_session_factory = None


def _worker_session():
    global _worker_session_factory
    if _worker_session_factory is None:
        from sqlalchemy.orm import sessionmaker
        engine = create_worker_engine(pool_size=2, max_overflow=1)
        _worker_session_factory = sessionmaker(autocommit=False, autoflush=False, bind=engine)
    return _worker_session_factory()


def process_batch(iso_ids: List[int], force: bool, session_factory=None) -> Tuple[int, int, list]:
    """Process a list of isoform ids. `session_factory` defaults to the worker
    engine (ProcessPool path); the serial path passes SessionLocal."""
    make = session_factory or _worker_session
    db = make()
    pae_service = PAEService()
    ok = err = 0
    errors: list = []
    try:
        for i, iso_id in enumerate(iso_ids):
            try:
                _process_one(db, pae_service, iso_id, force)
                ok += 1
                if ok % 10 == 0:
                    db.commit()
            except Exception as e:
                db.rollback()
                err += 1
                errors.append((iso_id, f"{type(e).__name__}: {e}"))
            if i % 50 == 49:
                db.commit(); db.close(); db = make()
        db.commit()
    finally:
        db.close()
    return ok, err, errors


def target_ids(gene_symbols=None, limit=None, force=False) -> List[int]:
    db = SessionLocal()
    query = db.query(bio_models.Isoform.id)
    if gene_symbols:
        query = query.join(bio_models.Gene, bio_models.Isoform.gene_id == bio_models.Gene.id).filter(
            bio_models.Gene.symbol.in_(gene_symbols))
    elif not force:
        # Only unfinished targets (resume-safe)
        query = query.outerjoin(bio_models.StructuralMetrics,
                                bio_models.Isoform.id == bio_models.StructuralMetrics.isoform_id).filter(
            (bio_models.StructuralMetrics.id == None) | (bio_models.StructuralMetrics.plddt_series == None))
    if limit:
        query = query.limit(limit)
    ids = [r.id for r in query.all()]      # Row (column query) or ORM object both expose .id
    db.close()
    return ids


def backfill_v4_metrics(gene_symbols=None, limit=None, force=False, workers: int = 1,
                        batch_size: int = 100):
    print("[*] Backfilling v4 Structural Metrics (Robust Targeted Refresh)...")
    isoform_ids = target_ids(gene_symbols, limit, force)
    print(f"[*] Found {len(isoform_ids)} target isoforms to analyze.  workers={workers}")
    t_start = time.time()
    processed = errors = 0
    err_list: list = []

    if workers <= 1:
        processed, errors, err_list = process_batch(isoform_ids, force, session_factory=SessionLocal)
        for iso_id, msg in err_list[:20]:
            print(f"  [ERROR] ID {iso_id}: {msg}")
    else:
        batches = [isoform_ids[i:i + batch_size] for i in range(0, len(isoform_ids), batch_size)]
        done = 0
        with ProcessPoolExecutor(max_workers=workers) as ex:
            futs = {ex.submit(process_batch, b, force): len(b) for b in batches}
            for fut in as_completed(futs):
                n = futs[fut]
                try:
                    ok, err, errs = fut.result()
                except Exception as e:
                    ok, err, errs = 0, n, [(-1, f"batch failure: {e}")]
                processed += ok; errors += err; err_list.extend(errs)
                done += n
                if done % 2000 < batch_size or done == len(isoform_ids):
                    rate = done / max(time.time() - t_start, 1)
                    print(f"  [{done}/{len(isoform_ids)}] ok={processed} err={errors} "
                          f"({rate:.1f} iso/s, ETA {(len(isoform_ids)-done)/max(rate,0.01)/60:.0f} min)",
                          flush=True)
        for iso_id, msg in err_list[:20]:
            print(f"  [ERROR] ID {iso_id}: {msg}")

    elapsed = time.time() - t_start
    print(f"\n[+] Done in {elapsed:.1f}s — {processed} processed, {errors} errors.")
    return processed, errors


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Backfill structural metrics v4")
    parser.add_argument("--genes", type=str, help="Comma-separated gene symbols")
    parser.add_argument("--limit", type=int, help="Max isoforms to process")
    parser.add_argument("--force", action="store_true", help="Re-process even if metrics exist")
    parser.add_argument("--workers", type=int, default=max(1, min((os.cpu_count() or 4) - 2, 10)),
                        help="Parallel worker processes (default: min(cpu-2, 10))")
    parser.add_argument("--batch-size", type=int, default=100)
    args = parser.parse_args()

    symbols = args.genes.split(",") if args.genes else None
    n_ok, n_err = backfill_v4_metrics(gene_symbols=symbols, limit=args.limit, force=args.force,
                                      workers=args.workers, batch_size=args.batch_size)
    # run_local_backfills.sh marks PASS on exit 0. If every target errored
    # nothing was done, so fail loudly instead of reporting PASS (STR-1 / ENG-2).
    from _backfill_io import exit_status
    sys.exit(exit_status(n_ok + n_err, n_err))
