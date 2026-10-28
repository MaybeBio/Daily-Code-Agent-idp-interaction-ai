#!/usr/bin/env python3
"""
backfill_pae_blocks.py
======================
Computes PAE summary statistics (``pae_mean``) and PAE rigid blocks
(``pae_blocks``) on ``structural_metrics`` for every isoform that has a
structure with a PAE file, from any source (GENCODE/AFDB, ccRCC, BC).

The PAE path is taken from ``structures.af3_pae_url`` — the column the
ingester fills for every structure — and written back to
``structural_metrics.pae_json_path`` once resolved.  (Review STR-1: the old
query required ``pae_json_path IS NOT NULL`` AND ``transcript_id LIKE 'PB.%'``,
so it selected 0 isoforms and every ``pae_mean``/``pae_blocks`` stayed NULL.)

``*_summary_confidences.json`` paths are upgraded to ``*_confidences.json``
(the full matrix) when the file is present.

Usage:
    .venv/bin/python backend/scripts/backfill_pae_blocks.py [--genes GENE1,GENE2] [--limit N] [--dry-run]
"""

import os
import sys
import json
import time
import argparse
import numpy as np

sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))

from app.database import SessionLocal
from app.models import bio as bio_models
from app.services.pae_service import PAEService


def summary_to_full_path(summary_path: str) -> str:
    """Convert a *_summary_confidences.json path to *_confidences.json."""
    return summary_path.replace("_summary_confidences.json", "_confidences.json")


def main(gene_symbols=None, limit=None, dry_run=False, force=False):
    db = SessionLocal()
    pae_service = PAEService()

    t_start = time.time()
    updated = 0
    already_done = 0
    no_file = 0
    wrong_format = 0
    errors = 0

    print("[*] Backfilling PAE blocks + pae_mean for all isoforms with a PAE file")
    print(f"    dry_run={dry_run}  force={force}")
    print()

    # Source of truth for the PAE file is structures.af3_pae_url (set at ingest
    # for every structure).  structural_metrics is outer-joined: a missing row
    # is created below rather than silently skipping the isoform.
    query = (
        db.query(bio_models.Isoform, bio_models.Structure.af3_pae_url)
        .join(bio_models.Structure, bio_models.Structure.isoform_id == bio_models.Isoform.id)
        .outerjoin(bio_models.StructuralMetrics,
                   bio_models.StructuralMetrics.isoform_id == bio_models.Isoform.id)
        .filter(bio_models.Structure.af3_pae_url.isnot(None),
                bio_models.Structure.af3_pae_url != "")
    )

    if gene_symbols:
        query = query.join(bio_models.Gene, bio_models.Isoform.gene_id == bio_models.Gene.id) \
                     .filter(bio_models.Gene.symbol.in_(gene_symbols))

    if not force:
        # Only process isoforms that still have no pae_blocks (incl. no metrics row)
        query = query.filter(
            (bio_models.StructuralMetrics.id.is_(None)) |
            (bio_models.StructuralMetrics.pae_blocks.is_(None)) |
            (bio_models.StructuralMetrics.pae_blocks == "[]")
        )

    query = query.order_by(bio_models.Isoform.id)
    if limit:
        query = query.limit(limit)

    # One row per isoform (an isoform may own several structure rows).
    seen_ids: set = set()
    rows = []
    for isoform, af3_pae_url in query.all():
        if isoform.id in seen_ids:
            continue
        seen_ids.add(isoform.id)
        rows.append((isoform, af3_pae_url))
    total = len(rows)
    print(f"    Found {total} isoforms to process\n")

    for idx, (isoform, af3_pae_url) in enumerate(rows):
        metrics = isoform.structural_metrics
        if metrics is None and not dry_run:
            metrics = bio_models.StructuralMetrics(isoform_id=isoform.id)
            isoform.structural_metrics = metrics  # attach via relationship (see STR-1)
            db.add(metrics)
        current_path = af3_pae_url

        # Derive the full confidences path
        if "_summary_confidences.json" in current_path:
            full_path = summary_to_full_path(current_path)
        elif "_confidences.json" in current_path:
            full_path = current_path  # already pointing to full file
        else:
            full_path = current_path  # unknown format, try as-is

        # Try to load the full PAE matrix; fall back to the ingested path if the
        # derived *_confidences.json is not present (summary files hold no matrix).
        matrix = pae_service.load_pae_matrix(full_path)
        if matrix is None and full_path != current_path:
            matrix = pae_service.load_pae_matrix(current_path)
            if matrix is not None:
                full_path = current_path

        if matrix is None:
            # Distinguish between truly missing and wrong-format files
            from app.services.storage import STATIC_DATA_DIR
            if full_path.startswith("/static/data/"):
                resolved = str(STATIC_DATA_DIR / full_path.replace("/static/data/", "", 1).lstrip("/"))
            elif full_path.startswith("static/data/"):
                resolved = str(STATIC_DATA_DIR / full_path.replace("static/data/", "", 1).lstrip("/"))
            else:
                resolved = full_path
            if os.path.exists(resolved):
                wrong_format += 1
                if wrong_format <= 3:
                    print(f"  [BAD FMT] {isoform.transcript_id}: file exists but has no pae key → {full_path}")
                    print(f"            Hint: run script.08.fix_pae_symlinks.v1.0.py on the cluster, then re-rsync.")
            else:
                no_file += 1
                if no_file <= 5:
                    print(f"  [MISSING] {isoform.transcript_id}: {full_path}")
            continue

        try:
            blocks = pae_service.find_rigid_blocks(matrix)
            pae_mean = float(np.mean(matrix))
        except Exception as e:  # noqa: BLE001 — per-isoform, counted and reported below
            errors += 1
            if errors <= 5:
                print(f"  [ERROR] {isoform.transcript_id}: {e}")
            continue

        if dry_run:
            print(f"  [DRY] {isoform.transcript_id}: {len(blocks)} rigid blocks, mean_pae={pae_mean:.2f} ({full_path})")
            updated += 1
            continue

        # Update the record
        metrics.pae_json_path = full_path
        metrics.pae_blocks = json.dumps(blocks)
        metrics.pae_mean = pae_mean

        updated += 1

        if updated % 100 == 0:
            db.commit()
            elapsed = time.time() - t_start
            rate = updated / elapsed
            eta = (total - updated) / rate if rate > 0 else 0
            print(f"  [{updated}/{total}] {isoform.transcript_id} "
                  f"— {len(blocks)} blocks | {rate:.1f}/s | ETA {eta/60:.1f}min")

    if not dry_run:
        db.commit()

    elapsed = time.time() - t_start
    print()
    print(f"[+] Done in {elapsed:.1f}s")
    print(f"    Updated:          {updated}")
    print(f"    Wrong format:     {wrong_format}  (file exists but contains summary data, not full PAE)")
    print(f"    File not found:   {no_file}")
    print(f"    Errors:           {errors}")
    if wrong_format > 0:
        print()
        print("  → To fix wrong-format files: a wrong-format PAE file is usually a")
        print("    symlink to a stub or to the summary JSON instead of the full")
        print("    *_confidences.json(.gz). Repoint it at the full confidences file")
        print("    on the prediction host, re-sync the structures tree, and re-run")
        print("    this backfill script.")
    if no_file > 0:
        print()
        print("  → If many files are missing, run script.rsync_pae_confidences.v1.0.sh first.")
    if updated > 0 and not dry_run:
        print()
        print("  → Next: re-run neoantigen backfill to fix Gold Medal tiers:")
        print("      .venv/bin/python backend/scripts/backfill_neoantigen.py --force")

    db.close()

    # PASS ≠ did-work (ENG-2): fail when there were targets and none was updated.
    if total > 0 and updated == 0:
        print(f"[FAIL] {total} targets, 0 updated "
              f"(missing={no_file}, wrong_format={wrong_format}, errors={errors})", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Backfill PAE rigid blocks + pae_mean for all isoforms with a PAE file")
    parser.add_argument("--genes", type=str, help="Comma-separated gene symbols (default: all)")
    parser.add_argument("--limit", type=int, help="Max isoforms to process")
    parser.add_argument("--dry-run", action="store_true", help="Show what would be done without writing")
    parser.add_argument("--force", action="store_true", help="Re-process even if pae_blocks already set")
    args = parser.parse_args()

    gene_symbols = [g.strip() for g in args.genes.split(",")] if args.genes else None
    main(gene_symbols=gene_symbols, limit=args.limit, dry_run=args.dry_run, force=args.force)
