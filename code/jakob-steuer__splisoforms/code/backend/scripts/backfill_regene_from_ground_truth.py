#!/usr/bin/env python3
"""
backfill_regene_from_ground_truth.py — STEP 2 of the ground-truth re-attribution.

Applies data/regene/targets.tsv (written by backfill_gene_registry.py): for every
novel isoform whose ground-truth gene (BC GTF ENSG / ccRCC SQANTI associated_gene)
differs from its current gene, move it and recompute delta_length vs the new gene's
canonical length (left NULL when the new gene has no canonical yet — the downstream
canonical-designation pass fills it).

Idempotent. Does NOT run the downstream canonical / analysis recompute — that is a
separate, deliberate step.

Usage:
    python backend/scripts/backfill_regene_from_ground_truth.py            # dry-run
    python backend/scripts/backfill_regene_from_ground_truth.py --apply
"""
import os
import sys
import csv
import argparse
from pathlib import Path

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from sqlalchemy import func
from app.database import SessionLocal
from app.models.bio import Gene, Isoform

# backend/scripts/ -> backend -> repo root
TARGETS = Path(__file__).resolve().parents[2] / "data" / "regene" / "targets.tsv"


def main():
    ap = argparse.ArgumentParser(description="Apply ground-truth gene re-attribution")
    ap.add_argument("--apply", action="store_true")
    args = ap.parse_args()
    if not TARGETS.exists():
        sys.exit(f"missing {TARGETS} — run backfill_gene_registry.py --apply first")

    db = SessionLocal()

    # rows where the target differs from current
    moves = []  # (iso_id, target_gid)
    with open(TARGETS) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            cur, tgt = int(r["current_gene_id"]), int(r["target_gene_id"])
            if cur != tgt:
                moves.append((int(r["isoform_id"]), tgt))
    print(f"re-attribution moves: {len(moves)}", flush=True)

    canon_len = dict(db.query(Gene.id, Gene.canonical_length).all())
    seqlen = dict(db.query(Isoform.id, func.length(Isoform.sequence))
                  .filter(Isoform.id.in_([m[0] for m in moves])).all()) if moves else {}

    # audit: distinct source + target genes touched
    from collections import Counter
    tgt_counter = Counter(t for _, t in moves)
    print(f"  distinct target genes       : {len(tgt_counter)}")
    print(f"  → to newly-created genes     : "
          f"{sum(n for g, n in tgt_counter.items() if canon_len.get(g) is None)}")

    if not args.apply:
        print("\n(dry-run — re-run with --apply)")
        db.close(); return

    updates = []
    for iso_id, tgt in moves:
        cl = canon_len.get(tgt)
        sl = seqlen.get(iso_id)
        updates.append({"id": iso_id, "gene_id": tgt,
                        "delta_length": (sl - cl) if (cl and sl) else None})
    for i in range(0, len(updates), 5000):
        db.bulk_update_mappings(Isoform, updates[i:i+5000]); db.commit()

    print(f"\napplied {len(updates)} re-attributions.")
    print("NEXT (separate step): re-run canonical designation, drop orphan stub genes,")
    print("  and recompute canonical-relative analyses for moved isoforms.")
    db.close()


if __name__ == "__main__":
    main()
