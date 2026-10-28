#!/usr/bin/env python3
"""
finalize_reattribution.py — STEP 4/5 of the ground-truth gene re-attribution.

After isoforms have been re-homed (backfill_regene_from_ground_truth.py) and one
canonical per gene re-chosen (backfill_canonical_designation.py), this closes the
loop with two set-based passes:

  4. Recompute `gene.canonical_sequence` + `gene.canonical_length` (from the
     designated canonical isoform) and `isoform.delta_length` (= len(isoform) −
     canonical_length) for the WHOLE DB — all went stale the moment re-attribution
     moved isoforms to new canonicals. `canonical_sequence` is THE sequence-level
     baseline every canonical-relative analysis (Pfam diff, NMD truncation ratio,
     neoantigen novelty) compares against, so it must equal the designated canonical
     isoform — not a UniProt sequence that may correspond to no isoform (review
     ING-3: 130 genes, 119 of them novel long-read sequences, had drifted).
  5. Drop orphan stub genes — gene rows left with ZERO isoforms after re-homing
     (their `canonical_transcript_id` is nulled first so the NO-ACTION self-FK
     doesn't block the delete).

Designed to run in the INGEST stage — after core ingest+finalize+regene, BEFORE any
analysis or drug ingest — so the only gene-referencing rows are isoforms themselves
and the orphan delete is clean. If a later table already references an orphan gene
(a catch-up run on an analysed DB), the delete will raise and name the blocker
rather than corrupt state; re-run after a clean rebuild in that case.

Idempotent. READ-ONLY unless --apply.

Usage:
    python backend/scripts/finalize_reattribution.py            # dry-run
    python backend/scripts/finalize_reattribution.py --apply
"""
import argparse
import os
import sys

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from sqlalchemy import text
from app.database import SessionLocal


def main():
    ap = argparse.ArgumentParser(description="Finalize gene re-attribution")
    ap.add_argument("--apply", action="store_true", help="write changes (default: dry-run)")
    args = ap.parse_args()

    db = SessionLocal()
    orphan_pred = (
        "id NOT IN (SELECT gene_id FROM isoforms WHERE gene_id IS NOT NULL)"
    )
    n_orphans = db.execute(text(f"SELECT count(*) FROM genes WHERE {orphan_pred}")).scalar()
    n_genes = db.execute(text("SELECT count(*) FROM genes")).scalar()
    print(f"genes total                 : {n_genes}")
    print(f"orphan genes (0 isoforms)   : {n_orphans}")

    if not args.apply:
        print("\n(dry-run — re-run with --apply)")
        db.close()
        return

    # ── 4. recompute canonical_sequence + canonical_length + delta_length ───────
    db.execute(text("""
        UPDATE genes g
           SET canonical_sequence = sub.seq,
               canonical_length   = sub.len
          FROM (SELECT gene.id AS gid, iso.sequence AS seq, length(iso.sequence) AS len
                  FROM genes gene
                  JOIN isoforms iso ON iso.id = gene.canonical_transcript_id) sub
         WHERE g.id = sub.gid
    """))
    db.execute(text("""
        UPDATE isoforms i
           SET delta_length = length(i.sequence) - g.canonical_length
          FROM genes g
         WHERE i.gene_id = g.id
           AND g.canonical_length IS NOT NULL
           AND i.sequence IS NOT NULL
    """))
    db.commit()
    print("recomputed canonical_sequence + canonical_length + delta_length")

    # ── 5. drop orphan stub genes ──────────────────────────────────────────────
    db.execute(text(f"UPDATE genes SET canonical_transcript_id = NULL WHERE {orphan_pred}"))
    deleted = db.execute(text(f"DELETE FROM genes WHERE {orphan_pred}")).rowcount
    db.commit()
    print(f"dropped {deleted} orphan genes → {db.execute(text('SELECT count(*) FROM genes')).scalar()} genes remain")
    db.close()


if __name__ == "__main__":
    main()
