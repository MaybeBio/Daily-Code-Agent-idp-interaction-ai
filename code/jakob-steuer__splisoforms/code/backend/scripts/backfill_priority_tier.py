#!/usr/bin/env python3
"""
backfill_priority_tier.py — (re)compute neoepitope ranking v2 over existing rows.

The loader (load_neoantigen_results.py) already sets presentation_tier / priority_tier /
is_immunogenic / is_expressed at load time. This script recomputes them IN PLACE from the
stored tool scores — so you can re-tier the whole table (e.g. after tuning
BIGMHC_IM_THRESH in app/analysis/neoantigen_tiers.py) WITHOUT a full HPC reload.

Set-based SQL (fast over millions of rows). Thresholds come from the tiers module so
the SQL and the Python loader can never drift. Supersedes backfill_gold_medal_universal.py
(the old splice/pLDDT-gated "Gold Medal" tier).

Usage:
    python backend/scripts/backfill_priority_tier.py            # apply
    python backend/scripts/backfill_priority_tier.py --dry-run  # counts only
"""
import argparse
import os
import sys

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from sqlalchemy import text
from app.database import SessionLocal
from app.analysis import neoantigen_tiers as T


def main():
    ap = argparse.ArgumentParser(description="Recompute neoepitope ranking v2")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()
    db = SessionLocal()

    # per-tool "presented" vote counts (NULL tools don't vote); thresholds from the module
    # NB: the mhcflurry percentile is stored in the `percentile_rank` column
    # (the loader maps percentile_rank = parquet's mhcflurry_percentile).
    n_strong = (
        f"((CASE WHEN percentile_rank IS NOT NULL AND percentile_rank <= {T.STRONG_RANK} THEN 1 ELSE 0 END)"
        f" + (CASE WHEN netmhcpan_rank_el IS NOT NULL AND netmhcpan_rank_el <= {T.STRONG_RANK} THEN 1 ELSE 0 END)"
        f" + (CASE WHEN bigmhc_el_score IS NOT NULL AND bigmhc_el_score >= {T.BIGMHC_EL_STRONG} THEN 1 ELSE 0 END))"
    )
    n_weak = (
        f"((CASE WHEN percentile_rank IS NOT NULL AND percentile_rank <= {T.WEAK_RANK} THEN 1 ELSE 0 END)"
        f" + (CASE WHEN netmhcpan_rank_el IS NOT NULL AND netmhcpan_rank_el <= {T.WEAK_RANK} THEN 1 ELSE 0 END)"
        f" + (CASE WHEN bigmhc_el_score IS NOT NULL AND bigmhc_el_score >= {T.BIGMHC_EL_WEAK} THEN 1 ELSE 0 END))"
    )
    pres = (
        f"CASE WHEN {n_strong} >= 2 THEN 'STRONG'"
        f" WHEN {n_strong} = 1 OR {n_weak} >= 2 THEN 'MODERATE'"
        f" WHEN {n_weak} = 1 THEN 'WEAK' ELSE 'NONE' END"
    )
    immuno = f"(bigmhc_im_score IS NOT NULL AND bigmhc_im_score >= {T.BIGMHC_IM_THRESH})"

    total = db.execute(text("SELECT count(*) FROM neoantigen_predictions")).scalar()
    if total == 0:
        print("neoantigen_predictions is empty — nothing to tier.")
        db.close(); return
    print(f"rows: {total:,}")

    if args.dry_run:
        dist = db.execute(text(f"""
            WITH t AS (SELECT {pres} AS pt, {immuno} AS im FROM neoantigen_predictions)
            SELECT pt, im, count(*) FROM t GROUP BY pt, im ORDER BY pt, im
        """)).fetchall()
        for pt, im, n in dist:
            print(f"  presentation={pt:9s} immunogenic={im!s:5s} → {n:,}")
        print("(dry-run — no writes)")
        db.close(); return

    # 1. presentation_tier + is_immunogenic
    db.execute(text(f"""
        UPDATE neoantigen_predictions
           SET presentation_tier = {pres},
               is_immunogenic     = {immuno}
    """))
    # 2. priority_tier from the two axes (P1–P4; NULL when not presented)
    db.execute(text("""
        UPDATE neoantigen_predictions
           SET priority_tier = CASE
               WHEN presentation_tier = 'STRONG'   AND is_immunogenic THEN 'P1'
               WHEN presentation_tier = 'STRONG'                      THEN 'P2'
               WHEN presentation_tier = 'MODERATE' AND is_immunogenic THEN 'P3'
               WHEN presentation_tier = 'MODERATE'                    THEN 'P4'
               WHEN presentation_tier = 'WEAK'                        THEN 'P4'
               ELSE NULL END
    """))
    # 3. is_expressed (decorating flag) — isoform detected in expression data
    db.execute(text("""
        UPDATE neoantigen_predictions n
           SET is_expressed = EXISTS (
               SELECT 1 FROM isoform_expression e
                WHERE e.isoform_id = n.isoform_id AND e.fl_count > 0)
    """))
    # 4. is_representative — materialises the atlas dedup: the single best row per
    #    (peptide_sequence, isoform_id) among PRESENTED predictions (best tier, then
    #    best consensus_rank). /neoantigens/top filters this flag instead of DISTINCT ON.
    db.execute(text("UPDATE neoantigen_predictions SET is_representative = FALSE WHERE is_representative"))
    db.execute(text("""
        WITH ranked AS (
            SELECT id, row_number() OVER (
                       PARTITION BY peptide_sequence, isoform_id
                       ORDER BY CASE priority_tier WHEN 'P1' THEN 1 WHEN 'P2' THEN 2
                                                   WHEN 'P3' THEN 3 WHEN 'P4' THEN 4 ELSE 5 END ASC,
                                consensus_rank ASC NULLS LAST, id ASC) AS rn
              FROM neoantigen_predictions
             WHERE priority_tier IS NOT NULL)
        UPDATE neoantigen_predictions n
           SET is_representative = TRUE
          FROM ranked r
         WHERE n.id = r.id AND r.rn = 1
    """))
    db.commit()

    dist = db.execute(text(
        "SELECT priority_tier, count(*) FROM neoantigen_predictions "
        "GROUP BY priority_tier ORDER BY priority_tier NULLS LAST")).fetchall()
    print("priority_tier distribution:")
    for pt, n in dist:
        print(f"  {pt or 'NONE (not presented)':24s} {n:,}")
    repr_n = db.execute(text(
        "SELECT count(*) FROM neoantigen_predictions WHERE is_representative")).scalar()
    print(f"representative (deduped atlas) rows: {repr_n:,}")
    db.close()


if __name__ == "__main__":
    main()
