"""
Populate the Isoform Discovery Index — the Structural Impact Score (SIS v2).

SIS v2 (2026-07 rework) is a **normalized 0–100 composite** of five bounded axes,
each mapped to [0, 1] and combined with fixed weights that sum to 100.  Every axis
is saturating, so no single event type (e.g. an isoform that loses 36 zinc fingers)
can dominate the ranking — the failure mode of the old unbounded additive score.

    axis                weight   sub-score s ∈ [0, 1]
    ─────────────────────────────────────────────────────────────────────────────
    Domain disruption     30     1 − 0.4^(lost + 0.5·gained)     (1 lost→0.60, 2→0.84)
    Length change         20     min(|Δlen| / 150, 1)            (150+ aa → full)
    Neoantigen priority   20     P1→1.0 P2→0.75 P3→0.5 P4→0.25 none→0   (ranking v2)
    NMD fate              15     target→1.0 suspect→0.5 else→0
    PTM accessibility     15     lost→1.0 else exposed/buried→0.6 else→0

    base = Σ weightᵢ · sᵢ            ∈ [0, 100]
    SIS  = base · (0.7 + 0.3·pLDDT/100)     — gentle confidence temper, factor ∈ [0.70, 1.00]

The neoantigen axis reads the ranking-v2 `priority_tier` (best/lowest-numbered tier per
isoform), NOT a raw peptide count — a P1 neoepitope contributes full weight, a P4 a
quarter.  `has_neoantigen` means "has a ranked/presented neoepitope" (priority_tier IS
NOT NULL), not "has any predicted peptide row".

The whole table is rebuilt with a single set-based UPSERT (seconds), replacing the old
per-isoform N+1 loop (~790k round-trips).  Construction is documented for users at
/documentation/discovery.
"""
import os
import sys
import argparse

from sqlalchemy import text

# Add backend directory to Python path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from app.database import SessionLocal

# ── SIS v2 weights (sum = 100) and sub-score parameters ───────────────────────
W_DOMAIN, W_LENGTH, W_NEO, W_NMD, W_PTM = 30, 20, 20, 15, 15
DOMAIN_DECAY = 0.4        # s_domain = 1 − DOMAIN_DECAY^(lost + 0.5·gained)
LENGTH_SAT_AA = 150.0     # s_length = min(|Δlen| / LENGTH_SAT_AA, 1)
CONF_FLOOR = 0.7          # confidence temper = CONF_FLOOR + (1−CONF_FLOOR)·pLDDT/100
PLDDT_DEFAULT = 70.0      # assumed pLDDT when structural_metrics is missing

# One set-based statement computes every component and UPSERTs the whole table.
# Sub-scores and the SIS are evaluated in SQL so it stays a single round-trip.
_SIS_SQL = f"""
WITH dom AS (
    SELECT isoform_id,
           count(*) FILTER (WHERE domain_status = 'lost')   AS lost,
           count(*) FILTER (WHERE domain_status = 'gained') AS gained
      FROM pfam_domains GROUP BY isoform_id
),
neo AS (
    SELECT isoform_id,
           min(CASE priority_tier WHEN 'P1' THEN 1 WHEN 'P2' THEN 2
                                  WHEN 'P3' THEN 3 WHEN 'P4' THEN 4 END) AS best_rank,
           count(*) FILTER (WHERE is_representative)                        AS neo_count,
           count(*) FILTER (WHERE is_representative AND priority_tier = 'P1') AS neo_p1_count
      FROM neoantigen_predictions
     WHERE priority_tier IS NOT NULL
     GROUP BY isoform_id
),
ptm AS (
    SELECT isoform_id,
           bool_or(accessibility_change = 'exposed') AS exposed,
           bool_or(accessibility_change = 'buried')  AS buried,
           bool_or(accessibility_change = 'lost')    AS lost
      FROM ptm_accessibility GROUP BY isoform_id
),
base AS (
    SELECT
        i.id                                     AS isoform_id,
        sm.plddt_mean                            AS plddt_mean,
        i.delta_length                           AS delta_length,
        COALESCE(dom.lost,   0)                  AS num_domains_lost,
        COALESCE(dom.gained, 0)                  AS num_domains_gained,
        (neo.best_rank IS NOT NULL)              AS has_neoantigen,
        COALESCE(neo.neo_count, 0)               AS neoantigen_count,
        COALESCE(neo.neo_p1_count, 0)            AS neoantigen_p1_count,
        COALESCE(ptm.exposed, false)             AS has_ptm_exposure_change,
        COALESCE(ptm.buried,  false)             AS has_ptm_buried_change,
        COALESCE(ptm.lost,    false)             AS has_ptm_loss,
        nmd.nmd_status                           AS nmd_status,
        -- ── sub-scores ∈ [0,1] ──
        (1 - power({DOMAIN_DECAY}, COALESCE(dom.lost,0) + 0.5*COALESCE(dom.gained,0)))
                                                 AS s_dom,
        least(abs(COALESCE(i.delta_length,0)) / {LENGTH_SAT_AA}, 1.0)
                                                 AS s_len,
        CASE neo.best_rank WHEN 1 THEN 1.0 WHEN 2 THEN 0.75
                           WHEN 3 THEN 0.5 WHEN 4 THEN 0.25 ELSE 0.0 END
                                                 AS s_neo,
        CASE nmd.nmd_status WHEN 'nmd_target' THEN 1.0
                            WHEN 'nmd_suspect' THEN 0.5 ELSE 0.0 END
                                                 AS s_nmd,
        CASE WHEN COALESCE(ptm.lost,false) THEN 1.0
             WHEN COALESCE(ptm.exposed,false) OR COALESCE(ptm.buried,false) THEN 0.6
             ELSE 0.0 END                        AS s_ptm
      FROM isoforms i
      LEFT JOIN structural_metrics sm ON sm.isoform_id = i.id
      LEFT JOIN dom ON dom.isoform_id = i.id
      LEFT JOIN neo ON neo.isoform_id = i.id
      LEFT JOIN ptm ON ptm.isoform_id = i.id
      LEFT JOIN nmd_predictions nmd ON nmd.isoform_id = i.id
      {{gene_filter}}
),
scored AS (
    SELECT base.*,
           round((( {W_DOMAIN}*s_dom + {W_LENGTH}*s_len + {W_NEO}*s_neo
                  + {W_NMD}*s_nmd + {W_PTM}*s_ptm )
                  * ({CONF_FLOOR} + {1 - CONF_FLOOR}*COALESCE(plddt_mean,{PLDDT_DEFAULT})/100.0)
                 )::numeric, 3) AS impact_score
      FROM base
)
INSERT INTO isoform_discovery_index
    (isoform_id, plddt_mean, delta_length, num_domains_lost, num_domains_gained,
     has_neoantigen, neoantigen_count, neoantigen_p1_count,
     has_ptm_exposure_change, has_ptm_buried_change, has_ptm_loss,
     nmd_status, impact_score, indexed_at)
SELECT isoform_id, plddt_mean, delta_length, num_domains_lost, num_domains_gained,
       has_neoantigen, neoantigen_count, neoantigen_p1_count,
       has_ptm_exposure_change, has_ptm_buried_change, has_ptm_loss,
       nmd_status, impact_score, now()
  FROM scored
ON CONFLICT (isoform_id) DO UPDATE SET
    plddt_mean              = EXCLUDED.plddt_mean,
    delta_length            = EXCLUDED.delta_length,
    num_domains_lost        = EXCLUDED.num_domains_lost,
    num_domains_gained      = EXCLUDED.num_domains_gained,
    has_neoantigen          = EXCLUDED.has_neoantigen,
    neoantigen_count        = EXCLUDED.neoantigen_count,
    neoantigen_p1_count     = EXCLUDED.neoantigen_p1_count,
    has_ptm_exposure_change = EXCLUDED.has_ptm_exposure_change,
    has_ptm_buried_change   = EXCLUDED.has_ptm_buried_change,
    has_ptm_loss            = EXCLUDED.has_ptm_loss,
    nmd_status              = EXCLUDED.nmd_status,
    impact_score            = EXCLUDED.impact_score,
    indexed_at              = now();
"""


def backfill_discovery_index(gene_symbols=None):
    db = SessionLocal()
    try:
        print("🚀 Discovery Index (SIS v2) backfill — set-based rebuild …")

        params = {}
        if gene_symbols:
            gene_filter = ("WHERE i.gene_id IN "
                           "(SELECT id FROM genes WHERE symbol = ANY(:symbols))")
            params["symbols"] = list(gene_symbols)
        else:
            gene_filter = ""

        result = db.execute(text(_SIS_SQL.format(gene_filter=gene_filter)), params)
        db.commit()
        print(f"✅ Upserted {result.rowcount:,} discovery-index rows.")

        # Audit summary — never exit 0 on a silently empty result.
        row = db.execute(text("""
            SELECT count(*)                                   AS n,
                   round(min(impact_score)::numeric, 2)       AS lo,
                   round(max(impact_score)::numeric, 2)       AS hi,
                   round(avg(impact_score)::numeric, 3)       AS mean,
                   count(*) FILTER (WHERE has_neoantigen)     AS neo,
                   count(*) FILTER (WHERE impact_score >= 60) AS severe,
                   count(*) FILTER (WHERE impact_score >= 40) AS high
              FROM isoform_discovery_index
        """)).one()
        print(f"   SIS range [{row.lo}, {row.hi}]  mean {row.mean}  "
              f"| has_neoantigen {row.neo:,}  | ≥40 {row.high:,}  ≥60 {row.severe:,}")
        if row.n == 0:
            print("❌ Discovery index is empty after backfill.")
            sys.exit(1)
        print("✅ Backfill Complete.")
    except Exception as e:
        db.rollback()
        print(f"❌ Error during backfill: {e}")
        sys.exit(1)
    finally:
        db.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Populate the Isoform Discovery Index (SIS v2)")
    parser.add_argument("--genes", type=str,
                        help="Comma-separated gene symbols to restrict the rebuild (default: all)")
    args = parser.parse_args()

    symbols = args.genes.split(",") if args.genes else None
    backfill_discovery_index(gene_symbols=symbols)
