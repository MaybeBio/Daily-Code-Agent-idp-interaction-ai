"""NMD redesign — orthogonal truncation_severity + nmd_confidence + smORF/uORF candidate flags; relabel heuristic rows

Revision ID: d7f3e2a8b4c5
Revises: c3e9a4b7d2f1
Create Date: 2026-05-19 12:30:00.000000

Splits the previous single nmd_status column into:

  nmd_status               nmd_target | nmd_escape | non_translatable |
                           predicted_translatable | unknown
                           (50-nt-rule rows keep nmd_target / nmd_escape;
                            heuristic rows are relabeled `predicted_translatable`
                            when they showed no PTC signal, `unknown` otherwise;
                            the old `nmd_suspect` label is retired — it conflated
                            protein truncation with NMD evidence.)

  truncation_severity      intact | minor | severe | non_translatable
                           (new column, orthogonal signal derived from
                            truncation_ratio.  Surfaces protein-length
                            divergence WITHOUT claiming NMD evidence.)

  nmd_confidence           float in [0, 1], only set when nmd_status='nmd_target'.
                           Linear ramp from 0 (distance=50, threshold) to 1.0
                           (distance>=200) using last_ejc_distance.

  is_candidate_smorf       Bool — short isoform (<100 aa) with a short or
                           absent canonical (<150 aa).  Reflects the TransCODE
                           Consortium's `peptidein` framework (Nature 2026) —
                           these isoforms may encode bona fide microproteins
                           rather than truncated artifacts.

  is_candidate_uorf        Bool — short complete-ORF isoform (<50 aa) with a
                           long canonical (>150 aa).  Heuristic flag; uORFs
                           are NMD-resistant via ribosome reinitiation and
                           the standard 50-nt rule does not apply.

The migration relabels existing rows in place — no nmd_prediction.py re-run
is strictly required, though running the analysis with TransDecoder data
afterwards will populate nmd_confidence and apply the 5'-proximal escape
rule.
"""
from typing import Sequence, Union

from alembic import op
import sqlalchemy as sa


revision: str = "d7f3e2a8b4c5"
down_revision: Union[str, Sequence[str], None] = "c3e9a4b7d2f1"
branch_labels: Union[str, Sequence[str], None] = None
depends_on: Union[str, Sequence[str], None] = None


def upgrade() -> None:
    # ── 0. Widen nmd_status from VARCHAR(20) to VARCHAR(30) so the longer
    #       `predicted_translatable` value fits.
    op.alter_column("nmd_predictions", "nmd_status",
                    existing_type=sa.String(20),
                    type_=sa.String(30),
                    existing_nullable=False)

    # ── 1. Add the four new columns ───────────────────────────────────────────
    op.add_column("nmd_predictions",
                  sa.Column("truncation_severity", sa.String(30), nullable=True))
    op.add_column("nmd_predictions",
                  sa.Column("nmd_confidence",      sa.Float(),     nullable=True))
    op.add_column("nmd_predictions",
                  sa.Column("is_candidate_smorf", sa.Boolean(), server_default=sa.false(), nullable=False))
    op.add_column("nmd_predictions",
                  sa.Column("is_candidate_uorf",  sa.Boolean(), server_default=sa.false(), nullable=False))
    op.create_index("ix_nmd_truncation_severity",  "nmd_predictions", ["truncation_severity"])
    op.create_index("ix_nmd_candidate_smorf",      "nmd_predictions", ["is_candidate_smorf"])

    # ── 2. Backfill truncation_severity from existing truncation_ratio ───────
    #
    #   intact            ratio >= 0.90   OR ratio IS NULL (no canonical to compare)
    #   minor             0.70 <= ratio < 0.90
    #   severe            ratio < 0.70
    #   non_translatable  current nmd_status = 'non_translatable'
    op.execute("""
        UPDATE nmd_predictions
        SET truncation_severity = CASE
            WHEN nmd_status = 'non_translatable' THEN 'non_translatable'
            WHEN truncation_ratio IS NULL        THEN 'intact'
            WHEN truncation_ratio < 0.70         THEN 'severe'
            WHEN truncation_ratio < 0.90         THEN 'minor'
            ELSE 'intact'
        END;
    """)

    # ── 3. Backfill nmd_confidence for 50-nt-rule rows ───────────────────────
    #   linear ramp:  conf = (last_ejc_distance - 50) / 150 clamped to [0, 1]
    op.execute("""
        UPDATE nmd_predictions
        SET nmd_confidence = LEAST(1.0, GREATEST(0.0, (last_ejc_distance - 50)::float / 150))
        WHERE nmd_status = 'nmd_target' AND last_ejc_distance IS NOT NULL;
    """)

    # ── 4. Tag smORF / uORF candidates via isoform + gene join ────────────────
    #   smORF candidate: novel length < 100 AND canonical_length IS NULL or < 150
    op.execute("""
        UPDATE nmd_predictions np
        SET is_candidate_smorf = TRUE
        FROM isoforms i
        LEFT JOIN genes g ON g.id = i.gene_id
        WHERE np.isoform_id = i.id
          AND LENGTH(i.sequence) < 100
          AND (g.canonical_length IS NULL OR g.canonical_length < 150);
    """)

    #   uORF candidate: novel length < 50 AND canonical_length >= 150
    op.execute("""
        UPDATE nmd_predictions np
        SET is_candidate_uorf = TRUE
        FROM isoforms i
        JOIN genes g ON g.id = i.gene_id
        WHERE np.isoform_id = i.id
          AND LENGTH(i.sequence) < 50
          AND g.canonical_length >= 150;
    """)

    # ── 5. Relabel nmd_status for heuristic rows ──────────────────────────────
    #
    # The evidence prefix tells us which method produced the row.  We only
    # touch the heuristic ones; 50-nt-rule and CDS-sequence rows keep their
    # labels (they're real evidence).
    #
    #   heuristic + nmd_suspect              → unknown  (truncation_severity carries the info)
    #   heuristic + nmd_escape (has_ptc=t)   → unknown  (we lied about has_ptc)
    #   heuristic + nmd_escape (has_ptc=f)   → predicted_translatable
    #   heuristic + normal                   → predicted_translatable
    op.execute("""
        UPDATE nmd_predictions
        SET nmd_status = CASE
            WHEN nmd_status = 'nmd_suspect'                              THEN 'unknown'
            WHEN nmd_status = 'nmd_escape'  AND has_ptc = TRUE           THEN 'unknown'
            WHEN nmd_status = 'nmd_escape'  AND has_ptc = FALSE          THEN 'predicted_translatable'
            WHEN nmd_status = 'normal'                                   THEN 'predicted_translatable'
            ELSE nmd_status
        END
        WHERE substring(evidence FROM 1 FOR 11) = '[heuristic]';
    """)

    # Clear the misleading has_ptc=TRUE on heuristic rows — we don't actually
    # know where the stop codon is from a protein-length comparison.
    op.execute("""
        UPDATE nmd_predictions
        SET has_ptc = FALSE
        WHERE substring(evidence FROM 1 FOR 11) = '[heuristic]' AND has_ptc = TRUE;
    """)

    # smORF candidates: NMD rule doesn't apply, mark unknown
    op.execute("""
        UPDATE nmd_predictions
        SET nmd_status = 'unknown'
        WHERE is_candidate_smorf = TRUE AND nmd_status NOT IN ('non_translatable');
    """)


def downgrade() -> None:
    op.drop_index("ix_nmd_candidate_smorf",     table_name="nmd_predictions")
    op.drop_index("ix_nmd_truncation_severity", table_name="nmd_predictions")
    op.drop_column("nmd_predictions", "is_candidate_uorf")
    op.drop_column("nmd_predictions", "is_candidate_smorf")
    op.drop_column("nmd_predictions", "nmd_confidence")
    op.drop_column("nmd_predictions", "truncation_severity")
