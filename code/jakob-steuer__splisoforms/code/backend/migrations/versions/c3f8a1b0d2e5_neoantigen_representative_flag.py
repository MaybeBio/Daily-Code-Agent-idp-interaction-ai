"""neoantigen is_representative flag — materialise the atlas dedup

Revision ID: c3f8a1b0d2e5
Revises: b7e4c1a95f2d
Create Date: 2026-07-31

Adds `is_representative` to neoantigen_predictions: the single best
(peptide_sequence, isoform_id) row among presented predictions. This
materialises the deduplication the /neoantigens/top atlas used to do with a
DISTINCT ON over ~7.6M rows on every request (page + count + stats), which cost
30 s+. The flag is populated by backfill_priority_tier.py (set-based) and lets
the endpoint filter an indexed ~2.6M-row subset instead.

A partial index on (priority_tier, consensus_rank) WHERE is_representative serves
the default priority sort and the count. NOTE: after this migration the atlas is
empty until backfill_priority_tier.py runs to set the flag (as with the ranking-v2
columns) — the coordinated rebuild runs that step.
"""
from alembic import op
import sqlalchemy as sa

revision = "c3f8a1b0d2e5"
down_revision = "b7e4c1a95f2d"
branch_labels = None
depends_on = None


def upgrade():
    op.add_column("neoantigen_predictions",
                  sa.Column("is_representative", sa.Boolean(), nullable=True))
    # Partial index: the deduped atlas set, ordered for the default priority sort.
    op.execute(
        "CREATE INDEX IF NOT EXISTS ix_neo_representative "
        "ON neoantigen_predictions (priority_tier, consensus_rank) WHERE is_representative"
    )
    # Sorted partial indexes so the atlas stats' COUNT(DISTINCT ...) run as index-only
    # scans (peptide-distinct 4.3s→0.1s, isoform-distinct 1.7s→0.1s over ~2.6M rows).
    op.execute(
        "CREATE INDEX IF NOT EXISTS ix_neo_repr_peptide "
        "ON neoantigen_predictions (peptide_sequence) WHERE is_representative"
    )
    op.execute(
        "CREATE INDEX IF NOT EXISTS ix_neo_repr_isoform "
        "ON neoantigen_predictions (isoform_id) WHERE is_representative"
    )
    # Clean up an ad-hoc index that may exist on dev DBs (never shipped in a migration).
    op.execute("DROP INDEX IF EXISTS ix_neo_dedup")


def downgrade():
    op.execute("DROP INDEX IF EXISTS ix_neo_repr_isoform")
    op.execute("DROP INDEX IF EXISTS ix_neo_repr_peptide")
    op.execute("DROP INDEX IF EXISTS ix_neo_representative")
    op.drop_column("neoantigen_predictions", "is_representative")
