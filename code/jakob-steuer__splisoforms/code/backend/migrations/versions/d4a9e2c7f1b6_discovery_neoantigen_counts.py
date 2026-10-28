"""discovery index neoantigen counts — materialise Total / P1 for the Impact Explorer

Revision ID: d4a9e2c7f1b6
Revises: c3f8a1b0d2e5
Create Date: 2026-08-08

Adds `neoantigen_count` (total representative neoepitopes) and
`neoantigen_p1_count` (P1 tier) to isoform_discovery_index so the Impact Explorer
can display + sort those columns by reading two indexed integers, instead of
aggregating the ~7.98M-row neoantigen_predictions table at request time (the cost
that makes the Neoepitope Atlas slow). The one-time aggregation happens in the
backfill (backfill_discovery_index.py), not per request.

NOTE: after this migration the columns are 0 until the discovery-index backfill
runs to populate them — the coordinated rebuild / DB sync carries that step.
"""
from alembic import op
import sqlalchemy as sa

revision = "d4a9e2c7f1b6"
down_revision = "c3f8a1b0d2e5"
branch_labels = None
depends_on = None


def upgrade():
    op.add_column("isoform_discovery_index",
                  sa.Column("neoantigen_count", sa.Integer(), nullable=True, server_default="0"))
    op.add_column("isoform_discovery_index",
                  sa.Column("neoantigen_p1_count", sa.Integer(), nullable=True, server_default="0"))
    op.execute("CREATE INDEX IF NOT EXISTS ix_idx_neoantigen_count ON isoform_discovery_index (neoantigen_count)")
    op.execute("CREATE INDEX IF NOT EXISTS ix_idx_neoantigen_p1_count ON isoform_discovery_index (neoantigen_p1_count)")


def downgrade():
    op.execute("DROP INDEX IF EXISTS ix_idx_neoantigen_p1_count")
    op.execute("DROP INDEX IF EXISTS ix_idx_neoantigen_count")
    op.drop_column("isoform_discovery_index", "neoantigen_p1_count")
    op.drop_column("isoform_discovery_index", "neoantigen_count")
