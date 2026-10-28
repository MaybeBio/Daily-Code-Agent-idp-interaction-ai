"""neoantigen ranking v2 — presentation_tier / priority_tier / is_immunogenic / is_expressed

Revision ID: b7e4c1a95f2d
Revises: a1f7c2e9d3b8
Create Date: 2026-07-24

Adds the ranking-v2 columns to neoantigen_predictions. Values are populated by
load_neoantigen_results.py at load time and can be recomputed by
backfill_priority_tier.py. See app/analysis/neoantigen_tiers.py.
"""
from alembic import op
import sqlalchemy as sa

revision = "b7e4c1a95f2d"
down_revision = "a1f7c2e9d3b8"
branch_labels = None
depends_on = None


def upgrade():
    op.add_column("neoantigen_predictions", sa.Column("presentation_tier", sa.String(length=10), nullable=True))
    op.add_column("neoantigen_predictions", sa.Column("priority_tier", sa.String(length=4), nullable=True))
    op.add_column("neoantigen_predictions", sa.Column("is_immunogenic", sa.Boolean(), nullable=True))
    op.add_column("neoantigen_predictions", sa.Column("is_expressed", sa.Boolean(), nullable=True))
    op.create_index("ix_neoantigen_predictions_presentation_tier", "neoantigen_predictions", ["presentation_tier"])
    op.create_index("ix_neoantigen_predictions_priority_tier", "neoantigen_predictions", ["priority_tier"])
    op.create_index("ix_neoantigen_predictions_is_immunogenic", "neoantigen_predictions", ["is_immunogenic"])


def downgrade():
    op.drop_index("ix_neoantigen_predictions_is_immunogenic", table_name="neoantigen_predictions")
    op.drop_index("ix_neoantigen_predictions_priority_tier", table_name="neoantigen_predictions")
    op.drop_index("ix_neoantigen_predictions_presentation_tier", table_name="neoantigen_predictions")
    op.drop_column("neoantigen_predictions", "is_expressed")
    op.drop_column("neoantigen_predictions", "is_immunogenic")
    op.drop_column("neoantigen_predictions", "priority_tier")
    op.drop_column("neoantigen_predictions", "presentation_tier")
