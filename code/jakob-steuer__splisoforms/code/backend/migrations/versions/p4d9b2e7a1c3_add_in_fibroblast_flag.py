"""add in_fibroblast tissue flag to isoforms + memberships

Revision ID: p4d9b2e7a1c3
Revises: n7c2b5d9a4e8
Create Date: 2026-06-26 18:40:00.000000

The breast-cancer dataset (Head 2026) detects isoforms in three tissue tracks:
healthy breast, tumor, and cultured fibroblast. The original BC ingest
repurposed `in_metastasis` for fibroblast — a hack the unified-ingest cleanup
removes. This adds a dedicated `in_fibroblast` boolean to both `isoforms`
(OR-merged fast-lookup view) and `isoform_source_memberships` (per-dataset
detection), so fibroblast is kept distinct from `in_normal` (healthy tissue)
and `in_metastasis` stays semantically correct (False for BC).

Note: a clean rebuild via reset_db.py recreates the schema from the models and
does not need this migration; it exists so environments that migrate forward
(rather than reset) stay in sync.
"""
from typing import Sequence, Union

from alembic import op
import sqlalchemy as sa


revision: str = "p4d9b2e7a1c3"
down_revision: Union[str, Sequence[str], None] = "n7c2b5d9a4e8"
branch_labels: Union[str, Sequence[str], None] = None
depends_on: Union[str, Sequence[str], None] = None


def upgrade() -> None:
    op.add_column(
        "isoforms",
        sa.Column("in_fibroblast", sa.Boolean(), nullable=True,
                  server_default=sa.false()),
    )
    op.add_column(
        "isoform_source_memberships",
        sa.Column("in_fibroblast", sa.Boolean(), nullable=False,
                  server_default=sa.false()),
    )


def downgrade() -> None:
    op.drop_column("isoform_source_memberships", "in_fibroblast")
    op.drop_column("isoforms", "in_fibroblast")
