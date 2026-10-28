"""add gene.protein_name and gene.function_summary

Revision ID: a7c4f9d2e1b3
Revises: 8073253a6f2b
Create Date: 2026-05-18 15:35:00.000000

Backfilled from UniProt REST proteinDescription.recommendedName.fullName
and the first entry of comments[type=FUNCTION].texts[0].value (see
scripts/backfill_gene_descriptions.py).  Both columns are nullable —
existing rows remain valid with NULL until the backfill runs.
"""
from typing import Sequence, Union

from alembic import op
import sqlalchemy as sa


revision: str = "a7c4f9d2e1b3"
down_revision: Union[str, Sequence[str], None] = "8073253a6f2b"
branch_labels: Union[str, Sequence[str], None] = None
depends_on: Union[str, Sequence[str], None] = None


def upgrade() -> None:
    op.add_column("genes", sa.Column("protein_name", sa.Text(), nullable=True))
    op.add_column("genes", sa.Column("function_summary", sa.Text(), nullable=True))


def downgrade() -> None:
    op.drop_column("genes", "function_summary")
    op.drop_column("genes", "protein_name")
