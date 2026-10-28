"""Deduplicate neoantigen_predictions and add unique constraint

Revision ID: c1d2e3f4a5b6
Revises: b9c3d7e2f1a0
Branch labels: None
Depends on: None

Removes duplicate (isoform_id, peptide_sequence, hla_allele, peptide_start) rows
that accumulated from multiple backfill runs, then adds a unique constraint to
prevent future duplicates.  The row with the highest ID (most recently written,
most likely to carry ensemble scores) is kept per group.
"""
from alembic import op
import sqlalchemy as sa

revision      = "c1d2e3f4a5b6"
down_revision = "b9c3d7e2f1a0"
branch_labels = None
depends_on    = None


def upgrade() -> None:
    # 1. Remove duplicates — keep the row with the highest id per unique key.
    #    Uses a single DELETE ... WHERE id NOT IN (SELECT max(id) ...) pattern
    #    that works on PostgreSQL without CTEs being required.
    op.execute("""
        DELETE FROM neoantigen_predictions
        WHERE id NOT IN (
            SELECT MAX(id)
            FROM neoantigen_predictions
            GROUP BY isoform_id, peptide_sequence, hla_allele, peptide_start
        )
    """)

    # 2. Add unique constraint to prevent future duplicates.
    op.create_unique_constraint(
        "uq_neoantigen_isoform_pep_allele_start",
        "neoantigen_predictions",
        ["isoform_id", "peptide_sequence", "hla_allele", "peptide_start"],
    )


def downgrade() -> None:
    op.drop_constraint(
        "uq_neoantigen_isoform_pep_allele_start",
        "neoantigen_predictions",
        type_="unique",
    )
