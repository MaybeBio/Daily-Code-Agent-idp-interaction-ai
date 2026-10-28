"""add exon_structural_impacts table

Revision ID: d5a1b3c72e89
Revises: c9d4e5f12a67
Create Date: 2026-03-18

"""
from alembic import op
import sqlalchemy as sa

# revision identifiers, used by Alembic.
revision = "d5a1b3c72e89"
down_revision = "c9d4e5f12a67"
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.create_table(
        "exon_structural_impacts",
        sa.Column("id", sa.Integer(), primary_key=True, autoincrement=True),
        sa.Column("isoform_id", sa.Integer(), sa.ForeignKey("isoforms.id", ondelete="CASCADE"), nullable=False, index=True),
        sa.Column("residue_number", sa.Integer(), nullable=False),
        sa.Column("dssp_assignment", sa.String(1), nullable=False),
        sa.Column("phi_angle", sa.Float(), nullable=True),
        sa.Column("psi_angle", sa.Float(), nullable=True),
        sa.Column("residue_plddt", sa.Float(), nullable=True),
        sa.Column("exon_number", sa.Integer(), nullable=True),
        sa.Column("canonical_dssp", sa.String(1), nullable=True),
        sa.Column("structural_change", sa.String(30), nullable=True),
        sa.Column("canonical_isoform_id", sa.Integer(), sa.ForeignKey("isoforms.id"), nullable=True),
    )
    # Composite index for per-isoform residue lookups
    op.create_index(
        "ix_exon_structural_impacts_isoform_residue",
        "exon_structural_impacts",
        ["isoform_id", "residue_number"],
    )


def downgrade() -> None:
    op.drop_index("ix_exon_structural_impacts_isoform_residue", table_name="exon_structural_impacts")
    op.drop_table("exon_structural_impacts")
