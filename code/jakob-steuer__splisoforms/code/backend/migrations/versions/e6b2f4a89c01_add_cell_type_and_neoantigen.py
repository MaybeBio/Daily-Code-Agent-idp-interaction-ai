"""add cell_type_expression and neoantigen_predictions tables

Revision ID: e6b2f4a89c01
Revises: d5a1b3c72e89
Create Date: 2026-03-18

"""
from alembic import op
import sqlalchemy as sa

# revision identifiers, used by Alembic.
revision = "e6b2f4a89c01"
down_revision = "d5a1b3c72e89"
branch_labels = None
depends_on = None


def upgrade() -> None:
    # Add ensembl_gene_id column to genes table
    op.add_column("genes", sa.Column("ensembl_gene_id", sa.String(30), nullable=True))
    op.create_index("ix_genes_ensembl_gene_id", "genes", ["ensembl_gene_id"])

    # Cell-type expression table (gene-level, from HPA)
    op.create_table(
        "cell_type_expression",
        sa.Column("id", sa.Integer(), primary_key=True, autoincrement=True),
        sa.Column("gene_id", sa.Integer(), sa.ForeignKey("genes.id", ondelete="CASCADE"), nullable=False, index=True),
        sa.Column("cell_type", sa.String(100), nullable=False),
        sa.Column("expression_ncpm", sa.Float(), nullable=False),
        sa.Column("source", sa.String(50), nullable=False),
        sa.Column("specificity_category", sa.String(50), nullable=True),
    )
    op.create_index(
        "ix_cell_type_expression_gene_cell",
        "cell_type_expression",
        ["gene_id", "cell_type"],
        unique=True,
    )

    # Neoantigen predictions table (isoform-level)
    op.create_table(
        "neoantigen_predictions",
        sa.Column("id", sa.Integer(), primary_key=True, autoincrement=True),
        sa.Column("isoform_id", sa.Integer(), sa.ForeignKey("isoforms.id", ondelete="CASCADE"), nullable=False, index=True),
        sa.Column("peptide_sequence", sa.String(30), nullable=False),
        sa.Column("peptide_start", sa.Integer(), nullable=False),
        sa.Column("peptide_end", sa.Integer(), nullable=False),
        sa.Column("hla_allele", sa.String(30), nullable=False),
        sa.Column("binding_affinity_nm", sa.Float(), nullable=False),
        sa.Column("binding_category", sa.String(20), nullable=False),
        sa.Column("percentile_rank", sa.Float(), nullable=True),
        sa.Column("junction_type", sa.String(30), nullable=False),
        sa.Column("canonical_isoform_id", sa.Integer(), sa.ForeignKey("isoforms.id"), nullable=True),
    )
    op.create_index(
        "ix_neoantigen_isoform_hla",
        "neoantigen_predictions",
        ["isoform_id", "hla_allele"],
    )


def downgrade() -> None:
    op.drop_index("ix_neoantigen_isoform_hla", table_name="neoantigen_predictions")
    op.drop_table("neoantigen_predictions")
    op.drop_index("ix_cell_type_expression_gene_cell", table_name="cell_type_expression")
    op.drop_table("cell_type_expression")
    op.drop_index("ix_genes_ensembl_gene_id", table_name="genes")
    op.drop_column("genes", "ensembl_gene_id")
