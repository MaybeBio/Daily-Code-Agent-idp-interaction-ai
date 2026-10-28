"""add isoform_expression table

Revision ID: b2c3d4e5f6a1
Revises: a1b2c3d4e5f6
Branch labels: None
Depends on: None

Stores per-sample FL read counts and derived expression metrics (CPM,
isoform fraction ψ) extracted from the per-sample SQANTI3 CustomFiltered
classification files.  Only isoforms present in those files (~40% of the DB)
will have rows here; the remainder were not retained through the custom filter
and have no quantified expression.
"""
from alembic import op
import sqlalchemy as sa

revision     = "b2c3d4e5f6a1"
down_revision = "a1b2c3d4e5f6"
branch_labels = None
depends_on    = None


def upgrade() -> None:
    op.create_table(
        "isoform_expression",
        sa.Column("id",           sa.Integer(),     primary_key=True),
        sa.Column("isoform_id",   sa.Integer(),     sa.ForeignKey("isoforms.id", ondelete="CASCADE"),
                  nullable=False, index=True),
        sa.Column("sample",       sa.String(50),    nullable=False,
                  comment="Sample name: Normal | ccRCC2 | ccRCC3 | ccRCC4 | ccRCC5"),
        sa.Column("condition",    sa.String(20),    nullable=False,
                  comment="'normal' or 'tumor'"),
        sa.Column("fl_count",     sa.Integer(),     nullable=False,
                  comment="Raw Full-Length read count for this isoform in this sample"),
        sa.Column("cpm",          sa.Float(),       nullable=True,
                  comment="Counts Per Million: fl_count / total_sample_fl * 1e6"),
        sa.Column("isoform_frac", sa.Float(),       nullable=True,
                  comment="Isoform fraction ψ: fl_count / sum(fl for all isoforms of same gene)"),
        sa.Column("n_cells",      sa.Integer(),     nullable=True,
                  comment="Number of cells with ≥1 UMI for this isoform (from sparse matrix)"),
        sa.Column("cell_frac",    sa.Float(),       nullable=True,
                  comment="Fraction of cells in sample expressing this isoform"),
    )
    # Unique constraint: one row per isoform per sample
    op.create_unique_constraint(
        "uq_isoform_expression_isoform_sample",
        "isoform_expression",
        ["isoform_id", "sample"],
    )
    # Index for fast per-sample queries (e.g. top-expressed in ccRCC2)
    op.create_index("ix_isoform_expression_sample",       "isoform_expression", ["sample"])
    op.create_index("ix_isoform_expression_cpm",          "isoform_expression", ["cpm"])
    op.create_index("ix_isoform_expression_isoform_frac", "isoform_expression", ["isoform_frac"])


def downgrade() -> None:
    op.drop_index("ix_isoform_expression_isoform_frac", table_name="isoform_expression")
    op.drop_index("ix_isoform_expression_cpm",          table_name="isoform_expression")
    op.drop_index("ix_isoform_expression_sample",       table_name="isoform_expression")
    op.drop_constraint("uq_isoform_expression_isoform_sample", "isoform_expression", type_="unique")
    op.drop_table("isoform_expression")
