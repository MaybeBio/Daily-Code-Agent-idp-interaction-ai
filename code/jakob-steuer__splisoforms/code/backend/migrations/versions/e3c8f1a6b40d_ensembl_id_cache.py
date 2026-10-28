"""ensembl_id_cache — resolve unknown Ensembl IDs to their gene (cached)

Revision ID: e3c8f1a6b40d
Revises: d4a9e2c7f1b6
Create Date: 2026-09-11

Backs the search "not in database — did you mean…" resolution: when a user
searches an Ensembl ID (ENST/ENSP/ENSG) we didn't ingest, the resolver maps it
to its gene via the Ensembl REST API and caches the answer here, including
negative results (status='not_found') so repeated misses don't re-hit Ensembl.
`ensembl_id` is stored version-stripped (base accession).
"""
from alembic import op
import sqlalchemy as sa

revision = "e3c8f1a6b40d"
down_revision = "d4a9e2c7f1b6"
branch_labels = None
depends_on = None


def upgrade():
    op.create_table(
        "ensembl_id_cache",
        sa.Column("id", sa.Integer(), nullable=False),
        sa.Column("ensembl_id", sa.String(length=40), nullable=False),
        sa.Column("id_type", sa.String(length=10), nullable=True),
        sa.Column("resolved_ensg", sa.String(length=40), nullable=True),
        sa.Column("gene_symbol", sa.String(length=100), nullable=True),
        sa.Column("resolved_gene_id", sa.Integer(), nullable=True),
        sa.Column("status", sa.String(length=20), nullable=False, server_default="found"),
        sa.Column("checked_at", sa.DateTime(), nullable=True),
        sa.ForeignKeyConstraint(["resolved_gene_id"], ["genes.id"], ondelete="SET NULL"),
        sa.PrimaryKeyConstraint("id"),
    )
    op.create_index("ix_ensembl_id_cache_ensembl_id", "ensembl_id_cache", ["ensembl_id"], unique=True)
    op.create_index("ix_ensembl_id_cache_resolved_ensg", "ensembl_id_cache", ["resolved_ensg"])
    op.create_index("ix_ensembl_id_cache_resolved_gene_id", "ensembl_id_cache", ["resolved_gene_id"])
    op.create_index("ix_ensembl_id_cache_status", "ensembl_id_cache", ["status"])


def downgrade():
    op.drop_index("ix_ensembl_id_cache_status", table_name="ensembl_id_cache")
    op.drop_index("ix_ensembl_id_cache_resolved_gene_id", table_name="ensembl_id_cache")
    op.drop_index("ix_ensembl_id_cache_resolved_ensg", table_name="ensembl_id_cache")
    op.drop_index("ix_ensembl_id_cache_ensembl_id", table_name="ensembl_id_cache")
    op.drop_table("ensembl_id_cache")
