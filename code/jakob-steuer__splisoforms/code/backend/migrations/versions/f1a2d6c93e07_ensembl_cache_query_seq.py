"""ensembl_id_cache — store the fetched query protein (for live closest-isoform)

Revision ID: f1a2d6c93e07
Revises: e3c8f1a6b40d
Create Date: 2026-09-11

Adds the query protein sequence + md5 to ensembl_id_cache so the "most similar
isoform" line can be computed WITHOUT re-fetching from Ensembl each request.

Deliberately we cache only this EXTERNAL, stable data (the Ensembl protein for a
stable id). The DB-relative facts — is the id now present, which of our isoforms
is closest, and the identity % — are recomputed live against the current DB on
every request, so a future dataset that ingests the isoform is reflected
immediately with nothing stale to invalidate.

`protein_seq_checked` distinguishes "not fetched yet" from "Ensembl has no
protein for this id" (non-coding / no translation) so we don't re-fetch misses.
"""
from alembic import op
import sqlalchemy as sa

revision = "f1a2d6c93e07"
down_revision = "e3c8f1a6b40d"
branch_labels = None
depends_on = None


def upgrade():
    op.add_column("ensembl_id_cache", sa.Column("query_protein_seq", sa.Text(), nullable=True))
    op.add_column("ensembl_id_cache", sa.Column("query_protein_md5", sa.String(length=32), nullable=True))
    op.add_column("ensembl_id_cache",
                  sa.Column("protein_seq_checked", sa.Boolean(), nullable=False, server_default=sa.false()))
    op.create_index("ix_ensembl_id_cache_query_protein_md5", "ensembl_id_cache", ["query_protein_md5"])


def downgrade():
    op.drop_index("ix_ensembl_id_cache_query_protein_md5", table_name="ensembl_id_cache")
    op.drop_column("ensembl_id_cache", "protein_seq_checked")
    op.drop_column("ensembl_id_cache", "query_protein_md5")
    op.drop_column("ensembl_id_cache", "query_protein_seq")
