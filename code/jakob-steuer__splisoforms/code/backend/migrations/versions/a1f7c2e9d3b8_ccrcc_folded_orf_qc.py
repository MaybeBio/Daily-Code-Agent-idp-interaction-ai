"""ccRCC folded-ORF QC fields (integrity + transparency)

Adds three transcript_metadata columns that flag ccRCC transcripts whose folded
ORF (= the AF3 structure's sequence) has a protein identity diverging from the
transcript's genomic-locus gene, and record the folded protein's best UniProt
match (provenance rules in backfill_ccrcc_orf_flag.py).

Revision ID: a1f7c2e9d3b8
Revises: d5b2e9a740c1
Create Date: 2026-07-17
"""
from alembic import op
import sqlalchemy as sa

revision = "a1f7c2e9d3b8"
down_revision = "d5b2e9a740c1"
branch_labels = None
depends_on = None


def upgrade():
    op.add_column("transcript_metadata",
                  sa.Column("folded_orf_mismatch", sa.Boolean(), server_default=sa.false(), nullable=True))
    op.add_column("transcript_metadata",
                  sa.Column("folded_best_ref", sa.String(length=100), nullable=True))
    op.add_column("transcript_metadata",
                  sa.Column("folded_best_ref_identity", sa.Float(), nullable=True))
    op.create_index("ix_transcript_metadata_folded_orf_mismatch",
                    "transcript_metadata", ["folded_orf_mismatch"])


def downgrade():
    op.drop_index("ix_transcript_metadata_folded_orf_mismatch", table_name="transcript_metadata")
    op.drop_column("transcript_metadata", "folded_best_ref_identity")
    op.drop_column("transcript_metadata", "folded_best_ref")
    op.drop_column("transcript_metadata", "folded_orf_mismatch")
