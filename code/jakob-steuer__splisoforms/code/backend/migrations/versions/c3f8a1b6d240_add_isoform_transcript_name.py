"""Add isoforms.transcript_name — a human-readable display name

`transcript_id` holds the SOURCE identifier (ENST00000397752, PB.1969727307.p1). The PacBio
IDs are opaque, so an earlier build (rename_to_gencode.py, now in legacy/) overwrote
transcript_id with the Ensembl transcript_name (MET-201) and gave novel isoforms a
`<closest>-NOV-<6char_md5>` name. The 2026-06 unified rebuild dropped that scheme and no
isoform carries such a name today.

Rather than rename transcript_id again, this adds a SEPARATE display column:
  * transcript_id   stays the stable source identifier — the API, exports, and the HPC
                    parquet/.done files all key on it, and it must not churn.
  * transcript_name is derived and recomputable after any re-ingest (see
                    backfill_transcript_names.py).

Nullable: novel isoforms with no canonical_match_id get no name, and the UI falls back to
transcript_id.

Revision ID: c3f8a1b6d240
Revises: b2e9f3c14d67
Create Date: 2026-07-14
"""
from typing import Sequence, Union

import sqlalchemy as sa
from alembic import op


revision: str = 'c3f8a1b6d240'
down_revision: Union[str, Sequence[str], None] = 'b2e9f3c14d67'
branch_labels: Union[str, Sequence[str], None] = None
depends_on: Union[str, Sequence[str], None] = None


def upgrade() -> None:
    op.add_column("isoforms", sa.Column("transcript_name", sa.String(255), nullable=True))
    op.create_index("ix_isoforms_transcript_name", "isoforms", ["transcript_name"])


def downgrade() -> None:
    op.drop_index("ix_isoforms_transcript_name", table_name="isoforms")
    op.drop_column("isoforms", "transcript_name")
