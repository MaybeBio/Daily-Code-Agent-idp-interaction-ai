"""transcript_metadata source-aware fields (sqanti3 vs gencode vs custom)

Revision ID: n7c2b5d9a4e8
Revises: m3f8a1c7e2b9
Create Date: 2026-06-02 14:30:00.000000

Adds nullable columns to transcript_metadata so the same table can cleanly
hold metadata from multiple annotation pipelines without overloading
`structural_category` with values that only make sense for SQANTI3
(historically the AF3 ingestion was dumping the predict-queue tier like
`tier1` into `structural_category`, and the ensembl-reference import was
dumping the literal string "ensembl_reference" there).

New columns (all nullable, so backfill is incremental):

  * metadata_kind  TEXT
      Pipeline this row came from. One of:
        sqanti3          — long-read SQANTI3 classification (FSM/ISM/NIC/NNC/...)
        gencode          — GENCODE/Ensembl reference annotation
        afdb_reference   — AlphaFold DB canonical / alternative imported on its own
        custom           — manual or other custom ingestion
        (NULL)           — legacy / pre-this-migration data

  * biotype  TEXT
      GFF biotype string. Universal field, populated whenever the source GFF
      provides it. Common values: protein_coding, lncRNA, processed_transcript,
      retained_intron, nonsense_mediated_decay.

  * appris_tag  TEXT
      GENCODE/APPRIS principal isoform tag. Values like 'principal:1',
      'alternative:1', etc. Only populated for metadata_kind == 'gencode'.

  * tsl  SMALLINT
      Transcript Support Level 1-5 (1 = best evidence, 5 = no single transcript
      supports the model). GENCODE-specific.

  * gencode_tags  TEXT[]
      Free-form tag set from a GENCODE / Ensembl GFF (e.g. ['basic', 'CCDS',
      'MANE_Select']). Stored as Postgres array so downstream queries can
      filter with ANY().

The old `structural_category` and `subcategory` columns are left as-is for
backwards-compatible reads of existing data. After backfill they should only
contain SQANTI3 categories (FSM / ISM / NIC / NNC / antisense / fusion /
genic / intergenic) and ingestion scripts have been updated to stop writing
non-SQANTI values to them.
"""
from typing import Sequence, Union

from alembic import op
import sqlalchemy as sa


revision: str = "n7c2b5d9a4e8"
down_revision: Union[str, Sequence[str], None] = "m3f8a1c7e2b9"
branch_labels: Union[str, Sequence[str], None] = None
depends_on: Union[str, Sequence[str], None] = None


def upgrade() -> None:
    op.add_column(
        "transcript_metadata",
        sa.Column("metadata_kind", sa.String(length=32), nullable=True),
    )
    op.create_index(
        "ix_transcript_metadata_metadata_kind",
        "transcript_metadata",
        ["metadata_kind"],
    )

    op.add_column(
        "transcript_metadata",
        sa.Column("biotype", sa.String(length=64), nullable=True),
    )
    op.create_index(
        "ix_transcript_metadata_biotype",
        "transcript_metadata",
        ["biotype"],
    )

    op.add_column(
        "transcript_metadata",
        sa.Column("appris_tag", sa.String(length=64), nullable=True),
    )

    op.add_column(
        "transcript_metadata",
        sa.Column("tsl", sa.SmallInteger(), nullable=True),
    )

    op.add_column(
        "transcript_metadata",
        sa.Column(
            "gencode_tags",
            sa.dialects.postgresql.ARRAY(sa.String(length=64)),
            nullable=True,
        ),
    )


def downgrade() -> None:
    op.drop_column("transcript_metadata", "gencode_tags")
    op.drop_column("transcript_metadata", "tsl")
    op.drop_column("transcript_metadata", "appris_tag")
    op.drop_index("ix_transcript_metadata_biotype", table_name="transcript_metadata")
    op.drop_column("transcript_metadata", "biotype")
    op.drop_index("ix_transcript_metadata_metadata_kind", table_name="transcript_metadata")
    op.drop_column("transcript_metadata", "metadata_kind")
