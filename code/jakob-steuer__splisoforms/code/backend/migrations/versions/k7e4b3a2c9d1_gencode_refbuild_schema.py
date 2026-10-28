"""gencode reference build — Phase 1: Isoform Ensembl-aligned columns + prediction_jobs table

Revision ID: k7e4b3a2c9d1
Revises: j9d0c4f6b8a1
Create Date: 2026-05-20 11:00:00.000000

Adds the columns needed for the GENCODE v47 reference build (design notes
in the chat around 2026-05-20).  Phase 1 is schema-only — Phases 2-5 populate
the new fields via Python scripts.

New columns on `isoforms`
-------------------------
  ensembl_transcript_id     ENST00000... (no version), nullable, indexed.
                            Set when the isoform matches a GENCODE transcript
                            (either canonical or via sequence-md5 dedup).
  ensembl_protein_id        ENSP00000... (no version), nullable.
  gencode_version           e.g. 'v47'; tags which release this row was
                            last reconciled against.
  legacy_transcript_id      The original transcript_id (PB.* or AFDB_*)
                            preserved verbatim for traceability.
  canonical_match_id        FK to isoforms.id; for novel rows, the closest
                            GENCODE-anchored isoform by sequence identity.
  protein_md5               md5 of the protein sequence — drives dedup
                            and cross-version matching.
  is_mane_select            Boolean.  Promoted from transcript_metadata.
  gencode_match_identity    Float in [0, 1]; the alignment identity of a
                            novel row to its canonical_match.
  structure_provenance      e.g. 'afdb_legacy', 'af3_local_v1',
                            'af3_local_2026Q2'.  Lets downstream consumers
                            know which model produced the stored structure.
  ms_evidence               Boolean; populated by Phase 4 from PeptideAtlas /
                            TransCODE.

New table: prediction_jobs
--------------------------
Tracks the lifecycle of every AF3 prediction we ask Leonardo to run.
States: queued → submitted → predicting → predicted → ingested → failed/skipped.

Backfill done in this migration (all cheap SQL)
-----------------------------------------------
  * legacy_transcript_id = transcript_id (every row)
  * protein_md5 = md5(sequence) where sequence is non-null
  * structure_provenance =
        'afdb_legacy'     when transcript_id LIKE 'AFDB_%'
        'af3_local_v1'    when transcript_id LIKE 'PB.%'
        NULL              otherwise
  * is_mane_select promoted from transcript_metadata.is_mane_select where set
"""
from typing import Sequence, Union

from alembic import op
import sqlalchemy as sa


revision: str = "k7e4b3a2c9d1"
down_revision: Union[str, Sequence[str], None] = "j9d0c4f6b8a1"
branch_labels: Union[str, Sequence[str], None] = None
depends_on: Union[str, Sequence[str], None] = None


def upgrade() -> None:
    # ── 1. New columns on isoforms ────────────────────────────────────────────
    op.add_column("isoforms", sa.Column("ensembl_transcript_id",  sa.String(20),  nullable=True))
    op.add_column("isoforms", sa.Column("ensembl_protein_id",     sa.String(20),  nullable=True))
    op.add_column("isoforms", sa.Column("gencode_version",        sa.String(10),  nullable=True))
    op.add_column("isoforms", sa.Column("legacy_transcript_id",   sa.String(100), nullable=True))
    op.add_column("isoforms", sa.Column(
        "canonical_match_id",
        sa.Integer(),
        sa.ForeignKey("isoforms.id", name="fk_isoforms_canonical_match", use_alter=True),
        nullable=True,
    ))
    op.add_column("isoforms", sa.Column("protein_md5",            sa.String(32),  nullable=True))
    op.add_column("isoforms", sa.Column("is_mane_select",         sa.Boolean(),   nullable=True))
    op.add_column("isoforms", sa.Column("gencode_match_identity", sa.Float(),     nullable=True))
    op.add_column("isoforms", sa.Column("structure_provenance",   sa.String(50),  nullable=True))
    op.add_column("isoforms", sa.Column("ms_evidence",            sa.Boolean(),   nullable=False, server_default=sa.false()))

    op.create_index("ix_isoforms_ensembl_transcript_id",  "isoforms", ["ensembl_transcript_id"])
    op.create_index("ix_isoforms_legacy_transcript_id",   "isoforms", ["legacy_transcript_id"])
    op.create_index("ix_isoforms_protein_md5",            "isoforms", ["protein_md5"])
    op.create_index("ix_isoforms_is_mane_select",         "isoforms", ["is_mane_select"])
    op.create_index("ix_isoforms_structure_provenance",   "isoforms", ["structure_provenance"])

    # ── 2. Cheap backfill of the new fields from existing data ──────────────
    # legacy_transcript_id mirrors the current transcript_id (preserved verbatim).
    op.execute("UPDATE isoforms SET legacy_transcript_id = transcript_id;")

    # protein_md5 — PostgreSQL has md5() native.  Apply to non-empty sequences.
    op.execute("""
        UPDATE isoforms
        SET protein_md5 = md5(sequence)
        WHERE sequence IS NOT NULL AND sequence <> '' AND sequence <> 'UNKNOWN';
    """)

    # structure_provenance — tag existing rows by their origin convention.
    op.execute("""
        UPDATE isoforms SET structure_provenance =
            CASE
                WHEN substring(transcript_id FROM 1 FOR 5) = 'AFDB_' THEN 'afdb_legacy'
                WHEN substring(transcript_id FROM 1 FOR 3) = 'PB.'  THEN 'af3_local_v1'
                ELSE NULL
            END;
    """)

    # is_mane_select promotion from transcript_metadata — true where the
    # metadata row says so; leave NULL otherwise (the GENCODE rename pass
    # will set TRUE/FALSE for every reconciled isoform).
    op.execute("""
        UPDATE isoforms i
        SET is_mane_select = TRUE
        FROM transcript_metadata tm
        WHERE tm.isoform_id = i.id AND tm.is_mane_select = TRUE;
    """)

    # ── 3. prediction_jobs ──────────────────────────────────────────────────
    op.create_table(
        "prediction_jobs",
        sa.Column("id", sa.Integer(), primary_key=True),
        # Target identity (this is what the eventual Isoform.transcript_id
        # will be once the prediction lands and is ingested).
        sa.Column("target_id", sa.String(100), nullable=False, unique=True, index=True),
        sa.Column("protein_md5", sa.String(32), nullable=False, index=True),
        sa.Column("protein_sequence", sa.Text(), nullable=False),
        sa.Column("protein_length", sa.Integer(), nullable=False),
        # Priority + scope tags
        sa.Column("priority_tier", sa.SmallInteger(), nullable=False, server_default="2"),
        sa.Column("queued_for", sa.String(50), nullable=True),
        # Workflow state
        sa.Column("status", sa.String(20), nullable=False, server_default="queued", index=True),
        # HPC bookkeeping
        sa.Column("hpc_job_id", sa.String(100), nullable=True),
        sa.Column("submitted_at", sa.DateTime(), nullable=True),
        sa.Column("completed_at", sa.DateTime(), nullable=True),
        sa.Column("output_path", sa.String(511), nullable=True),
        sa.Column("error_message", sa.Text(), nullable=True),
        sa.Column("gencode_version", sa.String(10), nullable=True),
        sa.Column("created_at", sa.DateTime(), server_default=sa.func.now(), nullable=False),
        sa.Column("updated_at", sa.DateTime(), server_default=sa.func.now(), nullable=False),
    )
    op.create_index("ix_prediction_jobs_priority_tier", "prediction_jobs", ["priority_tier"])
    op.create_index("ix_prediction_jobs_md5_status",    "prediction_jobs", ["protein_md5", "status"])


def downgrade() -> None:
    op.drop_table("prediction_jobs")

    op.drop_index("ix_isoforms_structure_provenance",  table_name="isoforms")
    op.drop_index("ix_isoforms_is_mane_select",        table_name="isoforms")
    op.drop_index("ix_isoforms_protein_md5",           table_name="isoforms")
    op.drop_index("ix_isoforms_legacy_transcript_id",  table_name="isoforms")
    op.drop_index("ix_isoforms_ensembl_transcript_id", table_name="isoforms")

    op.drop_column("isoforms", "ms_evidence")
    op.drop_column("isoforms", "structure_provenance")
    op.drop_column("isoforms", "gencode_match_identity")
    op.drop_column("isoforms", "is_mane_select")
    op.drop_column("isoforms", "protein_md5")
    op.drop_column("isoforms", "canonical_match_id")
    op.drop_column("isoforms", "legacy_transcript_id")
    op.drop_column("isoforms", "gencode_version")
    op.drop_column("isoforms", "ensembl_protein_id")
    op.drop_column("isoforms", "ensembl_transcript_id")
