"""Index the NO-ACTION FK columns that reference isoforms.id

Deleting isoforms validates every NO-ACTION FK that points at isoforms.id. Several
of those columns (the analysis tables' canonical_isoform_id, genes.canonical_transcript_id,
isoforms.canonical_match_id, ppi_networks.isoform_id_1, clinical_frequencies.isoform_id)
were unindexed, so each deleted row triggered full-table scans — an isoform cleanup
of ~5k rows took many minutes / effectively hung. These indexes make FK validation
(and any isoform delete / re-home) index lookups. Idempotent (IF NOT EXISTS).

Revision ID: b2e9f3c14d67
Revises: a1d4c7f0e2b9
Create Date: 2026-07-11
"""
from typing import Sequence, Union

from alembic import op


revision: str = 'b2e9f3c14d67'
down_revision: Union[str, Sequence[str], None] = 'a1d4c7f0e2b9'
branch_labels: Union[str, Sequence[str], None] = None
depends_on: Union[str, Sequence[str], None] = None

_INDEXES = [
    ("ix_isoforms_canonical_match_id", "isoforms", "canonical_match_id"),
    ("ix_genes_canonical_transcript_id", "genes", "canonical_transcript_id"),
    ("ix_clinfreq_isoform_id", "clinical_frequencies", "isoform_id"),
    ("ix_ppi_isoform_id_1", "ppi_networks", "isoform_id_1"),
    ("ix_pfam_canon_iso", "pfam_domains", "canonical_isoform_id"),
    ("ix_ptm_canon_iso", "ptm_accessibility", "canonical_isoform_id"),
    ("ix_nmd_canon_iso", "nmd_predictions", "canonical_isoform_id"),
    ("ix_exonimpact_canon_iso", "exon_structural_impacts", "canonical_isoform_id"),
    ("ix_neoantigen_canon_iso", "neoantigen_predictions", "canonical_isoform_id"),
]


def upgrade() -> None:
    for name, table, col in _INDEXES:
        op.execute(f"CREATE INDEX IF NOT EXISTS {name} ON {table} ({col})")


def downgrade() -> None:
    for name, _table, _col in _INDEXES:
        op.execute(f"DROP INDEX IF EXISTS {name}")
