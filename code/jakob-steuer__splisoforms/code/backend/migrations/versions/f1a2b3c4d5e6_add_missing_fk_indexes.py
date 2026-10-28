"""Add missing FK indexes for AF3 ingestion performance

Adds indexes on high-cardinality foreign key columns that were missing:
- isoforms.gene_id: used in all gene-to-isoform joins
- isoforms.source_id: used in source filtering
- structures.isoform_id: used in structure lookups
- genes.hgnc_symbol: used in search queries

These are critical for 1M+ row performance.

Revision ID: f1a2b3c4d5e6
Revises: bdbb59b74816
Create Date: 2026-04-08 14:00:00.000000

"""
from typing import Sequence, Union

from alembic import op

# revision identifiers, used by Alembic.
revision: str = 'f1a2b3c4d5e6'
down_revision: Union[str, Sequence[str], None] = 'bdbb59b74816'
branch_labels: Union[str, Sequence[str], None] = None
depends_on: Union[str, Sequence[str], None] = None


def upgrade() -> None:
    op.create_index('ix_isoforms_gene_id', 'isoforms', ['gene_id'])
    op.create_index('ix_isoforms_source_id', 'isoforms', ['source_id'])
    op.create_index('ix_structures_isoform_id', 'structures', ['isoform_id'])
    op.create_index('ix_genes_hgnc_symbol', 'genes', ['hgnc_symbol'])


def downgrade() -> None:
    op.drop_index('ix_genes_hgnc_symbol', table_name='genes')
    op.drop_index('ix_structures_isoform_id', table_name='structures')
    op.drop_index('ix_isoforms_source_id', table_name='isoforms')
    op.drop_index('ix_isoforms_gene_id', table_name='isoforms')
