"""add pfam_domains table

Revision ID: a3f8b2d91e4c
Revises: cf66127306a5
Create Date: 2026-03-18
"""
from typing import Sequence, Union

from alembic import op
import sqlalchemy as sa


revision: str = 'a3f8b2d91e4c'
down_revision: Union[str, Sequence[str], None] = 'cf66127306a5'
branch_labels: Union[str, Sequence[str], None] = None
depends_on: Union[str, Sequence[str], None] = None


def upgrade() -> None:
    """Create pfam_domains table for storing Pfam domain annotations."""
    op.create_table(
        'pfam_domains',
        sa.Column('id', sa.Integer(), primary_key=True, autoincrement=True),
        sa.Column('isoform_id', sa.Integer(), sa.ForeignKey('isoforms.id', ondelete='CASCADE'), nullable=False),
        sa.Column('pfam_accession', sa.String(20), nullable=False),
        sa.Column('pfam_name', sa.String(100), nullable=False),
        sa.Column('pfam_description', sa.Text(), nullable=True),
        sa.Column('pfam_clan', sa.String(20), nullable=True),
        sa.Column('seq_start', sa.Integer(), nullable=False),
        sa.Column('seq_end', sa.Integer(), nullable=False),
        sa.Column('hmm_start', sa.Integer(), nullable=False),
        sa.Column('hmm_end', sa.Integer(), nullable=False),
        sa.Column('hmm_length', sa.Integer(), nullable=False),
        sa.Column('domain_coverage', sa.Float(), nullable=False),
        sa.Column('e_value', sa.Float(), nullable=False),
        sa.Column('bit_score', sa.Float(), nullable=False),
        sa.Column('domain_status', sa.String(20), nullable=False),
        sa.Column('is_kinase', sa.Boolean(), server_default='false'),
        sa.Column('canonical_isoform_id', sa.Integer(), sa.ForeignKey('isoforms.id'), nullable=True),
    )
    op.create_index('ix_pfam_domains_id', 'pfam_domains', ['id'])
    op.create_index('ix_pfam_domains_isoform_id', 'pfam_domains', ['isoform_id'])
    op.create_index('ix_pfam_domains_pfam_accession', 'pfam_domains', ['pfam_accession'])


def downgrade() -> None:
    """Drop pfam_domains table."""
    op.drop_index('ix_pfam_domains_pfam_accession', table_name='pfam_domains')
    op.drop_index('ix_pfam_domains_isoform_id', table_name='pfam_domains')
    op.drop_index('ix_pfam_domains_id', table_name='pfam_domains')
    op.drop_table('pfam_domains')
