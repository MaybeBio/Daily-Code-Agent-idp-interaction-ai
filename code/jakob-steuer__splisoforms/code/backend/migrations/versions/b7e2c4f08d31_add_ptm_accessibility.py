"""add ptm_accessibility table

Revision ID: b7e2c4f08d31
Revises: a3f8b2d91e4c
Create Date: 2026-03-18
"""
from typing import Sequence, Union

from alembic import op
import sqlalchemy as sa


revision: str = 'b7e2c4f08d31'
down_revision: Union[str, Sequence[str], None] = 'a3f8b2d91e4c'
branch_labels: Union[str, Sequence[str], None] = None
depends_on: Union[str, Sequence[str], None] = None


def upgrade() -> None:
    """Create ptm_accessibility table for storing PTM site SASA comparisons."""
    op.create_table(
        'ptm_accessibility',
        sa.Column('id', sa.Integer(), primary_key=True, autoincrement=True),
        sa.Column('isoform_id', sa.Integer(), sa.ForeignKey('isoforms.id', ondelete='CASCADE'), nullable=False),
        sa.Column('uniprot_id', sa.String(50), nullable=False),
        sa.Column('residue_position', sa.Integer(), nullable=False),
        sa.Column('novel_residue_position', sa.Integer(), nullable=True),
        sa.Column('residue_type', sa.String(5), nullable=False),
        sa.Column('ptm_type', sa.String(50), nullable=False),
        sa.Column('source_db', sa.String(50), nullable=False),
        sa.Column('canonical_sasa', sa.Float(), nullable=True),
        sa.Column('novel_sasa', sa.Float(), nullable=True),
        sa.Column('delta_sasa', sa.Float(), nullable=True),
        sa.Column('accessibility_change', sa.String(20), nullable=False),
        sa.Column('canonical_isoform_id', sa.Integer(), sa.ForeignKey('isoforms.id'), nullable=True),
    )
    op.create_index('ix_ptm_accessibility_id', 'ptm_accessibility', ['id'])
    op.create_index('ix_ptm_accessibility_isoform_id', 'ptm_accessibility', ['isoform_id'])
    op.create_index('ix_ptm_accessibility_uniprot_pos', 'ptm_accessibility', ['uniprot_id', 'residue_position'])


def downgrade() -> None:
    """Drop ptm_accessibility table."""
    op.drop_index('ix_ptm_accessibility_uniprot_pos', table_name='ptm_accessibility')
    op.drop_index('ix_ptm_accessibility_isoform_id', table_name='ptm_accessibility')
    op.drop_index('ix_ptm_accessibility_id', table_name='ptm_accessibility')
    op.drop_table('ptm_accessibility')
