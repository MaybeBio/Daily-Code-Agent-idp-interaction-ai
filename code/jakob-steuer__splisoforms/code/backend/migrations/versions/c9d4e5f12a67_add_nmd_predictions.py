"""add nmd_predictions table and cds_sequence column

Revision ID: c9d4e5f12a67
Revises: b7e2c4f08d31
Create Date: 2026-03-18
"""
from typing import Sequence, Union

from alembic import op
import sqlalchemy as sa


revision: str = 'c9d4e5f12a67'
down_revision: Union[str, Sequence[str], None] = 'b7e2c4f08d31'
branch_labels: Union[str, Sequence[str], None] = None
depends_on: Union[str, Sequence[str], None] = None


def upgrade() -> None:
    """Create nmd_predictions table and add cds_sequence to isoforms."""
    op.create_table(
        'nmd_predictions',
        sa.Column('id', sa.Integer(), primary_key=True, autoincrement=True),
        sa.Column('isoform_id', sa.Integer(), sa.ForeignKey('isoforms.id', ondelete='CASCADE'), nullable=False, unique=True),
        sa.Column('nmd_status', sa.String(20), nullable=False),
        sa.Column('has_ptc', sa.Boolean(), nullable=False),
        sa.Column('truncation_ratio', sa.Float(), nullable=True),
        sa.Column('orf_type', sa.String(50), nullable=True),
        sa.Column('ptc_position', sa.Integer(), nullable=True),
        sa.Column('last_ejc_distance', sa.Integer(), nullable=True),
        sa.Column('evidence', sa.Text(), nullable=False),
        sa.Column('canonical_isoform_id', sa.Integer(), sa.ForeignKey('isoforms.id'), nullable=True),
    )
    op.create_index('ix_nmd_predictions_id', 'nmd_predictions', ['id'])
    op.create_index('ix_nmd_predictions_isoform_id', 'nmd_predictions', ['isoform_id'], unique=True)

    # Add cds_sequence column to isoforms for future NMD 50-nt rule support
    op.add_column('isoforms', sa.Column('cds_sequence', sa.Text(), nullable=True))


def downgrade() -> None:
    """Drop nmd_predictions table and cds_sequence column."""
    op.drop_column('isoforms', 'cds_sequence')
    op.drop_index('ix_nmd_predictions_isoform_id', table_name='nmd_predictions')
    op.drop_index('ix_nmd_predictions_id', table_name='nmd_predictions')
    op.drop_table('nmd_predictions')
