"""add_isoform_source_memberships

Revision ID: 8073253a6f2b
Revises: c1d2e3f4a5b6
Create Date: 2026-05-13 15:38:58.797688

"""
from typing import Sequence, Union

from alembic import op
import sqlalchemy as sa


# revision identifiers, used by Alembic.
revision: str = '8073253a6f2b'
down_revision: Union[str, Sequence[str], None] = 'c1d2e3f4a5b6'
branch_labels: Union[str, Sequence[str], None] = None
depends_on: Union[str, Sequence[str], None] = None


def upgrade() -> None:
    op.create_table(
        'isoform_source_memberships',
        sa.Column('id', sa.Integer(), nullable=False),
        sa.Column('isoform_id', sa.Integer(), nullable=False),
        sa.Column('source_id', sa.Integer(), nullable=False),
        sa.Column('in_normal', sa.Boolean(), nullable=False, server_default='false'),
        sa.Column('in_primary_tumor', sa.Boolean(), nullable=False, server_default='false'),
        sa.Column('in_metastasis', sa.Boolean(), nullable=False, server_default='false'),
        sa.ForeignKeyConstraint(['isoform_id'], ['isoforms.id'], ondelete='CASCADE'),
        sa.ForeignKeyConstraint(['source_id'], ['sources.id'], ondelete='CASCADE'),
        sa.PrimaryKeyConstraint('id'),
        sa.UniqueConstraint('isoform_id', 'source_id', name='uq_isoform_source'),
    )
    op.create_index('ix_isoform_source_memberships_id',       'isoform_source_memberships', ['id'])
    op.create_index('ix_isoform_source_memberships_isoform_id', 'isoform_source_memberships', ['isoform_id'])
    op.create_index('ix_isoform_source_memberships_source_id',  'isoform_source_memberships', ['source_id'])


def downgrade() -> None:
    op.drop_index('ix_isoform_source_memberships_source_id',  table_name='isoform_source_memberships')
    op.drop_index('ix_isoform_source_memberships_isoform_id', table_name='isoform_source_memberships')
    op.drop_index('ix_isoform_source_memberships_id',         table_name='isoform_source_memberships')
    op.drop_table('isoform_source_memberships')
