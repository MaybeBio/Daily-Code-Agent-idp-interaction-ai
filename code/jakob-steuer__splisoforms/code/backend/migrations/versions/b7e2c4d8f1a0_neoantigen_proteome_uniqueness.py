"""neoantigen_predictions — proteome-wide peptide uniqueness flags

Revision ID: b7e2c4d8f1a0
Revises: a4c2e7d9b1f3
Create Date: 2026-10-25

Adds `proteome_uniqueness` ('unique' | 'self_gene' | 'other_gene') and the derived
`is_proteome_unique` boolean. Both NULL until `backfill_neoantigen_uniqueness.py`
runs. A neoepitope is "novel relative to the canonical" by construction; these
columns record whether the peptide nevertheless occurs in some other GENCODE v47
reference protein (review F-4: 22 % of sampled P1 peptides did, within the same
gene) — the proteome-wide filter the preprint lists as planned.
"""
from alembic import op
import sqlalchemy as sa

revision = "b7e2c4d8f1a0"
down_revision = "a4c2e7d9b1f3"
branch_labels = None
depends_on = None


def upgrade():
    op.add_column("neoantigen_predictions",
                  sa.Column("proteome_uniqueness", sa.String(20), nullable=True))
    op.add_column("neoantigen_predictions",
                  sa.Column("is_proteome_unique", sa.Boolean(), nullable=True))
    op.create_index("ix_neoantigen_predictions_proteome_uniqueness",
                    "neoantigen_predictions", ["proteome_uniqueness"])
    op.create_index("ix_neoantigen_predictions_is_proteome_unique",
                    "neoantigen_predictions", ["is_proteome_unique"])


def downgrade():
    op.drop_index("ix_neoantigen_predictions_is_proteome_unique", table_name="neoantigen_predictions")
    op.drop_index("ix_neoantigen_predictions_proteome_uniqueness", table_name="neoantigen_predictions")
    op.drop_column("neoantigen_predictions", "is_proteome_unique")
    op.drop_column("neoantigen_predictions", "proteome_uniqueness")
