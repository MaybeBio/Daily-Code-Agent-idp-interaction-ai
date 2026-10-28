"""Add molecular formula / weight / IUPAC name to drugs

The drug layer was ingested from the drug-dataset manifest, which carries only pubchem_cid,
broad_id and name — `smiles`, `inchikey` and `chembl_id` were declared but never populated,
and there was no formula or mass at all. The drug page therefore had no chemistry to show.

These are filled from PubChem by backfill_drug_chemistry.py (batched: one request per ~100
CIDs). `smiles` and `inchikey` already exist as columns and are populated by the same pass.

Revision ID: d5b2e9a740c1
Revises: c3f8a1b6d240
Create Date: 2026-07-14
"""
from typing import Sequence, Union

import sqlalchemy as sa
from alembic import op


revision: str = 'd5b2e9a740c1'
down_revision: Union[str, Sequence[str], None] = 'c3f8a1b6d240'
branch_labels: Union[str, Sequence[str], None] = None
depends_on: Union[str, Sequence[str], None] = None


def upgrade() -> None:
    op.add_column("drugs", sa.Column("molecular_formula", sa.String(128), nullable=True))
    op.add_column("drugs", sa.Column("molecular_weight", sa.Float(), nullable=True))
    op.add_column("drugs", sa.Column("iupac_name", sa.Text(), nullable=True))


def downgrade() -> None:
    op.drop_column("drugs", "iupac_name")
    op.drop_column("drugs", "molecular_weight")
    op.drop_column("drugs", "molecular_formula")
