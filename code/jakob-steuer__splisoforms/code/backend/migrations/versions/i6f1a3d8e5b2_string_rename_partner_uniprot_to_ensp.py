"""string_interactions: rename partner_uniprot_id → partner_ensp_id

Revision ID: i6f1a3d8e5b2
Revises: h2a5e9c1f8b6
Create Date: 2026-05-19 19:00:00.000000

Audit-5 finding: the `partner_uniprot_id` column never contained a UniProt
accession — _run_string_analysis() in app/services/ppi_service.py reads
STRING's `stringId_B`/`stringId_A` (format `9606.ENSP00000…`) and splits on
'.' which yields the Ensembl protein ID (ENSP*), not the UniProt accession.
The data in the column is correct as-is (Ensembl protein IDs); only the
column name is misleading.  This migration is a pure rename — no data
transformation.

Downstream consumers expecting a UniProt accession in this column have
always been getting an ENSP*.  Renaming to `partner_ensp_id` makes the
mismatch impossible to ignore.
"""
from typing import Sequence, Union

from alembic import op
import sqlalchemy as sa


revision: str = "i6f1a3d8e5b2"
down_revision: Union[str, Sequence[str], None] = "h2a5e9c1f8b6"
branch_labels: Union[str, Sequence[str], None] = None
depends_on: Union[str, Sequence[str], None] = None


def upgrade() -> None:
    op.alter_column("string_interactions", "partner_uniprot_id",
                    new_column_name="partner_ensp_id")


def downgrade() -> None:
    op.alter_column("string_interactions", "partner_ensp_id",
                    new_column_name="partner_uniprot_id")
