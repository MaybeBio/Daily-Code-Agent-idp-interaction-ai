"""PTM accessibility redesign — add relative_delta_sasa, min_plddt, change_strength; introduce `undetermined` status

Revision ID: e4b8c1a7d9f2
Revises: d7f3e2a8b4c5
Create Date: 2026-05-19 14:00:00.000000

Adds three new columns on `ptm_accessibility` for the redesigned classifier:

  relative_delta_sasa    delta_sasa / canonical_sasa  (NULL when canonical_sasa
                         is 0 or undefined).  Better signal than the absolute
                         Δ alone because a 10 Å² shift on a buried (baseline
                         5 Å²) residue is enormous, on an exposed (baseline
                         150 Å²) residue is noise.

  min_plddt              min(canonical_plddt[site], novel_plddt[site]).  SASA
                         at low-confidence residues (pLDDT < 70) is unreliable;
                         the rewritten classifier marks those as `undetermined`
                         rather than `exposed` / `buried`.

  change_strength        [0, 1] magnitude of the relative change, capped at
                         1.0 (i.e. 100% loss of SASA = 1.0).  Drives a
                         continuous confidence bar in the UI.

Adds `undetermined` as a new accessibility_change value and relabels the
existing `unchanged-when-delta-is-null` cohort to surface the missing-data
case honestly.  Other relabels (pLDDT masking, relative-threshold flips)
require a recompute via backfill_ptm.py and are NOT applied here.
"""
from typing import Sequence, Union

from alembic import op
import sqlalchemy as sa


revision: str = "e4b8c1a7d9f2"
down_revision: Union[str, Sequence[str], None] = "d7f3e2a8b4c5"
branch_labels: Union[str, Sequence[str], None] = None
depends_on: Union[str, Sequence[str], None] = None


def upgrade() -> None:
    op.add_column("ptm_accessibility",
                  sa.Column("relative_delta_sasa", sa.Float(), nullable=True))
    op.add_column("ptm_accessibility",
                  sa.Column("min_plddt",            sa.Float(), nullable=True))
    op.add_column("ptm_accessibility",
                  sa.Column("change_strength",     sa.Float(), nullable=True))

    # Cheap relabel: rows where we have NO delta_sasa value were silently
    # labeled `unchanged`.  Move them to `undetermined` so the UI can
    # distinguish "we computed it and it didn't change" from "we couldn't
    # compute it".  No pLDDT-masking is applied here — that requires a
    # recompute via backfill_ptm.py --force.
    op.execute("""
        UPDATE ptm_accessibility
        SET accessibility_change = 'undetermined'
        WHERE accessibility_change = 'unchanged' AND delta_sasa IS NULL;
    """)


def downgrade() -> None:
    op.execute("""
        UPDATE ptm_accessibility
        SET accessibility_change = 'unchanged'
        WHERE accessibility_change = 'undetermined' AND delta_sasa IS NULL;
    """)
    op.drop_column("ptm_accessibility", "change_strength")
    op.drop_column("ptm_accessibility", "min_plddt")
    op.drop_column("ptm_accessibility", "relative_delta_sasa")
