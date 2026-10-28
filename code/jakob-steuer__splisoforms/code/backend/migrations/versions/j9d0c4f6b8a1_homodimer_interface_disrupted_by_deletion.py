"""homodimer_analysis: add interface_disrupted_by_deletion flag

Revision ID: j9d0c4f6b8a1
Revises: i6f1a3d8e5b2
Create Date: 2026-05-20 09:00:00.000000

Audit-9 finding: the homodimer disruption_score conflates two biologies —

  (a) Splicing physically removes the interface region (large `count_absent`,
      disruption_score → 1.0).  This is "the interface isn't there" — not the
      same biological signal as a conformationally-disrupted interface.

  (b) Splicing perturbs the interface conformationally (high count_sasa_changed,
      moderate count_absent).  This is what the panel is really trying to surface.

A consumer ranking by `disruption_score` alone treats (a) and (b) identically.
This migration adds an `interface_disrupted_by_deletion` boolean: TRUE when at
least 50 % of the canonical interface residues are absent in the novel isoform.
Downstream consumers (Impact Explorer, UI badges) can then de-prioritise the
deletion cohort or surface it separately.

No new analysis is required — count_absent / total_interface_residues are
already stored, so the backfill is a pure SQL UPDATE.
"""
from typing import Sequence, Union

from alembic import op
import sqlalchemy as sa


revision: str = "j9d0c4f6b8a1"
down_revision: Union[str, Sequence[str], None] = "i6f1a3d8e5b2"
branch_labels: Union[str, Sequence[str], None] = None
depends_on: Union[str, Sequence[str], None] = None


def upgrade() -> None:
    op.add_column("homodimer_analysis",
                  sa.Column("interface_disrupted_by_deletion",
                            sa.Boolean(),
                            server_default=sa.false(), nullable=False))
    op.create_index("ix_homodimer_disrupted_by_deletion",
                    "homodimer_analysis", ["interface_disrupted_by_deletion"])

    # Backfill: a row is "deletion-disrupted" when at least half of the
    # canonical interface residues are absent in the novel.  Skip rows where
    # total_interface_residues is NULL or zero (can't compute the ratio).
    op.execute("""
        UPDATE homodimer_analysis
        SET interface_disrupted_by_deletion = TRUE
        WHERE total_interface_residues IS NOT NULL
          AND total_interface_residues > 0
          AND count_absent IS NOT NULL
          AND count_absent * 2 >= total_interface_residues;
    """)


def downgrade() -> None:
    op.drop_index("ix_homodimer_disrupted_by_deletion", table_name="homodimer_analysis")
    op.drop_column("homodimer_analysis", "interface_disrupted_by_deletion")
