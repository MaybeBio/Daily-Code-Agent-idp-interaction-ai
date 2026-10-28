"""neoantigen: fix NetMHCpan affinity_nm / rank_ba column swap; re-derive ba_rescue

Revision ID: g8e1d4c7b2a9
Revises: f5c2b8d6a3e4
Create Date: 2026-05-19 16:30:00.000000

Audit of the neoantigen pipeline (Audit 3) surfaced a parser bug:
backend/app/analysis/netmhcpan_runner.py:147-148 reads
  affinity_nm = parts[14]
  rank_ba     = parts[15]
but NetMHCpan-4.2's combined -BA output orders columns as:
  Pos MHC Peptide Core Of Gp Gl Ip Il Icore Identity Score_EL %Rank_EL Score_BA %Rank_BA Affinity(nM) BindLevel
so the correct mapping is:
  parts[14] = %Rank_BA   (range 0-100, the percentile rank)
  parts[15] = Affinity(nM) (range ~1 to 50000, the IC50)

Live data confirms the swap — netmhcpan_affinity_nm currently ranges
0-99 (a percentile) and netmhcpan_rank_ba ranges 1-49,141 (an IC50).
This migration swaps the values in place and re-derives the ba_rescue
flag using the corrected %Rank_BA.

The parser is fixed separately in this commit; future writes will be
correct without further migration.
"""
from typing import Sequence, Union

from alembic import op
import sqlalchemy as sa


revision: str = "g8e1d4c7b2a9"
down_revision: Union[str, Sequence[str], None] = "f5c2b8d6a3e4"
branch_labels: Union[str, Sequence[str], None] = None
depends_on: Union[str, Sequence[str], None] = None


def upgrade() -> None:
    # ── Atomic swap of the two columns ────────────────────────────────────────
    # PostgreSQL evaluates the RHS before assigning, so `SET a = b, b = a`
    # swaps the values atomically per row.
    op.execute("""
        UPDATE neoantigen_predictions
        SET netmhcpan_affinity_nm = netmhcpan_rank_ba,
            netmhcpan_rank_ba     = netmhcpan_affinity_nm
        WHERE netmhcpan_affinity_nm IS NOT NULL OR netmhcpan_rank_ba IS NOT NULL;
    """)

    # ── Re-apply the ba_rescue rule with the now-correctly-labeled %Rank_BA ──
    # Old rule (in f5c2b8d6a3e4) fired 0 times because it was running against
    # the IC50 column.  Same logic, now against the real percentile rank:
    #   rank_ba < 1.0  AND  rank_el > 2 * rank_ba  AND  ≥ 2 cysteines.
    op.execute("UPDATE neoantigen_predictions SET is_ba_rescue = FALSE;")
    op.execute("""
        UPDATE neoantigen_predictions
        SET is_ba_rescue = TRUE
        WHERE netmhcpan_rank_ba IS NOT NULL
          AND netmhcpan_rank_ba < 1.0
          AND netmhcpan_rank_el IS NOT NULL
          AND netmhcpan_rank_el > 2.0 * netmhcpan_rank_ba
          AND char_length(peptide_sequence) - char_length(replace(peptide_sequence, 'C', '')) >= 2;
    """)


def downgrade() -> None:
    # Swap back
    op.execute("""
        UPDATE neoantigen_predictions
        SET netmhcpan_affinity_nm = netmhcpan_rank_ba,
            netmhcpan_rank_ba     = netmhcpan_affinity_nm
        WHERE netmhcpan_affinity_nm IS NOT NULL OR netmhcpan_rank_ba IS NOT NULL;
    """)
