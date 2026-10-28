"""neoantigen Option A — tools_agreeing + agreement_tier, ba_rescue rewrite, ptm_canonical_residue_lost split

Revision ID: f5c2b8d6a3e4
Revises: e4b8c1a7d9f2
Create Date: 2026-05-19 16:00:00.000000

Adds two new columns to expose **inter-tool agreement** as a first-class
signal (the audit found only 5.1% of predictions agree across MHCflurry,
NetMHCpan and BigMHC, but the panel previously surfaced no agreement
metric), rewrites `is_ba_rescue` (currently dead — 0 rows fire), and adds
a flag distinguishing canonical-PTM-residue-retained vs canonical-PTM-
residue-lost neoantigens.

New columns
-----------
  tools_agreeing            smallint (0-3) — count of in-agreement tools:
                              MHCflurry %rank < 0.5
                              NetMHCpan %rank_el < 0.5
                              BigMHC EL score > 0.5

  agreement_tier            'high' (3) | 'medium' (2) | 'low' (0-1)

  ptm_canonical_residue_lost  boolean — for Direct PTM-Neoantigen rows,
                              True when the canonical PTM residue is
                              recorded as `lost` in ptm_accessibility for
                              this same isoform.  Splits the conflated
                              "Direct PTM-Neoantigen" label into the two
                              biologies (PTM site retained in novel
                              context vs PTM site abolished).

Rule rewrites
-------------
  is_ba_rescue              old rule required NetMHCpan-BA rank < 2.0 AND
                              ≥2 cysteines AND BOTH EL tools >2.0 — the
                              both-EL-miss conjunction makes the rule
                              never fire on the live corpus.  Replaced
                              with: BA rank < 1.0 AND NetMHCpan EL >
                              2× BA rank AND ≥2 cysteines.  Signals
                              tool-disagreement on a cysteine-rich
                              peptide rather than total EL miss.

No HMMER / MHCflurry re-scan required — all backfill is from existing
columns plus a join to ptm_accessibility.
"""
from typing import Sequence, Union

from alembic import op
import sqlalchemy as sa


revision: str = "f5c2b8d6a3e4"
down_revision: Union[str, Sequence[str], None] = "e4b8c1a7d9f2"
branch_labels: Union[str, Sequence[str], None] = None
depends_on: Union[str, Sequence[str], None] = None


def upgrade() -> None:
    op.add_column("neoantigen_predictions",
                  sa.Column("tools_agreeing", sa.SmallInteger(), nullable=True))
    op.add_column("neoantigen_predictions",
                  sa.Column("agreement_tier", sa.String(10), nullable=True))
    op.add_column("neoantigen_predictions",
                  sa.Column("ptm_canonical_residue_lost", sa.Boolean(),
                            server_default=sa.false(), nullable=False))
    op.create_index("ix_neoantigen_agreement_tier",
                    "neoantigen_predictions", ["agreement_tier"])
    op.create_index("ix_neoantigen_tools_agreeing",
                    "neoantigen_predictions", ["tools_agreeing"])

    # ── Backfill tools_agreeing ───────────────────────────────────────────────
    # Count how many of (MHCflurry, NetMHCpan, BigMHC) call this peptide as a
    # strong binder under each tool's natural threshold.
    op.execute("""
        UPDATE neoantigen_predictions
        SET tools_agreeing =
          (CASE WHEN percentile_rank   IS NOT NULL AND percentile_rank   < 0.5 THEN 1 ELSE 0 END)
        + (CASE WHEN netmhcpan_rank_el IS NOT NULL AND netmhcpan_rank_el < 0.5 THEN 1 ELSE 0 END)
        + (CASE WHEN bigmhc_el_score   IS NOT NULL AND bigmhc_el_score   > 0.5 THEN 1 ELSE 0 END);
    """)
    op.execute("""
        UPDATE neoantigen_predictions
        SET agreement_tier = CASE
            WHEN tools_agreeing >= 3 THEN 'high'
            WHEN tools_agreeing  = 2 THEN 'medium'
            ELSE 'low'
        END;
    """)

    # ── Rewrite is_ba_rescue with the relaxed rule ────────────────────────────
    # Old rule had BOTH EL tools required to miss (>2.0) — never fires.
    # New rule: NetMHCpan-BA is strong (rank < 1.0) AND NetMHCpan-EL is
    # more than 2× looser than BA (signals genuine BA-EL disagreement),
    # AND the peptide is cysteine-rich.
    #
    # Cysteine count via length(peptide) - length(replace(peptide, 'C', '')).
    op.execute("""
        UPDATE neoantigen_predictions
        SET is_ba_rescue = (
            netmhcpan_rank_ba IS NOT NULL
            AND netmhcpan_rank_ba < 1.0
            AND netmhcpan_rank_el IS NOT NULL
            AND netmhcpan_rank_el > 2.0 * netmhcpan_rank_ba
            AND char_length(peptide_sequence) - char_length(replace(peptide_sequence, 'C', '')) >= 2
        );
    """)

    # ── Backfill ptm_canonical_residue_lost ───────────────────────────────────
    # For predictions tagged Direct PTM-Neoantigen, set the flag when at least
    # one canonical PTM site within the peptide window is recorded as `lost`
    # in ptm_accessibility for this isoform — i.e. the splice event physically
    # removed the modifiable residue.
    op.execute("""
        UPDATE neoantigen_predictions np
        SET ptm_canonical_residue_lost = TRUE
        WHERE np.ptm_impact_status = 'Direct PTM-Neoantigen'
          AND EXISTS (
              SELECT 1 FROM ptm_accessibility pa
              WHERE pa.isoform_id = np.isoform_id
                AND pa.residue_position BETWEEN np.peptide_start AND np.peptide_end
                AND pa.accessibility_change = 'lost'
          );
    """)


def downgrade() -> None:
    op.drop_index("ix_neoantigen_tools_agreeing", table_name="neoantigen_predictions")
    op.drop_index("ix_neoantigen_agreement_tier", table_name="neoantigen_predictions")
    op.drop_column("neoantigen_predictions", "ptm_canonical_residue_lost")
    op.drop_column("neoantigen_predictions", "agreement_tier")
    op.drop_column("neoantigen_predictions", "tools_agreeing")
