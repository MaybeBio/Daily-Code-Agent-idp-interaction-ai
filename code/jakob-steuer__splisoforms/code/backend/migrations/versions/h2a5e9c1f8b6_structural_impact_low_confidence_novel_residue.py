"""structural_impact: add `low_confidence` + `novel_residue` to structural_change vocabulary; relabel low-pLDDT rows

Revision ID: h2a5e9c1f8b6
Revises: g8e1d4c7b2a9
Create Date: 2026-05-19 18:00:00.000000

Audit-4 of the per-residue DSSP structural impact pipeline raised three
specific issues addressable without a full pyDSSP re-run:

  1. structural_change H/E/C transition calls are emitted regardless of
     pLDDT confidence.  27% of residues have pLDDT < 50 and 38% are
     below the conventional 70 cutoff — many of the ~1.7M "changes" are
     noise from low-confidence side-chain placement.

  2. There is no way to distinguish "novel residue with no canonical
     counterpart" from "secondary-structure rearrangement of an aligned
     residue".  The old pipeline force-aligned novel exon residues against
     unrelated canonical residues and reported coil_to_helix etc.

This migration introduces two new vocabulary values:

  low_confidence    structural_change is uninterpretable because at least
                    one of the canonical / novel pLDDT is < 70.

  novel_residue     residue has no canonical counterpart (no alignment
                    match) — splicing introduced a new sequence segment.
                    Set by the rewritten analyze_isoform on new writes.

The migration relabels existing rows in two cheap SQL passes:

  • rows with H/E/C transition labels AND residue_plddt < 70 → low_confidence
  • novel_residue is NOT backfilled here — requires the splice-aware
    aligner to re-run (user does this separately).

No new columns are added — we reuse structural_change.  H/E/C values
remain the same; `same` rows are not touched (a `same` call between two
low-pLDDT residues is rarely informative anyway, but we'd rather under-
flag than over-flag).
"""
from typing import Sequence, Union

from alembic import op
import sqlalchemy as sa


revision: str = "h2a5e9c1f8b6"
down_revision: Union[str, Sequence[str], None] = "g8e1d4c7b2a9"
branch_labels: Union[str, Sequence[str], None] = None
depends_on: Union[str, Sequence[str], None] = None


def upgrade() -> None:
    # Relabel: any H/E/C transition where residue_plddt < 70 → low_confidence.
    # canonical_plddt isn't stored per-row in this table, so we can only mask
    # using the novel-side pLDDT.  The full bilateral mask (min(can, nov) < 70)
    # is applied at recompute time by the updated analyze_isoform.
    op.execute("""
        UPDATE exon_structural_impacts
        SET structural_change = 'low_confidence'
        WHERE structural_change IN (
            'coil_to_helix', 'helix_to_coil',
            'coil_to_sheet', 'sheet_to_coil',
            'helix_to_sheet', 'sheet_to_helix'
        )
        AND residue_plddt IS NOT NULL
        AND residue_plddt < 70;
    """)


def downgrade() -> None:
    # Cannot reverse: we don't store the original transition direction once
    # we've labeled it low_confidence.  Downgrade is a no-op; downstream
    # code that expects only the H/E/C-pair labels should be updated.
    pass
