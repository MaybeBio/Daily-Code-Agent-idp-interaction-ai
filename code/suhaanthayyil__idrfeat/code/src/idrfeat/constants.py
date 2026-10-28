"""Amino acid alphabet and the residue groupings still computed locally.

Physicochemical scales (charge, hydropathy, disorder-promoting composition) are no longer
defined here; those values come from localCIDER. Only the amino acid alphabet and the aromatic
and polar groups remain, and the group fractions are summed from localCIDER's own amino acid
fractions.
"""

from __future__ import annotations

AA = "ACDEFGHIKLMNPQRSTVWY"
AA_SET = frozenset(AA)

AROMATIC = frozenset("FWY")
POLAR = frozenset("STNQCYH")  # polar uncharged plus histidine
