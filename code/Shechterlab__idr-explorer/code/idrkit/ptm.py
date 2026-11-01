"""How modifications change the charge of a sequence: the charge each modified residue carries, and the net charge and
charge patterning of a sequence with its modified residues.

One modification per residue, at pH 7.4, the modifying group monomeric:

  phospho      S, T, Y     -2   (phosphate, pKa2 about 5.6-6: mostly dianionic)
  acetyl       K            0   (also the neutral acyls: propionyl, butyryl, crotonyl, 2-hydroxyisobutyryl, lactyl,
                                 beta-hydroxybutyryl, benzoyl)
  acyl_acidic  K           -1   (succinyl, malonyl, glutaryl: the lysine's +1 replaced by a carboxylate)
  citrulline   R            0
  adp_ribosyl  S           -2   (the two phosphates of ADP-ribose)
               D, E        -2   (an ester of the carboxylate: its -1 replaced by the -2 of ADP-ribose)
                                (on R, K or H the charge of the modified side chain is not settled here: left out)
  glutamyl     E           -2   (one glutamate added to the side chain: each further one adds -1, so this is a floor)
  methyl       R, K        unchanged (methylarginines and methyllysines keep their charge)

Ubiquitin and SUMO are left out: they neutralise the lysine but attach a protein, which no sequence measure captures.
Charge patterning (delta) counts a residue as positive or negative by its sign, as localCIDER does, so a phosphoserine
joins the negative residues; net charge (NCPR) uses the values above.
"""
from __future__ import annotations

import numpy as np

from .charge import _delta_q, charge_vector

# charge of a modified residue, by modification and residue (None: the residue's own charge is kept)
MODIFIED_CHARGE: dict[str, dict[str, float | None]] = {
    "phospho": {"S": -2.0, "T": -2.0, "Y": -2.0},
    "acetyl": {"K": 0.0},
    "acyl_acidic": {"K": -1.0},
    "citrulline": {"R": 0.0},
    "adp_ribosyl": {"S": -2.0, "D": -2.0, "E": -2.0},
    "glutamyl": {"E": -2.0},
    "methyl": {"R": None, "K": None},
}


def modified_charge(seq: str, mods: dict[int, str]) -> np.ndarray:
    """Per-residue charge (charge_vector, histidine 0) with the residues in mods (0-based position -> modification of
    MODIFIED_CHARGE) set to their modified charge. A modification the residue cannot carry raises ValueError."""
    q = charge_vector(seq)
    for i, m in mods.items():
        table = MODIFIED_CHARGE.get(m)
        if table is None:
            raise ValueError(f"unknown modification {m!r}")
        res = seq[i]
        if res not in table:
            raise ValueError(f"{m} on {res}{i + 1}")
        if table[res] is not None:
            q[i] = table[res]
    return q


def charge_shift(seq: str, mods: dict[int, str]) -> dict[str, float]:
    """Net charge per residue, fraction of charged residues and delta of seq unmodified and with every modification in
    mods made. delta is Das & Pappu's (idrkit.delta), on the sign of each residue's charge."""
    q0, q1 = charge_vector(seq), modified_charge(seq, mods)
    return {"ncpr0": float(q0.mean()), "ncpr1": float(q1.mean()),
            "fcr0": float(np.count_nonzero(q0) / len(q0)), "fcr1": float(np.count_nonzero(q1) / len(q1)),
            "delta0": _delta_q(q0), "delta1": _delta_q(q1)}
