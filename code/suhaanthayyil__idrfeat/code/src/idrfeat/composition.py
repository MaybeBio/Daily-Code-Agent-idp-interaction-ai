"""Amino acid composition and hydropathy features from localCIDER.

Amino acid fractions, grouped fractions, the fraction of disorder-promoting residues, and mean
hydropathy are all read from localCIDER's SequenceParameters (Holehouse et al., 2017). Mean
hydropathy is on localCIDER's normalized Kyte-Doolittle scale (0 to 9). Grouped fractions that
localCIDER does not expose directly (aromatic, polar) are summed from its own amino acid
fractions, so every number still traces to a single localCIDER computation.
"""

from __future__ import annotations

from .constants import AA, AROMATIC, POLAR
from .physchem import params


def _empty_composition() -> dict[str, float]:
    feats = {f"aa_frac_{aa}": 0.0 for aa in AA}
    feats.update(
        frac_charged=0.0,
        frac_polar=0.0,
        frac_aromatic=0.0,
        frac_proline=0.0,
        frac_glycine=0.0,
        frac_disorder_promoting=0.0,
    )
    return feats


def composition_features(seq: str) -> dict[str, float]:
    """Per-residue and grouped composition fractions, all from localCIDER."""
    if not seq:
        return _empty_composition()
    sp = params(seq)
    fractions = sp.get_amino_acid_fractions()
    feats = {f"aa_frac_{aa}": float(fractions.get(aa, 0.0)) for aa in AA}
    feats["frac_charged"] = float(sp.get_FCR())
    feats["frac_polar"] = sum(float(fractions.get(aa, 0.0)) for aa in POLAR)
    feats["frac_aromatic"] = sum(float(fractions.get(aa, 0.0)) for aa in AROMATIC)
    feats["frac_proline"] = float(fractions.get("P", 0.0))
    feats["frac_glycine"] = float(fractions.get("G", 0.0))
    feats["frac_disorder_promoting"] = float(sp.get_fraction_disorder_promoting())
    return feats


def mean_hydropathy(seq: str) -> float:
    """localCIDER mean normalized Kyte-Doolittle hydropathy (0 to 9)."""
    if not seq:
        return 0.0
    return float(params(seq).get_mean_hydropathy())
