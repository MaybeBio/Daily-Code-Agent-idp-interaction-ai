"""Charge, charge patterning, and net-charge features from localCIDER.

Every value is read from localCIDER's SequenceParameters (Holehouse et al., 2017): the fraction
of charged residues (FCR), net charge per residue (NCPR), net charge, the fraction of positive
and negative residues, the sequence charge decoration (SCD, Sawle and Ghosh 2015), and the two
charge-patterning parameters kappa (Das and Pappu 2013) and Omega (Martin et al. 2016). kappa is
-1 when undefined, which happens when the segment carries fewer than two charge signs; the value
is meaningful only when both positive and negative residues are present, so that case is reported
as undefined rather than as localCIDER's degenerate 1.
"""

from __future__ import annotations

from .physchem import params

_EMPTY = {
    "fcr": 0.0,
    "ncpr": 0.0,
    "net_charge": 0,
    "kappa": -1.0,
    "omega": -1.0,
    "scd": 0.0,
    "frac_positive": 0.0,
    "frac_negative": 0.0,
}


def charge_features(seq: str) -> dict[str, float]:
    if not seq:
        return dict(_EMPTY)
    sp = params(seq)
    return {
        "fcr": float(sp.get_FCR()),
        "ncpr": float(sp.get_NCPR()),
        "net_charge": int(sp.get_countPos() - sp.get_countNeg()),
        "kappa": kappa(seq),
        "omega": float(sp.get_Omega()),
        "scd": float(sp.get_SCD()),
        "frac_positive": float(sp.get_fraction_positive()),
        "frac_negative": float(sp.get_fraction_negative()),
    }


def kappa(seq: str) -> float:
    """localCIDER kappa, or -1 when the segment lacks both charge signs (kappa undefined)."""
    if not seq:
        return -1.0
    sp = params(seq)
    if sp.get_countPos() == 0 or sp.get_countNeg() == 0:
        return -1.0
    return float(sp.get_kappa())


def scd(seq: str) -> float:
    """localCIDER sequence charge decoration (Sawle and Ghosh 2015)."""
    if not seq:
        return 0.0
    return float(params(seq).get_SCD())
