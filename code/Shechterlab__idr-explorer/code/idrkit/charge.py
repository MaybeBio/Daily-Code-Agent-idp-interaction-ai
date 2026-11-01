"""Residue charge, net and total charge, charge patterning, and basic patches.

One definition of each, for the Shechter lab's repositories (glu_proteomics, MicroDomains, PTM_Pipeline and the IDR
Explorer). Positions are 0-based half-open [start, end) intervals, as in glu_proteomics, unless a function says so.
"""
from __future__ import annotations

import numpy as np

# Histidine is the only residue the three repositories counted differently; each keeps its value under a name.
HIS = {
    "plain": 0.0,       # glu_proteomics (config sequence.charge_histidine), the IDR Explorer's IDR tables
    "md": 0.1,          # MicroDomains md_figure_style.CHARGE
    "emboss": 0.5,      # EMBOSS charge scale: the Methylarginine Site Explorer's charge track
}


def charge_vector(seq: str, his: float = 0.0) -> np.ndarray:
    """Per-residue charge: K and R +1, D and E -1, H `his` (see HIS), everything else 0."""
    lut = {"D": -1.0, "E": -1.0, "K": 1.0, "R": 1.0, "H": his}
    return np.fromiter((lut.get(c, 0.0) for c in seq), dtype=float, count=len(seq))


def ncpr(seq: str, his: float = 0.0) -> float:
    """Net charge per residue."""
    return float(charge_vector(seq, his).mean()) if seq else float("nan")


def fcr(seq: str) -> float:
    """Fraction of charged residues (D, E, K, R)."""
    return sum(c in "DEKR" for c in seq) / len(seq) if seq else float("nan")


# Charge type, in words, from two numbers: net charge per residue first (a segment with 20% D/E and 4% K/R is acidic
# however few its charges), then the fraction of charged residues (FCR 0.25 is the Das & Pappu boundary of weak
# polyampholytes, PNAS 110:13392, 2013).
NET_CUT = 0.10      # |NCPR| at or above: acidic or basic (net charge 5 or more per 50 residues)
FCR_CUT = 0.25      # FCR at or above, net charge below NET_CUT: mixed; below: few charges
CHARGE_TYPES = ("few", "mixed", "acidic", "basic")


def charge_type(fcr_value: float, ncpr_value: float) -> str:
    """'acidic' (NCPR <= -0.10), 'basic' (NCPR >= +0.10), else 'mixed' (FCR >= 0.25) or 'few' (fewer charges)."""
    if ncpr_value <= -NET_CUT:
        return "acidic"
    if ncpr_value >= NET_CUT:
        return "basic"
    return "mixed" if fcr_value >= FCR_CUT else "few"


def scd(seq: str, his: float = 0.0) -> float:
    """Sequence charge decoration: (1/N) sum over i < j of q_i q_j sqrt(j - i).

    Negative when opposite charges are segregated into blocks, near zero when they are well mixed, positive for a
    chain with net charge. (K8E8)3 gives -4.290. Memory is O(N); time O(N^2), a fraction of a second at N = 5,000.
    """
    return _scd_q(charge_vector(seq, his))


def _scd_q(q: np.ndarray) -> float:
    n = len(q)
    if n < 2:
        return 0.0
    nz = np.flatnonzero(q)
    if len(nz) < 2:
        return 0.0
    qs, pos = q[nz], nz.astype(float)
    total = 0.0
    for k in range(len(nz) - 1):
        total += qs[k] * float((qs[k + 1:] * np.sqrt(pos[k + 1:] - pos[k])).sum())
    return total / n


def scd_z(seq: str, n_shuffles: int = 40, seed: int = 0, his: float = 0.0) -> tuple[float, float]:
    """SCD and its z-score against shuffles of the same residues (length and composition kept): below 0, the charges
    are more segregated into blocks than the composition alone gives; above 0, better mixed."""
    q = charge_vector(seq, his)
    obs = _scd_q(q)
    if np.count_nonzero(q) < 2:
        return obs, 0.0
    rng = np.random.default_rng(seed)
    null = np.array([_scd_q(rng.permutation(q)) for _ in range(n_shuffles)])
    sd = null.std(ddof=1)
    return float(obs), float((obs - null.mean()) / sd) if sd > 0 else 0.0


def _delta_form(pos: np.ndarray, neg: np.ndarray, blob: int) -> float:
    n = len(pos)
    if n < blob:
        return 0.0
    fp, fn = pos.mean(), neg.mean()
    sigma = (fp - fn) ** 2 / (fp + fn) if fp + fn else 0.0
    cp, cn = np.concatenate([[0], np.cumsum(pos)]), np.concatenate([[0], np.cumsum(neg)])
    bp, bn = (cp[blob:] - cp[:-blob]) / blob, (cn[blob:] - cn[:-blob]) / blob
    bf = bp + bn
    bsig = np.divide((bp - bn) ** 2, bf, out=np.zeros_like(bf), where=bf > 0)
    return float(((sigma - bsig) ** 2).mean())


def _delta_q(q: np.ndarray) -> float:
    pos, neg = (q > 0).astype(float), (q < 0).astype(float)
    return (_delta_form(pos, neg, 5) + _delta_form(pos, neg, 6)) / 2


def delta(seq: str) -> float:
    """Das & Pappu's delta (PNAS 110:13392, 2013), as localCIDER computes it: the mean over blobs of 5 and of 6
    residues of (sigma - sigma_blob)^2, sigma = NCPR^2 / FCR (K, R positive; D, E negative). Higher: charges in local
    blocks. kappa divides it by the largest delta of the composition; delta_z compares it with shuffles instead."""
    return _delta_q(charge_vector(seq))


def delta_z(seq: str, n_shuffles: int = 40, seed: int = 0) -> tuple[float, float]:
    """delta and its z-score against shuffles of the same residues: above 0, charges sit in more local blocks than the
    composition gives (DDX4's alternating blocks); below 0, better mixed."""
    q = charge_vector(seq)
    obs = _delta_q(q)
    if np.count_nonzero(q) < 2:
        return obs, 0.0
    rng = np.random.default_rng(seed)
    null = np.array([_delta_q(rng.permutation(q)) for _ in range(n_shuffles)])
    sd = null.std(ddof=1)
    return float(obs), float((obs - null.mean()) / sd) if sd > 0 else 0.0


def runs(mask: np.ndarray) -> list[tuple[int, int]]:
    """[start, end) intervals where mask is True."""
    if mask.size == 0:
        return []
    d = np.diff(np.concatenate([[0], mask.astype(int), [0]]))
    return list(zip(np.flatnonzero(d == 1), np.flatnonzero(d == -1)))


def basic_patches(seq: str, window: int = 10, min_net: float = 3.0, his: float = 0.0) -> list[tuple[int, int]]:
    """Merged windows of `window` residues with net charge >= min_net (glu_proteomics ANALYSIS_NOTES 0.10: net >= +3
    in 10 residues)."""
    q = charge_vector(seq, his)
    if len(q) < window:
        return []
    net = np.convolve(q, np.ones(window), "valid")
    mark = np.zeros(len(seq), bool)
    for i in np.flatnonzero(net >= min_net):
        mark[i:i + window] = True
    return [(int(a), int(b)) for a, b in runs(mark)]
