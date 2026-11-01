"""Sequence descriptors: the network's inputs and the closed form's inputs, computed from the amino-acid sequence alone.

Residue parameters are those of the HPS1 model (Dignon et al., PLoS Comput. Biol. 14, e1005941, 2018): charge q (e;
histidine +0.5) and hydropathy lambda. The model was trained on descriptors rounded to float32; `descriptors` returns
them rounded the same way.

    scd_pos, scd_neg   s = asinh(5 SCD / N) split at zero: s+ = max(s, 0), s- = min(s, 0)
                       SCD = sum_{i<j} q_i q_j |i - j|^(1/2) / N  (sequence charge decoration, Sawle & Ghosh 2015)
    shd                sum_{i<j} (lambda_i + lambda_j) |i - j|^-1 / N  (sequence hydropathy decoration, Zheng et al. 2020)
    contact_energy     loop-entropy-weighted contact score of residue windows (see `contact_energy`)
    fcr                fraction of charged residues
    patch_neg/_neu     fraction of adjacent residue pairs that are both charge -1 / both neutral
    frac_sticky        fraction of residues with lambda >= 0.85 (I L F P W Y V)
    mean_hydro         mean lambda, <lambda>
"""
import numpy as np
import pandas as pd

#               q       lambda
RESIDUES = {
    "A": (0.0, 0.730), "R": (1.0, 0.000), "N": (0.0, 0.432), "D": (-1.0, 0.378), "C": (0.0, 0.595),
    "Q": (0.0, 0.514), "E": (-1.0, 0.459), "G": (0.0, 0.649), "H": (0.5, 0.514), "I": (0.0, 0.973),
    "L": (0.0, 0.973), "K": (1.0, 0.514), "M": (0.0, 0.838), "F": (0.0, 1.000), "P": (0.0, 1.000),
    "S": (0.0, 0.595), "T": (0.0, 0.676), "W": (0.0, 0.946), "Y": (0.0, 0.865), "V": (0.0, 0.892),
}
_INDEX = {a: i for i, a in enumerate(RESIDUES)}
_Q = np.array([v[0] for v in RESIDUES.values()], np.float32)
_LAM = np.array([v[1] for v in RESIDUES.values()], np.float32)

# The network's inputs, in model order, and the law covariate <lambda>.
INPUTS = ("log_seq_length", "scd_pos", "scd_neg", "shd", "contact_energy", "fcr", "patch_neg", "patch_neu", "frac_sticky")
COLUMNS = INPUTS[1:] + ("mean_hydro",)

SCD_K = 5.0
STICKY = 0.85

# Contact energy. Pair energy of two residue windows with mean charge q and mean hydropathy lambda: a linear fit
# c0 + alpha q_a q_b + beta (lambda_a + lambda_b) / 2 to the HPS1 pair potential's reduced second virial coefficients at
# 298 K and 150 mM over the 210 residue pairs (R^2 0.979); tau sets how strongly the strongest contacts dominate.
CE_C0, CE_ALPHA, CE_BETA = -0.018783755317377603, 0.45324834008573206, -0.17632799639803304
CE_TAU = 0.06
# Window widths w and the standardization constants (mu_w, sigma_w): the mean and sample standard deviation (ddof 1) of
# S_w over the training proteins that use width w (2w <= N), frozen so the descriptor is a fixed function of the
# sequence; "N/8" is round(N / 8), at least 3.
CE_WIDTHS = {5: (18.34902013609721, 5.978780467194333), 9: (26.102342172719517, 7.617873214997021),
             15: (37.972913187600774, 10.166576667199715), "N/8": (55.82912431720738, 31.535882252919006)}


def encode(sequence):
    """Residue charges and hydropathies (float32) of a one-letter sequence; raises on anything but the 20 standard
    amino acids."""
    s = "".join(sequence.split()).upper()
    bad = sorted({a for a in s if a not in _INDEX})
    if bad:
        raise ValueError(f"non-standard residue(s) {''.join(bad)!r}: only the 20 standard amino acids are supported")
    if len(s) < 2:
        raise ValueError("sequence shorter than 2 residues")
    i = np.array([_INDEX[a] for a in s])
    return _Q[i], _LAM[i]


def scd(q):
    n = len(q)
    d = np.sqrt(np.abs(np.subtract.outer(np.arange(n), np.arange(n))))
    return float(np.tril(np.outer(q, q) * d, -1).sum() / n)


def shd(lam):
    n = len(lam)
    d = np.abs(np.subtract.outer(np.arange(n), np.arange(n))).astype(float)
    return float(np.tril((lam[:, None] + lam[None, :]) / np.where(d == 0, 1.0, d), -1).sum() / n)


def patch(q, charge):
    """Fraction of the N - 1 adjacent pairs whose two residues both have exactly this charge."""
    return float(np.sum((q[:-1] == charge) & (q[1:] == charge)) / (len(q) - 1))


def contact_energy(q, lam):
    """For windows of w residues (mean charge q_a, mean hydropathy lambda_a), E_ab = w (c0 + alpha q_a q_b +
    beta (lambda_a + lambda_b) / 2) for every non-overlapping pair a < b, and
        S_w = ln sum_ab (b - a)^(-3/2) K_ab exp(-E_ab / tau),
    K_ab = the fractional drop of a Gaussian chain's Rg^2 when the window centres a, b (fractions of N) are joined.
    The descriptor is the mean over the widths of the standardized S_w (widths with 2w > N skipped)."""
    q, lam = np.asarray(q, float), np.asarray(lam, float)
    n, zs = len(q), []
    for wl, (mu, sd) in CE_WIDTHS.items():
        w = max(3, int(round(n / 8))) if wl == "N/8" else wl
        if 2 * w > n:
            continue
        box = np.ones(w) / w
        qa, la = np.convolve(q, box, "valid"), np.convolve(lam, box, "valid")
        i, j = np.triu_indices(len(qa), w)
        if not len(i):
            continue
        x, y = (i + w / 2) / n, (j + w / 2) / n
        loop, tails = y - x, 1 - (y - x)
        k = 1 - 6 * (loop ** 3 / 12 + tails ** 3 / 6 + loop * (x ** 2 + (1 - y) ** 2) / 2 + tails * loop ** 2 / 6)
        e = w * (CE_C0 + CE_ALPHA * qa[i] * qa[j] + CE_BETA * (la[i] + la[j]) / 2)
        t = -e / CE_TAU - 1.5 * np.log((j - i).astype(float)) + np.log(np.clip(k, 1e-12, None))
        s = t.max() + np.log(np.exp(t - t.max()).sum())
        zs.append((s - mu) / sd)
    return float(np.mean(zs)) if zs else 0.0


def descriptors(sequence):
    """{N, log_seq_length, scd_pos, scd_neg, shd, contact_energy, fcr, patch_neg, patch_neu, frac_sticky, mean_hydro}."""
    q, lam = encode(sequence)
    n = len(q)
    s = float(np.arcsinh(SCD_K * scd(q) / n))
    d = dict(scd_pos=max(s, 0.0), scd_neg=min(s, 0.0), shd=shd(lam), contact_energy=contact_energy(q, lam),
             fcr=float(np.sum(q != 0) / n), patch_neg=patch(q, -1.0), patch_neu=patch(q, 0.0),
             frac_sticky=float(np.mean(lam >= STICKY)), mean_hydro=float(np.mean(lam)))
    return dict(N=n, log_seq_length=float(np.log(n)), **{k: float(np.float32(v)) for k, v in d.items()})


def table(sequences):
    """One row of `descriptors` per sequence."""
    rows = []
    for k, s in enumerate(sequences):
        try:
            rows.append(descriptors(s))
        except ValueError as e:
            raise ValueError(f"sequence {k + 1}: {e}") from None
    return pd.DataFrame(rows)
