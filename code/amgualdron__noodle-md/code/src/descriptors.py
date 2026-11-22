"""Sequence descriptors: numbers computed from the amino-acid sequence alone, before any simulation.

./noodle launch puts them in the run's manifest.csv, one row per protein and physics set.
They use the residue charges q and hydropathies lambda of the model itself (HPS1, the
table in forcefield.py: histidine has charge +0.5), and the same formulas and column names
as the group's idp-rg project, so the tables can be joined directly.

    N                number of residues
    net_charge       sum of the charges, in e
    ncpr             net charge per residue: net_charge / N
    fcr              fraction of charged residues: residues with q != 0 (histidine included)
    f_plus, f_minus  positive and negative charge per residue (sum of q > 0, and of |q| for q < 0, over N)
    mean_hydro       mean hydropathy <lambda>
    kappa            how mixed (0) or segregated into blocks (1) the + and - charges are (Das & Pappu 2013);
                     empty for a sequence without charges. It is measured against the most segregated
                     arrangement of the same charges in blocks (see kappa below). For ~2% of the
                     shared sequences, mostly with few charges (fcr < 0.2), spreading the charges out
                     beats any block arrangement and kappa comes out above 1: there it says little.
    scd              sequence charge decoration, sum over i < j of q_i q_j |i - j|^(1/2), divided by N
                     (Sawle & Ghosh 2015): negative when opposite charges sit in blocks
    shd              sequence hydropathy decoration, sum over i < j of (lambda_i + lambda_j) |i - j|^(-1),
                     divided by N (Zheng et al. 2020)
"""

import math

import numpy as np

from forcefield import AA_PARAMS

KAPPA_BLOBS = (5, 6)  # blob sizes kappa is averaged over (Das & Pappu 2013, as localCIDER does)


def scd(q):
    n = len(q)
    separation = np.abs(np.subtract.outer(np.arange(n), np.arange(n)))
    return float(np.tril(np.outer(q, q) * np.sqrt(separation), -1).sum() / n)


def shd(lam):
    n = len(lam)
    separation = np.abs(np.subtract.outer(np.arange(n), np.arange(n))).astype(float)
    separation[separation == 0] = 1.0  # the diagonal, left out by np.tril below
    return float(np.tril((lam[:, None] + lam[None, :]) / separation, -1).sum() / n)


def charge_asymmetry(plus, minus):
    """sigma = (f+ - f-)^2 / (f+ + f-), elementwise; 0 where there is no charge."""
    total = plus + minus
    return np.divide((plus - minus) ** 2, total, out=np.zeros_like(total), where=total > 0)


def delta(q, blob):
    """Das & Pappu's delta of one or several sequences of charges (one per row of q): the mean
    over every window of `blob` residues of (sigma of the window - sigma of the whole
    sequence)^2. f+ and f- are the positive and negative charge per residue, so a
    histidine (+0.5) counts half."""
    q = np.atleast_2d(q)
    plus, minus = np.clip(q, 0, None), np.clip(-q, 0, None)

    def window_means(x):  # mean of every `blob` consecutive values, along each row
        sums = np.cumsum(np.pad(x, ((0, 0), (1, 0))), axis=1)
        return (sums[:, blob:] - sums[:, :-blob]) / blob

    sigma_blobs = charge_asymmetry(window_means(plus), window_means(minus))
    sigma = charge_asymmetry(plus.mean(axis=1), minus.mean(axis=1))
    return np.mean((sigma_blobs - sigma[:, None]) ** 2, axis=1)


def kappa(q):
    """delta / delta_max, each summed over the blob sizes 5 and 6 (localCIDER's average).

    delta_max is the largest delta of a sequence with the same charges arranged to be as
    segregated as possible: all positive residues in one block, all negative ones in
    another, the neutral ones split between the two ends (as in localCIDER), or between
    the two blocks. Both block orders are tried, each with its histidines (+0.5) at either
    end of the positive block."""
    q = np.asarray(q, float)
    if not q.any():
        return math.nan
    positive, negative, n_neutral = np.sort(q[q > 0])[::-1], q[q < 0], int(np.sum(q == 0))
    candidates = []
    for plus in (positive, positive[::-1]):
        for first, second in ((plus, negative), (negative, plus)):
            charged = np.concatenate([first, second])
            # `before` neutral residues, the charged blocks, then the other neutral ones. Once
            # the blocks are a blob or more away from both ends, moving them changes no
            # window, so only the splits near either end need trying.
            edge = max(KAPPA_BLOBS)
            splits = set(range(0, min(edge, n_neutral) + 1)) | set(range(max(0, n_neutral - edge), n_neutral + 1))
            for before in sorted(splits):
                candidates.append(np.concatenate([np.zeros(before), charged, np.zeros(n_neutral - before)]))
            candidates.append(np.concatenate([first, np.zeros(n_neutral), second]))
    candidates = np.array(candidates)
    value = sum(delta(q, blob)[0] for blob in KAPPA_BLOBS)
    maximum = np.max(sum(delta(candidates, blob) for blob in KAPPA_BLOBS))
    return float(value / maximum) if maximum > 0 else math.nan


def descriptors(sequence):
    """{column: value} of the descriptors listed at the top of this file."""
    q = np.array([AA_PARAMS[aa]["charge"] for aa in sequence], dtype=float)
    lam = np.array([AA_PARAMS[aa]["hydropathy"] for aa in sequence], dtype=float)
    n = len(sequence)
    return dict(N=n, net_charge=q.sum(), ncpr=q.sum() / n, fcr=np.count_nonzero(q) / n,
                f_plus=q[q > 0].sum() / n, f_minus=-q[q < 0].sum() / n, mean_hydro=lam.mean(),
                kappa=kappa(q), scd=scd(q), shd=shd(lam))
