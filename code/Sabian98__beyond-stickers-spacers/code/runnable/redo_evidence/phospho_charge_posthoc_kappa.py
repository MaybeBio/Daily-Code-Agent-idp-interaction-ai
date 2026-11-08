"""
phospho_charge_posthoc_kappa.py — POST HOC (2026-10-14, user request), written after Part C and the |dSCD| post hoc were
known: Part C repeated with the change in the charge-patterning parameter kappa (Das & Pappu 2013) when a site is
phosphorylated, in the same way as phospho_charge_posthoc_scd.py (within-protein Spearman with direct rewiring, and the
within-protein permutation test of each predictor beyond the other).

kappa: residue charges K, R +1, D, E -1 (others 0). For blob size g in {5, 6}, sigma = (f+ - f-)^2 / (f+ + f-) for every
window of g residues (0 if the window is uncharged), delta = mean over windows of (sigma_window - sigma_sequence)^2, and
kappa = delta / delta_max, where delta_max is the largest delta over the six block arrangements (positive, neutral and
negative residues as contiguous blocks in every order) of the same multiset of charges; kappa is the mean over g = 5 and
6. f+ and f- are the summed positive and absolute negative charge divided by the window length, so a charge of -2
counts twice. The phosphosite is given charge -1 (phosphomimetic convention, kappa assumes unit charges) or -2 (as in the
simulations); both are reported. |dkappa| = |kappa(phosphorylated) - kappa(unmodified)|.
  ~/.venvs/finches/bin/python phospho_charge_posthoc_kappa.py   ->   phospho_charge/posthoc_kappa.txt
"""
import itertools
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

DYN = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, DYN)
import phospho_charge_test as C  # noqa: E402
from phospho_direct_test import check_frozen  # noqa: E402

DS = os.path.join(DYN, 'phospho_charge', 'direct_sites.csv')
Q = {'K': 1.0, 'R': 1.0, 'D': -1.0, 'E': -1.0}


def sigma(q):
    fp, fm = q[q > 0].sum() / len(q), -q[q < 0].sum() / len(q)
    return 0.0 if fp + fm == 0 else (fp - fm) ** 2 / (fp + fm)


def delta(q, g):
    s = sigma(q)
    return float(np.mean([(sigma(q[i:i + g]) - s) ** 2 for i in range(len(q) - g + 1)]))


def kappa(q):
    out = []
    for g in (5, 6):
        pos, neu, neg = q[q > 0], q[q == 0], q[q < 0]
        dmax = max(delta(np.concatenate(o), g) for o in itertools.permutations((pos, neu, neg)))
        out.append(delta(q, g) / dmax if dmax > 0 else 0.0)
    return float(np.mean(out))


def charges(seq):
    return np.array([Q.get(a, 0.0) for a in seq])


def dkappa(seq, pos, qsite):
    q = charges(seq)
    k0 = kappa(q)
    q = q.copy()
    q[pos - 1] = qsite
    return kappa(q) - k0


def main():
    check_frozen(DS)
    d = pd.read_csv(DS)
    g = d.UniProt_ID.to_numpy()
    L = ['POST HOC (2026-10-14): Part C with |dkappa| (Das & Pappu kappa) as the sequence control', '']
    for qs in (-1.0, -2.0):
        d[f'abs_dk{int(-qs)}'] = [abs(dkappa(C.seq_of(u), p, qs)) for u, p in zip(d.UniProt_ID, d.pos)]
    c = {k: C.centre(d[k], g) for k in ('log_ess', 'abs_dk1', 'abs_dk2', 'seqpred', 'y_occ', 'y_deg')}
    L.append('Within-protein Spearman with direct rewiring (y = E / R(W, W), centred), 60 sites in 20 proteins:')
    for k, lab in (('abs_dk1', '|dkappa|, site -1'), ('abs_dk2', '|dkappa|, site -2')):
        L.append(f'  {lab:20s} {spearmanr(c[k], c["y_occ"]).correlation:+.3f} (occ), '
                 f'{spearmanr(c[k], c["y_deg"]).correlation:+.3f} (deg);  with -log ESS '
                 f'{spearmanr(c[k], -c["log_ess"]).correlation:+.3f}')
    L.append(f'  {"-log ESS":20s} {spearmanr(-c["log_ess"], c["y_occ"]).correlation:+.3f} (occ), '
             f'{spearmanr(-c["log_ess"], c["y_deg"]).correlation:+.3f} (deg)')
    L += ['', f'Permutation tests ({C.B} permutations within protein):']
    rng = np.random.default_rng(C.SEED)
    for k in ('abs_dk1', 'abs_dk2'):
        for name, y in (('E_occ', c['y_occ']), ('E_deg', c['y_deg'])):
            b, p = C.perm_coef(y, c['log_ess'], c[k], g, rng)
            L.append(f'  {k} {name}: log ESS beyond |dkappa|: b {b:+.4f}, P {p:.4g}')
            b, p = C.perm_coef(-y, c[k], c['log_ess'], g, rng)
            L.append(f'  {k} {name}: |dkappa| beyond log ESS: b {-b:+.4f}, P {p:.4g}')
    L += ['', 'kappa of the unmodified sequences: median {:.3f}, range {:.3f} to {:.3f}; |dkappa| (site -2) median {:.4f}'.format(
        *(lambda v: (np.median(v), v.min(), v.max()))(np.array([kappa(charges(C.seq_of(u))) for u in d.UniProt_ID.unique()])),
        d.abs_dk2.median())]
    txt = '\n'.join(L) + '\n'
    open(os.path.join(DYN, 'phospho_charge', 'posthoc_kappa.txt'), 'w').write(txt)
    print(txt)


if __name__ == '__main__':
    main()
