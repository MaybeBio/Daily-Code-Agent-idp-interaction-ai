"""
phospho_nmr_rw_diagnostic.py — POST HOC diagnostic, after the NMR primary (phospho_nmr_test.py a68098b9): why the charge-
model reweighting predicted the shift changes less well than the direct runs. Per protein and unmodified replicate:
effective sample fraction of the multisite weights (all sites, mutual terms, sample conditions; exactly as in the test
script), and the Spearman correlation of the reweighted contact-rewiring profile with the direct one (E_i), over all
residues. No new test; descriptive.

  ~/.venvs/finches/bin/python phospho_nmr_rw_diagnostic.py   ->   phospho_nmr/posthoc_rw_diagnostic.txt
"""
import csv
import os
import sys

DYN = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, DYN)
import numpy as np  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402

import phospho_nmr_test as T  # noqa: E402


def main():
    pairs = list(csv.DictReader(open(T.PAIRS)))
    L = ['POST HOC: charge-model reweighting diagnostic for the NMR case series', '',
         'protein    sites  total dq  ESS fraction per WT replicate   rho(reweighted, direct E)   per-replicate rho']
    for p in pairs:
        prot, Tk, ph, ionic = p['protein'], float(p['T_K']), float(p['pH']), float(p['ionic_M'])
        sites = [int(s[1:]) - 1 for s in p['sites'].split()]
        res = {}
        for form in ('WT', 'P'):
            for k in (1, 2, 3):
                name = f'NMR-{prot}-{form}-r{k}'
                d = os.path.join(T.RUNS, f'{p["T_K"]}K', name)
                xyz = T.load(os.path.join(d, f'{name}.dcd'), os.path.join(d, 'top.pdb'))
                res[name] = (xyz, T.contact_stack(xyz))
        W = [res[f'NMR-{prot}-WT-r{k}'] for k in (1, 2, 3)]
        P = [res[f'NMR-{prot}-P-r{k}'] for k in (1, 2, 3)]
        pw = [w[1].mean(0) for w in W]
        pp = [q[1].mean(0) for q in P]
        rr = lambda a, b: np.abs(a - b).sum(1)  # noqa: E731
        E = np.mean([rr(a, b) for a in pp for b in pw], 0) - np.mean([rr(pw[a], pw[b]) for a, b in ((0, 1), (0, 2), (1, 2))], 0)
        ess, rws = [], []
        for (xyz, C), p0 in zip(W, pw):
            phi, (kap, eps, rc, ec, beta) = T.phi_cond(xyz, p['wt_seq'], Tk, ionic, ph)
            dH = T.pn.DQ * phi[:, sites].sum(1)
            for a in range(len(sites)):
                for b in range(a + 1, len(sites)):
                    if abs(sites[a] - sites[b]) >= 2:
                        r = np.linalg.norm(xyz[:, sites[a]] - xyz[:, sites[b]], axis=1)
                        dH += T.pn.DQ * T.pn.DQ * eps * np.where(r < rc, np.exp(-kap * r) / r - ec / rc, 0.0)
            x = -beta * dH
            w = np.exp(x - x.max())
            w /= w.sum()
            ess.append(1.0 / (w ** 2).sum() / len(w))
            rws.append(np.abs(np.tensordot(w, C.astype(float), axes=1) - p0).sum(1))
        rho = spearmanr(np.mean(rws, 0), E).correlation
        per = [spearmanr(r, E).correlation for r in rws]
        L.append(f'{prot:10s} {len(sites):5d}  {-2 * len(sites):8d}  ' + ' '.join(f'{e:.3f}' for e in ess)
                 + f'        {rho:+.3f}                   ' + ' '.join(f'{v:+.2f}' for v in per))
    L += ['', 'ESS fraction = 1 / sum(w^2) / n_frames (1 = all 1,000 frames count equally). For comparison, the 60 single sites',
          'of phospho_direct had median ESS about 0.65 and reweighted Rg change agreed with direct runs at rho 0.82.']
    txt = '\n'.join(L) + '\n'
    open(os.path.join(T.OUT, 'posthoc_rw_diagnostic.txt'), 'w').write(txt)
    print(txt)


if __name__ == '__main__':
    main()
