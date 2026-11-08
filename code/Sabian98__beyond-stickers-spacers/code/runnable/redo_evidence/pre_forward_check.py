"""
pre_forward_check.py — diagnostic for pre_reweight.py: does the Calpha PRE forward model reproduce the PRE data with the
experts' own PRE-restrained PED ensembles (uniform weights)? If not, the forward model, not the reweighting, is at fault.
Also reports the distance range that dominates the predicted <r^-6>.
  ~/.venvs/bender/bin/python pre_forward_check.py
"""
import os
import sys
from multiprocessing import Pool

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ped_null as pn  # noqa: E402
import pre_reweight as pr  # noqa: E402


def main():
    jobs = []
    for name, prot, pool, ref in pr.VALIDATE:
        seq, _ = pn.load(ref)
        jobs.append((f'refpre_{prot}', 'ped', ref, prot, pr.offset_of(seq, prot)))
    with Pool(3) as p:
        list(p.imap_unordered(pr.build, jobs))
    print(f'near-label exclusion |i - s| <= {pr.EXCLUDE}')
    for name, prot, pool, ref in pr.VALIDATE:
        if prot in pr.TOTAL_R2:
            continue
        for tag in (f'refpre_{prot}', f'pool_{prot}'):
            c = pr.restrict(pr.load(tag))
            n = c['r6'].shape[0]
            u = np.full(n, 1.0 / n)
            fits = {t: pr.chi2(*pr.calc_matrix(c, t)[:3], u) for t in pr.TAUS}
            best = min(fits, key=fits.get)
            F, y, sd, lab = pr.calc_matrix(c, best)
            pred = u @ F
            # share of <r^-6> coming from frames closer than 8 A
            r = c['r6'] ** (-1 / 6)
            close = np.where(r < 8, c['r6'], 0).sum() / c['r6'].sum()
            print(f'{tag:14s} frames {n:5d} best tau_c {best * 1e9:4.0f} ns  chi2 {fits[best]:8.1f} | median pred/exp Gamma2 '
                  f'{np.median(pred / y):.2f} | share of <r^-6> from r < 8 A: {close:.2f}')


if __name__ == '__main__':
    main()
