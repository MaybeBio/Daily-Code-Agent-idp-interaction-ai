"""
figS_nmr_maps_prep.py — contact-probability maps for the six proteins of the NMR case series (phospho_nmr_plan.md), for
the Supplementary figure of Fig. 5f-style panels. Same runs, loader and contact rule as redo_evidence/phospho_nmr_test.py
(3 unmodified and 3 phosphorylated CALVADOS-2 replicates per protein under the conditions of the NMR sample, 1,000
frames each, C-alpha < 0.8 nm, |i-j| >= 3). Descriptive; no statistic is computed here.
  ~/.venvs/finches/bin/python figS_nmr_maps_prep.py  ->  cache/nmr_maps_<protein>.npz (pP_reps, pW_reps)
"""
import csv
import os
import sys
from multiprocessing import Pool

import numpy as np

EV = '/Volumes/User Homes/Nu Project/BENDER-DYN/redo_evidence'
sys.path.insert(0, EV)
import phospho_nmr_test as T  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))


def one(job):
    prot, tk, form, k = job
    name = f'NMR-{prot}-{form}-r{k}'
    d = os.path.join(T.RUNS, f'{tk}K', name)
    return job, T.contact_stack(T.load(os.path.join(d, f'{name}.dcd'), os.path.join(d, 'top.pdb'))).mean(0)


def main():
    pairs = list(csv.DictReader(open(os.path.join(EV, 'gaivi', 'phospho_nmr', 'pairs.csv'))))
    jobs = [(p['protein'], p['T_K'], f, k) for p in pairs for f in ('WT', 'P') for k in (1, 2, 3)
            if not os.path.exists(os.path.join(HERE, 'cache', f'nmr_maps_{p["protein"]}.npz'))]
    with Pool(6) as pool:
        res = pool.map(one, jobs)
    out = {}
    for (prot, tk, form, k), m in res:
        out.setdefault(prot, {}).setdefault(form, {})[k] = m
    for prot, d in out.items():
        np.savez(os.path.join(HERE, 'cache', f'nmr_maps_{prot}.npz'),
                 pW_reps=np.array([d['WT'][k] for k in (1, 2, 3)]), pP_reps=np.array([d['P'][k] for k in (1, 2, 3)]))
        print(prot, d['WT'][1].shape)


if __name__ == '__main__':
    main()
