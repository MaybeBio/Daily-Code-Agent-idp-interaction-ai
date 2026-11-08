"""
figS_nmr_maps_prep2.py — contact-probability maps for the four pairs of NMR Amendment 2 (phospho_nmr_amend2), CALVADOS-2
runs, same loader and contact rule as figS_nmr_maps_prep.py. Descriptive.
  ~/.venvs/finches/bin/python figS_nmr_maps_prep2.py  ->  cache/nmr_maps_<pair>.npz
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
RUNS = '/Volumes/User Homes/Nu Project/BENDER-DYN/phospho_nmr2'


def one(job):
    prot, tk, form, k = job
    name = f'NMR2-{prot}-{form}-r{k}'
    d = os.path.join(RUNS, f'{tk}K', name)
    return job, T.contact_stack(T.load(os.path.join(d, f'{name}.dcd'), os.path.join(d, 'top.pdb'))).mean(0)


def main():
    pairs = list(csv.DictReader(open(os.path.join(EV, 'phospho_nmr_amend2', 'pairs_nmr2.csv'))))
    jobs = [(p['protein'], p['T_K'], f, k) for p in pairs for f in ('WT', 'P') for k in (1, 2, 3)]
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
