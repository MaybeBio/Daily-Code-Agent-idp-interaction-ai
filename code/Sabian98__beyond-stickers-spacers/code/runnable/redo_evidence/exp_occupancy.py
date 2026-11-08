"""
exp_occupancy.py — hub occupancy (fraction of models in which a residue is among the top 10% by betweenness, the
same definition as hub_frames.py) for the PED ensembles of disordered proteins and their Rg-matched residue-neutral
nulls (ped_null.py), fixed 0.8 nm rule. Feeds the external test of occupancy-trained hub predictors.
Standard library only.  python3 exp_occupancy.py  ->  hub_predictor/experimental_occupancy.csv
"""
import csv
import math
import os
import sys
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import hub_frames as hf  # noqa: E402
import ped_null as pn  # noqa: E402

TOP = 0.10


def occupancy(seq, models):
    n = len(seq)
    n_top = max(1, round(TOP * n))
    occ = [0] * n
    used = 0
    rng = hf.tie_rng(seq)                 # unbiased, reproducible tie-breaking at the hub cutoff
    for m in models:
        if len(m[0]) != n:
            continue
        adj = hf.contact_lists([list(m[0]), list(m[1]), list(m[2])], n)
        bc = hf.betweenness(adj, n)
        for i in hf.top_hubs(bc, n_top, rng):
            occ[i] += 1
        used += 1
    return [o / used for o in occ], used


def one(job):
    k, (entry, path, prot, restraints) = job
    seq, models = pn.load(path)
    if len(models) > 1000:
        models = models[:: math.ceil(len(models) / 1000)]
    nm = pn.null_models(seq, models, 1000 + k)[0]
    oe, ue = occupancy(seq, models)
    on, un = occupancy(seq, nm)
    return entry, prot, restraints, seq, oe, on, ue


def main():
    jobs = list(enumerate(pn.ENSEMBLES))
    with Pool(9) as pool:
        res = list(pool.imap_unordered(one, jobs))
    out = os.path.join(HERE, 'hub_predictor', 'experimental_occupancy.csv')
    with open(out, 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['entry', 'protein', 'restraints', 'pos', 'aa', 'occ_exp', 'occ_null', 'models'])
        for entry, prot, restraints, seq, oe, on, ue in sorted(res):
            for i, a in enumerate(seq):
                w.writerow([entry, prot, restraints, i, a, f'{oe[i]:.4f}', f'{on[i]:.4f}', ue])
    print(f'{len(res)} ensembles -> {out}')


if __name__ == '__main__':
    main()
