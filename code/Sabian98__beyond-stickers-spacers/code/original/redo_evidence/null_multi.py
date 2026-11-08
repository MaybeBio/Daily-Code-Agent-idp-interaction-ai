"""
null_multi.py — the residue-neutral null baseline for the PED tests, averaged over several independent chain pools.

A single null draw (ped_null.null_models, one pool of 4,000 chains) is noisy at the level of individual residues
(independent draws correlate at Spearman 0.03-0.09 after removing chain position), which handicaps the null as a
baseline. The expected occupancy of a residue-neutral chain is better estimated by averaging K independent pools.
For each PED ensemble (the 18 of ped_null.ENSEMBLES and the 8 of ped_fresh_set.csv) this computes hub occupancy of the
Rg-matched null for K seeds and writes the mean (and the per-seed values) to hub_predictor/null_occupancy_multi.csv.
Same model subsampling, contact rule and hub definition as exp_occupancy.py.
  python3 null_multi.py [K=5]
"""
import csv
import math
import os
import sys
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import exp_occupancy as eo  # noqa: E402
import ped_null as pn  # noqa: E402

K = int(sys.argv[1]) if len(sys.argv) > 1 else 5


def one(job):
    entry, path, seed = job
    seq, models = pn.load(path)
    if len(models) > 1000:
        models = models[::math.ceil(len(models) / 1000)]
    nm = pn.null_models(seq, models, seed)[0]
    occ, _ = eo.occupancy(seq, nm)
    return entry, seed, seq, occ


def main():
    ens = [(e, f) for e, f, _, _ in pn.ENSEMBLES]
    ens += [(r['id'], r['file']) for r in csv.DictReader(open(os.path.join(HERE, 'ped_fresh_set.csv')))]
    jobs = [(e, f, 50000 + 97 * k + j) for j, (e, f) in enumerate(ens) for k in range(K)]
    res = {}
    with Pool(8) as p:
        for n, (e, seed, seq, occ) in enumerate(p.imap_unordered(one, jobs), 1):
            res.setdefault(e, (seq, []))[1].append(occ)
            if n % 13 == 0:
                print(f'{n}/{len(jobs)}', flush=True)
    out = os.path.join(HERE, 'hub_predictor', 'null_occupancy_multi.csv')
    with open(out, 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['entry', 'pos', 'aa', 'occ_null_mean', 'n_pools'] + [f'occ_null_{k}' for k in range(K)])
        for e, (seq, occs) in sorted(res.items()):
            for i, a in enumerate(seq):
                vals = [o[i] for o in occs]
                w.writerow([e, i, a, f'{sum(vals) / len(vals):.4f}', len(vals)] + [f'{v:.4f}' for v in vals])
    print(f'{len(res)} ensembles x {K} pools -> {out}')


if __name__ == '__main__':
    main()
