"""
hub_vs_degree_ped.py — Part B of hub_vs_degree_plan.md (hash 9212f51e).
Reproducibility of degree profiles and occupancy profiles between independent experimental ensembles, for the pairs of
hub_predictor/ped_ceiling_pairs.csv whose two ensembles are both on disk in ped_local/. Occupancy and degree are computed
here from the same models (every ceil(m/1000)-th model), with the contact rule of hub_frames.py.

  OMP_NUM_THREADS=1 ~/.venvs/bender/bin/python hub_vs_degree_ped.py  ->  hub_vs_degree/part_b.txt
"""
import csv
import math
import os
import sys
from multiprocessing import Pool

import numpy as np
from scipy.stats import rankdata, spearmanr

HERE = (__import__('os').environ['BSS_REPO'] + '/redo_evidence')
PED = os.path.join(HERE, 'ped_local')
sys.path.insert(0, HERE)
os.chdir(PED)                      # ped_null.load reads paths relative to ped_local
import exp_occupancy as EO  # noqa: E402
import hub_frames as hf  # noqa: E402
import ped_null as pn  # noqa: E402


def pdb_of(entry):
    """The one multi-model PDB of an entry directory (names differ: e001, e002, e003, e006)."""
    f = [x for x in os.listdir(entry) if x.endswith('.pdb') and not x.startswith('.')]
    return f'{entry}/{f[0]}' if len(f) == 1 else None


def profiles(entry):
    seq, models = pn.load(os.path.abspath(pdb_of(entry)))   # ped_null.PED points at an old location, an absolute path overrides it
    total = len(models)
    if total > 1000:
        models = models[:: math.ceil(total / 1000)]
    n = len(seq)
    deg, used = np.zeros(n), 0
    for m in models:
        if len(m[0]) != n:
            continue
        adj = hf.contact_lists([list(m[0]), list(m[1]), list(m[2])], n)
        deg += np.array([len(a) for a in adj])
        used += 1
    occ, _ = EO.occupancy(seq, models)
    return entry, seq, deg / used, np.array(occ), used


def partial(x, y):
    n = len(y)
    z = np.minimum(np.arange(n), n - 1 - np.arange(n)).astype(float)
    rx, ry, rz = rankdata(x), rankdata(y), rankdata(z)
    A = np.vstack([rz, np.ones_like(rz)]).T
    ex = rx - A @ np.linalg.lstsq(A, rx, rcond=None)[0]
    ey = ry - A @ np.linalg.lstsq(A, ry, rcond=None)[0]
    return float(np.corrcoef(ex, ey)[0, 1])


def main():
    pairs = list(csv.DictReader(open(os.path.join(HERE, 'hub_predictor/ped_ceiling_pairs.csv'))))
    local = {e for e in os.listdir(PED) if os.path.isdir(e) and pdb_of(e)}
    use = [p for p in pairs if p['entry_a'] in local and p['entry_b'] in local]
    entries = sorted({p['entry_a'] for p in use} | {p['entry_b'] for p in use})
    print(f'{len(pairs)} ceiling pairs, {len(use)} with both ensembles on disk, {len(entries)} ensembles', flush=True)
    with Pool(min(9, len(entries))) as pool:
        res = {r[0]: r for r in pool.map(profiles, entries, chunksize=1)}
    L = ['reproducibility between independent experimental ensembles, degree against occupancy (hub_vs_degree_plan.md, Part B)',
         f'ceiling pairs with both ensembles on disk: {len(use)} of {len(pairs)}; ensembles {len(entries)}', '',
         '  protein  pair                    n  | occupancy rho  partial | degree rho  partial']
    ro, po, rd, pd_ = [], [], [], []
    for p in use:
        a, b = res[p['entry_a']], res[p['entry_b']]
        if len(a[1]) != len(b[1]):
            continue
        oa, ob, da, db = a[3], b[3], a[2], b[2]
        row = (spearmanr(oa, ob).correlation, partial(oa, ob), spearmanr(da, db).correlation, partial(da, db))
        ro.append(row[0]); po.append(row[1]); rd.append(row[2]); pd_.append(row[3])
        L.append(f'  {p["protein"]:7s} {p["entry_a"]}-{p["entry_b"]}  {len(a[1]):4d}  |   {row[0]:+.3f}  {row[1]:+.3f}   |  {row[2]:+.3f}  {row[3]:+.3f}')
    n = len(ro)
    L.append('')
    L.append(f'medians over {n} pairs: occupancy {np.median(ro):+.3f} (partial {np.median(po):+.3f}), '
             f'degree {np.median(rd):+.3f} (partial {np.median(pd_):+.3f})')
    diff = np.median(rd) - np.median(ro)
    verdict = ('DEGREE more reproducible' if diff >= 0.05 else 'OCCUPANCY more reproducible' if diff <= -0.05 else 'TIED')
    L.append(f'Rule R2 (raw agreement, margin 0.05): {verdict} (degree minus occupancy {diff:+.3f})'
             + ('; descriptive only, fewer than 8 pairs' if n < 8 else ''))
    L.append(f'degree higher than occupancy in {sum(d > o for d, o in zip(rd, ro))} of {n} pairs')
    open(os.path.join(HERE, 'hub_vs_degree', 'part_b.txt'), 'w').write('\n'.join(L) + '\n')
    print('\n'.join(L))


if __name__ == '__main__':
    main()
