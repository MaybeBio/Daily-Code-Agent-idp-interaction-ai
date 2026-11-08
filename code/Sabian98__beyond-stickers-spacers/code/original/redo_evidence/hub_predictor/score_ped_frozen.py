"""
score_ped_frozen.py — per-ensemble scores of the frozen models (frozen_models_2026-10-02.pkl) on the 18 PED ensembles
used during development (experimental_occupancy.csv; NOT blind) and on the 8 fresh ones (fresh_ped_occupancy.csv; blind,
predictions frozen before labels). Raw Spearman and Spearman partial on distance from the chain ends, for the model
(consensus / pooled / CALVADOS-2-only), the residue-neutral null and position only.
  ~/.venvs/bender/bin/python score_ped_frozen.py  ->  ped_scores_frozen.csv
"""
import csv
import os
import pickle
from collections import defaultdict

import numpy as np
from scipy.stats import rankdata, spearmanr

import train_hub_model as T

HERE = os.path.dirname(os.path.abspath(__file__))


def partial(x, y, z):
    rx, ry, rz = (rankdata(v) for v in (x, y, z))
    A = np.vstack([rz, np.ones_like(rz)]).T
    ex = rx - A @ np.linalg.lstsq(A, rx, rcond=None)[0]
    ey = ry - A @ np.linalg.lstsq(A, ry, rcond=None)[0]
    return float(np.corrcoef(ex, ey)[0, 1])


def multi_null():
    """Null occupancy averaged over independent pools (null_multi.py); a single pool is noisy per residue."""
    f = os.path.join(HERE, 'null_occupancy_multi.csv')
    d = defaultdict(dict)
    if os.path.exists(f):
        for r in csv.DictReader(open(f)):
            d[r['entry']][int(r['pos'])] = float(r['occ_null_mean'])
    return {e: np.array([v[i] for i in range(len(v))]) for e, v in d.items()}


def main():
    models = pickle.load(open(os.path.join(HERE, 'frozen_models_2026-10-02.pkl'), 'rb'))
    MN = multi_null()
    print(f'null baseline: {"mean over independent pools" if MN else "single pool"} ({len(MN)} ensembles)')
    out = []
    # 18 used ensembles: predict now with the frozen models
    ex = defaultdict(dict)
    meta = {}
    for r in csv.DictReader(open(os.path.join(HERE, 'experimental_occupancy.csv'))):
        ex[r['entry']][int(r['pos'])] = (r['aa'], float(r['occ_exp']), float(r['occ_null']))
        meta[r['entry']] = (r['protein'], r['restraints'])
    sets = []
    for e, d in sorted(ex.items()):
        n = len(d)
        seq = ''.join(d[i][0] for i in range(n))
        X = T.features(seq, 10)
        pr = {k: m.predict(X) for k, m in models.items()}
        sets.append(('used18', e, meta[e][0], meta[e][1], n, np.array([d[i][1] for i in range(n)]),
                     np.array([d[i][2] for i in range(n)]), pr))
    # 8 fresh ensembles: frozen predictions stored with the labels
    fr = defaultdict(list)
    for r in csv.DictReader(open(os.path.join(HERE, 'fresh_ped_occupancy.csv'))):
        fr[r['entry']].append(r)
    for e, rows in sorted(fr.items()):
        n = len(rows)
        pr = {k: np.array([float(r[f'pred_{k}']) for r in rows]) for k in ('consensus', 'pooled', 'cal')}
        sets.append(('fresh8', e, rows[0]['protein'], rows[0]['restraints'], n,
                     np.array([float(r['occ_exp']) for r in rows]), np.array([float(r['occ_null']) for r in rows]), pr))
    for s, e, prot, rest, n, oe, on, pr in sets:
        if e in MN and len(MN[e]) == n:
            on = MN[e]
        end = np.minimum(np.arange(n), n - 1 - np.arange(n)).astype(float)
        row = [s, e, prot, rest, n]
        for name, v in (('consensus', pr['consensus']), ('pooled', pr['pooled']), ('cal', pr['cal']), ('null', on),
                        ('position', end)):
            row += [f'{spearmanr(v, oe).correlation:.4f}', f'{partial(v, oe, end):.4f}' if name != 'position' else '']
        out.append(row)
    hdr = ['set', 'entry', 'protein', 'restraints', 'N']
    for name in ('consensus', 'pooled', 'cal', 'null', 'position'):
        hdr += [f'rho_{name}', f'partial_{name}']
    with open(os.path.join(HERE, 'ped_scores_frozen.csv'), 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(hdr)
        w.writerows(out)
    for s in ('used18', 'fresh8'):
        R = [r for r in out if r[0] == s]
        med = lambda j: np.nanmedian([float(r[j]) for r in R])  # noqa: E731
        print(f'{s}: n={len(R)}  rho consensus {med(5):+.2f} null {med(11):+.2f} position {med(13):+.2f} | '
              f'partial consensus {med(6):+.2f} null {med(12):+.2f}')


if __name__ == '__main__':
    main()
