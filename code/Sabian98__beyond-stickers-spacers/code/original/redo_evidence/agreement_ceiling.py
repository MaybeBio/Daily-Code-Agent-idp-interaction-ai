"""
agreement_ceiling.py — the force-field agreement ceiling (methods_paper/RECOMB_PLAN.md, item 3 of "What the full
Mpipi-Recharged run changes"): per-protein Spearman correlation between CALVADOS-2 and Mpipi-Recharged hub occupancy
for every sequence labelled in both, raw and after removing the distance from the chain ends, stratified by organism
group, chain length and compaction (CALVADOS-2 scaling exponent nu from BENDER-BIO/BENDER_BIO.csv).
This bounds any force-field-agnostic predictor: a model cannot be expected to agree with one force field better than
the other force field does.
Strata fixed here: length < 100, 100-199, 200-399, >= 400 residues; nu < 0.40 (compact), 0.40-0.50, >= 0.50 (expanded).
  python3 agreement_ceiling.py  ->  agreement_ceiling_summary.txt, agreement_ceiling_per_protein.csv
"""
import csv
import math
import os
from collections import defaultdict

import numpy as np
from scipy.stats import rankdata, spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
BIO = '/Volumes/User Homes/Nu Project/BENDER-BIO/BENDER_BIO.csv'


def occ(d):
    o, g = defaultdict(dict), {}
    for r in csv.DictReader(open(os.path.join(HERE, d, 'per_residue.csv'))):
        o[r['UniProt_ID']][int(r['pos'])] = float(r['occupancy'])
        g[r['UniProt_ID']] = r['group']
    return {k: np.array([v[i] for i in range(len(v))]) for k, v in o.items()}, g


def partial(x, y):
    n = len(x)
    z = np.minimum(np.arange(n), n - 1 - np.arange(n)).astype(float)
    rx, ry, rz = (rankdata(v) for v in (x, y, z))
    A = np.vstack([rz, np.ones_like(rz)]).T
    ex = rx - A @ np.linalg.lstsq(A, rx, rcond=None)[0]
    ey = ry - A @ np.linalg.lstsq(A, ry, rcond=None)[0]
    return float(np.corrcoef(ex, ey)[0, 1])


def length_bin(n):
    return '<100' if n < 100 else '100-199' if n < 200 else '200-399' if n < 400 else '>=400'


def nu_bin(v):
    if v is None or math.isnan(v):
        return 'nu n/a'
    return 'nu<0.40' if v < 0.40 else 'nu 0.40-0.50' if v < 0.50 else 'nu>=0.50'


def main():
    cal, grp = occ('hub_frames_v3_cal')
    rec, _ = occ('hub_frames_v3_recharged')
    nu = {}
    for r in csv.DictReader(open(BIO)):
        try:
            nu[r['UniProt_ID']] = float(r['nu'])
        except (ValueError, KeyError):
            nu[r['UniProt_ID']] = float('nan')
    rows = []
    for u in sorted(set(cal) & set(rec)):
        a, b = cal[u], rec[u]
        if len(a) != len(b) or a.std() == 0 or b.std() == 0:
            continue
        rows.append(dict(id=u, group=grp.get(u, '?'), N=len(a), nu=nu.get(u, float('nan')),
                         raw=float(spearmanr(a, b).correlation), part=partial(a, b)))
    out = []

    def p(s=''):
        out.append(s)
        print(s, flush=True)
    p(f'Force-field agreement ceiling: CALVADOS-2 vs Mpipi-Recharged hub occupancy, {len(rows)} proteins labelled in both')
    p(f'  ({sum(math.isnan(r["part"]) for r in rows)} proteins with an undefined position-removed correlation are left out of those medians)')
    allr = np.array([r['raw'] for r in rows])
    allp = np.array([r['part'] for r in rows])
    p(f'  all: median Spearman {np.nanmedian(allr):.2f} (IQR {np.nanpercentile(allr, 25):.2f}-{np.nanpercentile(allr, 75):.2f}), '
      f'chain position removed {np.nanmedian(allp):.2f} (IQR {np.nanpercentile(allp, 25):.2f}-{np.nanpercentile(allp, 75):.2f})')
    for name, key in (('group', lambda r: r['group']), ('length', lambda r: length_bin(r['N'])),
                      ('compaction', lambda r: nu_bin(r['nu']))):
        p(f'\nby {name}:')
        strata = defaultdict(list)
        for r in rows:
            strata[key(r)].append(r)
        for k in sorted(strata):
            v = strata[k]
            p(f'  {k:16s} n = {len(v):5d}  median {np.nanmedian([r["raw"] for r in v]):.2f}  '
              f'position removed {np.nanmedian([r["part"] for r in v]):.2f}')
    with open(os.path.join(HERE, 'agreement_ceiling_per_protein.csv'), 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['protein_id', 'group', 'N', 'nu', 'spearman_cal_rech', 'partial_cal_rech'])
        w.writerows([r['id'], r['group'], r['N'], f'{r["nu"]:.4f}', f'{r["raw"]:.4f}', f'{r["part"]:.4f}'] for r in rows)
    open(os.path.join(HERE, 'agreement_ceiling_summary.txt'), 'w').write('\n'.join(out) + '\n')


if __name__ == '__main__':
    main()
