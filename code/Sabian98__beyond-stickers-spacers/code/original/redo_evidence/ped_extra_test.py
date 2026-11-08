"""
ped_extra_test.py — (1) the experimental ceiling: how well do independent PED ensembles of the SAME protein agree in hub
occupancy (drkN x 4, Sic1/pSic1 x 5)? (2) blind scores for the three PED ensembles added on 2026-10-03 (PED00158 drkN
unfolded state, PED00489 5-phospho 4E-BP2, PED00001 pSic1), whose predictions were frozen before these labels
(gaivi/predictions_frozen_2026-10-02_ped_extra.csv).
Labels use the code-review tie-breaking (hub_frames.top_hubs); the null is averaged over five independent pools.
Local copies of the PED files are in ped_local/.
  ~/.venvs/bender/bin/python ped_extra_test.py  ->  ped_extra_summary.txt, hub_predictor/ped_extra_occupancy.csv
"""
import csv
import hashlib
import itertools
import math
import os
import pickle
import sys
from multiprocessing import Pool

import numpy as np
from scipy.stats import rankdata, spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [HERE, os.path.join(HERE, 'hub_predictor')]
import exp_occupancy as eo  # noqa: E402
import ped_null as pn  # noqa: E402

LOCAL = os.path.join(HERE, 'ped_local')
FROZEN = os.path.join(HERE, 'gaivi', 'predictions_frozen_2026-10-02_ped_extra.csv')
GROUPS = {'drkN': [('PED00022', 'PED00022/PED00022e001.pdb', 'NOE/PRE'),
                   ('PED00156', 'PED00156/PED00156e002.pdb', 'NOE/PRE/smFRET'),
                   ('PED00157', 'PED00157/PED00157e002.pdb', 'NOE/PRE/RDC/SAXS'),
                   ('PED00158', 'PED00158/PED00158e006.pdb', 'NOE/PRE/smFRET, unfolded state')],
          'Sic1': [('PED00023', 'PED00023/PED00023e001.pdb', 'PRE, 16 models'),
                   ('PED00159', 'PED00159/PED00159e002.pdb', 'PRE/SAXS'),
                   ('PED00160', 'PED00160/PED00160e002.pdb', 'PRE/SAXS/CS'),
                   ('PED00161', 'PED00161/PED00161e002.pdb', 'PRE/SAXS, pSic1'),
                   ('PED00001', 'PED00001/PED00001e003.pdb', 'PRE, pSic1, 11 models')]}
NEW = {'PED00158', 'PED00489', 'PED00001'}
K_NULL = 5
out = []


def p(s=''):
    out.append(s)
    print(s, flush=True)


def partial(x, y):
    n = len(y)
    z = np.minimum(np.arange(n), n - 1 - np.arange(n)).astype(float)
    rx, ry, rz = (rankdata(v) for v in (x, y, z))
    A = np.vstack([rz, np.ones_like(rz)]).T
    ex = rx - A @ np.linalg.lstsq(A, rx, rcond=None)[0]
    ey = ry - A @ np.linalg.lstsq(A, ry, rcond=None)[0]
    return float(np.corrcoef(ex, ey)[0, 1])


def label(job):
    entry, path, with_null = job
    seq, models = pn.load(os.path.join(LOCAL, path))
    if len(models) > 1000:
        models = models[::math.ceil(len(models) / 1000)]
    oe, used = eo.occupancy(seq, models)
    nulls = []
    if with_null:
        for k in range(K_NULL):
            nm = pn.null_models(seq, models, 70000 + 97 * k + sum(map(ord, entry)))[0]
            nulls.append(eo.occupancy(seq, nm)[0])
    return entry, seq, np.array(oe), used, (np.mean(nulls, axis=0) if nulls else None)


def main():
    h = hashlib.sha256(open(FROZEN, 'rb').read()).hexdigest()
    if h not in open(os.path.join(HERE, 'gaivi', 'FROZEN.sha256')).read():
        raise SystemExit('frozen predictions for the extra ensembles are not registered; refusing to compute labels')
    jobs = [(e, f, True) for g in GROUPS.values() for e, f, _ in g]
    jobs.append(('PED00489', 'PED00489/pdbfile.pdb', True))
    with Pool(8) as pool:
        res = {r[0]: r for r in pool.imap_unordered(label, jobs)}
    models = pickle.load(open(os.path.join(HERE, 'hub_predictor', 'frozen_models_2026-10-02.pkl'), 'rb'))
    import train_hub_model as T
    frozen = {}
    for r in csv.DictReader(open(FROZEN)):
        frozen.setdefault(r['id'], []).append(float(r['consensus']))

    p('Experimental ceiling: hub occupancy agreement (Spearman) between independent PED ensembles of the same protein')
    p('(code-review tie-breaking; PED00158 and PED00001 added on 2026-10-03)\n')
    ceiling = []
    for name, ents in GROUPS.items():
        p(f'  {name}')
        for (a, _, ra), (b, _, rb) in itertools.combinations(ents, 2):
            r = spearmanr(res[a][2], res[b][2]).correlation
            ceiling.append(r)
            p(f'    {a} ({ra}) vs {b} ({rb}): {r:+.2f}')
    p(f'\n  median over {len(ceiling)} pairs: {np.median(ceiling):+.2f} (range {min(ceiling):+.2f} to {max(ceiling):+.2f})')

    p('\nModel (frozen consensus) against each ensemble: Spearman / partial on end distance; null and position for reference')
    rows = []
    for e in [x for g in GROUPS.values() for x, _, _ in g] + ['PED00489']:
        _, seq, oe, used, on = res[e]
        n = len(seq)
        pr = np.array(frozen[e]) if e in frozen else models['consensus'].predict(T.features(seq, 10))
        end = np.minimum(np.arange(n), n - 1 - np.arange(n)).astype(float)
        r = dict(model=spearmanr(pr, oe).correlation, model_p=partial(pr, oe), null=spearmanr(on, oe).correlation,
                 null_p=partial(on, oe), position=spearmanr(end, oe).correlation)
        tag = 'blind' if e in NEW else 'development'
        p(f'  {e} [{tag:11s}] N {n:3d} models {used:4d}: model {r["model"]:+.2f} (partial {r["model_p"]:+.2f}) | '
          f'null {r["null"]:+.2f} (partial {r["null_p"]:+.2f}) | position {r["position"]:+.2f}')
        for i in range(n):
            rows.append([e, tag, i, seq[i], f'{oe[i]:.4f}', f'{on[i]:.4f}', f'{pr[i]:.5f}'])
    with open(os.path.join(HERE, 'hub_predictor', 'ped_extra_occupancy.csv'), 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['entry', 'set', 'pos', 'aa', 'occ_exp', 'occ_null_mean5', 'pred_consensus'])
        w.writerows(rows)
    open(os.path.join(HERE, 'ped_extra_summary.txt'), 'w').write('\n'.join(out) + '\n')


if __name__ == '__main__':
    main()
