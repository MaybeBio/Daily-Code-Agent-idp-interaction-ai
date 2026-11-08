"""
ped_v2_score.py — scores for the PED v2 extension as pre-registered in ped_v2_plan.md (with Amendment 1).
Inputs: gaivi/predictions_frozen_2026-10-02_ped_v2.csv (frozen consensus predictions), hub_predictor/ped_v2_occupancy.csv
(labels and five-pool null), ped_v2_set.csv. Set C is compared with the earlier coarse-grained labels where available.
  ~/.venvs/bender/bin/python ped_v2_score.py  ->  ped_v2_scores.csv, ped_v2_summary.txt
"""
import csv
import itertools
import math
import os
from collections import defaultdict

import numpy as np
from scipy.stats import binomtest, rankdata, spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
NOVEL = {'Sendai NTAIL', 'KISS-1', 'MKK7', 'E1A 36-146', 'N-WASP domain V', 'Tif2 NRID', 'AR Tau-5 R2-R3', 'Tir 388-550',
         'SARS-CoV-2 N linker'}
EXCLUDED = {'IAPP'}
out = []


def p(s=''):
    out.append(s)
    print(s, flush=True)


def partial(x, y, z):
    rx, ry, rz = (rankdata(v) for v in (x, y, z))
    A = np.vstack([rz, np.ones_like(rz)]).T
    ex = rx - A @ np.linalg.lstsq(A, rx, rcond=None)[0]
    ey = ry - A @ np.linalg.lstsq(A, ry, rcond=None)[0]
    return float(np.corrcoef(ex, ey)[0, 1])


def main():
    pred = defaultdict(list)
    for r in csv.DictReader(open(os.path.join(HERE, 'gaivi', 'predictions_frozen_2026-10-02_ped_v2.csv'))):
        pred[r['id']].append(float(r['consensus']))
    lab = defaultdict(list)
    for r in csv.DictReader(open(os.path.join(HERE, 'hub_predictor', 'ped_v2_occupancy.csv'))):
        lab[r['entry']].append(r)
    meta = {r['id']: r for r in csv.DictReader(open(os.path.join(HERE, 'ped_v2_set.csv')))}
    rows = []
    for e, rs in lab.items():
        n = len(rs)
        oe = np.array([float(r['occ_exp']) for r in rs])
        on = np.array([float(r['occ_null_mean']) for r in rs])
        pr = np.array(pred[e])
        assert len(pr) == n, e
        end = np.minimum(np.arange(n), n - 1 - np.arange(n)).astype(float)
        rows.append(dict(set=meta[e]['set'], entry=e, protein=meta[e]['protein'], N=n, models=int(rs[0]['models']),
                         rho_model=spearmanr(pr, oe).correlation, rho_null=spearmanr(on, oe).correlation,
                         rho_position=spearmanr(end, oe).correlation, partial_model=partial(pr, oe, end),
                         partial_null=partial(on, oe, end), occ=oe))
    with open(os.path.join(HERE, 'ped_v2_scores.csv'), 'w', newline='') as fh:
        keys = ['set', 'entry', 'protein', 'N', 'models', 'rho_model', 'rho_null', 'rho_position', 'partial_model',
                'partial_null']
        w = csv.writer(fh)
        w.writerow(keys)
        for r in sorted(rows, key=lambda r: (r['set'], r['entry'])):
            w.writerow([r[k] if not isinstance(r[k], float) else f'{r[k]:.4f}' for k in keys])

    p('PED v2 extension: frozen consensus predictions (frozen_models_2026-10-02.pkl) against hub occupancy')
    for st in 'AB':
        p(f'\nSet {st}')
        p(f'  {"entry":9s} {"protein":22s} {"N":>4s} {"models":>6s} {"rho mod":>8s} {"null":>6s} {"pos":>6s} {"partial mod":>11s} {"null":>6s}')
        for r in sorted((r for r in rows if r['set'] == st), key=lambda r: -r['partial_model']):
            tag = ' [excluded]' if r['protein'] in EXCLUDED else ('' if st != 'A' or r['protein'] in NOVEL else ' [ortholog]')
            p(f'  {r["entry"]:9s} {(r["protein"] + tag)[:30]:30s} {r["N"]:4d} {r["models"]:6d} {r["rho_model"]:+8.2f} '
              f'{r["rho_null"]:+6.2f} {r["rho_position"]:+6.2f} {r["partial_model"]:+11.2f} {r["partial_null"]:+6.2f}')

    def per_protein(sel):
        g = defaultdict(list)
        for r in sel:
            g[r['protein']].append(r)
        return {k: {m: float(np.mean([r[m] for r in v])) for m in
                    ('rho_model', 'rho_null', 'rho_position', 'partial_model', 'partial_null')} for k, v in g.items()}

    for label, sel in (('PRIMARY: set A without IAPP', [r for r in rows if r['set'] == 'A' and r['protein'] not in EXCLUDED]),
                       ('SENSITIVITY: 9 strictly novel proteins', [r for r in rows if r['set'] == 'A' and r['protein'] in NOVEL]),
                       ('IAPP (excluded by AIUPred)', [r for r in rows if r['protein'] in EXCLUDED])):
        pp = per_protein(sel)
        if not pp:
            continue
        v = list(pp.values())
        med = {m: float(np.median([x[m] for x in v])) for m in v[0]}
        k_pos = sum(x['partial_model'] > 0 for x in v)
        k_null = sum(x['partial_model'] > x['partial_null'] for x in v)
        n = len(v)
        p(f'\n{label}: {n} proteins')
        p(f'  median Spearman: model {med["rho_model"]:+.2f}, null {med["rho_null"]:+.2f}, position {med["rho_position"]:+.2f}')
        p(f'  median partial:  model {med["partial_model"]:+.2f}, null {med["partial_null"]:+.2f}')
        p(f'  partial > 0 in {k_pos}/{n} (one-sided sign test P = {binomtest(k_pos, n, 0.5, "greater").pvalue:.4f}); '
          f'model > null in {k_null}/{n} (P = {binomtest(k_null, n, 0.5, "greater").pvalue:.4f})')
        if label.startswith('PRIMARY'):
            ok = med['partial_model'] > med['partial_null'] and binomtest(k_null, n, 0.5, 'greater').pvalue < 0.05
            p(f'  pre-registered decision rule met: {ok}')

    # set C: all-atom force-field concordance
    p('\nSet C: hub occupancy agreement between all-atom force fields (unweighted deposited frames)')
    byprot = defaultdict(dict)
    for r in rows:
        if r['set'] == 'C':
            prot, ff = r['protein'].rsplit(' ', 1)
            byprot[prot][ff] = r
    allp, allr = [], []
    for prot, d in sorted(byprot.items()):
        n = len(next(iter(d.values()))['occ'])
        end = np.minimum(np.arange(n), n - 1 - np.arange(n)).astype(float)
        for a, b in itertools.combinations(sorted(d), 2):
            x, y = d[a]['occ'], d[b]['occ']
            rho, pa = spearmanr(x, y).correlation, partial(x, y, end)
            allr.append(rho)
            allp.append(pa)
            p(f'  {prot:12s} {a:10s} vs {b:10s}  Spearman {rho:+.2f}  partial {pa:+.2f}')
        for ff, r in sorted(d.items()):
            p(f'  {prot:12s} {ff:10s} model {r["rho_model"]:+.2f} (partial {r["partial_model"]:+.2f}), null {r["rho_null"]:+.2f} '
              f'(partial {r["partial_null"]:+.2f})')
    if allr:
        p(f'  median over {len(allr)} pairs: Spearman {np.median(allr):+.2f}, partial {np.median(allp):+.2f}')
    open(os.path.join(HERE, 'ped_v2_summary.txt'), 'w').write('\n'.join(out) + '\n')


if __name__ == '__main__':
    main()
