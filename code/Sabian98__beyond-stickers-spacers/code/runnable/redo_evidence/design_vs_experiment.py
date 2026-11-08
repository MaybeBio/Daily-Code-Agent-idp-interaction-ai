"""
design_vs_experiment.py — DESCRIPTIVE (tier-1 head-to-head, no pre-registered metric): the simulation of each PED
construct in the design batch, the frozen prediction and the residue-neutral null, each against the hub occupancy of
the experimental PED ensembles of the same construct (phosphoresidues read as S/T/Y, the phosphomimetic Sic1 construct
mapped to pSic1; residues matched by difflib; values averaged over all matching ensembles). Same procedure as the
inline run of 2026-10-09 that produced design_response/design_vs_experiment.txt (Mpipi-Recharged).
  OMP_NUM_THREADS=1 ~/.venvs/bender/bin/python design_vs_experiment.py --labels labels_design_calvados --name CALVADOS-2
      ->  design_response/design_vs_experiment_<name>.txt
"""
import argparse
import csv
import difflib
import os
import pickle
import sys
from collections import defaultdict

import numpy as np
from scipy.stats import rankdata, spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, 'hub_predictor'))
import train_hub_model as Tm  # noqa: E402

MD = {'PED_Ab40', 'PED_ACTR'}
PH = {'x': 'S', 'y': 'T', 'z': 'Y'}


def partial(x, y):
    n = len(x)
    z = np.minimum(np.arange(n), n - 1 - np.arange(n)).astype(float)
    rx, ry, rz = (rankdata(v) for v in (x, y, z))
    A = np.vstack([rz, np.ones_like(rz)]).T
    ex = rx - A @ np.linalg.lstsq(A, rx, rcond=None)[0]
    ey = ry - A @ np.linalg.lstsq(A, ry, rcond=None)[0]
    return float(np.corrcoef(ex, ey)[0, 1])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--labels', required=True)
    ap.add_argument('--name', required=True)
    a = ap.parse_args()
    out = []

    def p(s=''):
        out.append(s)
        print(s, flush=True)
    model = pickle.load(open(os.path.join(HERE, 'hub_predictor', 'frozen_models_2026-10-02.pkl'), 'rb'))['consensus']
    occ, seq, nm = defaultdict(list), defaultdict(str), defaultdict(dict)
    for r in csv.DictReader(open(os.path.join(HERE, 'hub_predictor', 'null_occupancy_multi.csv'))):
        nm[r['entry']][int(r['pos'])] = float(r['occ_null_mean'])
    for f in ('experimental_occupancy.csv', 'fresh_ped_occupancy.csv', 'ped_extra_occupancy.csv'):
        for r in csv.DictReader(open(os.path.join(HERE, 'hub_predictor', f))):
            if r['entry'] in occ and len(occ[r['entry']]) > int(r['pos']):
                continue
            occ[r['entry']].append(float(r['occ_exp']))
            seq[r['entry']] += r['aa']
    plain = {e: ''.join(PH.get(c, c) for c in s) for e, s in seq.items()}
    sim = defaultdict(dict)
    for r in csv.DictReader(open(os.path.join(HERE, a.labels, 'per_residue.csv'))):
        sim[r['UniProt_ID']][int(r['pos'])] = float(r['occupancy'])
    pred = defaultdict(list)
    for r in csv.DictReader(open(os.path.join(HERE, 'gaivi', 'predictions_frozen_2026-10-02.csv'))):
        pred[r['id']].append(float(r['consensus']))
    T = {r['id']: r for r in csv.DictReader(open(os.path.join(HERE, 'gaivi', 'hub_design_table.csv')))
         if r['category'] == 'ped'}
    rows = []
    p(f'Tier-1 head-to-head (descriptive), {a.name}: Spearman with the experimental PED ensembles, raw / chain position')
    p('removed, averaged over all PED ensembles of the same construct')
    p(f'{"construct":11s} {"PED":9s} {"n":>4s} {"sim~exp":>13s} {"model~exp":>13s} {"null~exp":>13s} {"sim~model":>13s}  ens')
    for d, t in sorted(T.items()):
        s = t['sequence']
        S = np.array([sim[d][i] for i in range(len(s))])
        P = np.array(pred[d]) if len(pred[d]) == len(s) else model.predict(Tm.features(s, 10))
        target = plain[t['parent']]
        res = []
        for e in [e for e in occ if plain[e] == target]:
            q = ''.join('E' if (c in 'xyz' and d == 'PED_pSic1E') else PH.get(c, c) for c in seq[e])
            sm = difflib.SequenceMatcher(None, q, s, autojunk=False)
            pr = [(b.a + k, b.b + k) for b in sm.get_matching_blocks() for k in range(b.size)]
            if len(pr) < 0.8 * len(s):
                continue
            ia, ib = np.array([x for x, _ in pr]), np.array([y for _, y in pr])
            X = np.array(occ[e])[ia]
            N = np.array([nm[e][i] for i in ia]) if e in nm else None
            res.append([spearmanr(S[ib], X).correlation, partial(S[ib], X), spearmanr(P[ib], X).correlation,
                        partial(P[ib], X), spearmanr(N, X).correlation if N is not None else np.nan,
                        partial(N, X) if N is not None else np.nan])
        if not res:
            p(f'{d} no matching ensemble')
            continue
        m = np.nanmean(np.array(res), axis=0)
        rows.append((d, m))
        p(f'{d:11s} {t["parent"]:9s} {len(s):4d} {m[0]:+6.2f}/{m[1]:+5.2f} {m[2]:+6.2f}/{m[3]:+5.2f} '
          f'{m[4]:+6.2f}/{m[5]:+5.2f} {spearmanr(S, P).correlation:+6.2f}/{partial(S, P):+5.2f}  {len(res)}'
          + ('  (MD)' if d in MD else ''))
    E = np.array([m for d, m in rows if d not in MD])
    ok = ~np.isnan(E[:, 5])
    p(f'\nmedian over {len(E)} constructs with ensembles built from experimental data (MD excluded):')
    p(f'  {a.name} simulation vs experiment  raw {np.median(E[:, 0]):+.2f}, partial {np.median(E[:, 1]):+.2f}')
    p(f'  frozen model vs experiment          raw {np.median(E[:, 2]):+.2f}, partial {np.median(E[:, 3]):+.2f}')
    p(f'  residue-neutral null vs experiment  raw {np.nanmedian(E[:, 4]):+.2f}, partial {np.nanmedian(E[:, 5]):+.2f}'
      f'  (n = {int(ok.sum())})')
    p(f'  simulation partial > null partial in {int((E[ok, 1] > E[ok, 5]).sum())}/{int(ok.sum())}; model partial > '
      f'simulation partial in {int((E[:, 3] > E[:, 1]).sum())}/{len(E)}; simulation partial > 0 in '
      f'{int((E[:, 1] > 0).sum())}/{len(E)}')
    open(os.path.join(HERE, 'design_response', f'design_vs_experiment_{a.name}.txt'), 'w').write('\n'.join(out) + '\n')


if __name__ == '__main__':
    main()
