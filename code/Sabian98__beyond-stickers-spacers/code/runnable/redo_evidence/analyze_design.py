"""
analyze_design.py — prospective test on the GAIVI design batch (gaivi/hub_design_table.csv): simulated hub occupancy
(labels from gaivi/hub_labels_gaivi.py, one folder per force field) against predictions frozen before any design label
existed. Metrics declared before the labels (methods_paper/RECOMB_PLAN.md): raw Spearman and Spearman partial on
distance from the chain ends, always next to the position-only baseline.

  sequences    ped, pre, variant, host, replica: per-sequence raw and position-partial Spearman
  mutations    variants vs their parent construct: per-residue change in occupancy, predicted vs simulated (Spearman
               over residues), and the mean change at the mutated sites (sign agreement, correlation across variants)
  guest-host   change at the guest position (guest - host) across the 19 guests of each host, predicted vs simulated;
               context: the host-specific part (guest x host interaction after removing guest and host means)
  series       valence (mean occupancy at the inserted residues), pairs (occupancy of Y and partner), kappa and
               homopolymers (per-sequence Spearman)
  replicas     occupancy in the new run vs the original BENDER run of the same sequence and force field (between-run
               noise ceiling)
Labels for the force-field mean are used where all given force fields have the sequence.

  ~/.venvs/bender/bin/python analyze_design.py --labels CALVADOS-2=labels_design_calvados \
        --labels Mpipi-Recharged=labels_design_recharged --labels MPIPI-GG=labels_design_mpipi \
        [--pred gaivi/predictions_frozen_2026-10-02.csv] [--column consensus]
  --mock makes noisy fake labels from the predictions, to test the code path only.
"""
import argparse
import csv
import os
from collections import defaultdict

import numpy as np
from scipy.stats import rankdata, spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
GAIVI = os.path.join(HERE, 'gaivi')
FRESH = {'PED00501': 'PED_V2Rpm', 'PED00502': 'PED_GHSRpm', 'PED00503': 'PED_B2ARpm', 'PED00154': 'PED_MKK4',
         'PED00017': 'PED_TauF', 'PED00533': 'PED_Ab40', 'PED00216': 'PED_CBPID4', 'PED00536': 'PED_ACTR'}
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


def rho(x, y):
    return float(spearmanr(x, y).correlation)


def end_dist(n):
    return np.minimum(np.arange(n), n - 1 - np.arange(n)).astype(float)


def load_preds(paths, column):
    pred = {}
    for path in paths:
        for r in csv.DictReader(open(path)):
            uid = FRESH.get(r['id'], r['id'])
            pred.setdefault(uid, {})[int(r['pos'])] = float(r[column])
    return {u: np.array([d[i] for i in range(len(d))]) for u, d in pred.items()}


def load_labels(d):
    lab = defaultdict(dict)
    for r in csv.DictReader(open(os.path.join(d, 'per_residue.csv'))):
        lab[r['UniProt_ID']][int(r['pos'])] = float(r['occupancy'])
    return {u: np.array([v[i] for i in range(len(v))]) for u, v in lab.items()}


def bender_labels():
    lab = defaultdict(dict)
    for r in csv.DictReader(open(os.path.join(HERE, 'hub_predictor', 'data_full', 'residues.csv'))):
        lab[(r['force_field'], r['protein_id'])][int(r['pos'])] = float(r['occupancy'])
    return {k: np.array([v[i] for i in range(len(v))]) for k, v in lab.items()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--labels', action='append', default=[], help='FF=DIR (per_residue.csv from hub_labels_gaivi.py)')
    ap.add_argument('--pred', action='append', default=[])
    ap.add_argument('--column', default='consensus')
    ap.add_argument('--mock', action='store_true')
    a = ap.parse_args()
    table = {r['id']: r for r in csv.DictReader(open(os.path.join(GAIVI, 'hub_design_table.csv')))}
    preds = load_preds(a.pred or [os.path.join(GAIVI, f) for f in ('predictions_frozen_2026-10-02.csv',
                                                                   'predictions_frozen_2026-10-02_ped_fresh.csv',
                                                                   'predictions_frozen_2026-10-02_peptone_pre.csv')],
                       a.column)
    for u, r in table.items():                      # replica runs share the parent's (identical) sequence and prediction
        if r['category'] == 'gh_replica' and r['parent'] in preds:
            preds[u] = preds[r['parent']]
    missing = [u for u in table if u not in preds]
    p(f'Design batch: {len(table)} sequences; predictions ({a.column}) for {len(table) - len(missing)}'
      + (f'; missing {missing[:5]}...' if missing else ''))
    labels = {}
    if a.mock:
        rng = np.random.default_rng(0)
        labels['MOCK'] = {u: np.clip(v + rng.normal(0, 0.05, len(v)), 0, 1) for u, v in preds.items()}
        p('MOCK labels (prediction + noise): tests the code path only, numbers are meaningless')
    for spec in a.labels:
        ff, d = spec.split('=', 1)
        labels[ff] = load_labels(d if os.path.isabs(d) else os.path.join(HERE, d))
        p(f'  {ff}: labels for {sum(u in labels[ff] for u in table)} of {len(table)} design sequences')
    if len(labels) > 1:
        common = set.intersection(*[set(v) for v in labels.values()])
        labels['FF mean'] = {u: np.mean([labels[f][u] for f in labels if f != 'FF mean'], axis=0) for u in common}
    rows = []
    for ff, L in labels.items():
        p(f'\n== {ff} ==')
        # per-sequence accuracy
        cats = defaultdict(list)
        for u, r in table.items():
            if u in L and u in preds and len(L[u]) == len(preds[u]) and np.std(L[u]) > 0:
                n = len(L[u])
                e = end_dist(n)
                rr = (rho(preds[u], L[u]), partial(preds[u], L[u]), rho(e, L[u]))
                cats[r['category']].append(rr)
                rows.append([ff, 'sequence', u, r['category'], n, *[f'{x:.4f}' for x in rr]])
        p('  per sequence: median raw rho (model / position only), position-partial rho (model); * = in the training set')
        for c, v in sorted(cats.items()):
            v = np.array(v)
            if c == 'replica':     # BENDER sequences: in the training set, so their accuracy is in-sample
                c = 'replica*'
            p(f'    {c:12s} n={len(v):3d}  raw {np.nanmedian(v[:, 0]):+.2f} / {np.nanmedian(v[:, 2]):+.2f}   '
              f'partial {np.nanmedian(v[:, 1]):+.2f}  (> 0 in {int((v[:, 1] > 0).sum())}/{len(v)})')
        # mutation effects
        site, glob = [], []
        for u, r in table.items():
            if r['category'] != 'variant' or u not in L or r['parent'] not in L:
                continue
            ds, dp = L[u] - L[r['parent']], preds[u] - preds[r['parent']]
            pos = [int(x) for x in r['positions'].split()] if r['positions'] else []
            if np.std(ds) > 0:
                rows.append([ff, 'mutation', u, 'variant', len(ds), f'{rho(dp, ds):.4f}', '', ''])
                glob.append(rho(dp, ds))
            if pos:
                site.append((float(np.mean(dp[pos])), float(np.mean(ds[pos]))))
        if glob:
            p(f'  mutations: per-variant Spearman of predicted vs simulated per-residue change: median {np.nanmedian(glob):+.2f} '
              f'(> 0 in {sum(g > 0 for g in glob)}/{len(glob)})')
        if len(site) > 2:
            s = np.array(site)
            sign = np.mean(np.sign(s[:, 0]) == np.sign(s[:, 1]))
            p(f'    mean change at mutated sites: Spearman across {len(s)} variants {rho(s[:, 0], s[:, 1]):+.2f}; sign agreement {sign:.0%}')
        # guest-host
        def runs(uid):                              # the original run plus any independent replicas (tier 4)
            return [L[x][50] for x in (uid, f'{uid}_r2', f'{uid}_r3') if x in L]
        M, se = {}, []
        for h in ('GS', 'AVG', 'CHG', 'PLD'):
            hid = f'HOST_{h}'
            if hid not in L:
                continue
            hr = runs(hid)
            for u, r in table.items():
                if r['category'] == 'guest' and r['host'] == h and r['positions'] == '50' and u in L:
                    gr = runs(u)
                    M[(h, r['guest'])] = (preds[u][50] - preds[hid][50], np.mean(gr) - np.mean(hr))
                    if len(gr) > 1 and len(hr) > 1:
                        se.append(np.sqrt(np.var(gr, ddof=1) / len(gr) + np.var(hr, ddof=1) / len(hr)))
        if M:
            v = np.array(list(M.values()))
            p(f'  guest effects, all hosts pooled: Spearman across {len(v)} host-guest pairs {rho(v[:, 0], v[:, 1]):+.2f}'
              + (f'; between-run SE of the simulated change: median {np.median(se):.4f} '
                 f'(|change| > 2 SE for {int(np.sum(np.abs(v[:, 1]) > 2 * np.median(se)))}/{len(v)})' if se else
                 '; single runs only (no between-run SE)'))
            for h in ('GS', 'AVG', 'CHG', 'PLD'):
                v = np.array([M[k] for k in M if k[0] == h])
                if len(v) > 3:
                    p(f'  guest effects in {h:3s}: Spearman across {len(v)} guests, predicted vs simulated {rho(v[:, 0], v[:, 1]):+.2f}')
            hosts = sorted({k[0] for k in M})
            aas = sorted({k[1] for k in M})
            if len(hosts) == 4:
                P = np.array([[M.get((h, g), (np.nan, np.nan))[0] for h in hosts] for g in aas])
                S = np.array([[M.get((h, g), (np.nan, np.nan))[1] for h in hosts] for g in aas])
                ok = ~np.isnan(P).any(1) & ~np.isnan(S).any(1)
                P, S = P[ok], S[ok]
                Pi = P - P.mean(1, keepdims=True) - P.mean(0, keepdims=True) + P.mean()
                Si = S - S.mean(1, keepdims=True) - S.mean(0, keepdims=True) + S.mean()
                p(f'    context (guest x host interaction): Spearman predicted vs simulated {rho(Pi.ravel(), Si.ravel()):+.2f}; '
                  f'share of simulated variance that is context-specific {np.var(Si) / np.var(S - S.mean()):.0%}')
        # series
        for cat in ('valence', 'pair'):
            v = []
            for u, r in table.items():
                if r['category'] == cat and u in L:
                    pos = [int(x) for x in r['positions'].split()]
                    v.append((float(np.mean(preds[u][pos])), float(np.mean(L[u][pos]))))
            if len(v) > 3:
                v = np.array(v)
                p(f'  {cat}: mean occupancy at the inserted residues, Spearman across {len(v)} designs {rho(v[:, 0], v[:, 1]):+.2f}')
        # replicas
        if ff in ('CALVADOS-2', 'MPIPI-GG', 'Mpipi-Recharged'):
            B = bender_labels()
            rr = []
            for u, r in table.items():
                if r['category'] == 'replica' and u in L and (ff, r['parent']) in B and len(B[(ff, r['parent'])]) == len(L[u]):
                    rr.append(rho(L[u], B[(ff, r['parent'])]))
            if rr:
                p(f'  replicas: new run vs original run, median Spearman {np.median(rr):.2f} ({len(rr)} sequences) '
                  f'= between-run noise ceiling')
    with open(os.path.join(HERE, 'design_per_item.csv'), 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['force_field', 'test', 'id', 'category', 'N', 'rho_model', 'partial_model', 'rho_position'])
        w.writerows(rows)
    open(os.path.join(HERE, 'design_analysis_summary.txt'), 'w').write('\n'.join(out) + '\n')


if __name__ == '__main__':
    main()
