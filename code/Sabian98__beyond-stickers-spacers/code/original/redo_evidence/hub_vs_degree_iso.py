"""
hub_vs_degree_iso.py — POST-HOC robustness of Rule R1 of hub_vs_degree_plan.md (hash 9212f51e), added after Part A showed
a residual that replicates between force fields (median rho +0.264 from a LINEAR regression on degree).
Concern: if occupancy is a nonlinear, monotone function of degree, a linear residual keeps the nonlinear part of degree,
which would replicate between force fields without any bridging component. Here the within-protein dependence on degree
is removed with isotonic (monotone, nonparametric) regression instead, so the residual is orthogonal to every monotone
function of the position-adjusted degree.

  ~/.venvs/bender/bin/python -I hub_vs_degree_iso.py  ->  hub_vs_degree/part_a_iso.txt
Exploratory, labelled as such in the paper if reported.
"""
import csv
import os
from collections import defaultdict

import numpy as np
from scipy.stats import spearmanr
from sklearn.isotonic import IsotonicRegression

HERE = '/Volumes/User Homes/Nu Project/BENDER-DYN/redo_evidence'
OUT = os.path.join(HERE, 'hub_vs_degree')
data = {'CALVADOS-2': {}, 'Mpipi-Recharged': {}}
with open(os.path.join(HERE, 'hub_predictor/data_v3/residues.csv')) as fh:
    for r in csv.DictReader(fh):
        if r['occupancy'] == '' or r['force_field'] not in data:
            continue
        d = data[r['force_field']].setdefault(r['protein_id'], ([], []))
        d[0].append(float(r['occupancy']))
        d[1].append(float(r['deg_mean']))
data = {ff: {p: (np.array(o), np.array(g)) for p, (o, g) in v.items()} for ff, v in data.items()}


def end_bin(n):
    return np.minimum(np.minimum(np.arange(n), n - 1 - np.arange(n)), 30)


def profile(ff, which):
    acc = defaultdict(list)
    for o, g in data[ff].values():
        x = (o, g)[which]
        x = x - x.mean()
        for b, v in zip(end_bin(len(x)), x):
            acc[b].append(v)
    return {b: float(np.mean(v)) for b, v in acc.items()}


def adjust(x, prof):
    y = x - x.mean() - np.array([prof[i] for i in end_bin(len(x))])
    return y - y.mean()


res, adj = {}, {}
for ff in data:
    po, pg = profile(ff, 0), profile(ff, 1)
    adj[ff], res[ff] = {}, {}
    for p, (o, g) in data[ff].items():
        a, b = adjust(o, po), adjust(g, pg)
        adj[ff][p] = (a, b)
        iso = IsotonicRegression(out_of_bounds='clip').fit(b, a)      # occupancy as a monotone function of degree
        res[ff][p] = a - iso.predict(b)
paired = [p for p in sorted(set(adj['CALVADOS-2']) & set(adj['Mpipi-Recharged']))
          if len(adj['CALVADOS-2'][p][0]) == len(adj['Mpipi-Recharged'][p][0]) >= 30]
r_iso = np.array([spearmanr(res['CALVADOS-2'][p], res['Mpipi-Recharged'][p]).correlation for p in paired])
r_iso = r_iso[~np.isnan(r_iso)]
by_len = defaultdict(list)
for p in paired:
    by_len[len(adj['CALVADOS-2'][p][0])].append(p)
ctrl = []
for n, ps in by_len.items():
    if len(ps) < 2:
        continue
    for i, p in enumerate(ps):
        ctrl.append(spearmanr(res['CALVADOS-2'][p], res['Mpipi-Recharged'][ps[(i + 1) % len(ps)]]).correlation)
ctrl = np.array([x for x in ctrl if x == x])
# the isotonic fit also absorbs noise, so report how much variance of occupancy it leaves
leave = np.mean([res['CALVADOS-2'][p].var() / adj['CALVADOS-2'][p][0].var() for p in paired])
boot = np.random.RandomState(0)
ci = np.percentile([np.median(boot.choice(r_iso, len(r_iso))) for _ in range(500)], [2.5, 97.5])
L = ['POST-HOC robustness of R1: residual of occupancy after an isotonic fit on degree (hub_vs_degree_plan.md)', '',
     f'paired proteins {len(paired)}; the isotonic fit leaves {leave * 100:.1f}% of the variance of adjusted occupancy',
     f'cross-force-field Spearman of the isotonic residual: median {np.median(r_iso):+.3f} '
     f'(95% bootstrap interval {ci[0]:+.3f} to {ci[1]:+.3f}); control, same length different protein {np.median(ctrl):+.3f}',
     f'for reference, the linear residual of the plan gave +0.264 and the control +0.042']
open(os.path.join(OUT, 'part_a_iso.txt'), 'w').write('\n'.join(L) + '\n')
print('\n'.join(L))
