"""
hub_vs_degree_sim.py — Part A of hub_vs_degree_plan.md (hash 9212f51e). Labels only, no model.
  A1  within-protein position-adjusted Spearman(occupancy, degree), three force fields, by group and length tercile
  A2  the part of occupancy that degree does not explain, and its replication between CALVADOS-2 and Mpipi-Recharged
  A3  residue-type composition of occupancy-only hubs against degree-only hubs

  ~/.venvs/bender/bin/python -I hub_vs_degree_sim.py  ->  hub_vs_degree/part_a.txt
"""
import csv
import os
from collections import defaultdict

import numpy as np
from scipy.stats import spearmanr

HERE = '/Volumes/User Homes/Nu Project/BENDER-DYN/redo_evidence'
OUT = os.path.join(HERE, 'hub_vs_degree')
os.makedirs(OUT, exist_ok=True)
FFS = ('CALVADOS-2', 'Mpipi-Recharged', 'MPIPI-GG')

meta = {}
for r in csv.DictReader(open(os.path.join(HERE, 'hub_predictor/data_v3/proteins.csv'))):
    meta[(r['force_field'], r['protein_id'])] = (r['group'], int(r['N']))

data = {ff: {} for ff in FFS}
with open(os.path.join(HERE, 'hub_predictor/data_v3/residues.csv')) as fh:
    for r in csv.DictReader(fh):
        if r['occupancy'] == '' or r['force_field'] not in data:
            continue
        d = data[r['force_field']].setdefault(r['protein_id'], ([], [], []))
        d[0].append(float(r['occupancy']))
        d[1].append(float(r['deg_mean']))
        d[2].append(r['aa'])
data = {ff: {p: (np.array(o), np.array(g), ''.join(a)) for p, (o, g, a) in v.items()} for ff, v in data.items()}


def end_bin(n):
    return np.minimum(np.minimum(np.arange(n), n - 1 - np.arange(n)), 30)


def profile(ff, which):
    """Pooled mean of the within-protein-centred quantity at each end-distance bin."""
    acc = defaultdict(list)
    for o, g, a in data[ff].values():
        x = (o, g)[which]
        x = x - x.mean()
        for b, v in zip(end_bin(len(x)), x):
            acc[b].append(v)
    return {b: float(np.mean(v)) for b, v in acc.items()}


def adjust(x, prof):
    b = end_bin(len(x))
    y = x - x.mean() - np.array([prof[i] for i in b])
    return y - y.mean()


L = ['betweenness hub occupancy against contact degree, Part A (hub_vs_degree_plan.md)', '']
adj = {}
for ff in FFS:
    po, pg = profile(ff, 0), profile(ff, 1)
    adj[ff] = {p: (adjust(o, po), adjust(g, pg), a) for p, (o, g, a) in data[ff].items()}

# A1
L.append('A1  within-protein position-adjusted Spearman(occupancy, degree), median over proteins with n >= 30')
for ff in FFS:
    rho = {p: spearmanr(o, g).correlation for p, (o, g, a) in adj[ff].items() if len(o) >= 30}
    v = np.array(list(rho.values()))
    L.append(f'  {ff:16s} n {len(v):5d}  median {np.median(v):+.3f}  IQR {np.percentile(v, 25):+.3f} to {np.percentile(v, 75):+.3f}')
    byg = defaultdict(list)
    for p, r in rho.items():
        byg[meta[(ff, p)][0]].append(r)
    L.append('     by group: ' + '  '.join(f'{g} {np.median(v):+.3f}' for g, v in sorted(byg.items())))
    ns = np.array([meta[(ff, p)][1] for p in rho])
    q = np.percentile(ns, [33.3, 66.7])
    rr = np.array(list(rho.values()))
    L.append(f'     by length tercile (cut at {q[0]:.0f} and {q[1]:.0f} residues): '
             f'{np.median(rr[ns <= q[0]]):+.3f}  {np.median(rr[(ns > q[0]) & (ns <= q[1])]):+.3f}  {np.median(rr[ns > q[1]]):+.3f}')
L.append('')

# A2
L.append('A2  residual of occupancy after regressing on degree, within protein; replication between CALVADOS-2 and Recharged')
resid = {}
for ff in ('CALVADOS-2', 'Mpipi-Recharged'):
    resid[ff] = {}
    for p, (o, g, a) in adj[ff].items():
        beta = float(np.dot(o, g) / np.dot(g, g)) if np.dot(g, g) > 0 else 0.0
        resid[ff][p] = o - beta * g
paired = sorted(set(adj['CALVADOS-2']) & set(adj['Mpipi-Recharged']))
paired = [p for p in paired if len(adj['CALVADOS-2'][p][0]) == len(adj['Mpipi-Recharged'][p][0]) >= 30]


def med_pair(fx):
    v = [spearmanr(*fx(p)).correlation for p in paired]
    v = np.array([x for x in v if x == x])
    return v


r_res = med_pair(lambda p: (resid['CALVADOS-2'][p], resid['Mpipi-Recharged'][p]))
r_occ = med_pair(lambda p: (adj['CALVADOS-2'][p][0], adj['Mpipi-Recharged'][p][0]))
r_deg = med_pair(lambda p: (adj['CALVADOS-2'][p][1], adj['Mpipi-Recharged'][p][1]))
# control: pair each protein's CALVADOS-2 residual with the Recharged residual of another protein of the same length
by_len = defaultdict(list)
for p in paired:
    by_len[len(adj['CALVADOS-2'][p][0])].append(p)
ctrl = []
for n, ps in by_len.items():
    if len(ps) < 2:
        continue
    for i, p in enumerate(ps):
        q = ps[(i + 1) % len(ps)]
        ctrl.append(spearmanr(resid['CALVADOS-2'][p], resid['Mpipi-Recharged'][q]).correlation)
ctrl = np.array([x for x in ctrl if x == x])
boot = np.random.RandomState(0)
ci = np.percentile([np.median(boot.choice(r_res, len(r_res))) for _ in range(500)], [2.5, 97.5])
L.append(f'  paired proteins {len(paired)}')
L.append(f'  cross-force-field Spearman, median: position-adjusted occupancy {np.median(r_occ):+.3f}, '
         f'degree {np.median(r_deg):+.3f}, RESIDUAL {np.median(r_res):+.3f} (95% bootstrap interval {ci[0]:+.3f} to {ci[1]:+.3f})')
L.append(f'  control, same-length different protein (n {len(ctrl)}): median {np.median(ctrl):+.3f}')
m = np.median(r_res)
verdict = ('REPLICATED bridging component' if m >= 0.20 and m > 2 * np.median(ctrl)
           else 'ABSENT' if m < 0.10 else 'AMBIGUOUS')
L.append(f'  Rule R1 (replicated if median >= 0.20 and above twice the control; absent if < 0.10): {verdict}')
L.append('')

# A3
L.append('A3  residue types among occupancy-only hubs and degree-only hubs, enrichment over all residues')
for ff in ('CALVADOS-2', 'Mpipi-Recharged'):
    cnt = {k: defaultdict(int) for k in ('all', 'occ_only', 'deg_only')}
    for p, (o, g, a) in data[ff].items():
        n = len(o)
        k = max(1, int(round(0.1 * n)))
        so = set(np.lexsort((np.arange(n), -o))[:k])
        sg = set(np.lexsort((np.arange(n), -g))[:k])
        for i, aa in enumerate(a):
            cnt['all'][aa] += 1
            if i in so and i not in sg:
                cnt['occ_only'][aa] += 1
            if i in sg and i not in so:
                cnt['deg_only'][aa] += 1
    tot = {k: sum(v.values()) for k, v in cnt.items()}
    L.append(f'  {ff}: occupancy-only hubs {tot["occ_only"]}, degree-only hubs {tot["deg_only"]}, '
             f'of {tot["all"]} residues')
    for key in ('occ_only', 'deg_only'):
        e = {aa: (cnt[key][aa] / tot[key]) / (cnt['all'][aa] / tot['all']) for aa in cnt['all']}
        top = sorted(e.items(), key=lambda t: -t[1])
        L.append(f'     {key:9s} most enriched ' + ' '.join(f'{a} {v:.2f}' for a, v in top[:5]) +
                 '   most depleted ' + ' '.join(f'{a} {v:.2f}' for a, v in top[-5:]))
open(os.path.join(OUT, 'part_a.txt'), 'w').write('\n'.join(L) + '\n')
print('\n'.join(L))
