"""
fig2_crossff_export.py — per-protein values behind hub_vs_degree/part_a.txt A2 (cross-force-field Spearman between
CALVADOS-2 and Mpipi-Recharged of position-adjusted hub occupancy and degree), for plotting in Fig. 2. The definitions
(end-distance profile, adjustment, pairing, n >= 30) are copied verbatim from redo_evidence/hub_vs_degree_sim.py
(hash 9212f51e plan); no new statistics. Writes cache/crossff_per_protein.csv (A2, medians +0.616, +0.691) and
cache/within_per_protein.csv (A1, within-protein Spearman of the two readouts, medians +0.935 / +0.953 / +0.964) and
cache/a3_enrichment.csv (A3, residue-type enrichment among occupancy-only and degree-only top-10% residues, all 20 types;
check CALVADOS-2 Y 2.08 occupancy-only, G 4.04 degree-only) and cache/a2_residual.csv (A2 residual of occupancy
after degree, cross-force-field Spearman per protein and the same-length control; check medians +0.264 / +0.042).
  ~/.venvs/finches/bin/python fig2_crossff_export.py
"""
import csv
import os
from collections import defaultdict

import numpy as np
from scipy.stats import spearmanr

HERE = '/Volumes/User Homes/Nu Project/BENDER-DYN/redo_evidence'
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'cache', 'crossff_per_protein.csv')
FFS = ('CALVADOS-2', 'Mpipi-Recharged', 'MPIPI-GG')

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


adj = {}
for ff in FFS:
    po, pg = profile(ff, 0), profile(ff, 1)
    adj[ff] = {p: (adjust(o, po), adjust(g, pg), a) for p, (o, g, a) in data[ff].items()}
paired = sorted(set(adj['CALVADOS-2']) & set(adj['Mpipi-Recharged']))
paired = [p for p in paired if len(adj['CALVADOS-2'][p][0]) == len(adj['Mpipi-Recharged'][p][0]) >= 30]
rows = []
for p in paired:
    ro = spearmanr(adj['CALVADOS-2'][p][0], adj['Mpipi-Recharged'][p][0]).correlation
    rd = spearmanr(adj['CALVADOS-2'][p][1], adj['Mpipi-Recharged'][p][1]).correlation
    rows.append((p, len(adj['CALVADOS-2'][p][0]), ro, rd))
with open(OUT, 'w') as fh:
    w = csv.writer(fh)
    w.writerow(['protein_id', 'N', 'rho_occupancy', 'rho_degree'])
    w.writerows(rows)
ro = np.array([r[2] for r in rows if r[2] == r[2]])
rd = np.array([r[3] for r in rows if r[3] == r[3]])
print(f'paired {len(rows)}  median occupancy {np.median(ro):+.3f}  degree {np.median(rd):+.3f}  '
      f'degree higher in {np.mean([r[3] > r[2] for r in rows]) * 100:.0f}%')
with open(os.path.join(os.path.dirname(OUT), 'within_per_protein.csv'), 'w') as fh:
    w = csv.writer(fh)
    w.writerow(['force_field', 'protein_id', 'N', 'rho_occ_deg'])
    for ff in FFS:
        v = []
        for p, (o, g, a) in adj[ff].items():
            if len(o) >= 30:
                r = spearmanr(o, g).correlation
                w.writerow([ff, p, len(o), r])
                v.append(r)
        print(f'A1 {ff}: n {len(v)} median {np.median(v):+.3f}')
with open(os.path.join(os.path.dirname(OUT), 'a3_enrichment.csv'), 'w') as fh:
    w = csv.writer(fh)
    w.writerow(['force_field', 'set', 'aa', 'enrichment', 'count'])
    for ff in ('CALVADOS-2', 'Mpipi-Recharged'):
        cnt = {k: defaultdict(int) for k in ('all', 'occ_only', 'deg_only')}
        for p, (o, g, a) in data[ff].items():          # raw labels, as A3
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
        for key in ('occ_only', 'deg_only'):
            for aa in sorted(cnt['all']):
                e = (cnt[key][aa] / tot[key]) / (cnt['all'][aa] / tot['all'])
                w.writerow([ff, key, aa, e, cnt[key][aa]])
            top = sorted(((aa, (cnt[key][aa] / tot[key]) / (cnt['all'][aa] / tot['all'])) for aa in cnt['all']), key=lambda t: -t[1])
            print(f'A3 {ff} {key}: ' + ' '.join(f'{a} {v:.2f}' for a, v in top[:5]))
# A2 residual, verbatim from hub_vs_degree_sim.py
resid = {}
for ff in ('CALVADOS-2', 'Mpipi-Recharged'):
    resid[ff] = {}
    for p, (o, g, a) in adj[ff].items():
        beta = float(np.dot(o, g) / np.dot(g, g)) if np.dot(g, g) > 0 else 0.0
        resid[ff][p] = o - beta * g
r_res = [(p, spearmanr(resid['CALVADOS-2'][p], resid['Mpipi-Recharged'][p]).correlation) for p in paired]
by_len = defaultdict(list)
for p in paired:
    by_len[len(adj['CALVADOS-2'][p][0])].append(p)
ctrl = []
for n, ps in by_len.items():
    if len(ps) < 2:
        continue
    for i, p in enumerate(ps):
        q = ps[(i + 1) % len(ps)]
        ctrl.append((p, spearmanr(resid['CALVADOS-2'][p], resid['Mpipi-Recharged'][q]).correlation))
with open(os.path.join(os.path.dirname(OUT), 'a2_residual.csv'), 'w') as fh:
    w = csv.writer(fh)
    w.writerow(['kind', 'protein_id', 'rho'])
    w.writerows([('same protein', p, r) for p, r in r_res] + [('control', p, r) for p, r in ctrl])
v = np.array([r for _, r in r_res if r == r]); c = np.array([r for _, r in ctrl if r == r])
print(f'A2 residual: same protein n {len(v)} median {np.median(v):+.3f}; control n {len(c)} median {np.median(c):+.3f}')
