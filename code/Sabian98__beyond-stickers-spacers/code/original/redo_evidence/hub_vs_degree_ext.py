"""
hub_vs_degree_ext.py — Parts D and E of hub_vs_degree_plan.md (hash 9212f51e). Needs hub_vs_degree/degree_model_cal.pkl,
hashed into gaivi/FROZEN.sha256 before this is run.
  D  activation-domain tiles with an absolute quantity (second look at the same tiles)
  E  responses to substitutions in the CALVADOS-2 design labels (not blind)

  OMP_NUM_THREADS=8 ~/.venvs/bender/bin/python -I hub_vs_degree_ext.py  ->  hub_vs_degree/part_de.txt
"""
import csv
import os
import pickle
import sys
from collections import defaultdict

import numpy as np
from scipy import stats

HERE = '/Volumes/User Homes/Nu Project/BENDER-DYN/redo_evidence'
sys.path.insert(0, os.path.join(HERE, 'hub_predictor'))
sys.path.insert(0, HERE)
import train_hub_model as T  # noqa: E402
import clinvar_hub as CH  # noqa: E402
import sticker_spacer_test as SS  # noqa: E402

OUT = os.path.join(HERE, 'hub_vs_degree')
TILES = '/Volumes/User Homes/Nu Project/BENDER-DYN/external_data/ad_sanborn_tiles/sanborn_raw_no_reps.csv'
deg_model = pickle.load(open(os.path.join(OUT, 'degree_model_cal.pkl'), 'rb'))
occ_model = pickle.load(open(os.path.join(HERE, 'hub_predictor/models_v3_2026-10-04.pkl'), 'rb'))[('CALVADOS-2', 'full')]
L = []

# ------------------------------------------------------------------ D
score = {}
with open(TILES, newline='') as fh:
    rd = csv.reader(fh)
    next(rd)
    for row in rd:
        try:
            score[row[0]] = float(row[-1])
        except Exception:
            pass
P = [r for r in csv.DictReader(open(os.path.join(HERE, 'finches_hub/tiles_predictions.csv'))) if r['name'] in score]
f = lambda r, k: float(r[k]) if r[k] != '' else np.nan

# positional baseline of degree: trees on length and position only, CALVADOS-2 degree labels
seqs, ys = defaultdict(list), defaultdict(list)
for r in csv.DictReader(open(os.path.join(HERE, 'hub_predictor/data_v3/residues.csv'))):
    if r['force_field'] == 'CALVADOS-2' and r['deg_mean'] != '':
        seqs[r['protein_id']].append(r['aa'])
        ys[r['protein_id']].append(float(r['deg_mean']))
Xb = np.vstack([CH.pos_features(len(v)) for v in seqs.values()])
yb = np.concatenate([ys[k] for k in seqs])
from sklearn.ensemble import HistGradientBoostingRegressor  # noqa: E402
posm = HistGradientBoostingRegressor(max_iter=300, max_leaf_nodes=63, learning_rate=0.08, l2_regularization=1.0,
                                     random_state=0).fit(Xb, yb)

D, Dexc = [], []
for r in P:
    s = r['seq']
    p = deg_model.predict(T.features(s, 10))
    b = posm.predict(CH.pos_features(len(s)))
    idx = [i for i, a in enumerate(s) if a in 'FWYL']
    D.append(float(p[idx].mean()) if len(idx) >= 2 else np.nan)
    Dexc.append(float((p - b)[idx].mean()) if len(idx) >= 2 else np.nan)
D, Dexc = np.array(D), np.array(Dexc)
sc = np.array([score[r['name']] for r in P])
att, rep, slf = (np.array([f(r, f'{c}_cal') for r in P]) for c in ('attr', 'rep', 'self'))
comp = np.column_stack([np.array([f(r, c) for r in P]) for c in ['DE', 'FWY', 'L', 'KR', 'P', 'GS', 'net']])
clus = np.array([r['name'].rsplit(':', 1)[0] for r in P])
inS = (att <= np.percentile(att, 100 / 3)) & (rep <= np.median(rep))


def resid(y, Z):
    A = np.column_stack([np.ones(len(y)), Z])
    return y - A @ np.linalg.lstsq(A, y, rcond=None)[0]


def prho(x, idx):
    Z = np.column_stack([slf[idx], comp[idx]])
    Zr = np.column_stack([stats.rankdata(Z[:, j]) for j in range(Z.shape[1])])
    return float(np.corrcoef(resid(stats.rankdata(sc[idx]), Zr), resid(stats.rankdata(x[idx]), Zr))[0, 1]), len(idx)


def test(x, mask, label):
    idx = np.where(mask & ~np.isnan(x))[0]
    r, n = prho(x, idx)
    k = comp.shape[1] + 1
    df = n - k - 2
    t = r * np.sqrt(df / (1 - r * r))
    p1, p2 = stats.t.cdf(t, df), 2 * min(stats.t.cdf(t, df), stats.t.sf(t, df))
    L.append(f'{label:62s} n {n:5d}  partial rho {r:+.3f}  P one-sided(-) {p1:.4f}  P two-sided {p2:.4f}')
    return idx, r


L.append('D  activation-domain tiles with an absolute quantity, second look (hub_vs_degree_plan.md, Part D)')
L.append('   compare with the occupancy result of finches_hub/part1_results.txt: primary partial rho -0.075, one-sided P 0.0002')
L.append('')
idx, r0 = test(Dexc * 0 + D, inS, 'PRIMARY: D (mean predicted degree of F/W/Y/L, raw) in S')
test(Dexc, inS, 'secondary: D position-excess (minus the degree positional baseline) in S')
test(D, np.ones(len(P), bool), 'secondary: D over all tiles')
test(Dexc, np.ones(len(P), bool), 'secondary: D position-excess over all tiles')
rng = np.random.RandomState(1)
cl = {c: idx[clus[idx] == c] for c in set(clus[idx])}
keys = list(cl)
boots = []
for _ in range(1000):
    pick = rng.choice(len(keys), len(keys), replace=True)
    i2 = np.concatenate([cl[keys[k]] for k in pick])
    Z = np.column_stack([slf[i2], comp[i2]])
    Zr = np.column_stack([stats.rankdata(Z[:, j]) for j in range(Z.shape[1])])
    boots.append(np.corrcoef(resid(stats.rankdata(sc[i2]), Zr), resid(stats.rankdata(D[i2]), Zr))[0, 1])
boots = np.array(boots)
L.append(f'   cluster bootstrap over {len(keys)} source proteins in S: 95% interval {np.percentile(boots, 2.5):+.3f} to '
         f'{np.percentile(boots, 97.5):+.3f}; resamples negative {np.mean(boots < 0) * 100:.1f}%')
ok = ~np.isnan(D)
L.append(f'   Spearman(D, epsilon_hom) = {stats.spearmanr(D[ok], slf[ok]).correlation:+.3f};  '
         f'Spearman(D, score) in S = {stats.spearmanr(D[ok & inS], sc[ok & inS]).correlation:+.3f}')
L.append('')

# ------------------------------------------------------------------ E
table = {r['id']: r for r in csv.DictReader(open(os.path.join(HERE, 'gaivi/hub_design_table.csv')))}


def load_labels(path, key):
    o = defaultdict(lambda: ([], []))
    for r in csv.DictReader(open(path)):
        d = o[r['UniProt_ID']]
        d[0].append(float(r['occupancy']))
        d[1].append(float(r['deg_mean']))
    return {k: (np.array(a), np.array(b)) for k, (a, b) in o.items()}


lab = load_labels(os.path.join(HERE, 'labels_design_calvados/per_residue.csv'), 'UniProt_ID')
res = {'occupancy': [], 'degree': []}
used = []
for uid, r in table.items():
    if r['category'] != 'variant' or r['parent'] not in table or uid not in lab or r['parent'] not in lab:
        continue
    sv, sp = r['sequence'], table[r['parent']]['sequence']
    if len(sv) != len(sp) or sv == sp or len(lab[uid][0]) != len(sv) or len(lab[r['parent']][0]) != len(sp):
        continue
    used.append(uid)
    Fv, Fp = T.features(sv, 10), T.features(sp, 10)
    for k, m, j in (('occupancy', occ_model, 0), ('degree', deg_model, 1)):
        dp = m.predict(Fv) - m.predict(Fp)
        ds = lab[uid][j] - lab[r['parent']][j]
        res[k].append(stats.spearmanr(dp, ds).correlation if np.std(dp) > 0 and np.std(ds) > 0 else np.nan)
# reference: between-run Spearman of the level on the replica pairs (design replicas against the BENDER run)
bender = defaultdict(lambda: ([], []))
for r in csv.DictReader(open(os.path.join(HERE, 'hub_frames_v3_cal/per_residue.csv'))):
    d = bender[r['UniProt_ID']]
    d[0].append(float(r['occupancy']))
    d[1].append(float(r['deg_mean']))
ref = {'occupancy': [], 'degree': []}
for uid, r in table.items():
    if r['category'] in ('replica', 'gh_replica') and uid in lab and r['parent'] in bender:
        a, b = lab[uid], bender[r['parent']]
        if len(a[0]) == len(b[0]) >= 20:
            ref['occupancy'].append(stats.spearmanr(a[0], b[0]).correlation)
            ref['degree'].append(stats.spearmanr(a[1], b[1]).correlation)
L.append('E  responses to substitutions, CALVADOS-2 design labels (hub_vs_degree_plan.md, Part E; not blind)')
L.append(f'   variants with the same length as a labelled parent and at least one substitution: {len(used)}')
for k in ('occupancy', 'degree'):
    v = np.array([x for x in res[k] if x == x])
    c = np.median(ref[k]) if ref[k] else np.nan
    L.append(f'   {k:9s} median per-variant Spearman of predicted vs simulated change {np.median(v):+.3f}, '
             f'above zero in {int((v > 0).sum())} of {len(v)};  replica-run level agreement {c:.3f} '
             f'({len(ref[k])} pairs), ratio {np.median(v) / c:.2f}')
L.append('   the ceiling is the between-run agreement of the LEVEL (change-level replicates do not exist), as a reference only')
open(os.path.join(OUT, 'part_de.txt'), 'w').write('\n'.join(L) + '\n')
print('\n'.join(L))
