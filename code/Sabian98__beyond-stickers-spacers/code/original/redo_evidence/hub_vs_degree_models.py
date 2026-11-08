"""
hub_vs_degree_models.py — Part C of hub_vs_degree_plan.md (hash 9212f51e).
Degree model trained exactly like the occupancy models: CALVADOS-2 labels of data_v3, features train_hub_model.features
(seq, 10), booster sticker_spacer_test.hgb(), hash folds of proteins.csv.
  1  matched comparison: occupancy and degree models fitted on folds 0-3, scored on fold 4
  2  final degree model fitted on all CALVADOS-2 proteins, pickled (hash into gaivi/FROZEN.sha256 afterwards)

  OMP_NUM_THREADS=8 ~/.venvs/bender/bin/python hub_vs_degree_models.py
      ->  hub_vs_degree/part_c.txt, hub_vs_degree/degree_model_cal.pkl
"""
import os
import pickle
import sys
import time
from collections import defaultdict

import numpy as np
from scipy.stats import spearmanr

HERE = '/Volumes/User Homes/Nu Project/BENDER-DYN/redo_evidence'
sys.path.insert(0, os.path.join(HERE, 'hub_predictor'))
import sticker_spacer_test as SS  # noqa: E402
import train_hub_model as T  # noqa: E402

OUT = os.path.join(HERE, 'hub_vs_degree')
os.makedirs(OUT, exist_ok=True)
t0 = time.time()
prots, y_occ = T.load(os.path.join(HERE, 'hub_predictor', 'data_v3'), 'occupancy')
_, y_deg = T.load(os.path.join(HERE, 'hub_predictor', 'data_v3'), 'deg_mean')
keys = [k for k in prots if k[0] == 'CALVADOS-2' and len(y_occ.get(k, {})) == len(prots[k]['sequence'])
        and len(y_deg.get(k, {})) == len(prots[k]['sequence'])]
seqs = [prots[k]['sequence'] for k in keys]
fold = np.array([int(prots[k]['fold']) for k in keys])
X = np.vstack([T.features(s, 10).astype(np.float32) for s in seqs])
Yo = np.concatenate([[y_occ[k][i] for i in range(len(s))] for k, s in zip(keys, seqs)]).astype(np.float32)
Yd = np.concatenate([[y_deg[k][i] for i in range(len(s))] for k, s in zip(keys, seqs)]).astype(np.float32)
G = np.concatenate([np.full(len(s), j) for j, s in enumerate(seqs)])
F = fold[G]
print(f'{len(keys)} proteins, {len(Yo)} residues, features {X.shape}, {time.time() - t0:.0f} s', flush=True)

tr, te = F != 4, F == 4
mo = SS.hgb().fit(X[tr], Yo[tr])
print(f'occupancy model fitted, {time.time() - t0:.0f} s', flush=True)
md = SS.hgb().fit(X[tr], Yd[tr])
print(f'degree model fitted, {time.time() - t0:.0f} s', flush=True)
po, pd_ = mo.predict(X[te]), md.predict(X[te])


def end_bin(n):
    return np.minimum(np.minimum(np.arange(n), n - 1 - np.arange(n)), 30)


def per_protein(pred, true, g):
    """Raw and position-adjusted per-protein Spearman, adjusting each vector by its own pooled end-distance profile."""
    ids = np.unique(g)
    accp, acct = defaultdict(list), defaultdict(list)
    for j in ids:
        m = g == j
        b = end_bin(m.sum())
        for arr, acc in ((pred[m], accp), (true[m], acct)):
            x = arr - arr.mean()
            for bi, v in zip(b, x):
                acc[bi].append(v)
    pp = {b: np.mean(v) for b, v in accp.items()}
    pt = {b: np.mean(v) for b, v in acct.items()}
    raw, adj = [], []
    for j in ids:
        m = g == j
        if m.sum() < 30:
            continue
        b = end_bin(m.sum())
        raw.append(spearmanr(pred[m], true[m]).correlation)
        a = pred[m] - pred[m].mean() - np.array([pp[i] for i in b])
        t = true[m] - true[m].mean() - np.array([pt[i] for i in b])
        adj.append(spearmanr(a, t).correlation)
    return np.array(raw), np.array(adj)


g_te = G[te]
ro, ao = per_protein(po, Yo[te], g_te)
rd, ad = per_protein(pd_, Yd[te], g_te)
r2 = 1 - ((Yd[te] - pd_) ** 2).sum() / ((Yd[te] - Yd[te].mean()) ** 2).sum()
# between-protein level: per-protein mean predicted vs true degree
ids = np.unique(g_te)
mp = np.array([pd_[g_te == j].mean() for j in ids])
mt = np.array([Yd[te][g_te == j].mean() for j in ids])
L = ['degree model against the occupancy model, matched fits on folds 0-3, scored on fold 4 (hub_vs_degree_plan.md, Part C)',
     f'held-out proteins {len(ids)} (scored with >= 30 residues: {len(ro)}), CALVADOS-2', '',
     f'  occupancy model   per-protein Spearman raw {np.median(ro):.3f}   position-adjusted {np.median(ao):.3f}',
     f'  degree model      per-protein Spearman raw {np.median(rd):.3f}   position-adjusted {np.median(ad):.3f}',
     f'  degree model      pooled R2 over residues {r2:.3f};  Spearman of per-protein mean degree, predicted vs true '
     f'{spearmanr(mp, mt).correlation:.3f}']
open(os.path.join(OUT, 'part_c.txt'), 'w').write('\n'.join(L) + '\n')
print('\n'.join(L), flush=True)

md_all = SS.hgb().fit(X, Yd)
pickle.dump(md_all, open(os.path.join(OUT, 'degree_model_cal.pkl'), 'wb'))
print(f'final degree model fitted on all {len(keys)} CALVADOS-2 proteins and pickled, {time.time() - t0:.0f} s', flush=True)
