"""
train_hub_model.py — predict per-residue network hubness of disordered proteins from sequence alone.
Needs numpy + scikit-learn (GAIVI, or ~/.venvs/bender on the lab Mac). Data: proteins.csv / residues.csv /
experimental.csv from hub_dataset.py (same folder, or --data DIR).

Target (--target): z_btw (position-adjusted z of 1,000-frame time-averaged betweenness; default), occupancy
(fraction of frames as a top-10% hub) or excess_fixed (simulation minus residue-neutral null; chemistry only).
Features (sequence only): one-hot window ±W, local composition in windows ±5 and ±15, global composition, FCR, NCPR,
SCD, log length, relative and log absolute distance to the nearer chain end.
Models: ridge (alpha chosen on the training folds) and HistGradientBoosting.
Evaluation, 5-fold cross-validation by protein (folds fixed per protein across force fields):
  per-protein Spearman between predicted and simulated residue values (median over test proteins), and its noise
  ceiling (median split-half Spearman-Brown reliability of the labels, as a correlation ceiling sqrt(rel));
  top-10% hub precision (fraction of predicted top 10% in the true top 10%; chance 0.10);
  cross-force-field transfer: train on one force field, test on another on the held-out proteins;
  external test: models trained on all CALVADOS-2 proteins predict the experimental PED ensembles (per-protein Spearman
  vs experimental centrality, alongside the residue-neutral null's Spearman).
  python3 train_hub_model.py [--ff CALVADOS-2] [--target z_btw] [--window 10] [--max-proteins 0] [--data DIR]
"""
import argparse
import csv
import math
import os
import warnings
from collections import defaultdict

import numpy as np
from scipy.stats import spearmanr
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.linear_model import Ridge
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

# numpy 2.0 + Apple Accelerate raises spurious divide/overflow flags in matmul; checked: inputs and predictions are
# finite and match an einsum computation to 1e-15. Harmless elsewhere (e.g. on GAIVI).
warnings.filterwarnings('ignore', message='.*encountered in matmul', category=RuntimeWarning)

AAS = 'ACDEFGHIKLMNPQRSTVWY'
AI = {a: i for i, a in enumerate(AAS)}
PHOS = {'x': 'S', 'y': 'T', 'z': 'Y'}           # phosphoresidues in experimental.csv


def norm_aa(a):
    return PHOS.get(a, a)


def scd(seq):
    q = np.array([1 if a in 'KR' else -1 if a in 'DE' else 0 for a in seq], float)
    idx = np.nonzero(q)[0]
    if len(idx) < 2:
        return 0.0
    qi = q[idx]
    d = np.sqrt(np.abs(idx[:, None] - idx[None, :]))
    return float(np.triu(np.outer(qi, qi) * d, 1).sum() / len(seq))


def features(seq, window):
    """(N, F) feature matrix for one sequence."""
    seq = ''.join(norm_aa(a) for a in seq)
    n = len(seq)
    oh = np.zeros((n, 20))
    for i, a in enumerate(seq):
        if a in AI:
            oh[i, AI[a]] = 1.0
    pad = np.zeros((window, 20))
    ohp = np.vstack([pad, oh, pad])
    win = np.hstack([ohp[window + k: window + k + n] for k in range(-window, window + 1)])     # (n, 20 * (2W + 1))
    cs = np.vstack([np.zeros((1, 20)), np.cumsum(oh, axis=0)])

    def local(w):
        lo = np.clip(np.arange(n) - w, 0, n)
        hi = np.clip(np.arange(n) + w + 1, 0, n)
        cnt = (cs[hi] - cs[lo]) - oh
        return cnt / np.maximum(1, (hi - lo - 1))[:, None]
    comp = oh.mean(axis=0)
    fcr = comp[[AI['D'], AI['E'], AI['K'], AI['R']]].sum()
    ncpr = comp[[AI['K'], AI['R']]].sum() - comp[[AI['D'], AI['E']]].sum()
    glob = np.concatenate([comp, [fcr, ncpr, scd(seq), math.log(n)]])
    pos = np.arange(n)
    end = np.minimum(pos, n - 1 - pos)
    posf = np.stack([end / n, np.log1p(end), pos / max(1, n - 1)], axis=1)
    return np.hstack([win, local(5), local(15), np.tile(glob, (n, 1)), posf])


def load(data_dir, target, ff_filter=None, max_proteins=0):
    prots = {}
    for r in csv.DictReader(open(os.path.join(data_dir, 'proteins.csv'))):
        if ff_filter and r['force_field'] not in ff_filter:
            continue
        prots[(r['force_field'], r['protein_id'])] = r
    if max_proteins:
        keep = set(list(prots)[:max_proteins])
        prots = {k: v for k, v in prots.items() if k in keep}
    y = defaultdict(dict)
    for r in csv.DictReader(open(os.path.join(data_dir, 'residues.csv'))):
        k = (r['force_field'], r['protein_id'])
        if k in prots and r[target] != '':
            y[k][int(r['pos'])] = float(r[target])
    return prots, y


def per_protein_scores(pred, true, groups):
    rho, prec = [], []
    for g in np.unique(groups):
        m = groups == g
        if m.sum() < 10 or np.std(true[m]) == 0:
            continue
        rho.append(spearmanr(pred[m], true[m]).correlation)
        k = max(1, int(round(0.1 * m.sum())))
        tp = set(np.argsort(-pred[m])[:k])
        tt = set(np.argsort(-true[m])[:k])
        prec.append(len(tp & tt) / k)
    return float(np.nanmedian(rho)), float(np.mean(prec)), len(rho)


def build(prots, y, window, keys):
    X, Y, G, F = [], [], [], []
    for gi, k in enumerate(keys):
        seq = prots[k]['sequence']
        lab = y.get(k, {})
        if len(lab) != len(seq):
            continue
        X.append(features(seq, window))
        Y.append(np.array([lab[i] for i in range(len(seq))]))
        G.append(np.full(len(seq), gi))
        F.append(np.full(len(seq), int(prots[k]['fold'])))
    return np.vstack(X), np.concatenate(Y), np.concatenate(G), np.concatenate(F)


def fit(model, X, Y):
    if model == 'ridge':
        best = None
        n = len(Y)
        idx = np.random.default_rng(0).permutation(n)
        tr, va = idx[: int(0.8 * n)], idx[int(0.8 * n):]
        for a in (10.0, 100.0, 1000.0, 10000.0):
            m = make_pipeline(StandardScaler(), Ridge(alpha=a)).fit(X[tr], Y[tr])
            r = np.corrcoef(m.predict(X[va]), Y[va])[0, 1]
            if best is None or r > best[0]:
                best = (r, a)
        return make_pipeline(StandardScaler(), Ridge(alpha=best[1])).fit(X, Y)
    return HistGradientBoostingRegressor(max_iter=300, learning_rate=0.08, max_leaf_nodes=63, l2_regularization=1.0,
                                         random_state=0).fit(X, Y)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--data', default=os.path.dirname(os.path.abspath(__file__)))
    ap.add_argument('--ff', default='CALVADOS-2')
    ap.add_argument('--target', default='z_btw')
    ap.add_argument('--window', type=int, default=10)
    ap.add_argument('--max-proteins', type=int, default=0)
    ap.add_argument('--models', default='ridge,hgb')
    a = ap.parse_args()
    lines = []

    def p(s=''):
        lines.append(s)
        print(s, flush=True)

    prots, y = load(a.data, a.target, max_proteins=a.max_proteins)
    ffs = sorted({k[0] for k in prots})
    p(f'Sequence -> per-residue hub predictor; target {a.target}; window ±{a.window}; force fields in data: {ffs}')
    rel = [float(v['splithalf_sb']) for k, v in prots.items() if k[0] == a.ff and v['splithalf_sb'] not in ('', 'nan')]
    if rel:
        rel = np.array(rel)
        p(f'  label reliability ({a.ff}): median split-half Spearman-Brown {np.median(rel):.2f}; '
          f'correlation ceiling ~ {np.median(np.sqrt(np.clip(rel, 0, 1))):.2f}')
    keys = [k for k in prots if k[0] == a.ff]
    X, Y, G, F = build(prots, y, a.window, keys)
    p(f'  {a.ff}: {len(np.unique(G))} proteins, {len(Y)} residues, {X.shape[1]} features')
    fold_models = {}
    for model in a.models.split(','):
        pred = np.zeros_like(Y)
        for f in range(5):
            tr, te = F != f, F == f
            m = fit(model, X[tr], Y[tr])
            fold_models[(model, f)] = m
            pred[te] = m.predict(X[te])
        rho, prec, n = per_protein_scores(pred, Y, G)
        r2 = 1 - np.sum((Y - pred) ** 2) / np.sum((Y - Y.mean()) ** 2)
        p(f'  [{model}] 5-fold CV by protein: per-protein Spearman median {rho:.2f} ({n} proteins); top-10% precision '
          f'{prec:.2f} (chance 0.10); pooled R^2 {r2:.2f}')
    # cross-force-field transfer on the held-out folds
    for other in ffs:
        if other == a.ff:
            continue
        okeys = [k for k in prots if k[0] == other]
        if not okeys:
            continue
        Xo, Yo, Go, Fo = build(prots, y, a.window, okeys)
        for model in a.models.split(','):
            pred = np.zeros_like(Yo)
            for f in range(5):
                m = fold_models[(model, f)]          # same held-out proteins as in the cross-validation
                te = Fo == f
                pred[te] = m.predict(Xo[te])
            rho, prec, n = per_protein_scores(pred, Yo, Go)
            p(f'  [{model}] trained on {a.ff}, tested on {other} (held-out proteins): per-protein Spearman {rho:.2f} '
              f'({n} proteins); top-10% precision {prec:.2f}')
    # external test on experimental ensembles
    ef = os.path.join(a.data, 'experimental.csv')
    if os.path.exists(ef) and a.target in ('z_btw', 'excess_fixed'):
        ex = defaultdict(dict)
        names = {}
        for r in csv.DictReader(open(ef)):
            ex[r['entry']][int(r['pos'])] = (r['aa'], float(r['z_exp']), float(r['z_null']))
            names[r['entry']] = (r['protein'], r['restraints'])
        for model in a.models.split(','):
            m = fit(model, X, Y)
            p(f'\n  [{model}] external test, trained on all {a.ff} proteins: per-protein Spearman with the experimental '
              f'ensemble (model / residue-neutral null)')
            rows = []
            for e, d in sorted(ex.items()):
                seq = ''.join(d[i][0] for i in range(len(d)))
                ze = np.array([d[i][1] for i in range(len(d))])
                zn = np.array([d[i][2] for i in range(len(d))])
                target = ze - zn if a.target == 'excess_fixed' else ze
                pr = m.predict(features(seq, a.window))
                r_m = spearmanr(pr, target).correlation
                r_n = spearmanr(zn, ze).correlation
                rows.append((e, names[e][0], names[e][1], len(seq), r_m, r_n))
                p(f'    {e} {names[e][0][:24]:24s} {names[e][1][:22]:22s} N {len(seq):4d}: model {r_m:+.2f}   null {r_n:+.2f}')
            rm = [r[4] for r in rows if not math.isnan(r[4])]
            p(f'    median over {len(rm)} ensembles: model {np.median(rm):+.2f}; model > 0 in {sum(v > 0 for v in rm)}/{len(rm)}')
    of = os.path.join(a.data, 'experimental_occupancy.csv')
    if os.path.exists(of) and a.target == 'occupancy':
        ex = defaultdict(dict)
        names = {}
        for r in csv.DictReader(open(of)):
            ex[r['entry']][int(r['pos'])] = (r['aa'], float(r['occ_exp']), float(r['occ_null']))
            names[r['entry']] = (r['protein'], r['restraints'])
        for model in a.models.split(','):
            m = fit(model, X, Y)
            p(f'\n  [{model}] external test (occupancy), trained on all {a.ff} proteins: per-protein Spearman with the '
              f'experimental ensemble (model / residue-neutral null / position only)')
            rows = []
            for e, d in sorted(ex.items()):
                seq = ''.join(d[i][0] for i in range(len(d)))
                oe = np.array([d[i][1] for i in range(len(d))])
                on = np.array([d[i][2] for i in range(len(d))])
                n = len(seq)
                end = np.minimum(np.arange(n), n - 1 - np.arange(n))
                r_m = spearmanr(m.predict(features(seq, a.window)), oe).correlation
                r_n = spearmanr(on, oe).correlation
                r_p = spearmanr(end, oe).correlation
                rows.append((r_m, r_n, r_p))
                p(f'    {e} {names[e][0][:24]:24s} {names[e][1][:22]:22s} N {n:4d}: model {r_m:+.2f}   null {r_n:+.2f}   position {r_p:+.2f}')
            R = np.array(rows)
            p(f'    medians over {len(R)} ensembles: model {np.nanmedian(R[:, 0]):+.2f}; null {np.nanmedian(R[:, 1]):+.2f}; '
              f'position only {np.nanmedian(R[:, 2]):+.2f}; model > null in {int((R[:, 0] > R[:, 1]).sum())}/{len(R)}, '
              f'model > position in {int((R[:, 0] > R[:, 2]).sum())}/{len(R)}')
    out = os.path.join(a.data, f'results_{a.ff}_{a.target}_{a.models.replace(",", "+")}.txt')
    open(out, 'w').write('\n'.join(lines) + '\n')
    print(f'-> {out}')


if __name__ == '__main__':
    main()
