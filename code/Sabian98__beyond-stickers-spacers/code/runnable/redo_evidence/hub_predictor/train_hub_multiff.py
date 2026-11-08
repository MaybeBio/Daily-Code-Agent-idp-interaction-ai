"""
train_hub_multiff.py — a force-field-agnostic hub-propensity predictor.

Two ways of removing the force-field signature (same features as train_hub_model.py; no force-field label is ever
given to the model):
  pooled     train on all force fields together, residues weighted so that each force field contributes equally
  consensus  for proteins simulated in >= 2 force fields, train on the mean per-residue label across force fields
Evaluation:
  leave-one-force-field-out (pooled): train on the other force fields' proteins outside fold f, test on the held-out
      force field's proteins in fold f. This is the direct measure of force-field agnosticism, to compare with the
      CALVADOS-2-only transfer (0.42 / 0.46).
  within-force-field CV for each force field (pooled model), and CV for the consensus model against the consensus
      label and against each force field's own labels
  external test on the PED experimental ensembles of disordered proteins (target z_btw only)
  python3 train_hub_multiff.py [--target z_btw] [--window 10] [--model hgb|ridge] [--data DIR]
"""
import argparse
import csv
import os
from collections import defaultdict

import numpy as np
from scipy.stats import spearmanr
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.linear_model import Ridge
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

import train_hub_model as T                     # features(), load(), per_protein_scores(); silences Accelerate warnings


def fit_w(model, X, Y, w):
    if model == 'ridge':
        m = make_pipeline(StandardScaler(), Ridge(alpha=1000.0))
        return m.fit(X, Y, ridge__sample_weight=w)
    m = HistGradientBoostingRegressor(max_iter=300, learning_rate=0.08, max_leaf_nodes=63, l2_regularization=1.0,
                                      random_state=0)
    return m.fit(X, Y, sample_weight=w)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--data', default=os.path.dirname(os.path.abspath(__file__)))
    ap.add_argument('--target', default='z_btw')
    ap.add_argument('--window', type=int, default=10)
    ap.add_argument('--model', default='hgb')
    a = ap.parse_args()
    lines = []

    def p(s=''):
        lines.append(s)
        print(s, flush=True)

    prots, y = T.load(a.data, a.target)
    ffs = sorted({k[0] for k in prots})
    seq_of, fold_of = {}, {}
    for (ff, uid), r in prots.items():
        if uid in seq_of and seq_of[uid] != r['sequence']:
            continue                     # same ID, different sequence across force fields: keep the first
        seq_of[uid] = r['sequence']
        fold_of[uid] = int(r['fold'])
    feats = {}

    def F(uid):
        if uid not in feats:
            feats[uid] = T.features(seq_of[uid], a.window)
        return feats[uid]

    # pooled rows
    Xs, Ys, FFs, Gs, Fo = [], [], [], [], []
    for gi, ((ff, uid), r) in enumerate(prots.items()):
        lab = y.get((ff, uid), {})
        if uid not in seq_of or r['sequence'] != seq_of[uid] or len(lab) != len(seq_of[uid]):
            continue
        Xs.append(F(uid))
        Ys.append(np.array([lab[i] for i in range(len(lab))]))
        FFs.append(np.full(len(lab), ffs.index(ff)))
        Gs.append(np.full(len(lab), gi))
        Fo.append(np.full(len(lab), fold_of[uid]))
    X, Y, FF, G, FO = np.vstack(Xs), np.concatenate(Ys), np.concatenate(FFs), np.concatenate(Gs), np.concatenate(Fo)
    counts = {i: int((FF == i).sum()) for i in range(len(ffs))}
    w = np.array([1.0 / counts[i] for i in FF]) * len(Y) / len(ffs)
    p(f'Force-field-agnostic hub predictor ({a.model}); target {a.target}; window ±{a.window}')
    p('  data: ' + ', '.join(f'{ffs[i]} {len(np.unique(G[FF == i]))} proteins / {counts[i]} residues' for i in counts))

    p('\n  Pooled model, leave-one-force-field-out (train on the other force fields, held-out proteins of the unseen one):')
    for h in range(len(ffs)):
        pred = np.full(len(Y), np.nan)
        for f in range(5):
            tr = (FF != h) & (FO != f)
            te = (FF == h) & (FO == f)
            if te.sum() == 0 or tr.sum() == 0:
                continue
            m = fit_w(a.model, X[tr], Y[tr], w[tr])
            pred[te] = m.predict(X[te])
        sel = FF == h
        rho, prec, n = T.per_protein_scores(pred[sel], Y[sel], G[sel])
        p(f'    unseen {ffs[h]:16s}: per-protein Spearman {rho:.2f} ({n} proteins); top-10% precision {prec:.2f}')

    p('\n  Pooled model, all force fields in training, 5-fold CV by protein, scored per force field:')
    pred = np.full(len(Y), np.nan)
    for f in range(5):
        m = fit_w(a.model, X[FO != f], Y[FO != f], w[FO != f])
        pred[FO == f] = m.predict(X[FO == f])
    for h in range(len(ffs)):
        sel = FF == h
        rho, prec, n = T.per_protein_scores(pred[sel], Y[sel], G[sel])
        p(f'    {ffs[h]:16s}: per-protein Spearman {rho:.2f} ({n} proteins); top-10% precision {prec:.2f}')

    # consensus rows: proteins with labels in >= 2 force fields
    by_uid = defaultdict(dict)
    for (ff, uid), lab in y.items():
        if uid in seq_of and len(lab) == len(seq_of[uid]) and prots[(ff, uid)]['sequence'] == seq_of[uid]:
            by_uid[uid][ff] = np.array([lab[i] for i in range(len(lab))])
    multi = [u for u, d in by_uid.items() if len(d) >= 2]
    p(f'\n  Consensus model: {len(multi)} proteins with labels in >= 2 force fields '
      f'({sum(len(by_uid[u]) == 3 for u in multi)} in all three); target = mean label across force fields')
    if multi:
        Xc = np.vstack([F(u) for u in multi])
        Yc = np.concatenate([np.mean(list(by_uid[u].values()), axis=0) for u in multi])
        Gc = np.concatenate([np.full(len(seq_of[u]), k) for k, u in enumerate(multi)])
        Fc = np.concatenate([np.full(len(seq_of[u]), fold_of[u]) for u in multi])
        per_ff = {ff: np.concatenate([by_uid[u].get(ff, np.full(len(seq_of[u]), np.nan)) for u in multi]) for ff in ffs}
        pc = np.zeros(len(Yc))
        for f in range(5):
            m = fit_w(a.model, Xc[Fc != f], Yc[Fc != f], np.ones(int((Fc != f).sum())))
            pc[Fc == f] = m.predict(Xc[Fc == f])
        rho, prec, n = T.per_protein_scores(pc, Yc, Gc)
        p(f'    vs consensus label: per-protein Spearman {rho:.2f} ({n} proteins); top-10% precision {prec:.2f}')
        for ff in ffs:
            ok = ~np.isnan(per_ff[ff])
            if ok.sum():
                rho, prec, n = T.per_protein_scores(pc[ok], per_ff[ff][ok], Gc[ok])
                p(f'    vs {ff:16s} labels: per-protein Spearman {rho:.2f} ({n} proteins); top-10% precision {prec:.2f}')

    ef = os.path.join(a.data, 'experimental.csv')
    if os.path.exists(ef) and a.target == 'z_btw':
        ex = defaultdict(dict)
        names = {}
        for r in csv.DictReader(open(ef)):
            ex[r['entry']][int(r['pos'])] = (r['aa'], float(r['z_exp']), float(r['z_null']))
            names[r['entry']] = r['protein']
        models = {'pooled': fit_w(a.model, X, Y, w)}
        if multi:
            models['consensus'] = fit_w(a.model, Xc, Yc, np.ones(len(Yc)))
        for name, m in models.items():
            rs = []
            for e, d in sorted(ex.items()):
                seq = ''.join(d[i][0] for i in range(len(d)))
                ze = np.array([d[i][1] for i in range(len(d))])
                rs.append(spearmanr(m.predict(T.features(seq, a.window)), ze).correlation)
            rs = np.array(rs)
            p(f'\n  External PED test, {name} model: median per-protein Spearman {np.nanmedian(rs):+.2f}; '
              f'> 0 in {int((rs > 0).sum())}/{len(rs)}   (CALVADOS-2-only gradient boosting: +0.27)')
    of = os.path.join(a.data, 'experimental_occupancy.csv')
    if os.path.exists(of) and a.target == 'occupancy':
        ex = defaultdict(dict)
        for r in csv.DictReader(open(of)):
            ex[r['entry']][int(r['pos'])] = (r['aa'], float(r['occ_exp']), float(r['occ_null']))
        models = {'pooled': fit_w(a.model, X, Y, w)}
        if multi:
            models['consensus'] = fit_w(a.model, Xc, Yc, np.ones(len(Yc)))
        base_n, base_p = [], []
        for e, d in sorted(ex.items()):
            n = len(d)
            oe = np.array([d[i][1] for i in range(n)])
            base_n.append(spearmanr(np.array([d[i][2] for i in range(n)]), oe).correlation)
            base_p.append(spearmanr(np.minimum(np.arange(n), n - 1 - np.arange(n)), oe).correlation)
        for name, m in models.items():
            rs = []
            for e, d in sorted(ex.items()):
                seq = ''.join(d[i][0] for i in range(len(d)))
                oe = np.array([d[i][1] for i in range(len(d))])
                rs.append(spearmanr(m.predict(T.features(seq, a.window)), oe).correlation)
            rs = np.array(rs)
            p(f'\n  External PED test (occupancy), {name} model: median per-protein Spearman {np.nanmedian(rs):+.2f}; '
              f'beats the null in {int((rs > np.array(base_n)).sum())}/{len(rs)} and position-only in '
              f'{int((rs > np.array(base_p)).sum())}/{len(rs)}   (CALVADOS-2-only occupancy model: +0.39; null +0.30; position +0.31)')
    out = os.path.join(a.data, f'results_multiff_{a.target}_{a.model}.txt')
    open(out, 'w').write('\n'.join(lines) + '\n')
    print(f'-> {out}')


if __name__ == '__main__':
    main()
