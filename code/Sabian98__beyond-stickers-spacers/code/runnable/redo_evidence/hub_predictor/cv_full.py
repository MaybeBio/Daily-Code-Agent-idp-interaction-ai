"""
cv_full.py — 5-fold CV (folds fixed by protein-ID hash) on all CALVADOS-2 proteins in data_full/, writing per-protein
scores for Figure 3 instead of only medians.
Conditions (target occupancy unless noted):
  position     gradient boosting on position + length features only (the trivial baseline)
  ridge        ridge on the full 487 sequence features
  hgb          gradient boosting on the full features
  hgb_z        gradient boosting, target z_btw (position-adjusted z of time-averaged betweenness)
  transfer     the hgb fold models applied to the same held-out proteins in MPIPI-GG and Mpipi-Recharged
Also per protein: the label noise ceiling, sqrt of the Spearman-Brown split-half reliability of occupancy.
Writes cv_full_per_protein.csv (protein_id, force_field, N, condition, spearman, top10_precision, ceiling) and
cv_full_example.csv (per-residue truth and prediction for the held-out protein closest to the median score).
  ~/.venvs/bender/bin/python cv_full.py [--data data_full]
"""
import argparse
import csv
import os

import numpy as np
from scipy.stats import spearmanr
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.linear_model import Ridge
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

import train_hub_model as T

HERE = os.path.dirname(os.path.abspath(__file__))
POS_COLS = slice(483, 487)                   # log length (last global feature) + the 3 position features


def hgb():
    return HistGradientBoostingRegressor(max_iter=300, learning_rate=0.08, max_leaf_nodes=63, l2_regularization=1.0,
                                         random_state=0)


def scores(pred, true, n_top=None):
    rho = spearmanr(pred, true).correlation
    k = max(1, int(round(0.1 * len(true))))
    prec = len(set(np.argsort(-pred)[:k]) & set(np.argsort(-true)[:k])) / k
    return rho, prec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--data', default=os.path.join(HERE, 'data_full'))
    a = ap.parse_args()
    prots, y_occ = T.load(a.data, 'occupancy')
    _, y_z = T.load(a.data, 'z_btw')
    cal = [k for k in prots if k[0] == 'CALVADOS-2' and len(y_occ.get(k, {})) == len(prots[k]['sequence'])]
    feats = {}

    def F(seq):
        if seq not in feats:
            feats[seq] = T.features(seq, 10).astype(np.float32)
        return feats[seq]
    X = np.vstack([F(prots[k]['sequence']) for k in cal])
    assert X.shape[1] == 487, X.shape
    Yo = np.concatenate([[y_occ[k][i] for i in range(len(prots[k]['sequence']))] for k in cal]).astype(np.float32)
    Yz = np.concatenate([[y_z[k].get(i, np.nan) for i in range(len(prots[k]['sequence']))] for k in cal]).astype(np.float32)
    G = np.concatenate([np.full(len(prots[k]['sequence']), j) for j, k in enumerate(cal)])
    FO = np.concatenate([np.full(len(prots[k]['sequence']), int(prots[k]['fold'])) for k in cal])
    print(f'{len(cal)} CALVADOS-2 proteins, {len(Yo)} residues', flush=True)

    preds = {c: np.zeros(len(Yo), np.float32) for c in ('position', 'ridge', 'hgb', 'hgb_z')}
    fold_models = {}
    for f in range(5):
        tr, te = FO != f, FO == f
        m = hgb().fit(X[tr][:, POS_COLS], Yo[tr])
        preds['position'][te] = m.predict(X[te][:, POS_COLS])
        m = make_pipeline(StandardScaler(), Ridge(alpha=1000.0)).fit(X[tr], Yo[tr])
        preds['ridge'][te] = m.predict(X[te])
        m = hgb().fit(X[tr], Yo[tr])
        fold_models[f] = m
        preds['hgb'][te] = m.predict(X[te])
        ok = tr & ~np.isnan(Yz)
        m = hgb().fit(X[ok], Yz[ok])
        preds['hgb_z'][te] = m.predict(X[te])
        print(f'fold {f} done', flush=True)

    rows = []
    for j, k in enumerate(cal):
        sel = G == j
        r = prots[k]
        sb = float(r['splithalf_sb_occ']) if r.get('splithalf_sb_occ', '') not in ('', 'nan') else np.nan   # occupancy, not betweenness
        ceil = float(np.sqrt(max(sb, 0))) if sb == sb else np.nan
        for c in ('position', 'ridge', 'hgb'):
            rho, prec = scores(preds[c][sel], Yo[sel])
            rows.append([k[1], 'CALVADOS-2', sel.sum(), c, rho, prec, ceil])
        zt = Yz[sel]
        if not np.isnan(zt).any():
            rho, prec = scores(preds['hgb_z'][sel], zt)
            rows.append([k[1], 'CALVADOS-2', sel.sum(), 'hgb_z', rho, prec, ceil])
    # transfer to the other force fields, same held-out folds
    for ff in ('MPIPI-GG', 'Mpipi-Recharged'):
        for k in [k for k in prots if k[0] == ff]:
            seq = prots[k]['sequence']
            lab = y_occ.get(k, {})
            if len(lab) != len(seq):
                continue
            m = fold_models[int(prots[k]['fold'])]
            rho, prec = scores(m.predict(F(seq)), np.array([lab[i] for i in range(len(seq))]))
            rows.append([k[1], ff, len(seq), 'transfer', rho, prec, np.nan])
    with open(os.path.join(HERE, 'cv_full_per_protein.csv'), 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['protein_id', 'force_field', 'N', 'condition', 'spearman', 'top10_precision', 'ceiling'])
        w.writerows(rows)
    import collections
    by = collections.defaultdict(list)
    for r in rows:
        by[(r[1], r[3])].append(r[4])
    for key, v in sorted(by.items()):
        print(f'{key[0]:16s} {key[1]:9s} n={len(v):5d} median Spearman {np.nanmedian(v):.3f}', flush=True)

    # example: the held-out protein closest to the median hgb score, length 80-200
    hg = [(r[0], r[4]) for r in rows if r[3] == 'hgb' and r[1] == 'CALVADOS-2' and 80 <= r[2] <= 200]
    med = np.nanmedian([v for _, v in hg])
    ex = min(hg, key=lambda t: abs(t[1] - med))[0]
    j = cal.index(('CALVADOS-2', ex))
    sel = G == j
    with open(os.path.join(HERE, 'cv_full_example.csv'), 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['protein_id', 'pos', 'aa', 'occupancy', 'pred_hgb', 'pred_position'])
        seq = prots[cal[j]]['sequence']
        for i, (t, p1, p2) in enumerate(zip(Yo[sel], preds['hgb'][sel], preds['position'][sel])):
            w.writerow([ex, i, seq[i], f'{t:.4f}', f'{p1:.4f}', f'{p2:.4f}'])
    print('example protein:', ex)


if __name__ == '__main__':
    main()
