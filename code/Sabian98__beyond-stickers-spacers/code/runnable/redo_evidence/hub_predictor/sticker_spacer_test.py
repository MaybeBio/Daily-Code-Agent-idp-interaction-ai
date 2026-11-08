"""
sticker_spacer_test.py — is hub occupancy just a restatement of the sticker-spacer picture?

On the v3 CALVADOS-2 labels (10,281 proteins, same folds and booster settings as cv_full.py) this asks:
  A  overlap: which residues are hubs (top 10% occupancy of their protein)? Enrichment of stickers (aromatic F/W/Y, and
     aromatic plus R) and of the classic spacers G/S among hubs; share of stickers that are hubs; within-protein Spearman
     between occupancy and the residue's own stickiness (CALVADOS-2 lambda).
  B  variance: share of the within-protein variance of position-adjusted occupancy explained by residue identity alone
     (between-type variance / total), i.e. how much of hubness a per-type sticker/spacer label could ever capture.
  C  prediction (5-fold CV by protein, folds of data_v3): gradient boosting on
       position          position + length (as cv_full.py)
       identity          position + one-hot of the residue itself
       sticker_spacer    position + the residue's own stickiness and charge + sticker density (mean lambda, aromatic
                         fraction), spacer content (G/S fraction) and net charge in windows of +-2, +-5, +-10 residues,
                         + whole-sequence mean lambda, aromatic fraction and net charge per residue
     compared with the full sequence model (cv_full_per_protein.csv, condition hgb).
  ~/.venvs/bender/bin/python sticker_spacer_test.py  ->  sticker_spacer_summary.txt, sticker_spacer_per_protein.csv
"""
import csv
import os

import numpy as np
from scipy.stats import spearmanr
from sklearn.ensemble import HistGradientBoostingRegressor

import train_hub_model as T

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, 'data_v3')
CAL_TABLE = os.path.expanduser((__import__('os').environ['BSS_DATA'] + '/Scripts/CALVADOS IDP/input/residues_CALVADOS2.csv'))
AAS = 'ACDEFGHIKLMNPQRSTVWY'
POS_COLS = slice(483, 487)                   # as in cv_full.py: log length + the 3 position features
AROM, STICK, SPACER = set('FWY'), set('FWYR'), set('GS')
out = []


def p(s=''):
    out.append(s)
    print(s, flush=True)


def hgb():
    return HistGradientBoostingRegressor(max_iter=300, learning_rate=0.08, max_leaf_nodes=63, l2_regularization=1.0,
                                         random_state=0)


def window_mean(x, w):
    """Mean of x over i-w..i+w (truncated at the chain ends)."""
    c = np.concatenate([[0.0], np.cumsum(x)])
    n = len(x)
    lo = np.clip(np.arange(n) - w, 0, n)
    hi = np.clip(np.arange(n) + w + 1, 0, n)
    return (c[hi] - c[lo]) / (hi - lo)


def sticker_features(seq, lam, chg):
    l = np.array([lam[a] for a in seq])
    q = np.array([chg[a] for a in seq], float)
    aro = np.array([a in AROM for a in seq], float)
    gs = np.array([a in SPACER for a in seq], float)
    cols = [l, q]
    for w in (2, 5, 10):
        cols += [window_mean(l, w), window_mean(aro, w), window_mean(gs, w), window_mean(q, w)]
    n = len(seq)
    cols += [np.full(n, l.mean()), np.full(n, aro.mean()), np.full(n, q.mean())]
    return np.vstack(cols).T


def main():
    tab = list(csv.DictReader(open(CAL_TABLE)))
    lam = {r['one']: float(r['lambdas']) for r in tab}
    chg = {r['one']: float(r['q']) for r in tab}
    chg['H'] = 0.0
    prots, y = T.load(DATA, 'occupancy')
    cal = [k for k in prots if k[0] == 'CALVADOS-2' and len(y.get(k, {})) == len(prots[k]['sequence'])]
    seqs = [prots[k]['sequence'] for k in cal]
    Y = [np.array([y[k][i] for i in range(len(s))]) for k, s in zip(cal, seqs)]
    p(f'{len(cal)} CALVADOS-2 proteins, {sum(map(len, seqs))} residues')

    # A  overlap of hubs with stickers and spacers
    rng = np.random.default_rng(0)
    tot = {'res': 0, 'hub': 0}
    cnt = {name: {'res': 0, 'hub': 0} for name in ('aromatic FWY', 'stickers FWYR', 'spacers GS')}
    sets = {'aromatic FWY': AROM, 'stickers FWYR': STICK, 'spacers GS': SPACER}
    rho_lam = []
    for s, occ in zip(seqs, Y):
        n = len(s)
        k = max(1, int(round(0.1 * n)))
        hubs = set(np.argsort(-(occ + rng.uniform(0, 1e-9, n)))[:k])
        tot['res'] += n
        tot['hub'] += k
        for name, st in sets.items():
            idx = [i for i, a in enumerate(s) if a in st]
            cnt[name]['res'] += len(idx)
            cnt[name]['hub'] += sum(i in hubs for i in idx)
        lv = np.array([lam[a] for a in s])
        if lv.std() > 0 and occ.std() > 0:
            rho_lam.append(spearmanr(lv, occ).correlation)
    p('\nA  hubs (top 10% occupancy) versus stickers and spacers')
    for name, c in cnt.items():
        share_hub = c['hub'] / tot['hub']
        share_seq = c['res'] / tot['res']
        p(f'  {name:14s}: {100 * share_hub:5.1f}% of hubs vs {100 * share_seq:5.1f}% of residues '
          f'(enrichment {share_hub / share_seq:.2f}); {100 * c["hub"] / c["res"]:.1f}% of them are hubs (10% expected)')
    p(f'  within-protein Spearman between occupancy and own stickiness (lambda): median {np.median(rho_lam):+.2f}')

    # B  variance of position-adjusted occupancy explained by residue identity
    resid, aa_all = [], []
    for s, occ in zip(seqs, Y):
        n = len(s)
        x = np.arange(n) / max(n - 1, 1)
        A = np.vstack([np.ones(n), x, x ** 2]).T
        r = occ - A @ np.linalg.lstsq(A, occ, rcond=None)[0]
        if r.std() > 0:
            resid.append(r / r.std())
            aa_all.append(np.array(list(s)))
    r = np.concatenate(resid)
    a = np.concatenate(aa_all)
    means = {t: r[a == t].mean() for t in AAS if (a == t).any()}
    between = sum((a == t).sum() * m ** 2 for t, m in means.items())
    p(f'\nB  residue identity explains {100 * between / (r ** 2).sum():.1f}% of the within-protein variance of '
      f'position-adjusted occupancy')

    # C  cross-validated prediction
    feats = {}

    def F(seq):
        if seq not in feats:
            feats[seq] = T.features(seq, 10).astype(np.float32)[:, POS_COLS]
        return feats[seq]
    P = np.vstack([F(s) for s in seqs])
    I = np.vstack([np.eye(20, dtype=np.float32)[[AAS.index(c) if c in AAS else 0 for c in s]] for s in seqs])
    S = np.vstack([sticker_features(s, lam, chg) for s in seqs]).astype(np.float32)
    Yc = np.concatenate(Y).astype(np.float32)
    G = np.concatenate([np.full(len(s), j) for j, s in enumerate(seqs)])
    FO = np.concatenate([np.full(len(s), int(prots[k]['fold'])) for k, s in zip(cal, seqs)])
    X = {'position': P, 'identity': np.hstack([P, I]), 'sticker_spacer': np.hstack([P, S])}
    pred = {c: np.zeros(len(Yc), np.float32) for c in X}
    for f in range(5):
        tr, te = FO != f, FO == f
        for c, M in X.items():
            pred[c][te] = hgb().fit(M[tr], Yc[tr]).predict(M[te])
        print(f'fold {f} done', flush=True)
    full = {r['protein_id']: float(r['spearman']) for r in csv.DictReader(open(os.path.join(HERE, 'cv_full_per_protein.csv')))
            if r['force_field'] == 'CALVADOS-2' and r['condition'] == 'hgb' and r['spearman'] not in ('', 'nan')}
    rows, med = [], {c: [] for c in X}
    for j, k in enumerate(cal):
        sel = G == j
        row = [k[1], int(sel.sum())]
        for c in X:
            v = spearmanr(pred[c][sel], Yc[sel]).correlation
            med[c].append(v)
            row.append(f'{v:.4f}')
        row.append(f'{full.get(k[1], np.nan):.4f}')
        rows.append(row)
    with open(os.path.join(HERE, 'sticker_spacer_per_protein.csv'), 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['protein_id', 'N'] + list(X) + ['full_model'])
        w.writerows(rows)
    p('\nC  held-out per-protein Spearman (median over proteins)')
    for c in X:
        p(f'  {c:15s} {np.nanmedian(med[c]):.3f}')
    p(f'  {"full sequence":15s} {np.nanmedian(list(full.values())):.3f}  (cv_full.py, 487 features)')
    open(os.path.join(HERE, 'sticker_spacer_summary.txt'), 'w').write('\n'.join(out) + '\n')


if __name__ == '__main__':
    main()
