"""
phospho_charge_test.py — sequence charge model and Part C of phospho_charge_function_plan.md (hash 28befc80).
Written and hashed before it was run. Refuses to run unless its inputs match gaivi/FROZEN.sha256.

1. Fits the frozen sequence charge model: OLS of log ESS on 18 sequence features (K/R and D/E counts in six separation
   bins, log distances to the termini, log N, FCR, NCPR, Thr indicator) over the Ser/Thr of the 63 phosphosite-test
   proteins not among the 20 phospho_direct proteins. Also fits the literal competitor (K/R and D/E counts within
   |i-j| <= 5, plus intercept) on the same sites. Coefficients -> phospho_charge/seq_model.json (used unchanged by Part F).
2. Part C on the 60 direct sites: does log ESS predict E_occ / R_occ(W, W) beyond the sequence-predicted log ESS,
   within protein (Freedman-Lane style permutation of the residual of log ESS on the prediction, 20,000 permutations).
   Secondary: E_deg, the reverse question, the literal competitor in place of the sequence model.

  ~/.venvs/finches/bin/python phospho_charge_test.py   ->   phospho_charge/{seq_model.json, summary.txt, direct_sites.csv}
"""
import json
import math
import os
import random
import sys

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

DYN = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, DYN)
from phospho_direct_test import check_frozen  # noqa: E402

SN = os.path.join(DYN, 'phospho_network', 'sites_network.csv')
PROF = os.path.join(DYN, 'phospho_network', 'profiles')
SEL = os.path.join(DYN, 'gaivi', 'phospho_direct', 'selection_amend1.csv')
SD = os.path.join(DYN, 'phospho_direct', 'sites_direct.csv')
WS = os.path.join(DYN, 'phospho_direct', 'wt_vs_stored.csv')
OUT = os.path.join(DYN, 'phospho_charge')
SEED, B = 20261013, 20000
BINS = [(1, 2), (3, 5), (6, 10), (11, 20), (21, 40), (41, 10 ** 9)]
NAMES = ([f'KR_{a}_{b}' for a, b in BINS] + [f'DE_{a}_{b}' for a, b in BINS]
         + ['log_dist_N', 'log_dist_C', 'log_N', 'FCR', 'NCPR', 'is_T'])


def features(seq, pos):
    """The 18 features of the plan for the residue at 1-based position pos."""
    i, n = pos - 1, len(seq)
    kr, de = [0] * len(BINS), [0] * len(BINS)
    for j, c in enumerate(seq):
        d = abs(j - i)
        if d == 0 or c not in 'KRDE':
            continue
        k = next(b for b, (lo, hi) in enumerate(BINS) if lo <= d <= hi)
        if c in 'KR':
            kr[k] += 1
        else:
            de[k] += 1
    npos, nneg = sum(seq.count(c) for c in 'KR'), sum(seq.count(c) for c in 'DE')
    return kr + de + [math.log(pos), math.log(n - pos + 1), math.log(n), (npos + nneg) / n, (npos - nneg) / n,
                      1.0 if seq[i] == 'T' else 0.0]


def local_counts(seq, pos):
    i = pos - 1
    w = [seq[j] for j in range(max(0, i - 5), min(len(seq), i + 6)) if j != i]
    return [sum(c in 'KR' for c in w), sum(c in 'DE' for c in w)]


def ols(X, y):
    A = np.column_stack([np.ones(len(X)), X])
    return np.linalg.lstsq(A, y, rcond=None)[0]


def predict(coef, X):
    return coef[0] + np.asarray(X) @ coef[1:]


def seq_of(uid):
    return str(np.load(os.path.join(PROF, f'{uid}.npz'))['seq'])


def centre(v, g):
    v = pd.Series(np.asarray(v, float))
    return (v - v.groupby(g).transform('mean')).to_numpy()


def perm_coef(y, x1, x2, groups, rng):
    """b1 of y ~ x1 + x2 (centred, no intercept) and its permutation null (residual of x1 on x2, permuted in group)."""
    def fit(a, b):
        M = np.array([[a @ a, a @ b], [a @ b, b @ b]])
        return np.linalg.solve(M, np.array([a @ y, b @ y]))[0]
    obs = fit(x1, x2)
    r = x1 - (x1 @ x2) / (x2 @ x2) * x2
    g = pd.factorize(groups)[0]
    base = np.lexsort((np.arange(len(r)), g))                  # positions grouped by protein, original order
    null = np.empty(B)
    for k in range(B):
        order = np.lexsort((rng.random(len(r)), g))            # the same groups, random order inside each
        rp = np.empty_like(r)
        rp[base] = r[order]                                    # each protein's positions refilled by its own values
        null[k] = fit(rp, x2)
    return obs, (1 + int((null <= obs).sum())) / (B + 1)


def main():
    for f in (SN, SEL, SD, WS):
        check_frozen(f)
    os.makedirs(OUT, exist_ok=True)
    sn = pd.read_csv(SN)
    sel = pd.read_csv(SEL)
    direct = set(sel.UniProt_ID)
    seqs = {u: seq_of(u) for u in sn.UniProt_ID.unique()}
    sn['log_ess'] = np.log(sn.ess_frac_net)
    X = np.array([features(seqs[u], p) for u, p in zip(sn.UniProt_ID, sn.pos)])
    Xl = np.array([local_counts(seqs[u], p) for u, p in zip(sn.UniProt_ID, sn.pos)])
    tr = ~sn.UniProt_ID.isin(direct).to_numpy()
    coef = ols(X[tr], sn.log_ess.to_numpy()[tr])
    coef_l = ols(Xl[tr], sn.log_ess.to_numpy()[tr])
    json.dump({'features': NAMES, 'bins': BINS, 'intercept': coef[0], 'coef': dict(zip(NAMES, coef[1:])),
               'literal_features': ['KR_within_5', 'DE_within_5'], 'literal_intercept': coef_l[0],
               'literal_coef': list(coef_l[1:]), 'n_train_sites': int(tr.sum()),
               'n_train_proteins': int(sn[tr].UniProt_ID.nunique()), 'target': 'log(ess_frac_net), CALVADOS-2'},
              open(os.path.join(OUT, 'seq_model.json'), 'w'), indent=1)

    # protein-grouped 5-fold CV on the training proteins
    prots = sorted(sn[tr].UniProt_ID.unique())
    random.Random(SEED).shuffle(prots)
    fold = {u: k % 5 for k, u in enumerate(prots)}
    T = sn[tr].reset_index(drop=True)
    Xt, Xlt, yt = X[tr], Xl[tr], T.log_ess.to_numpy()
    cv, cvl = np.empty(len(T)), np.empty(len(T))
    fk = T.UniProt_ID.map(fold).to_numpy()
    for k in range(5):
        m = fk == k
        cv[m] = predict(ols(Xt[~m], yt[~m]), Xt[m])
        cvl[m] = predict(ols(Xlt[~m], yt[~m]), Xlt[m])
    r2 = lambda p, y: 1 - ((y - p) ** 2).sum() / ((y - y.mean()) ** 2).sum()  # noqa: E731
    wcv = lambda p, y, g: spearmanr(centre(p, g), centre(y, g)).correlation  # noqa: E731

    # Part C
    sd = pd.read_csv(SD)
    ws = pd.read_csv(WS).set_index('UniProt_ID')
    sd['seqpred'] = [predict(coef, features(seqs[u], p)) for u, p in zip(sd.UniProt_ID, sd.pos)]
    sd['litpred'] = [predict(coef_l, local_counts(seqs[u], p)) for u, p in zip(sd.UniProt_ID, sd.pos)]
    sd['y_occ'] = sd.E_occ / sd.UniProt_ID.map(ws.R_occ_WW)
    sd['y_deg'] = sd.E_deg / sd.UniProt_ID.map(ws.R_deg_WW)
    sd.to_csv(os.path.join(OUT, 'direct_sites.csv'), index=False)
    g = sd.UniProt_ID.to_numpy()
    c = {k: centre(sd[k], g) for k in ('log_ess', 'seqpred', 'litpred', 'y_occ', 'y_deg')}
    rng = np.random.default_rng(SEED)
    res = {}
    for name, y, a, b in (('C   E_occ: log ESS beyond sequence model', c['y_occ'], c['log_ess'], c['seqpred']),
                          ('    E_deg: log ESS beyond sequence model', c['y_deg'], c['log_ess'], c['seqpred']),
                          ('    E_occ: sequence model beyond log ESS', c['y_occ'], c['seqpred'], c['log_ess']),
                          ('    E_deg: sequence model beyond log ESS', c['y_deg'], c['seqpred'], c['log_ess']),
                          ('    E_occ: log ESS beyond literal counts', c['y_occ'], c['log_ess'], c['litpred']),
                          ('    E_deg: log ESS beyond literal counts', c['y_deg'], c['log_ess'], c['litpred'])):
        res[name] = perm_coef(y, a, b, g, rng)

    L = ['phospho_charge Part C (plan 28befc80)', '',
         f'Sequence charge model: OLS, 18 features, trained on {int(tr.sum())} Ser/Thr of '
         f'{sn[tr].UniProt_ID.nunique()} proteins (the 20 direct proteins excluded)',
         f'  protein-grouped 5-fold CV R^2 for log ESS: sequence model {r2(cv, yt):.3f}, literal counts (|i-j| <= 5) '
         f'{r2(cvl, yt):.3f}',
         f'  within-protein Spearman (CV predictions): sequence model {wcv(cv, yt, T.UniProt_ID):+.3f}, literal '
         f'{wcv(cvl, yt, T.UniProt_ID):+.3f}',
         f'  60 direct sites (out of sample): R^2 {r2(sd.seqpred.to_numpy(), sd.log_ess.to_numpy()):.3f}, Spearman '
         f'{spearmanr(sd.seqpred, sd.log_ess).correlation:+.3f}, within-protein Spearman '
         f'{spearmanr(c["seqpred"], c["log_ess"]).correlation:+.3f}', '',
         'Direct rewiring, within protein (y = E / R(W, W), centred), Spearman with:',
         f'  -log ESS {spearmanr(-c["log_ess"], c["y_occ"]).correlation:+.3f} (occ), '
         f'{spearmanr(-c["log_ess"], c["y_deg"]).correlation:+.3f} (deg); '
         f'-sequence prediction {spearmanr(-c["seqpred"], c["y_occ"]).correlation:+.3f} (occ), '
         f'{spearmanr(-c["seqpred"], c["y_deg"]).correlation:+.3f} (deg); '
         f'-literal {spearmanr(-c["litpred"], c["y_occ"]).correlation:+.3f} (occ)', '',
         f'Permutation tests, coefficient of the first predictor given the second, one-sided negative, {B} permutations']
    for name, (b1, p) in res.items():
        L.append(f'  {name}: b {b1:+.4f}, P {p:.4g}')
    pc = res['C   E_occ: log ESS beyond sequence model'][1]
    L += ['', f'C DECISION: {"POSITIVE, the ensemble predicts direct rewiring beyond the sequence charge model" if pc < 0.05 else "NOT POSITIVE, sequence charge is not shown to be insufficient"} (P {pc:.4g})']
    txt = '\n'.join(L) + '\n'
    open(os.path.join(OUT, 'summary.txt'), 'w').write(txt)
    print(txt)


if __name__ == '__main__':
    main()
