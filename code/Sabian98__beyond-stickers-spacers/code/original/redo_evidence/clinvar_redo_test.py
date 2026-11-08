"""
clinvar_redo_test.py — step 2 of clinvar_redo_plan.md (hash 3ee66781). Run only after clinvar_redo/predictions.csv is
hashed. Classes are joined here with the code of clinvar_hub.analyse.

  ~/.venvs/bender/bin/python -I clinvar_redo_test.py  ->  clinvar_redo/summary.txt
"""
import csv
import gzip
import hashlib
import os
import sys
from collections import defaultdict

import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

HERE = '/Volumes/User Homes/Nu Project/BENDER-DYN/redo_evidence'
sys.path.insert(0, HERE)
import clinvar_hub as CH  # noqa: E402

OUT = os.path.join(HERE, 'clinvar_redo')
AAS = 'ACDEFGHIKLMNPQRSTVWY'
L = []


def p(s=''):
    L.append(s)
    print(s, flush=True)


# ---------------------------------------------------------------- classes, exactly as clinvar_hub.analyse
cls = {}
rd = csv.reader(gzip.open(os.path.join(CH.EXT, 'missense_variants.csv.gz'), 'rt'))
h = next(rd)
ix = {c: i for i, c in enumerate(h)}
for row in rd:
    s = row[ix['significance']]
    cls[row[0]] = 1 if 'athogenic' in s else 0 if 'enign' in s else None
rep_cls = {}
with gzip.open(os.path.join(CH.EXT, 'variant_summary_2026-10-05.txt.gz'), 'rt') as fh:
    r2 = csv.reader(fh, delimiter='\t')
    hh = next(r2)
    jx = {c.lstrip('#'): i for i, c in enumerate(hh)}
    for row in r2:
        if row[jx['Assembly']] == 'GRCh38':
            cs = row[jx['ClinicalSignificance']]
            if cs in CH.PATH_RC or cs in CH.BEN_RC:
                rep_cls[row[jx['VariationID']]] = 1 if cs in CH.PATH_RC else 0
R = []
for r in csv.DictReader(open(os.path.join(OUT, 'predictions.csv'))):
    y = cls.get(r['variant']) if r['set'] == 'primary' else rep_cls.get(r['variation_id'])
    if y is None:
        continue
    r['y'] = y
    R.append(r)
num = lambda r, k: float(r[k]) if r[k] != '' else np.nan


def subset(st, collagen=False, notch3=True):
    return [r for r in R if r['set'] == st and (collagen or not r['gene'].startswith('COL'))
            and (notch3 or r['gene'] != 'NOTCH3')]


def onehot(rows, key):
    M = np.zeros((len(rows), 20))
    for i, r in enumerate(rows):
        M[i, AAS.index(r[key])] = 1
    return M


def covars(rows):
    return np.column_stack([onehot(rows, 'wt'), onehot(rows, 'mut'), [num(r, 'cons') for r in rows],
                            [np.log1p(num(r, 'dist_end')) for r in rows], [np.log(num(r, 'N')) for r in rows]])


# ---------------------------------------------------------------- Primary A
def within_idr(rows, feat, nperm=10000, seed=0, sign=None):
    rows = [r for r in rows if r['cons'] != '']
    Z = covars(rows)
    A = np.column_stack([np.ones(len(rows)), Z])
    x = np.array([num(r, feat) for r in rows])
    res = x - A @ np.linalg.lstsq(A, x, rcond=None)[0]
    sd = res.std()
    g = defaultdict(list)
    for i, r in enumerate(rows):
        g[r['idr']].append(i)
    keep = [v for v in g.values() if 0 < sum(rows[i]['y'] for i in v) < len(v)]
    idx = np.concatenate([np.array(v) for v in keep])
    grp = np.concatenate([np.full(len(v), j) for j, v in enumerate(keep)])
    yy = np.array([rows[i]['y'] for i in idx])
    rr = res[idx]
    G = len(keep)
    n_g = np.bincount(grp, minlength=G)
    np_g = np.bincount(grp, weights=yy, minlength=G)
    nb_g = n_g - np_g
    S_g = np.bincount(grp, weights=rr, minlength=G)
    w = np_g * nb_g / n_g

    def stat(sp):
        d = sp / np_g - (S_g - sp) / nb_g
        return float((w * d).sum() / w.sum() / sd)
    t0 = stat(np.bincount(grp, weights=rr * yy, minlength=G))
    rng = np.random.RandomState(seed)
    order = np.lexsort((np.arange(len(grp)), grp))          # rows sorted by group
    starts = np.concatenate([[0], np.cumsum(n_g)[:-1]])
    rank_in_g = np.arange(len(grp)) - starts[grp[order]]
    perm = np.empty(nperm)
    for b in range(nperm):
        key = rng.random_sample(len(grp))
        o = np.lexsort((key, grp))                         # random order within each group
        pos = np.empty(len(grp), int)
        pos[o] = np.arange(len(grp)) - starts[grp[o]]
        lab = pos < np_g[grp]                              # the first n_p of each group are 'pathogenic'
        perm[b] = stat(np.bincount(grp, weights=rr * lab, minlength=G))
    p2 = (np.sum(np.abs(perm) >= abs(t0)) + 1) / (nperm + 1)
    p1 = None
    if sign is not None:
        p1 = (np.sum(perm * sign >= t0 * sign) + 1) / (nperm + 1)
    return t0, p2, p1, G, len(idx), int(yy.sum())


p('ClinVar second look (clinvar_redo_plan.md): predicted degree within IDRs, and gain over AlphaMissense')
p(f'variants with a class: primary {sum(r["set"] == "primary" for r in R)}, replication {sum(r["set"] == "replication" for r in R)}')
p('')
p('PRIMARY A  within-IDR effect of dsite (mutant minus wild-type predicted degree at the site), residualised on wild-type and')
p('           mutant amino acid, conservation, distance from the IDR end and IDR length; collagens excluded')
t0, p2, _, G, n, npat = within_idr(subset('primary'), 'dsite')
p(f'  primary:     standardised effect {t0:+.4f}  two-sided permutation P {p2:.4f}  ({G} IDRs, {n} variants, {npat} pathogenic)')
s = 1 if t0 > 0 else -1
t1, q2, q1, G1, n1, np1 = within_idr(subset('replication'), 'dsite', sign=s)
p(f'  replication: standardised effect {t1:+.4f}  one-sided P (same sign) {q1:.4f}  two-sided {q2:.4f}  ({G1} IDRs, {n1} variants, {np1} pathogenic)')
okA = p2 < 0.05 and np.sign(t1) == s and q1 < 0.05
p(f'  decision A: {"SUCCESS" if okA else "NOT met"}'
  + (f' (pathogenic substitutions {"RAISE" if s > 0 else "LOWER"} predicted degree at the site more than benign ones in the same IDR)' if okA else ''))
p('  secondary (no decision):')
for feat, lab in (('dglob', 'dglob, size of the perturbation over the IDR'), ('wsite', 'wsite, wild-type degree at the site')):
    a, b, _, G2, n2, _ = within_idr(subset('primary'), feat, nperm=5000)
    c, d, _, _, _, _ = within_idr(subset('replication'), feat, nperm=5000)
    p(f'    {lab:46s} primary {a:+.4f} (P {b:.4f});  replication {c:+.4f} (P {d:.4f})')
p('')


# ---------------------------------------------------------------- Primary B
def feats(rows, extra):
    am = np.clip(np.array([num(r, 'am') for r in rows]), 1e-4, 1 - 1e-4)
    base = np.column_stack([np.log(am / (1 - am)), covars(rows)])
    if not extra:
        return base
    return np.column_stack([base] + [[num(r, k) for r in rows] for k in extra])


def model():
    return make_pipeline(StandardScaler(), LogisticRegression(C=1.0, max_iter=5000))


def eligible(rows, need=('am', 'cons')):
    return [r for r in rows if all(r[k] != '' for k in need)]


def gene_boot(y, pa, pb, genes, nboot=1000, seed=1):
    idx = defaultdict(list)
    for i, g in enumerate(genes):
        idx[g].append(i)
    keys = list(idx)
    rng = np.random.RandomState(seed)
    ds = []
    for _ in range(nboot):
        pick = np.concatenate([idx[keys[j]] for j in rng.choice(len(keys), len(keys), replace=True)])
        if len(set(y[pick])) < 2:
            continue
        ds.append(roc_auc_score(y[pick], pa[pick]) - roc_auc_score(y[pick], pb[pick]))
    return np.percentile(ds, [2.5, 97.5])


def oof(rows, extra):
    X = feats(rows, extra)
    y = np.array([r['y'] for r in rows])
    fold = np.array([int(hashlib.md5(r['gene'].encode()).hexdigest(), 16) % 5 for r in rows])
    pr = np.zeros(len(rows))
    for k in range(5):
        tr, te = fold != k, fold == k
        pr[te] = model().fit(X[tr], y[tr]).predict_proba(X[te])[:, 1]
    return y, pr


def primary_b(rows_p, rows_r, extra, label):
    y, pb = oof(rows_p, [])
    _, pf = oof(rows_p, extra)
    g = [r['gene'] for r in rows_p]
    ab, af = roc_auc_score(y, pb), roc_auc_score(y, pf)
    ci = gene_boot(y, pf, pb, g)
    Xb, Xf = feats(rows_p, []), feats(rows_p, extra)
    mb, mf = model().fit(Xb, y), model().fit(Xf, y)
    yr = np.array([r['y'] for r in rows_r])
    rb, rf = mb.predict_proba(feats(rows_r, []))[:, 1], mf.predict_proba(feats(rows_r, extra))[:, 1]
    arb, arf = roc_auc_score(yr, rb), roc_auc_score(yr, rf)
    cir = gene_boot(yr, rf, rb, [r['gene'] for r in rows_r])
    p(f'  {label}')
    p(f'    primary (genes held out, n {len(rows_p)}): AUC base {ab:.4f}, full {af:.4f}, delta {af - ab:+.4f} '
      f'(95% gene bootstrap {ci[0]:+.4f} to {ci[1]:+.4f})')
    p(f'    replication (n {len(rows_r)}):             AUC base {arb:.4f}, full {arf:.4f}, delta {arf - arb:+.4f} '
      f'({cir[0]:+.4f} to {cir[1]:+.4f})')
    return af - ab, ci, arf - arb, y, pb, pf, rows_p


EP, ER = eligible(subset('primary')), eligible(subset('replication'))
p('PRIMARY B  gain over AlphaMissense + conservation + substitution identity + position (collagens excluded)')
p(f'  AlphaMissense matched: primary {sum(r["am"] != "" for r in subset("primary"))} of {len(subset("primary"))}, '
  f'replication {sum(r["am"] != "" for r in subset("replication"))} of {len(subset("replication"))}')
d, ci, dr, y, pb, pf, rows_p = primary_b(EP, ER, ['dsite', 'dglob', 'wsite'], 'base vs base + dsite + dglob + wsite')
okB = ci[0] > 0 and dr > 0
p(f'  decision B: {"SUCCESS" if okB else "NOT met"}'
  + (f'; {"meaningful (>= 0.01)" if d >= 0.01 else "detectable but negligible (< 0.01)"}' if okB else ''))
p('  secondary:')
both = defaultdict(set)
for r in rows_p:
    both[r['idr']].add(r['y'])
m = np.array([len(both[r['idr']]) == 2 for r in rows_p])
p(f'    B2 within IDRs carrying both classes ({int(m.sum())} variants): AUC base {roc_auc_score(y[m], pb[m]):.4f}, full {roc_auc_score(y[m], pf[m]):.4f}')
p(f'    B3 AlphaMissense alone: primary AUC {roc_auc_score(y, [num(r, "am") for r in rows_p]):.4f}; '
  f'replication {roc_auc_score([r["y"] for r in ER], [num(r, "am") for r in ER]):.4f}')
EPn = eligible(subset('primary'), need=('am', 'cons', 'delta_nu'))
yb4, pb4 = oof(EPn, [])
_, pn4 = oof(EPn, ['delta_nu'])
_, pd4 = oof(EPn, ['dsite', 'dglob', 'wsite'])
_, pa4 = oof(EPn, ['delta_nu', 'dsite', 'dglob', 'wsite'])
p(f'    B4 primary set with IDRome Delta_nu (n {len(EPn)}): AUC base {roc_auc_score(yb4, pb4):.4f}, + Delta_nu {roc_auc_score(yb4, pn4):.4f}, '
  f'+ degree {roc_auc_score(yb4, pd4):.4f}, + both {roc_auc_score(yb4, pa4):.4f}')
p('')
p('SENSITIVITY (no decision)')
t, q, _, _, _, _ = within_idr(subset('primary', notch3=False), 'dsite', nperm=5000)
p(f'  A without NOTCH3: effect {t:+.4f} (P {q:.4f})')
t, q, _, _, _, _ = within_idr(subset('primary', collagen=True), 'dsite', nperm=5000)
p(f'  A with collagens: effect {t:+.4f} (P {q:.4f})')
primary_b(eligible(subset('primary', notch3=False)), eligible(subset('replication', notch3=False)),
          ['dsite', 'dglob', 'wsite'], 'B without NOTCH3')
primary_b(eligible(subset('primary', collagen=True)), eligible(subset('replication', collagen=True)),
          ['dsite', 'dglob', 'wsite'], 'B with collagens')
open(os.path.join(OUT, 'summary.txt'), 'w').write('\n'.join(L) + '\n')
