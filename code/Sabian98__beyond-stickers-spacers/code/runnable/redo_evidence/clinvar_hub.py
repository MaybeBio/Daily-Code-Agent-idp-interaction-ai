"""
clinvar_hub.py — do pathogenic missense variants in human IDRs fall on predicted hubs? (clinvar_hub_plan.md, hashed
2026-10-09 before any variant was examined).

Stages (run in order; each stage hashes its outputs into gaivi/FROZEN.sha256):
  predict   frozen consensus model on the wild-type IDRs and on every mutant; positional baseline (position and length
            features only, trained on the CALVADOS-2 labels exactly as cv_full.py does) on the wild-type IDRs.
            Uses variant positions only, never the clinical classes.
  analyse   joins the classes and runs the pre-registered tests.
  ~/.venvs/bender/bin/python clinvar_hub.py predict|analyse
"""
import csv
import datetime
import gzip
import hashlib
import math
import os
import pickle
import re
import sys
from collections import defaultdict
from multiprocessing import Pool

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, 'hub_predictor'))
import train_hub_model as T  # noqa: E402

EXT = os.path.join(os.path.dirname(HERE), 'external_data', 'idrome_clinvar')
OUT = os.path.join(HERE, 'clinvar_hub')
FROZEN = os.path.join(HERE, 'gaivi', 'FROZEN.sha256')
MODELS = os.path.join(HERE, 'hub_predictor', 'frozen_models_2026-10-02.pkl')
AA3 = {'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C', 'Gln': 'Q', 'Glu': 'E', 'Gly': 'G', 'His': 'H',
       'Ile': 'I', 'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P', 'Ser': 'S', 'Thr': 'T', 'Trp': 'W',
       'Tyr': 'Y', 'Val': 'V'}
AAS = 'ACDEFGHIKLMNPQRSTVWY'
PATH_RC = {'Pathogenic', 'Likely pathogenic', 'Pathogenic/Likely pathogenic'}
BEN_RC = {'Benign', 'Likely benign', 'Benign/Likely benign'}


def register(path, note):
    h = hashlib.sha256(open(path, 'rb').read()).hexdigest()
    with open(FROZEN, 'a') as fh:
        fh.write(f'{h}  {os.path.relpath(path, os.path.dirname(FROZEN))}   # '
                 f'{datetime.datetime.now().isoformat(timespec="seconds")}; {note}\n')


def idrome():
    rows = {}
    for r in csv.DictReader(open(os.path.join(EXT, 'IDRome_DB.csv'))):
        rows[r['seq_name']] = dict(uniprot=r['UniProt_ID'], gene=r['gene_name'], seq=r['fasta'],
                                   first=int(float(r['first'])), last=int(float(r['last'])))
    return rows


def primary_variants():
    """(variant_id, idr, pos_in_idr, wt, mut) from IDRome's mapped ClinVar set; the class column is not read here."""
    out = []
    r = csv.reader(gzip.open(os.path.join(EXT, 'missense_variants.csv.gz'), 'rt'))
    h = next(r)
    ix = {c: i for i, c in enumerate(h)}
    for row in r:
        vid, idr = row[0], row[ix['WT']]
        m = re.match(r'(.+)_([A-Z])(\d+)([A-Z])$', vid)
        assert m and m.group(1) == idr, vid
        wt, pos, mut = m.group(2), int(m.group(3)), m.group(4)
        first = int(float(row[ix['first']]))
        k = pos - first
        seq = row[ix['fasta_WT']]
        out.append((vid, idr, k, wt, mut, seq[k] == wt))
    return out


def replication_variants(D, primary_keys):
    """Newer ClinVar missense variants (P/LP or B/LB, the two pooled) mapped to IDRome IDRs by gene and position."""
    by_gene = defaultdict(list)
    for name, d in D.items():
        by_gene[d['gene']].append(name)
    seen, out, stats = set(), [], defaultdict(int)
    with gzip.open(os.path.join(EXT, 'variant_summary_2026-10-05.txt.gz'), 'rt') as fh:
        rd = csv.reader(fh, delimiter='\t')
        h = next(rd)
        ix = {c.lstrip('#'): i for i, c in enumerate(h)}
        for row in rd:
            if row[ix['Assembly']] != 'GRCh38' or row[ix['Type']] != 'single nucleotide variant':
                continue
            cs = row[ix['ClinicalSignificance']]
            if cs not in PATH_RC and cs not in BEN_RC:
                continue
            m = re.search(r'\(p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})\)', row[ix['Name']])
            if not m or m.group(1) not in AA3 or m.group(3) not in AA3 or m.group(1) == m.group(3):
                continue
            stats['missense P/LP or B/LB'] += 1
            wt, pos, mut, gene = AA3[m.group(1)], int(m.group(2)), AA3[m.group(3)], row[ix['GeneSymbol']]
            hit = [n for n in by_gene.get(gene, []) if D[n]['first'] <= pos <= D[n]['last']]
            if not hit:
                continue
            stats['in an IDRome IDR'] += 1
            idr = hit[0]
            k = pos - D[idr]['first']
            if D[idr]['seq'][k] != wt:
                stats['reference residue mismatch'] += 1
                continue
            key = (idr, pos, mut)
            if key in seen:
                continue
            seen.add(key)
            if key in primary_keys:
                stats['already in primary set'] += 1
                continue
            out.append((f'{idr}_{wt}{pos}{mut}', idr, k, wt, mut, True, row[ix['VariationID']]))
    stats['kept'] = len(out)
    return out, dict(stats)


def positional_baseline():
    """HGB on log length + 3 position features, CALVADOS-2 occupancy labels (as cv_full.py 'position', all proteins)."""
    from sklearn.ensemble import HistGradientBoostingRegressor
    seqs, ys = defaultdict(list), defaultdict(list)
    for r in csv.DictReader(open(os.path.join(HERE, 'hub_frames_v3_cal', 'per_residue.csv'))):
        seqs[r['UniProt_ID']].append(r['aa'])
        ys[r['UniProt_ID']].append(float(r['occupancy']))
    X, Y = [], []
    for u in seqs:
        n = len(seqs[u])
        X.append(pos_features(n))
        Y.append(ys[u])
    m = HistGradientBoostingRegressor(max_iter=300, max_leaf_nodes=63, learning_rate=0.08, l2_regularization=1.0,
                                      random_state=0)
    m.fit(np.vstack(X), np.concatenate(Y))
    return m


def pos_features(n):
    pos = np.arange(n)
    end = np.minimum(pos, n - 1 - pos)
    return np.stack([np.full(n, math.log(n)), end / n, np.log1p(end), pos / max(1, n - 1)], axis=1)


_M = None


def _init():
    global _M
    _M = pickle.load(open(MODELS, 'rb'))['consensus']


def _wt(job):
    name, seq = job
    return name, _M.predict(T.features(seq, 10))


def _mut(job):
    vid, seq, k, mut, wtp = job
    s = seq[:k] + mut + seq[k + 1:]
    p = _M.predict(T.features(s, 10))
    return vid, float(np.abs(p - wtp).sum()), float(p[k] - wtp[k])


def predict():
    os.makedirs(OUT, exist_ok=True)
    h = hashlib.sha256(open(MODELS, 'rb').read()).hexdigest()
    assert h in open(FROZEN).read(), 'frozen models not registered'
    D = idrome()
    prim = primary_variants()
    print(f'primary: {len(prim)} variants, reference mismatch {sum(not v[5] for v in prim)}', flush=True)
    rep, rstats = replication_variants(D, {(v[1], int(v[0].rsplit("_", 1)[1][1:-1]), v[4]) for v in prim})
    print('replication:', rstats, flush=True)
    idrs = sorted({v[1] for v in prim if v[5]} | {v[1] for v in rep})
    with Pool(10, initializer=_init) as pool:
        wt = dict(pool.imap_unordered(_wt, [(n, D[n]['seq']) for n in idrs], chunksize=8))
    posm = positional_baseline()
    with open(os.path.join(OUT, 'wt_profiles.csv'), 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['idr', 'pos', 'aa', 'pred_consensus', 'pred_position'])
        for n in idrs:
            pp = posm.predict(pos_features(len(D[n]['seq'])))
            w.writerows([n, i, a, f'{wt[n][i]:.5f}', f'{pp[i]:.5f}'] for i, a in enumerate(D[n]['seq']))
    register(os.path.join(OUT, 'wt_profiles.csv'), f'ClinVar hub test: wild-type predictions for {len(idrs)} IDRs, '
                                                   'before any class was joined')
    print('wild-type profiles written and hashed', flush=True)
    jobs = [(v[0], D[v[1]]['seq'], v[2], v[4], wt[v[1]]) for v in prim if v[5]] + \
           [(v[0], D[v[1]]['seq'], v[2], v[4], wt[v[1]]) for v in rep]
    with Pool(10, initializer=_init) as pool:
        eff = {vid: (a, b) for vid, a, b in pool.imap_unordered(_mut, jobs, chunksize=32)}
    with open(os.path.join(OUT, 'variants.csv'), 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['set', 'variant', 'idr', 'gene', 'k', 'wt', 'mut', 'N', 'delta_sum', 'delta_site', 'variation_id'])
        for v in prim:
            if v[5]:
                w.writerow(['primary', v[0], v[1], D[v[1]]['gene'], v[2], v[3], v[4], len(D[v[1]]['seq']),
                            f'{eff[v[0]][0]:.5f}', f'{eff[v[0]][1]:.5f}', ''])
        for v in rep:
            w.writerow(['replication', v[0], v[1], D[v[1]]['gene'], v[2], v[3], v[4], len(D[v[1]]['seq']),
                        f'{eff[v[0]][0]:.5f}', f'{eff[v[0]][1]:.5f}', v[6]])
    register(os.path.join(OUT, 'variants.csv'), 'ClinVar hub test: variant positions and predicted effects, before any '
                                                'class was joined')
    print(f'variants written and hashed: primary {sum(v[5] for v in prim)}, replication {len(rep)}', flush=True)


# ---------------------------------------------------------------------------------------------------------- analysis
def logit(X, y, groups, ridge=1e-4):
    """Logistic regression by Newton steps; cluster-robust (by group) covariance with the G/(G-1) correction."""
    b = np.zeros(X.shape[1])
    for _ in range(100):
        p = 1 / (1 + np.exp(-(X @ b)))
        W = p * (1 - p)
        H = X.T @ (X * W[:, None]) + ridge * np.eye(X.shape[1])
        step = np.linalg.solve(H, X.T @ (y - p) - ridge * b)
        b += step
        if np.max(np.abs(step)) < 1e-9:
            break
    p = 1 / (1 + np.exp(-(X @ b)))
    Hinv = np.linalg.inv(X.T @ (X * (p * (1 - p))[:, None]) + ridge * np.eye(X.shape[1]))
    u = X * (y - p)[:, None]
    G = defaultdict(lambda: np.zeros(X.shape[1]))
    for g, ui in zip(groups, u):
        G[g] += ui
    S = np.zeros((X.shape[1], X.shape[1]))
    for v in G.values():
        S += np.outer(v, v)
    nG = len(G)
    V = Hinv @ S @ Hinv * nG / max(1, nG - 1)
    return b, V


def design(rows, score_key):
    cols = ['score', 'logend', 'loglen'] + [f'wt_{a}' for a in AAS[1:]] + [f'mut_{a}' for a in AAS[1:]] + ['const']
    X = []
    for r in rows:
        x = [r[score_key], math.log1p(r['end']), math.log(r['N'])]
        x += [1.0 if r['wt'] == a else 0.0 for a in AAS[1:]] + [1.0 if r['mut'] == a else 0.0 for a in AAS[1:]] + [1.0]
        X.append(x)
    X = np.array(X)
    sd = X[:, 1:3].std(0)
    X[:, 1:3] = (X[:, 1:3] - X[:, 1:3].mean(0)) / np.where(sd > 0, sd, 1)
    return X, cols


def fit_or(rows, score_key, per=0.25):
    X, _ = design(rows, score_key)
    y = np.array([r['y'] for r in rows], float)
    b, V = logit(X, y, [r['gene'] for r in rows])
    se = math.sqrt(V[0, 0])
    from scipy.stats import norm
    z = b[0] / se
    return dict(OR=math.exp(b[0] * per), lo=math.exp((b[0] - 1.96 * se) * per), hi=math.exp((b[0] + 1.96 * se) * per),
                P=2 * norm.sf(abs(z)), n=len(rows), n_path=int(y.sum()), beta=b[0])


def auc(pos, neg):
    from scipy.stats import mannwhitneyu
    if not pos or not neg:
        return float('nan')
    return mannwhitneyu(pos, neg).statistic / (len(pos) * len(neg))


def analyse():
    lines = []

    def p(s=''):
        lines.append(s)
        print(s, flush=True)

    # classes: IDRome's significance column (primary) and ClinVar's ClinicalSignificance (replication)
    cls = {}
    r = csv.reader(gzip.open(os.path.join(EXT, 'missense_variants.csv.gz'), 'rt'))
    h = next(r)
    ix = {c: i for i, c in enumerate(h)}
    sig_values = defaultdict(int)
    for row in r:
        s = row[ix['significance']]
        sig_values[s] += 1
        cls[row[0]] = 1 if 'athogenic' in s else 0 if 'enign' in s else None
    vids = set()
    with gzip.open(os.path.join(EXT, 'variant_summary_2026-10-05.txt.gz'), 'rt') as fh:
        rd = csv.reader(fh, delimiter='\t')
        hh = next(rd)
        jx = {c.lstrip('#'): i for i, c in enumerate(hh)}
        rep_cls = {}
        for row in rd:
            if row[jx['Assembly']] == 'GRCh38':
                cs = row[jx['ClinicalSignificance']]
                if cs in PATH_RC or cs in BEN_RC:
                    rep_cls[row[jx['VariationID']]] = 1 if cs in PATH_RC else 0
    wt = defaultdict(dict)
    for r in csv.DictReader(open(os.path.join(OUT, 'wt_profiles.csv'))):
        wt[r['idr']][int(r['pos'])] = (float(r['pred_consensus']), float(r['pred_position']))
    pct = {}
    for idr, d in wt.items():
        n = len(d)
        for j, key in enumerate(('cons', 'posb')):
            vals = np.array([d[i][j] for i in range(n)])
            ranks = vals.argsort().argsort()
            pct[(idr, key)] = ranks / max(1, n - 1)
    rows = []
    for r in csv.DictReader(open(os.path.join(OUT, 'variants.csv'))):
        y = cls.get(r['variant']) if r['set'] == 'primary' else rep_cls.get(r['variation_id'])
        if y is None:
            continue
        k, N = int(r['k']), int(r['N'])
        rows.append(dict(set=r['set'], variant=r['variant'], idr=r['idr'], gene=r['gene'], wt=r['wt'], mut=r['mut'],
                         N=N, end=min(k, N - 1 - k), y=y, score=pct[(r['idr'], 'cons')][k],
                         score_pos=pct[(r['idr'], 'posb')][k], dsum=float(r['delta_sum']),
                         dsite=float(r['delta_site'])))
    P = [r for r in rows if r['set'] == 'primary']
    R = [r for r in rows if r['set'] == 'replication']
    p('ClinVar hub-position test (clinvar_hub_plan.md, pre-registered 2026-10-09)')
    p(f'IDRome significance values: {dict(sig_values)}')
    p(f'primary set: {len(P)} variants ({sum(r["y"] for r in P)} pathogenic) in {len({r["idr"] for r in P})} IDRs, '
      f'{len({r["gene"] for r in P})} genes; replication set: {len(R)} variants ({sum(r["y"] for r in R)} pathogenic)')

    p('\nPRIMARY (odds ratio per 0.25 of within-IDR percentile of predicted occupancy; covariates wild-type and mutant '
      'amino acid, distance from the IDR ends, IDR length; SEs clustered by gene)')
    main = fit_or(P, 'score')
    p(f'  primary set:      OR {main["OR"]:.3f} [{main["lo"]:.3f}, {main["hi"]:.3f}], P = {main["P"]:.2g}, '
      f'n = {main["n"]} ({main["n_path"]} pathogenic)')
    noN = fit_or([r for r in P if r['gene'] != 'NOTCH3'], 'score')
    p(f'  without NOTCH3:   OR {noN["OR"]:.3f} [{noN["lo"]:.3f}, {noN["hi"]:.3f}], P = {noN["P"]:.2g}, n = {noN["n"]}')
    genes = defaultdict(int)
    for r in P:
        genes[r['gene']] += 1
    infl = []
    for g in genes:
        sub = [r for r in P if r['gene'] != g]
        infl.append((abs(fit_or(sub, 'score')['beta'] - main['beta']), g))
    g_top = max(infl)[1]
    noG = fit_or([r for r in P if r['gene'] != g_top], 'score')
    p(f'  without {g_top} (most influential gene): OR {noG["OR"]:.3f} [{noG["lo"]:.3f}, {noG["hi"]:.3f}], '
      f'P = {noG["P"]:.2g}')
    rep = fit_or(R, 'score') if R and 0 < sum(r['y'] for r in R) < len(R) else None
    if rep:
        p(f'  replication set:  OR {rep["OR"]:.3f} [{rep["lo"]:.3f}, {rep["hi"]:.3f}], P = {rep["P"]:.2g}, '
          f'n = {rep["n"]} ({rep["n_path"]} pathogenic)')
    both = defaultdict(set)
    for r in P:
        both[r['idr']].add(r['y'])
    sub = [r for r in P if both[r['idr']] == {0, 1}]
    p(f'  AUC of the score alone in IDRs with both classes: {auc([r["score"] for r in sub if r["y"]], [r["score"] for r in sub if not r["y"]]):.3f} '
      f'({len(sub)} variants, {len({r["idr"] for r in sub})} IDRs)')
    ok = (main['OR'] > 1 and main['P'] < 0.05 and noN['OR'] > 1 and noG['OR'] > 1 and rep is not None and rep['OR'] > 1)
    p(f'  pre-registered decision rule met: {ok}')

    p('\nSECONDARY')
    for key, lab in (('score_pos', 'positional baseline percentile'),):
        o = fit_or(P, key)
        p(f'  {lab}: OR {o["OR"]:.3f} [{o["lo"]:.3f}, {o["hi"]:.3f}], P = {o["P"]:.2g}')
    for key, lab in (('dsum', 'summed |change| of the predicted profile'), ('dsite', 'change at the variant residue')):
        vals = np.array([r[key] for r in P])
        mu, sd = vals.mean(), vals.std()
        for r in P:
            r[key + '_z'] = (r[key] - mu) / sd
        o = fit_or(P, key + '_z', per=1.0)
        p(f'  effect test, {lab} (per SD): OR {o["OR"]:.3f} [{o["lo"]:.3f}, {o["hi"]:.3f}], P = {o["P"]:.2g}')
    classes = {'charged': 'DEKR', 'aromatic': 'FWY', 'G/S/P': 'GSP', 'other': 'ACHILMNQTV'}
    for lab, aas in classes.items():
        sub = [r for r in P if r['wt'] in aas]
        if 0 < sum(r['y'] for r in sub) < len(sub):
            X, _ = design(sub, 'score')
            keep = [0, 1, 2] + [j for j in range(3, X.shape[1] - 1) if X[:, j].any()] + [X.shape[1] - 1]
            y = np.array([r['y'] for r in sub], float)
            b, V = logit(X[:, keep], y, [r['gene'] for r in sub])
            se = math.sqrt(V[0, 0])
            p(f'  wild-type {lab:8s}: OR {math.exp(b[0] * 0.25):.3f} [{math.exp((b[0] - 1.96 * se) * 0.25):.3f}, '
              f'{math.exp((b[0] + 1.96 * se) * 0.25):.3f}], n = {len(sub)} ({int(y.sum())} pathogenic)')
    p('\n  genes with the most variants (n, pathogenic, mean percentile pathogenic vs benign):')
    for g, n in sorted(genes.items(), key=lambda kv: -kv[1])[:12]:
        sp = [r['score'] for r in P if r['gene'] == g and r['y']]
        sb = [r['score'] for r in P if r['gene'] == g and not r['y']]
        p(f'    {g:10s} {n:4d} {len(sp):4d}  {np.mean(sp) if sp else float("nan"):.2f} vs {np.mean(sb) if sb else float("nan"):.2f}')
    open(os.path.join(OUT, 'clinvar_hub_summary.txt'), 'w').write('\n'.join(lines) + '\n')


if __name__ == '__main__':
    {'predict': predict, 'analyse': analyse}[sys.argv[1]]()
