"""
clinvar_hub_explore.py — EXPLORATORY follow-up of the pre-registered ClinVar test (clinvar_hub.py), written after its
result was known (2026-10-09). Asks whether the depletion of pathogenic variants at predicted hubs survives control for
local sequence composition and complexity, i.e. whether hub occupancy carries information beyond "this stretch is
glycine/serine-rich or low-complexity".

Added covariates (all from the wild-type IDR sequence, window of +-10 residues around the variant, site excluded):
  fGS    fraction of glycine + serine        fG   fraction of glycine
  H      Shannon entropy of the window composition (bits), a standard local-complexity measure
  idrGS  glycine + serine fraction of the whole IDR
Same regression, clustering and score definition as the pre-registered analysis.
  ~/.venvs/bender/bin/python clinvar_hub_explore.py  ->  clinvar_hub/clinvar_hub_explore.txt
"""
import csv
import gzip
import math
import os
import sys
from collections import Counter

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import clinvar_hub as C  # noqa: E402

W = 10


def local_feats(seq, k):
    win = seq[max(0, k - W):k] + seq[k + 1:k + W + 1]
    n = max(1, len(win))
    c = Counter(win)
    H = -sum(v / n * math.log2(v / n) for v in c.values())
    return dict(fGS=(c['G'] + c['S']) / n, fG=c['G'] / n, H=H,
                idrGS=(seq.count('G') + seq.count('S')) / len(seq))


def fit(rows, score_key, extra, per=0.25):
    X, _ = C.design(rows, score_key)
    if extra:
        E = np.array([[r[e] for e in extra] for r in rows], float)
        E = (E - E.mean(0)) / np.where(E.std(0) > 0, E.std(0), 1)
        X = np.hstack([X[:, :-1], E, X[:, -1:]])
    y = np.array([r['y'] for r in rows], float)
    b, V = C.logit(X, y, [r['gene'] for r in rows])
    se = math.sqrt(V[0, 0])
    from scipy.stats import norm
    return math.exp(b[0] * per), math.exp((b[0] - 1.96 * se) * per), math.exp((b[0] + 1.96 * se) * per), \
        2 * norm.sf(abs(b[0] / se))


def main():
    lines = []

    def p(s=''):
        lines.append(s)
        print(s, flush=True)

    D = C.idrome()
    cls = {}
    r = csv.reader(gzip.open(os.path.join(C.EXT, 'missense_variants.csv.gz'), 'rt'))
    h = next(r)
    ix = {c: i for i, c in enumerate(h)}
    for row in r:
        s = row[ix['significance']]
        cls[row[0]] = 1 if 'athogenic' in s else 0 if 'enign' in s else None
    rep_cls = {}
    with gzip.open(os.path.join(C.EXT, 'variant_summary_2026-10-05.txt.gz'), 'rt') as fh:
        rd = csv.reader(fh, delimiter='\t')
        hh = next(rd)
        jx = {c.lstrip('#'): i for i, c in enumerate(hh)}
        for row in rd:
            if row[jx['Assembly']] == 'GRCh38':
                cs = row[jx['ClinicalSignificance']]
                if cs in C.PATH_RC or cs in C.BEN_RC:
                    rep_cls[row[jx['VariationID']]] = 1 if cs in C.PATH_RC else 0
    wt = {}
    cur = None
    for rr in csv.DictReader(open(os.path.join(C.OUT, 'wt_profiles.csv'))):
        wt.setdefault(rr['idr'], []).append(float(rr['pred_consensus']))
    pct = {k: np.array(v).argsort().argsort() / max(1, len(v) - 1) for k, v in wt.items()}
    rows = []
    for rr in csv.DictReader(open(os.path.join(C.OUT, 'variants.csv'))):
        y = cls.get(rr['variant']) if rr['set'] == 'primary' else rep_cls.get(rr['variation_id'])
        if y is None:
            continue
        k, N = int(rr['k']), int(rr['N'])
        seq = D[rr['idr']]['seq']
        row = dict(set=rr['set'], idr=rr['idr'], gene=rr['gene'], wt=rr['wt'], mut=rr['mut'], N=N,
                   end=min(k, N - 1 - k), y=y, score=pct[rr['idr']][k], dsite=float(rr['delta_site']))
        row.update(local_feats(seq, k))
        rows.append(row)
    for s in ('primary', 'replication'):
        S = [r for r in rows if r['set'] == s]
        mu, sd = np.mean([r['dsite'] for r in S]), np.std([r['dsite'] for r in S])
        for r in S:
            r['dsite_z'] = (r['dsite'] - mu) / sd

    p('EXPLORATORY (written after the pre-registered result was known): does the hub signal survive control for local '
      'composition and complexity?')
    p('Odds ratio per 0.25 of within-IDR percentile of predicted occupancy (position test) and per SD of the predicted '
      'change at the variant residue (effect test); base covariates as pre-registered')
    sets = [('none (pre-registered model)', []), ('+ local G+S fraction', ['fGS']), ('+ local G fraction', ['fG']),
            ('+ local complexity (entropy)', ['H']), ('+ local G+S, G, entropy, IDR G+S', ['fGS', 'fG', 'H', 'idrGS'])]
    for s in ('primary', 'replication'):
        S = [r for r in rows if r['set'] == s]
        p(f'\n{s} set: {len(S)} variants ({sum(r["y"] for r in S)} pathogenic)')
        p(f'  {"added covariates":38s} {"position test OR [95%], P":34s} {"effect test OR [95%], P":30s}')
        for lab, extra in sets:
            a = fit(S, 'score', extra)
            e = fit(S, 'dsite_z', extra, per=1.0)
            p(f'  {lab:38s} {a[0]:.3f} [{a[1]:.3f}, {a[2]:.3f}], {a[3]:.1e}   {e[0]:.3f} [{e[1]:.3f}, {e[2]:.3f}], {e[3]:.1e}')
        NC = [r for r in S if not r['gene'].startswith('COL')]
        a = fit(NC, 'score', ['fGS', 'fG', 'H', 'idrGS'])
        e = fit(NC, 'dsite_z', ['fGS', 'fG', 'H', 'idrGS'], per=1.0)
        p(f'  {"full controls, collagens (COL*) removed":38s} {a[0]:.3f} [{a[1]:.3f}, {a[2]:.3f}], {a[3]:.1e}   '
          f'{e[0]:.3f} [{e[1]:.3f}, {e[2]:.3f}], {e[3]:.1e}   (n = {len(NC)})')
        # how strongly the score itself tracks local composition
        sc = np.array([r['score'] for r in S])
        for f in ('fGS', 'fG', 'H'):
            from scipy.stats import spearmanr
            p(f'  Spearman(score, {f}) = {spearmanr(sc, [r[f] for r in S]).correlation:+.2f}')
    open(os.path.join(C.OUT, 'clinvar_hub_explore.txt'), 'w').write('\n'.join(lines) + '\n')


if __name__ == '__main__':
    main()
