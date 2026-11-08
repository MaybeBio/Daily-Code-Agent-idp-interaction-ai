"""
clinvar_hub_idrome_style.py — the ClinVar comparison done the way IDRome (Tesei et al. 2024, Fig. 5e-g and Extended Data
Fig. 9) did it: unadjusted two-sided Brunner-Munzel tests of pathogenic versus benign variants, with Cohen's d (pathogenic
minus benign, pooled SD; standard error from 1,000 bootstrap resamples) and the Brunner-Munzel degrees of freedom, plus
IDRome's robustness check of removing the NOTCH3 cysteine variants. Collagens are also shown separately.
Properties: predicted hub occupancy of the wild-type residue (raw and as within-IDR percentile), the predicted change of
occupancy at the variant residue, and the summed absolute change over the IDR (as IDRome's Delta-nu and Delta-Sconf/N).
Written after the pre-registered result (clinvar_hub.py) was known, so exploratory.
  ~/.venvs/bender/bin/python clinvar_hub_idrome_style.py  ->  clinvar_hub/clinvar_hub_idrome_style.txt
"""
import csv
import gzip
import os
import sys

import numpy as np
from scipy.stats import brunnermunzel, rankdata

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import clinvar_hub as C  # noqa: E402


def bm_df(x, y):
    """Satterthwaite-type degrees of freedom of the Brunner-Munzel test."""
    nx, ny = len(x), len(y)
    r = rankdata(np.concatenate([x, y]))
    rx, ry = r[:nx], r[nx:]
    rxi, ryi = rankdata(x), rankdata(y)
    sx = np.sum((rx - rxi - rx.mean() + (nx + 1) / 2) ** 2) / (nx - 1)
    sy = np.sum((ry - ryi - ry.mean() + (ny + 1) / 2) ** 2) / (ny - 1)
    a, b = sx / ny, sy / nx
    return (nx * a + ny * b) ** 2 / ((nx * a) ** 2 / (nx - 1) + (ny * b) ** 2 / (ny - 1))


def cohen(x, y, rng, nboot=1000):
    def d(a, b):
        s = np.sqrt(((len(a) - 1) * a.var(ddof=1) + (len(b) - 1) * b.var(ddof=1)) / (len(a) + len(b) - 2))
        return (a.mean() - b.mean()) / s
    boots = [d(rng.choice(x, len(x)), rng.choice(y, len(y))) for _ in range(nboot)]
    return d(x, y), float(np.std(boots))


def main():
    lines = []

    def p(s=''):
        lines.append(s)
        print(s, flush=True)

    cls = {}
    rd = csv.reader(gzip.open(os.path.join(C.EXT, 'missense_variants.csv.gz'), 'rt'))
    h = next(rd)
    ix = {c: i for i, c in enumerate(h)}
    for row in rd:
        s = row[ix['significance']]
        cls[row[0]] = 1 if 'athogenic' in s else 0 if 'enign' in s else None
    rep_cls = {}
    with gzip.open(os.path.join(C.EXT, 'variant_summary_2026-10-05.txt.gz'), 'rt') as fh:
        r2 = csv.reader(fh, delimiter='\t')
        hh = next(r2)
        jx = {c.lstrip('#'): i for i, c in enumerate(hh)}
        for row in r2:
            if row[jx['Assembly']] == 'GRCh38':
                cs = row[jx['ClinicalSignificance']]
                if cs in C.PATH_RC or cs in C.BEN_RC:
                    rep_cls[row[jx['VariationID']]] = 1 if cs in C.PATH_RC else 0
    wt = {}
    for rr in csv.DictReader(open(os.path.join(C.OUT, 'wt_profiles.csv'))):
        wt.setdefault(rr['idr'], []).append(float(rr['pred_consensus']))
    pct = {k: np.array(v).argsort().argsort() / max(1, len(v) - 1) for k, v in wt.items()}
    rows = []
    for rr in csv.DictReader(open(os.path.join(C.OUT, 'variants.csv'))):
        y = cls.get(rr['variant']) if rr['set'] == 'primary' else rep_cls.get(rr['variation_id'])
        if y is None:
            continue
        k = int(rr['k'])
        rows.append(dict(set=rr['set'], gene=rr['gene'], wt=rr['wt'], y=y, occ=wt[rr['idr']][k], pct=pct[rr['idr']][k],
                         dsite=float(rr['delta_site']), dsum=float(rr['delta_sum'])))
    props = [('occ', 'predicted occupancy of the wild-type residue'), ('pct', 'within-IDR percentile of that occupancy'),
             ('dsite', 'predicted change at the variant residue'), ('dsum', 'summed |change| over the IDR')]
    rng = np.random.default_rng(0)
    p('ClinVar comparison in the style of IDRome (unadjusted two-sided Brunner-Munzel tests, pathogenic vs benign; '
      'Cohen\'s d = pathogenic minus benign). EXPLORATORY.')
    for s, sname in (('primary', 'IDRome ClinVar set'), ('replication', 'newer ClinVar (5 Oct 2026)')):
        S = [r for r in rows if r['set'] == s]
        subsets = [('all variants', S),
                   ('without NOTCH3 cysteine variants (IDRome check)', [r for r in S if not (r['gene'] == 'NOTCH3' and r['wt'] == 'C')]),
                   ('without collagens (COL*)', [r for r in S if not r['gene'].startswith('COL')])]
        p(f'\n{sname}: {len(S)} variants ({sum(r["y"] for r in S)} pathogenic)')
        for lab, sub in subsets:
            p(f'  {lab} (n = {len(sub)}, {sum(r["y"] for r in sub)} pathogenic)')
            for key, desc in props:
                x = np.array([r[key] for r in sub if r['y']])
                yb = np.array([r[key] for r in sub if not r['y']])
                bm = brunnermunzel(x, yb)
                d, se = cohen(x, yb, rng)
                p(f'    {desc:46s} median {np.median(x):+.4f} vs {np.median(yb):+.4f}   Cohen\'s d {d:+.2f} +- {se:.2f}   '
                  f'DoF {bm_df(x, yb):,.0f}   P = {bm.pvalue:.1e}')
    open(os.path.join(C.OUT, 'clinvar_hub_idrome_style.txt'), 'w').write('\n'.join(lines) + '\n')


if __name__ == '__main__':
    main()
