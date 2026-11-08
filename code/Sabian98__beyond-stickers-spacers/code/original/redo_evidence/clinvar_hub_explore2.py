"""
clinvar_hub_explore2.py — EXPLORATORY check 1 (2026-10-09): does the effect-test signal (substitutions that raise the
predicted occupancy of their own residue are more often pathogenic) survive control for how drastic the substitution is?
Adds BLOSUM62 score and Grantham distance of each wild-type -> mutant pair to the full composition-controlled model of
clinvar_hub_explore.py, with and without collagens.
Grantham distances are computed from Grantham's (1974) composition, polarity and volume values with his weights.
  ~/.venvs/bender/bin/python clinvar_hub_explore2.py  ->  clinvar_hub/clinvar_hub_explore2.txt
"""
import csv
import gzip
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import clinvar_hub as C  # noqa: E402
import clinvar_hub_explore as E  # noqa: E402

ORDER = 'ARNDCQEGHILKMFPSTWYV'
B62 = """
 4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0
-1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3
-2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3
-2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3
 0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1
-1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2
-1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2
 0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3
-2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3
-1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3
-1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1
-1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2
-1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1
-2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1
-1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2
 1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2
 0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0
-3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3
-2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1
 0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4"""
M = np.array([[int(x) for x in line.split()] for line in B62.strip().split('\n')])
assert M.shape == (20, 20) and (M == M.T).all()
BLOSUM = {(a, b): int(M[i, j]) for i, a in enumerate(ORDER) for j, b in enumerate(ORDER)}
# Grantham (1974): composition c, polarity p, volume v
GP = {'S': (1.42, 9.2, 32), 'R': (0.65, 10.5, 124), 'L': (0, 4.9, 111), 'P': (0.39, 8.0, 32.5), 'T': (0.71, 8.6, 61),
      'A': (0, 8.1, 31), 'V': (0, 5.9, 84), 'G': (0.74, 9.0, 3), 'I': (0, 5.2, 111), 'F': (0, 5.2, 132),
      'Y': (0.20, 6.2, 136), 'C': (2.75, 5.5, 55), 'H': (0.58, 10.4, 96), 'Q': (0.89, 10.5, 85), 'N': (1.33, 11.6, 56),
      'K': (0.33, 11.3, 119), 'D': (1.38, 13.0, 54), 'E': (0.92, 12.3, 83), 'M': (0, 5.7, 105), 'W': (0.13, 5.4, 170)}


def grantham(a, b):
    (c1, p1, v1), (c2, p2, v2) = GP[a], GP[b]
    return 50.723 * math.sqrt(1.833 * (c1 - c2) ** 2 + 0.1018 * (p1 - p2) ** 2 + 0.000399 * (v1 - v2) ** 2)


def main():
    assert abs(grantham('R', 'S') - 110) < 3 and abs(grantham('C', 'W') - 215) < 3, 'Grantham check'
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
    for rr in csv.DictReader(open(os.path.join(C.OUT, 'wt_profiles.csv'))):
        wt.setdefault(rr['idr'], []).append(float(rr['pred_consensus']))
    pct = {k: np.array(v).argsort().argsort() / max(1, len(v) - 1) for k, v in wt.items()}
    rows = []
    for rr in csv.DictReader(open(os.path.join(C.OUT, 'variants.csv'))):
        y = cls.get(rr['variant']) if rr['set'] == 'primary' else rep_cls.get(rr['variation_id'])
        if y is None:
            continue
        k, N = int(rr['k']), int(rr['N'])
        row = dict(set=rr['set'], idr=rr['idr'], gene=rr['gene'], wt=rr['wt'], mut=rr['mut'], N=N,
                   end=min(k, N - 1 - k), y=y, score=pct[rr['idr']][k], dsite=float(rr['delta_site']),
                   blosum=BLOSUM[(rr['wt'], rr['mut'])], grantham=grantham(rr['wt'], rr['mut']))
        row.update(E.local_feats(D[rr['idr']]['seq'], k))
        rows.append(row)
    comp = ['fGS', 'fG', 'H', 'idrGS']
    p('EXPLORATORY check 1: substitution severity (BLOSUM62 score, Grantham distance) added to the composition-controlled '
      'model. Effect test = OR per SD of the predicted change in occupancy at the variant residue.')
    for s in ('primary', 'replication'):
        S = [r for r in rows if r['set'] == s]
        mu, sd = np.mean([r['dsite'] for r in S]), np.std([r['dsite'] for r in S])
        for r in S:
            r['dsite_z'] = (r['dsite'] - mu) / sd
        from scipy.stats import spearmanr
        p(f'\n{s} set: {len(S)} variants ({sum(r["y"] for r in S)} pathogenic)')
        p(f'  Spearman(predicted change, BLOSUM62) = {spearmanr([r["dsite"] for r in S], [r["blosum"] for r in S]).correlation:+.2f}, '
          f'(predicted change, Grantham) = {spearmanr([r["dsite"] for r in S], [r["grantham"] for r in S]).correlation:+.2f}')
        for lab, extra, sub in (('composition (as before)', comp, S),
                                ('composition + BLOSUM62', comp + ['blosum'], S),
                                ('composition + Grantham', comp + ['grantham'], S),
                                ('composition + BLOSUM62 + Grantham', comp + ['blosum', 'grantham'], S),
                                ('same, collagens removed', comp + ['blosum', 'grantham'],
                                 [r for r in S if not r['gene'].startswith('COL')])):
            e = E.fit(sub, 'dsite_z', extra, per=1.0)
            a = E.fit(sub, 'score', extra)
            p(f'  {lab:36s} effect {e[0]:.3f} [{e[1]:.3f}, {e[2]:.3f}], P = {e[3]:.1e}   '
              f'position {a[0]:.3f} [{a[1]:.3f}, {a[2]:.3f}], P = {a[3]:.1e}   (n = {len(sub)})')
        # the severity terms themselves, for scale
        Xs = [r for r in S]
        for key in ('blosum', 'grantham'):
            vals = np.array([r[key] for r in Xs], float)
            for r, v in zip(Xs, (vals - vals.mean()) / vals.std()):
                r[key + '_z'] = v
            o = E.fit(Xs, key + '_z', comp, per=1.0)
            p(f'  for scale: {key} alone (per SD, composition-controlled) OR {o[0]:.3f} [{o[1]:.3f}, {o[2]:.3f}]')
    open(os.path.join(C.OUT, 'clinvar_hub_explore2.txt'), 'w').write('\n'.join(lines) + '\n')


if __name__ == '__main__':
    main()
