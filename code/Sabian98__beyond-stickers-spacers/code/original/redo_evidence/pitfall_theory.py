"""
pitfall_theory.py — the three analysis pitfalls as derived statements, with the numbers quoted in the paper.

1  Near-neighbour contacts. With |i-j| >= 2, residues i and i+2 (C-alpha atoms at most 2 x 0.38 nm apart) are in contact
   in every snapshot. If these were the only edges, the contact graph would be two paths (even and odd residues), and the
   residue at index j of a path of m residues lies on exactly j(m-1-j) shortest paths. Tested on the 1,198 CALVADOS-2
   proteins analysed with both rules (as in Fig. 1): variance of time-averaged betweenness explained by this closed form,
   and the residue-type ranking it predicts (mean within-protein z over the positions each type occupies) against the
   observed one.
2  Sampling. For independent snapshots the reliability of an n-snapshot average follows Spearman-Brown,
   R(n) = n r1 / (1 + (n - 1) r1). One single-snapshot reliability r1 is fitted to the split-half medians of
   frames_reliability_summary.txt, and n* = R*(1 - r1) / (r1 (1 - R*)) snapshots give reliability R*.
3  Size-scaled cutoffs. In the residue-neutral null ensembles matched to the 908 paired CALVADOS-2 ensembles
   (sim_null_per_residue_all.csv), protein-level mean position-adjusted z per residue type, tested against zero with
   Wilcoxon signed-rank tests over proteins (Bonferroni over 20 types), for the fixed 0.8 nm rule and the size-scaled rule.
  ~/.venvs/bender/bin/python pitfall_theory.py  ->  pitfall_theory_summary.txt
"""
import csv
import os
import re
from collections import defaultdict

import numpy as np
from scipy.stats import pearsonr, spearmanr, wilcoxon

HERE = os.path.dirname(os.path.abspath(__file__))
AAS = 'ACDEFGHIKLMNPQRSTVWY'
out = []


def p(s=''):
    out.append(s)
    print(s, flush=True)


def two_path(n):
    b = np.zeros(n)
    for par in (0, 1):
        idx = np.arange(par, n, 2)
        j = np.arange(len(idx))
        b[idx] = j * (len(idx) - 1 - j)
    return b


def load(d):
    by = defaultdict(list)
    for r in csv.DictReader(open(os.path.join(HERE, d, 'per_residue.csv'))):
        by[r['UniProt_ID']].append((int(r['pos']), r['aa'], float(r['mean_btw'])))
    return {u: sorted(v) for u, v in by.items()}


def part1():
    d2, d3 = load('hub_frames_sep2'), load('hub_frames')
    ids = sorted(set(d2) & set(d3))
    p(f'1  Near-neighbour contacts ({len(ids)} CALVADOS-2 proteins analysed with both rules)')
    for d, label in ((d2, '|i-j| >= 2'), (d3, '|i-j| >= 3')):
        r2, obs, pred = [], defaultdict(list), defaultdict(list)
        for u in ids:
            rows = d[u]
            b = np.array([x[2] for x in rows])
            t = two_path(len(b))
            r2.append(pearsonr(t, b)[0] ** 2 if b.std() > 0 else np.nan)
            zb = (b - b.mean()) / (b.std() or 1)
            zt = (t - t.mean()) / (t.std() or 1)
            for (_, a, _), x, y in zip(rows, zb, zt):
                obs[a].append(x)
                pred[a].append(y)
        o = np.array([np.mean(obs[a]) for a in AAS])
        q = np.array([np.mean(pred[a]) for a in AAS])
        p(f'  {label}: closed form explains a median {100 * np.nanmedian(r2):.0f}% of the variance of mean betweenness; '
          f'predicted vs observed residue-type ranking r = {pearsonr(o, q)[0]:.2f} (Spearman {spearmanr(o, q).correlation:.2f})')


def part2():
    txt = open(os.path.join(HERE, 'frames_reliability_summary.txt')).read()
    pts = [(int(a), float(b)) for a, b in re.findall(r'^\s+(\d+)\s+([0-9.]+)\s+[0-9.]+\s*$', txt, flags=re.M)]
    r1s = [R / (n - (n - 1) * R) for n, R in pts]
    r1 = float(np.median(r1s))
    p(f'\n2  Sampling: single-snapshot reliability r1 = {r1:.4f} (implied by each point: '
      + ', '.join(f'{v:.4f}' for v in r1s) + ')')
    p('  Spearman-Brown with this r1 vs observed split-half medians: '
      + ', '.join(f'{n}: {n * r1 / (1 + (n - 1) * r1):.2f}/{R:.2f}' for n, R in pts))
    p('  snapshots needed: ' + ', '.join(f'R* = {Rs}: {Rs * (1 - r1) / (r1 * (1 - Rs)):.0f}' for Rs in (0.8, 0.9, 0.95))
      + f'; reliability of 1,000 snapshots {1000 * r1 / (1 + 999 * r1):.2f}')


def part3():
    per = {rule: defaultdict(lambda: defaultdict(list)) for rule in ('fixed', 'size')}
    for r in csv.DictReader(open(os.path.join(HERE, 'sim_null_per_residue_all.csv'))):
        if r['force_field'] != 'CALVADOS-2':
            continue
        for rule in ('fixed', 'size'):
            if r[f'z_null_{rule}']:
                per[rule][r['aa']][r['UniProt_ID']].append(float(r[f'z_null_{rule}']))
    n_prot = len({u for a in per['fixed'] for u in per['fixed'][a]})
    p(f'\n3  Size-scaled cutoffs: residue-neutral null ensembles matched to {n_prot} CALVADOS-2 ensembles')
    for rule in ('fixed', 'size'):
        rows = []
        for a in AAS:
            pm = np.array([np.mean(v) for v in per[rule][a].values()])
            rows.append((a, pm.mean(), wilcoxon(pm).pvalue))
        rows.sort(key=lambda x: -x[1])
        sig = sum(pv < 0.05 / 20 for _, _, pv in rows)
        big = [a for a, m, _ in rows if abs(m) >= 0.2]
        p(f'  {rule:5s}: {sig}/20 types differ from zero (Bonferroni); {len(big)} with |mean z| >= 0.2 ({"".join(big)}); '
          f'range {rows[-1][1]:+.2f} to {rows[0][1]:+.2f}')
        p('         most central ' + ', '.join(f'{a} {m:+.2f}' for a, m, _ in rows[:4])
          + ' | least central ' + ', '.join(f'{a} {m:+.2f}' for a, m, _ in rows[-4:]))


if __name__ == '__main__':
    part1()
    part2()
    part3()
    open(os.path.join(HERE, 'pitfall_theory_summary.txt'), 'w').write('\n'.join(out) + '\n')
