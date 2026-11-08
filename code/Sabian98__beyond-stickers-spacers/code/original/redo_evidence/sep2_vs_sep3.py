"""
sep2_vs_sep3.py — does the original contact rule (|i-j| >= 2) reproduce the original manuscript's residue ranking on
the same trajectories, and does it disappear with |i-j| >= 3? (Gap G1 of methods_paper/PLAN.md.)

Inputs: hub_frames_sep2/ (|i-j| >= 2, 200 frames) and hub_frames/ (|i-j| >= 3, 1,000 frames), same 1,198 CALVADOS-2
proteins. Paper method: within-protein z of time-averaged betweenness over all residues (no position adjustment),
pooled by residue type and group. Also: z of the first residue, and the manuscript's published values.
Standard library only.  python3 sep2_vs_sep3.py  ->  sep2_vs_sep3_summary.txt
"""
import csv
import math
import os
import sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from hub_z_v2_stdlib import zs, quad_resid  # noqa: E402

GROUPS = ['Bacteria', 'Fungi', 'Mammals', 'Plants', 'Protists', 'Viruses']
SHOW = 'EDRPKMCFYWSHGATLIVNQ'
PAPER = {'E': '+0.036…+0.144', 'D': '+0.036…+0.110', 'R': '+0.026…+0.121', 'P': '+0.034…+0.075', 'K': '−0.003 (mean)',
         'M': '−0.616…−1.067', 'C': '−0.083…−0.530', 'F': 'anti 5/6 (Pl −0.199, Vi −0.153)', 'Y': 'Ba +0.176 … Vi −0.169',
         'W': 'Vi −0.230', 'S': '≈ −0.046 (Ba, Pl)'}
out = []


def p(s=''):
    out.append(s)
    print(s, flush=True)


def load(dirname):
    prot = defaultdict(dict)
    for r in csv.DictReader(open(os.path.join(HERE, dirname, 'per_residue.csv'))):
        prot[r['UniProt_ID']][int(r['pos'])] = (r['aa'], float(r['mean_btw']), r['group'])
    res = {}
    for u, d in prot.items():
        n = len(d)
        seq = ''.join(d[i][0] for i in range(n))
        b = [d[i][1] for i in range(n)]
        z = zs(b)
        if z is None:
            continue
        res[u] = (d[0][2], seq, z, quad_resid(z))
    return res


def table(data, which):
    acc = defaultdict(lambda: [0.0, 0.0, 0])
    for u, (g, seq, z, zr) in data.items():
        v = z if which == 'raw' else zr
        for a, x in zip(seq, v):
            for gg in (g, 'All'):
                t = acc[(gg, a)]
                t[0] += x
                t[1] += x * x
                t[2] += 1
    return {k: (s / n, math.sqrt(max(ss / n - (s / n) ** 2, 0)) / math.sqrt(n)) for k, (s, ss, n) in acc.items()}


def main():
    d2, d3 = load('hub_frames_sep2'), load('hub_frames')
    common = sorted(set(d2) & set(d3))
    d2 = {u: d2[u] for u in common}
    d3 = {u: d3[u] for u in common}
    p(f'Original rule (|i-j| >= 2, 200 frames) vs corrected rule (|i-j| >= 3, 1,000 frames), same {len(common)} CALVADOS-2 proteins')
    t2r, t3r, t2p, t3p = table(d2, 'raw'), table(d3, 'raw'), table(d2, 'pos'), table(d3, 'pos')
    p('\nPaper method (within-protein z, all residues, no position adjustment), all groups; per-group range for |i-j| >= 2')
    p(f'  {"aa":3s} {"|i-j|>=2":>9s} {"range over groups (>=2)":>26s} {"|i-j|>=3":>9s} {">=2 pos-adj":>12s} {">=3 pos-adj":>12s}   manuscript')
    for aa in SHOW:
        rng = [t2r[(g, aa)][0] for g in GROUPS if (g, aa) in t2r]
        p(f'  {aa:3s} {t2r[("All", aa)][0]:+9.3f} {min(rng):+12.3f} … {max(rng):+.3f} {t3r[("All", aa)][0]:+9.3f} '
          f'{t2p[("All", aa)][0]:+12.3f} {t3p[("All", aa)][0]:+12.3f}   {PAPER.get(aa, "")}')
    for name, d in (('|i-j| >= 2', d2), ('|i-j| >= 3', d3)):
        first = [z[0] for g, seq, z, zr in d.values()]
        last = [z[-1] for g, seq, z, zr in d.values()]
        second = [z[1] for g, seq, z, zr in d.values()]
        p(f'\n  {name}: mean z of residue 1 {sum(first) / len(first):+.2f}, residue 2 {sum(second) / len(second):+.2f}, '
          f'last residue {sum(last) / len(last):+.2f}')
    # which charged/aromatic contrasts flip
    p('\nContrasts (all groups, paper method):')
    for name, t in (('|i-j| >= 2', t2r), ('|i-j| >= 3', t3r)):
        chg = sum(t[('All', a)][0] for a in 'EDR') / 3
        aro = sum(t[('All', a)][0] for a in 'FWY') / 3
        p(f'  {name}: mean(E, D, R) {chg:+.3f}; mean(F, W, Y) {aro:+.3f}; E − Y {t[("All", "E")][0] - t[("All", "Y")][0]:+.3f}; '
          f'R − K {t[("All", "R")][0] - t[("All", "K")][0]:+.3f}; M {t[("All", "M")][0]:+.3f}')
    open(os.path.join(HERE, 'sep2_vs_sep3_summary.txt'), 'w').write('\n'.join(out) + '\n')


if __name__ == '__main__':
    main()
