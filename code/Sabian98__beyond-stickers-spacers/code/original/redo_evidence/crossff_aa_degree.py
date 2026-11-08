"""Per-amino-acid within-sequence z of contact degree (k, |i-j|>=3) and long-range degree (k_lr, |i-j|>=6),
CALVADOS-2 vs MPIPI-GG, 957 sequences (Pilot_hub_context_crossff/per_residue.csv, corrected rules).
Pearson r across the 20 amino acids per group, as in the paper's cross-FF claim (which used betweenness)."""
import csv, math
from collections import defaultdict
PR = '/Volumes/User Homes/Nu Project/BENDER-MT/Traj Level Analysis/Pilot_hub_context_crossff/per_residue.csv'
BIO = '/Volumes/User Homes/Nu Project/BENDER-BIO/BENDER_BIO.csv'
grp = {r['UniProt_ID']: r['kingdom'] for r in csv.DictReader(open(BIO))}
seqs = defaultdict(lambda: defaultdict(list))   # (uid, ff) -> col -> list of (aa, val)
for r in csv.DictReader(open(PR)):
    for col in ('k', 'k_lr'):
        seqs[(r['UniProt_ID'], r['ff'])][col].append((r['aa'], float(r[col])))
acc = defaultdict(lambda: [0.0, 0])   # (ff, col, group, aa)
for (uid, ff), cols in seqs.items():
    g = grp.get(uid, '?')
    for col, L in cols.items():
        v = [x for _, x in L]
        m = sum(v) / len(v); sd = math.sqrt(sum((x - m) ** 2 for x in v) / len(v))
        if sd < 1e-12: continue
        for aa, x in L:
            for gg in (g, 'All'):
                a = acc[(ff, col, gg, aa)]; a[0] += (x - m) / sd; a[1] += 1
def pear(a, b):
    n = len(a); ma, mb = sum(a) / n, sum(b) / n
    sa = math.sqrt(sum((x - ma) ** 2 for x in a)); sb = math.sqrt(sum((y - mb) ** 2 for y in b))
    return sum((x - ma) * (y - mb) for x, y in zip(a, b)) / (sa * sb)
AAS = 'ACDEFGHIKLMNPQRSTVWY'
groups = sorted({k[2] for k in acc})
for col in ('k', 'k_lr'):
    print(f'\n== {col}: Pearson r across 20 AAs (CALVADOS-2 vs MPIPI-GG mean within-sequence z)')
    for g in groups:
        c = [acc[('calvados', col, g, aa)][0] / max(acc[('calvados', col, g, aa)][1], 1) for aa in AAS]
        m = [acc[('mpipi', col, g, aa)][0] / max(acc[('mpipi', col, g, aa)][1], 1) for aa in AAS]
        print(f'  {g:10s} r = {pear(c, m):+.3f}')
    print('  All, per AA (CAL / MPIPI):', '  '.join(f"{aa} {acc[('calvados', col, 'All', aa)][0] / acc[('calvados', col, 'All', aa)][1]:+.2f}/{acc[('mpipi', col, 'All', aa)][0] / acc[('mpipi', col, 'All', aa)][1]:+.2f}" for aa in AAS))
