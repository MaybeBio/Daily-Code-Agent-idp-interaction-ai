"""
hub_frames_crossff.py — CALVADOS-2 vs MPIPI-GG on the same sequences (hub_frames_paired_calvados/ vs hub_frames_mpipi/).

Per protein: within-protein Pearson r of time-averaged betweenness and of hub occupancy between the force fields,
against the ceiling set by each force field's own noise (sqrt of the two Spearman-Brown split-half reliabilities);
overlap (Jaccard) of the two force fields' hub sets (top 10% by time-averaged betweenness) vs random sets.
Per residue type: hub occupancy and long-range share of contacts in hub frames, each force field.
Standard library only.  python3 hub_frames_crossff.py  ->  prints and writes hub_frames_mpipi/crossff_summary.txt
"""
import csv
import math
import os
import statistics as st
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
CAL, MPI = os.path.join(HERE, 'hub_frames_paired_calvados'), os.path.join(HERE, 'hub_frames_mpipi')
AAS = 'GSACYFWNTPQHDRELIVMK'
out = []


def p(s=''):
    out.append(s)
    print(s)


def pearson(a, b):
    n = len(a)
    ma, mb = sum(a) / n, sum(b) / n
    sa = math.sqrt(sum((x - ma) ** 2 for x in a))
    sb = math.sqrt(sum((y - mb) ** 2 for y in b))
    return sum((x - ma) * (y - mb) for x, y in zip(a, b)) / (sa * sb) if sa and sb else None


def q(v, a):
    v = sorted(x for x in v if x is not None)
    return v[int(a * (len(v) - 1))]


def load(d):
    res = defaultdict(dict)
    for r in csv.DictReader(open(f'{d}/per_residue.csv')):
        res[r['UniProt_ID']][int(r['pos'])] = r
    prot = {r['UniProt_ID']: r for r in csv.DictReader(open(f'{d}/per_protein.csv'))}
    return res, prot


rc, pc = load(CAL)
rm, pm = load(MPI)
ids = sorted(set(rc) & set(rm))
p(f'CALVADOS-2 vs MPIPI-GG, same sequences: {len(ids)} proteins with both '
  f'({len(rc)} CALVADOS-2, {len(rm)} MPIPI-GG analysed)')

rows = []
for u in ids:
    a, b = rc[u], rm[u]
    if len(a) != len(b):
        continue
    n = len(a)
    ba = [float(a[i]['mean_btw']) for i in range(n)]
    bb = [float(b[i]['mean_btw']) for i in range(n)]
    oa = [float(a[i]['occupancy']) for i in range(n)]
    ob = [float(b[i]['occupancy']) for i in range(n)]
    r_b, r_o = pearson(ba, bb), pearson(oa, ob)
    sa, sb = pc[u]['splithalf_spearman_brown_btw'], pm[u]['splithalf_spearman_brown_btw']
    ceil = math.sqrt(float(sa) * float(sb)) if sa and sb and float(sa) > 0 and float(sb) > 0 else None
    ha = {i for i in range(n) if a[i]['class'] == 'hub'}
    hb = {i for i in range(n) if b[i]['class'] == 'hub'}
    m = len(ha)
    j = len(ha & hb) / len(ha | hb) if ha | hb else None
    rows.append((u, n, r_b, r_o, ceil, j, m / (2 * n - m) if n > m else None, a[0]['group']))

rb = [r[2] for r in rows if r[2] is not None]
ro = [r[3] for r in rows if r[3] is not None]
ce = [r[4] for r in rows if r[4] is not None]
sh = [r[2] / r[4] for r in rows if r[2] is not None and r[4]]
jj = [r[5] for r in rows if r[5] is not None]
jr = [r[6] for r in rows if r[6] is not None]
p(f'\n1. Which residues are central: within-protein agreement between the force fields (medians over proteins)')
p(f'   time-averaged betweenness r = {st.median(rb):.2f} (IQR {q(rb, .25):.2f}-{q(rb, .75):.2f}); noise ceiling '
  f'{st.median(ce):.2f}; share of ceiling {st.median(sh):.2f}')
p(f'   hub occupancy r = {st.median(ro):.2f} (IQR {q(ro, .25):.2f}-{q(ro, .75):.2f})')
p(f'   overlap of the two hub sets (Jaccard) {st.median(jj):.3f} (IQR {q(jj, .25):.3f}-{q(jj, .75):.3f}); '
  f'random sets {st.median(jr):.3f}')
p('   by group (median r of betweenness / Jaccard):  ' + '   '.join(
    f'{g} {st.median([r[2] for r in rows if r[7] == g and r[2] is not None]):.2f}/'
    f'{st.median([r[5] for r in rows if r[7] == g and r[5] is not None]):.3f}'
    for g in sorted({r[7] for r in rows})))

p('\n2. By residue type: hub occupancy (mean; null 0.10) and long-range share of contacts in hub frames (median)')
p('   aa     occupancy CAL / MPIPI      LR share in hub frames CAL / MPIPI     n residues')
occ = {}
for aa in AAS:
    va = [r for u in ids for r in rc[u].values() if r['aa'] == aa]
    vb = [r for u in ids for r in rm[u].values() if r['aa'] == aa]
    if not va or not vb:
        continue
    oa = sum(float(r['occupancy']) for r in va) / len(va)
    ob = sum(float(r['occupancy']) for r in vb) / len(vb)
    la = st.median([float(r['lr_share_hubframes']) for r in va if r['lr_share_hubframes'] and int(r['n_hubframes']) >= 5])
    lb = st.median([float(r['lr_share_hubframes']) for r in vb if r['lr_share_hubframes'] and int(r['n_hubframes']) >= 5])
    occ[aa] = (oa, ob)
    p(f'   {aa}      {oa:.3f} / {ob:.3f}                 {la:.2f} / {lb:.2f}                        {len(va)}')
a = [occ[x][0] for x in occ]
b = [occ[x][1] for x in occ]
p(f'   Pearson r across residue types (occupancy): {pearson(a, b):+.2f}')
open(f'{MPI}/crossff_summary.txt', 'w').write('\n'.join(out) + '\n')
