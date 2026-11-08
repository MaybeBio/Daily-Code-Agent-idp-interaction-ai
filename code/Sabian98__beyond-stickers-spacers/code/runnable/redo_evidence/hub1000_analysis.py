"""
hub1000_analysis.py — the residue-level hub analyses redone on 1,000-frame exact betweenness
(hub_frames*/per_residue.csv) instead of the stored 50-frame, k = 30-pivot values (npz_residue), which agree with
the 1,000-frame values at r = 0.56 only.

Data sets:
  cal     hub_frames/                   CALVADOS-2, 1,198 proteins (200 per group, <= 300 residues)
  calp    hub_frames_paired_calvados/   CALVADOS-2, the 908 sequences also simulated with MPIPI-GG
  mpipi   hub_frames_mpipi/             MPIPI-GG, the same 908 sequences
For each: (1) residue-type x group betweenness z (within-protein z, raw and position-adjusted), (2) R vs K paired
contrast, (3) copy-number test (fraction of copies in the protein's top 10%), and for cal (4) the same table from the
stored 50-frame values on the same proteins and (5) the CALVADOS-2 lambda/sigma regression; across force fields
(6) agreement of the residue-type profiles per group (the paper's "per-taxon r > 0.84" claim).
Standard library only.  python3 hub1000_analysis.py  ->  prints and writes hub1000_summary.txt + hub1000_by_group.csv
"""
import csv
import math
import os
import sys
import zipfile
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from hub_z_v2_stdlib import read_npy, zs, quad_resid  # noqa: E402

SETS = {'cal': 'hub_frames', 'calp': 'hub_frames_paired_calvados', 'mpipi': 'hub_frames_mpipi'}
NPZ = (__import__('os').environ['BSS_DATA'] + '/BENDER-MT/Traj Level Analysis/Analysis_300K_v2_graph/npz_residue')
PARAMS = (__import__('os').environ['BSS_DATA'] + '/Scripts/CALVADOS-Slab/input/residues_CALVADOS2.csv')
AAS = 'ACDEFGHIKLMNPQRSTVWY'
SHOW = 'GSACYFWNTPQHDRELIVMK'
GROUPS = ['Bacteria', 'Fungi', 'Mammals', 'Plants', 'Protists', 'Viruses']
BINS = [(1, 1), (2, 2), (3, 4), (5, 8), (9, 16), (17, 10 ** 9)]


def binlabel(lo, hi):
    return str(lo) if lo == hi else (f'{lo}-{hi}' if hi < 10 ** 9 else f'{lo}+')
out = []


def p(s=''):
    out.append(s)
    print(s, flush=True)


def pearson(a, b):
    n = len(a)
    ma, mb = sum(a) / n, sum(b) / n
    sa = math.sqrt(sum((x - ma) ** 2 for x in a))
    sb = math.sqrt(sum((y - mb) ** 2 for y in b))
    return sum((x - ma) * (y - mb) for x, y in zip(a, b)) / (sa * sb) if sa and sb else float('nan')


def proteins(dirname):
    """UniProt_ID -> (group, seq, mean betweenness list) from a hub_frames output folder."""
    prot = {}
    for r in csv.DictReader(open(os.path.join(HERE, dirname, 'per_residue.csv'))):
        d = prot.setdefault(r['UniProt_ID'], [r['group'], {}, {}])
        i = int(r['pos'])
        d[1][i] = r['aa']
        d[2][i] = float(r['mean_btw'])
    return {u: (g, ''.join(s[i] for i in range(len(s))), [b[i] for i in range(len(b))]) for u, (g, s, b) in prot.items()}


def per_protein_z(prot):
    """UniProt_ID -> (group, seq, raw z, position-adjusted z)."""
    res = {}
    for u, (g, seq, b) in prot.items():
        z = zs(b)
        if z is None:
            continue
        zr = quad_resid(z)
        if zr is None:
            continue
        res[u] = (g, seq, z, zr)
    return res


def type_table(pz, which):
    """(group, aa) -> (mean z, sem, n residues); which = 2 raw, 3 position-adjusted."""
    acc = defaultdict(lambda: [0.0, 0.0, 0])
    for u, v in pz.items():
        g, seq = v[0], v[1]
        for aa, x in zip(seq, v[which]):
            for gg in (g, 'All'):
                a = acc[(gg, aa)]
                a[0] += x
                a[1] += x * x
                a[2] += 1
    tab = {}
    for k, (s, ss, c) in acc.items():
        m = s / c
        tab[k] = (m, math.sqrt(max(ss / c - m * m, 0.0)) / math.sqrt(c), c)
    return tab


def show_table(tab, title):
    p(f'\n{title}  (mean z; * = |z| > 3 SEM)')
    p('AA ' + ''.join(f'{g[:8]:>10s}' for g in GROUPS + ['All']))
    for aa in SHOW:
        line = f'{aa}  '
        for g in GROUPS + ['All']:
            m = tab.get((g, aa))
            line += f'{m[0]:+9.3f}' + ('*' if abs(m[0]) > 3 * m[1] else ' ') if m else f'{"-":>10s}'
        p(line)


def rk(pz, which):
    by = defaultdict(list)
    for u, v in pz.items():
        g, seq, z = v[0], v[1], v[which]
        zr_ = [x for a, x in zip(seq, z) if a == 'R']
        zk_ = [x for a, x in zip(seq, z) if a == 'K']
        if zr_ and zk_:
            d = sum(zr_) / len(zr_) - sum(zk_) / len(zk_)
            by[g].append(d)
            by['All'].append(d)
    rows = []
    for g in GROUPS + ['All']:
        d = by.get(g, [])
        if len(d) > 2:
            m = sum(d) / len(d)
            sd = math.sqrt(sum((x - m) ** 2 for x in d) / (len(d) - 1))
            rows.append((g, len(d), m, m / sd, sum(x > 0 for x in d) / len(d)))
    return rows


def copy_number(pz):
    acc = defaultdict(lambda: [0, 0])     # (aa, bin) -> copies in top 10%, copies
    for u, v in pz.items():
        seq, z = v[1], v[3]
        n = len(seq)
        n_top = max(1, round(0.1 * n))
        top = set(sorted(range(n), key=lambda i: -z[i])[:n_top])
        cnt = defaultdict(list)
        for i, a in enumerate(seq):
            cnt[a].append(i)
        for a, idx in cnt.items():
            k = len(idx)
            b = next(binlabel(lo, hi) for lo, hi in BINS if lo <= k <= hi)
            c = acc[(a, b)]
            c[0] += sum(1 for i in idx if i in top)
            c[1] += k
    return acc


def stored_50frame(pz):
    """The same proteins' stored v2 values (50 frames, k = 30 pivots), same z-scoring."""
    prot = {}
    for u, v in pz.items():
        try:
            zf = zipfile.ZipFile(f'{NPZ}/{u}_residue_features.npz')
            b = read_npy(zf.read('betweenness_centrality_mean.npy'))
            seq = ''.join(read_npy(zf.read('sequence.npy')))
        except Exception:
            continue
        n = min(len(b), len(seq))
        if seq[:n] == v[1][:n] and n == len(v[1]):
            prot[u] = (v[0], seq[:n], b[:n])
    return per_protein_z(prot)


def main():
    data = {k: per_protein_z(proteins(d)) for k, d in SETS.items()}
    p('Residue-level hub analyses on 1,000-frame exact betweenness (CALVADOS-2 and MPIPI-GG, v2 contact rule)')
    for k, pz in data.items():
        p(f'  {k:6s} {len(pz)} proteins, {sum(len(v[1]) for v in pz.values())} residues')
    tabs = {}
    rows_out = []
    for k in SETS:
        for which, name in ((2, 'raw'), (3, 'posadj')):
            tabs[(k, name)] = type_table(data[k], which)
            for (g, aa), (m, s, c) in sorted(tabs[(k, name)].items()):
                rows_out.append([k, name, g, aa, f'{m:.5f}', f'{s:.5f}', c])
    with open(os.path.join(HERE, 'hub1000_by_group.csv'), 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['set', 'variant', 'group', 'aa', 'mean_z', 'sem', 'n_residues'])
        w.writerows(rows_out)

    show_table(tabs[('cal', 'raw')], '1a. CALVADOS-2 (1,198 proteins), 1,000 frames, within-protein z (paper method)')
    show_table(tabs[('cal', 'posadj')], '1b. CALVADOS-2 (1,198 proteins), 1,000 frames, position-adjusted z')
    show_table(tabs[('mpipi', 'posadj')], '1c. MPIPI-GG (908 proteins), 1,000 frames, position-adjusted z')

    s50 = stored_50frame(data['cal'])
    t50 = type_table(s50, 3)
    t1k = type_table({u: data['cal'][u] for u in s50}, 3)
    a = [t50[('All', aa)][0] for aa in AAS]
    b = [t1k[('All', aa)][0] for aa in AAS]
    p(f'\n2. Stored 50-frame values vs 1,000 frames on the same {len(s50)} CALVADOS-2 proteins (position-adjusted z, all groups)')
    p('   ' + '  '.join(f'{aa} {x:+.2f}/{y:+.2f}' for aa, x, y in zip(AAS, a, b)))
    p(f'   residue-type profile agreement r = {pearson(a, b):+.3f}; mean |difference| {sum(abs(x - y) for x, y in zip(a, b)) / 20:.3f}')

    p('\n3. R vs K, paired within protein (mean z(R) - mean z(K); d_z = mean / sd; fraction of proteins with R > K)')
    for k in SETS:
        for which, name in ((2, 'raw'), (3, 'pos-adj')):
            rr = rk(data[k], which)
            allr = [r for r in rr if r[0] == 'All'][0]
            p(f'   {k:6s} {name:8s} all: n {allr[1]}, diff {allr[2]:+.3f}, d_z {allr[3]:+.2f}, R > K {allr[4]:.0%};  by group d_z '
              + ' '.join(f'{r[0][:4]} {r[3]:+.2f}' for r in rr if r[0] != 'All'))

    p('\n4. Copy-number test: fraction of a type\'s copies in the protein\'s top 10% (position-adjusted; null 10%)')
    for k in ('cal', 'mpipi'):
        acc = copy_number(data[k])
        p(f'   {k}')
        p('   AA ' + ''.join(f'{binlabel(lo, hi):>8s}' for lo, hi in BINS))
        for aa in 'EDRKYFWSGP':
            line = f'   {aa}  '
            for lo, hi in BINS:
                b_ = binlabel(lo, hi)
                t, c = acc.get((aa, b_), (0, 0))
                line += f'{t / c:8.3f}' if c >= 30 else f'{"-":>8s}'
            p(line)

    prm = {r['one']: r for r in csv.DictReader(open(PARAMS))}
    lam = [float(prm[aa]['lambdas']) for aa in AAS]
    sig = [float(prm[aa]['sigmas']) for aa in AAS]
    z = [tabs[('cal', 'posadj')][('All', aa)][0] for aa in AAS]
    p(f'\n5. CALVADOS-2 parameters vs the 1,000-frame residue profile: r(z, lambda) {pearson(z, lam):+.2f}, '
      f'r(z, sigma) {pearson(z, sig):+.2f}, R^2(lambda + sigma) {r2_two(z, lam, sig):.2f}')

    p('\n6. Residue-type profiles, CALVADOS-2 vs MPIPI-GG on the same 908 sequences (Pearson r across 20 types; '
      'the paper claims r > 0.84 in every taxon)')
    for name in ('raw', 'posadj'):
        tc, tm = tabs[('calp', name)], tabs[('mpipi', name)]
        line = f'   {name:7s}'
        for g in GROUPS + ['All']:
            ks = [aa for aa in AAS if (g, aa) in tc and (g, aa) in tm and tc[(g, aa)][2] >= 20 and tm[(g, aa)][2] >= 20]
            r_ = pearson([tc[(g, aa)][0] for aa in ks], [tm[(g, aa)][0] for aa in ks])
            ks2 = [aa for aa in ks if aa != 'P']
            r2_ = pearson([tc[(g, aa)][0] for aa in ks2], [tm[(g, aa)][0] for aa in ks2])
            line += f'  {g[:4]} {r_:+.2f} ({r2_:+.2f} w/o P)'
        p(line)
    tc, tm = tabs[('calp', 'posadj')], tabs[('mpipi', 'posadj')]
    p('   all groups, position-adjusted z CAL / MPIPI: ' + '  '.join(
        f'{aa} {tc[("All", aa)][0]:+.2f}/{tm[("All", aa)][0]:+.2f}' for aa in SHOW))
    open(os.path.join(HERE, 'hub1000_summary.txt'), 'w').write('\n'.join(out) + '\n')


def r2_two(y, x1, x2):
    n = len(y)
    X = [[1.0, a, b] for a, b in zip(x1, x2)]
    S = [[sum(X[i][r] * X[i][c] for i in range(n)) for c in range(3)] for r in range(3)]
    t = [sum(X[i][r] * y[i] for i in range(n)) for r in range(3)]
    A = [S[r][:] + [t[r]] for r in range(3)]
    for c in range(3):
        piv = max(range(c, 3), key=lambda r: abs(A[r][c]))
        A[c], A[piv] = A[piv], A[c]
        for r in range(3):
            if r != c:
                f = A[r][c] / A[c][c]
                for k in range(c, 4):
                    A[r][k] -= f * A[c][k]
    b = [A[r][3] / A[r][r] for r in range(3)]
    yh = [b[0] + b[1] * a + b[2] * c for a, c in zip(x1, x2)]
    my = sum(y) / n
    return 1 - sum((u - v) ** 2 for u, v in zip(y, yh)) / sum((u - my) ** 2 for u in y)


if __name__ == '__main__':
    main()
