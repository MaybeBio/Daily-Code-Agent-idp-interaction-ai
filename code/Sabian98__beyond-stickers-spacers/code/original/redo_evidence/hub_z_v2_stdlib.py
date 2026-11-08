"""
hub_z_v2_stdlib.py — per-amino-acid, per-group betweenness (and degree) z-scores on the CORRECTED
contact networks (Analysis_300K_v2_graph/npz_residue: Calpha-Calpha < 0.8 nm, |i-j| >= 3, production half,
50 graph frames, k = 30 pivot betweenness), using only the Python standard library (no numpy on this machine).

Method as in the paper / Scripts/hub_zscores.py: per protein, z = (x_i - mean) / sd over the protein's residues
(population sd), pooled by amino acid and group (residue-weighted mean, SEM over residues).
Variants:
  all      all residues (paper's method)
  noends   first and last residue dropped before z-scoring (initiator-Met / chain-end artifact check)
  posadj   z replaced by its residual from a within-protein quadratic in normalised position (as arg_lys_hub_effect.py)
Also per protein: mean z per amino acid (for protein-level fractions and the paired R - K contrast).

Read-only on the NAS; writes to the scratchpad only.
"""
import array
import ast
import csv
import math
import os
import sys
import zipfile
from collections import defaultdict
from multiprocessing import Pool

BIO = '/Volumes/User Homes/Nu Project/BENDER-BIO/BENDER_BIO.csv'
NPZ = '/Volumes/User Homes/Nu Project/BENDER-MT/Traj Level Analysis/Analysis_300K_v2_graph/npz_residue'
OUT = os.path.dirname(os.path.abspath(__file__))  # tables written next to this script
AAS = 'ACDEFGHIKLMNPQRSTVWY'
METRICS = {'betweenness': 'betweenness_centrality_mean', 'degree': 'degree_centrality_mean'}
VARIANTS = ('all', 'noends', 'posadj')


def read_npy(b):
    assert b[:6] == b'\x93NUMPY'
    major = b[6]
    hl = int.from_bytes(b[8:10] if major == 1 else b[8:12], 'little')
    off = 10 if major == 1 else 12
    h = ast.literal_eval(b[off:off + hl].decode('latin1'))
    data = b[off + hl:]
    d = h['descr']
    if d == '<f8':
        a = array.array('d')
        a.frombytes(data)
        return list(a)
    if d.startswith('<U'):
        w = int(d[2:]) * 4
        return [data[i:i + w].decode('utf-32-le').rstrip('\x00') for i in range(0, len(data), w)]
    raise ValueError(d)


def zs(v):
    n = len(v)
    m = sum(v) / n
    sd = math.sqrt(sum((x - m) ** 2 for x in v) / n)
    if sd < 1e-12:
        return None
    return [(x - m) / sd for x in v]


def quad_resid(z):
    n = len(z)
    xs = [-1 + 2 * i / (n - 1) for i in range(n)]
    # least squares z ~ a + b x + c x^2 via normal equations
    S = [[0.0] * 3 for _ in range(3)]
    t = [0.0] * 3
    for x, y in zip(xs, z):
        p = (1.0, x, x * x)
        for r in range(3):
            t[r] += p[r] * y
            for c in range(3):
                S[r][c] += p[r] * p[c]
    # solve 3x3 by Gaussian elimination
    A = [S[r][:] + [t[r]] for r in range(3)]
    for c in range(3):
        piv = max(range(c, 3), key=lambda r: abs(A[r][c]))
        A[c], A[piv] = A[piv], A[c]
        if abs(A[c][c]) < 1e-15:
            return None
        for r in range(3):
            if r != c:
                f = A[r][c] / A[c][c]
                for k in range(c, 4):
                    A[r][k] -= f * A[c][k]
    coef = [A[r][3] / A[r][r] for r in range(3)]
    return [y - (coef[0] + coef[1] * x + coef[2] * x * x) for x, y in zip(xs, z)]


def one(args):
    uid, group = args
    path = f'{NPZ}/{uid}_residue_features.npz'
    try:
        zf = zipfile.ZipFile(path)
        seq = read_npy(zf.read('sequence.npy'))
        vals = {m: read_npy(zf.read(f'{key}.npy')) for m, key in METRICS.items()}
    except Exception as e:  # missing or unreadable
        return (uid, group, None, repr(e))
    n = min(len(seq), *(len(v) for v in vals.values()))
    if n < 8:
        return (uid, group, None, 'short')
    seq = seq[:n]
    out = {}
    for m, v in vals.items():
        v = v[:n]
        res = {}
        z = zs(v)
        if z is None:
            continue
        res['all'] = (seq, z)
        z2 = zs(v[1:-1])
        if z2 is not None:
            res['noends'] = (seq[1:-1], z2)
        zr = quad_resid(z)
        if zr is not None:
            res['posadj'] = (seq, zr)
        out[m] = res
    # compact: per variant/metric, per aa: (sum, sumsq, count); first residue identity and its z
    comp = {}
    for m, res in out.items():
        for var, (s, z) in res.items():
            acc = {}
            for aa, x in zip(s, z):
                if aa not in AAS:
                    continue
                a = acc.setdefault(aa, [0.0, 0.0, 0])
                a[0] += x
                a[1] += x * x
                a[2] += 1
            comp[(m, var)] = acc
    first = (seq[0], out.get('betweenness', {}).get('all', (None, [None]))[1][0])
    return (uid, group, comp, first)


def main():
    rows = list(csv.DictReader(open(BIO)))
    use_all = '--all' in sys.argv
    jobs = [(r['UniProt_ID'], r['kingdom']) for r in rows if use_all or r['analysis_ok'] == 'True']
    os.makedirs(OUT, exist_ok=True)
    pooled = defaultdict(lambda: [0.0, 0.0, 0])      # (metric, var, group, aa) -> sum, sumsq, n_res
    prot = defaultdict(list)                          # (metric, var, group, aa) -> per-protein mean z
    rk = defaultdict(list)                            # (metric, var, group) -> per-protein mean z(R) - mean z(K)
    first_res = defaultdict(list)                     # group -> (aa, z) of residue 0 (betweenness, all)
    failed = []
    done = 0
    with Pool(int(os.environ.get('WORKERS', '12'))) as pool:
        for uid, group, comp, info in pool.imap_unordered(one, jobs, chunksize=16):
            done += 1
            if done % 500 == 0:
                print(f'{done}/{len(jobs)}', flush=True)
            if comp is None:
                failed.append((uid, group, info))
                continue
            first_res[group].append(info)
            for (m, var), acc in comp.items():
                for g in (group, 'All'):
                    for aa, (s, ss, c) in acc.items():
                        p = pooled[(m, var, g, aa)]
                        p[0] += s
                        p[1] += ss
                        p[2] += c
                        prot[(m, var, g, aa)].append(s / c)
                    if 'R' in acc and 'K' in acc:
                        rk[(m, var, g)].append(acc['R'][0] / acc['R'][2] - acc['K'][0] / acc['K'][2])

    with open(f'{OUT}/hub_z_by_group.csv', 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['metric', 'variant', 'group', 'aa', 'mean_z', 'sem', 'n_residues', 'n_proteins',
                    'frac_proteins_z_pos', 'protein_mean_z'])
        for (m, var, g, aa), (s, ss, c) in sorted(pooled.items()):
            mean = s / c
            sd = math.sqrt(max(ss / c - mean * mean, 0.0))
            pm = prot[(m, var, g, aa)]
            w.writerow([m, var, g, aa, f'{mean:.5f}', f'{sd / math.sqrt(c):.5f}', c, len(pm),
                        f'{sum(x > 0 for x in pm) / len(pm):.4f}', f'{sum(pm) / len(pm):.5f}'])

    with open(f'{OUT}/rk_paired.csv', 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['metric', 'variant', 'group', 'n_proteins', 'mean_diff', 'sd_diff', 'cohens_dz', 'frac_R_gt_K'])
        for (m, var, g), d in sorted(rk.items()):
            n = len(d)
            mu = sum(d) / n
            sd = math.sqrt(sum((x - mu) ** 2 for x in d) / (n - 1))
            w.writerow([m, var, g, n, f'{mu:.5f}', f'{sd:.5f}', f'{mu / sd:.4f}', f'{sum(x > 0 for x in d) / n:.4f}'])

    with open(f'{OUT}/first_residue.csv', 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['group', 'n', 'frac_M_first', 'mean_z_first_residue_betweenness'])
        for g, L in sorted(first_res.items()):
            zz = [z for _, z in L if z is not None]
            w.writerow([g, len(L), f'{sum(a == "M" for a, _ in L) / len(L):.4f}', f'{sum(zz) / len(zz):.4f}'])

    with open(f'{OUT}/failed.csv', 'w', newline='') as fh:
        csv.writer(fh).writerows([('UniProt_ID', 'group', 'reason')] + failed)
    print(f'done: {done} proteins, {len(failed)} failed -> {OUT}')


if __name__ == '__main__':
    main()
