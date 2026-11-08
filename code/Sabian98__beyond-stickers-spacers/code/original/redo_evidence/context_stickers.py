"""
context_stickers.py — in what sequence context does an aromatic residue act as a network hub beyond geometry?

Target: per aromatic residue (F, W, Y), excess = simulation z - residue-neutral null z (sim_null_per_residue_all.csv;
position-adjusted betweenness z from 1,000 frames vs an Rg-matched self-avoiding null with the force field's bead
sizes). Fixed 0.8 nm rule (near-neutral) is primary; the bead-size rule is reported too.
Features of each aromatic residue (sequence only):
  local (excluding the residue), windows ±3 and ±7: fractions G, S, P, N+Q, D+E, K+R, R, other aromatics F/W/Y,
      aliphatic I/L/V/M, local net charge (K+R - D-E), mean CALVADOS-2 lambda
  spacing: log sequence distance to the nearest other aromatic; number of other aromatics within ±15
  protein: FCR, NCPR, aromatic fraction, P, G, S fractions, log length, SCD (charge patterning)
  position: relative and log absolute distance to the nearer chain end; identity W, Y (F = reference)
Models (per force field): univariate Spearman of each feature with the excess; ridge regression on standardised
features, 5-fold cross-validated by protein (CV R^2). Context effects shared by the three force fields are then
tested on experiment: the per-residue excess of aromatics in the PRE ensembles over their own null
(ped_null_per_residue.csv) vs the excess predicted from sequence context by each force field's model.
Standard library only.  python3 context_stickers.py  ->  context_stickers_summary.txt
"""
import csv
import math
import os
import random
import sys
import zlib
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import hub_frames_tests as ht  # noqa: E402

SIM = os.path.join(HERE, 'sim_null_per_residue_all.csv')
EXP = os.path.join(HERE, 'ped_null_per_residue.csv')
AROM = set('FWY')
PRIMARY = ['PED00024', 'PED00003', 'PED00006', 'PED00157', 'PED00160', 'PED00161', 'PED00225', 'PED00498',
           'PED00499', 'PED00500', 'PED00511']
LAM = {r['one']: float(r['lambdas']) for r in csv.DictReader(open(ht.CAL_PARAMS))}
LAM.update({'x': LAM['S'], 'y': LAM['T'], 'z': LAM['Y']})
ACID, BASE = set('DExyz'), set('KR')          # phosphoresidues (x, y, z) count as acidic
out = []


def p(s=''):
    out.append(s)
    print(s, flush=True)


def scd(seq):
    q = [1 if a in BASE else -1 if a in ACID else 0 for a in seq]
    idx = [i for i, v in enumerate(q) if v]
    n = len(seq)
    return sum(q[i] * q[j] * math.sqrt(j - i) for a, i in enumerate(idx) for j in idx[a + 1:]) / n if n else 0.0


def features(seq, i, prot_cache):
    n = len(seq)
    f = {}
    for w in (3, 7):
        win = [seq[j] for j in range(max(0, i - w), min(n, i + w + 1)) if j != i]
        L = max(1, len(win))
        f[f'G_w{w}'] = sum(a == 'G' for a in win) / L
        f[f'S_w{w}'] = sum(a == 'S' for a in win) / L
        f[f'P_w{w}'] = sum(a == 'P' for a in win) / L
        f[f'NQ_w{w}'] = sum(a in 'NQ' for a in win) / L
        f[f'acid_w{w}'] = sum(a in ACID for a in win) / L
        f[f'basic_w{w}'] = sum(a in BASE for a in win) / L
        f[f'R_w{w}'] = sum(a == 'R' for a in win) / L
        f[f'aro_w{w}'] = sum(a in AROM for a in win) / L
        f[f'aliph_w{w}'] = sum(a in 'ILVM' for a in win) / L
        f[f'ncpr_w{w}'] = (sum(a in BASE for a in win) - sum(a in ACID for a in win)) / L
        f[f'lam_w{w}'] = sum(LAM.get(a, 0.5) for a in win) / L
    others = [j for j, a in enumerate(seq) if a in AROM and j != i]
    f['log_dist_aro'] = math.log(1 + min([abs(j - i) for j in others] or [60]))
    f['n_aro_w15'] = sum(1 for j in others if abs(j - i) <= 15)
    if 'prot' not in prot_cache:
        fcr = sum(a in ACID or a in BASE for a in seq) / n
        prot_cache['prot'] = {'FCR': fcr, 'NCPR': (sum(a in BASE for a in seq) - sum(a in ACID for a in seq)) / n,
                              'f_aro': sum(a in AROM for a in seq) / n, 'f_P': seq.count('P') / n,
                              'f_G': seq.count('G') / n, 'f_S': seq.count('S') / n, 'logN': math.log(n),
                              'SCD': scd(seq)}
    f.update(prot_cache['prot'])
    f['end_rel'] = min(i, n - 1 - i) / n
    f['end_log'] = math.log(1 + min(i, n - 1 - i))
    f['is_W'] = 1.0 if seq[i] == 'W' else 0.0
    f['is_Y'] = 1.0 if seq[i] == 'Y' else 0.0
    return f


def ranks(v):
    o = sorted(range(len(v)), key=lambda i: v[i])
    r = [0.0] * len(v)
    i = 0
    while i < len(o):
        j = i
        while j + 1 < len(o) and v[o[j + 1]] == v[o[i]]:
            j += 1
        for k in range(i, j + 1):
            r[o[k]] = (i + j) / 2.0
        i = j + 1
    return r


def pearson(a, b):
    n = len(a)
    ma, mb = sum(a) / n, sum(b) / n
    sa = math.sqrt(sum((x - ma) ** 2 for x in a))
    sb = math.sqrt(sum((y - mb) ** 2 for y in b))
    return sum((x - ma) * (y - mb) for x, y in zip(a, b)) / (sa * sb) if sa and sb else float('nan')


def spearman(a, b):
    return pearson(ranks(a), ranks(b))


def solve(A, b):
    n = len(b)
    M = [A[r][:] + [b[r]] for r in range(n)]
    for c in range(n):
        piv = max(range(c, n), key=lambda r: abs(M[r][c]))
        M[c], M[piv] = M[piv], M[c]
        for r in range(n):
            if r != c and M[c][c]:
                f = M[r][c] / M[c][c]
                for k in range(c, n + 1):
                    M[r][k] -= f * M[c][k]
    return [M[r][n] / M[r][r] if M[r][r] else 0.0 for r in range(n)]


def ridge_fit(X, y, alpha):
    k = len(X[0])
    my = sum(y) / len(y)
    A = [[0.0] * k for _ in range(k)]
    b = [0.0] * k
    for row, yi in zip(X, y):
        for a in range(k):
            ra = row[a]
            b[a] += ra * (yi - my)
            Aa = A[a]
            for c in range(a, k):
                Aa[c] += ra * row[c]
    for a in range(k):
        for c in range(a):
            A[a][c] = A[c][a]
        A[a][a] += alpha
    return solve(A, b), my


def standardise(rows, keys, stats=None):
    if stats is None:
        stats = {}
        for k in keys:
            v = [r[k] for r in rows]
            m = sum(v) / len(v)
            sd = math.sqrt(sum((x - m) ** 2 for x in v) / len(v)) or 1.0
            stats[k] = (m, sd)
    return [[(r[k] - stats[k][0]) / stats[k][1] for k in keys] for r in rows], stats


def cv_r2(rows, y, groups, keys, alpha, folds=5):
    pred = [0.0] * len(y)
    fold = [zlib.crc32(g.encode()) % folds for g in groups]
    for f in range(folds):
        tr = [i for i in range(len(y)) if fold[i] != f]
        te = [i for i in range(len(y)) if fold[i] == f]
        Xtr, st = standardise([rows[i] for i in tr], keys)
        beta, my = ridge_fit(Xtr, [y[i] for i in tr], alpha)
        Xte, _ = standardise([rows[i] for i in te], keys, st)
        for i, row in zip(te, Xte):
            pred[i] = my + sum(bb * x for bb, x in zip(beta, row))
    ssr = sum((a - b) ** 2 for a, b in zip(y, pred))
    my = sum(y) / len(y)
    sst = sum((a - my) ** 2 for a in y)
    return 1 - ssr / sst


def load_sim():
    seqs = defaultdict(dict)
    vals = defaultdict(dict)
    for r in csv.DictReader(open(SIM)):
        key = (r['force_field'], r['UniProt_ID'])
        i = int(r['pos'])
        seqs[key][i] = r['aa']
        if r['z_null_fixed'] and r['z_sim_fixed']:
            vals[key][i] = (float(r['z_sim_fixed']) - float(r['z_null_fixed']),
                            float(r['z_sim_size']) - float(r['z_null_size']) if r['z_sim_size'] and r['z_null_size'] else None)
    data = defaultdict(list)
    for key, d in seqs.items():
        seq = ''.join(d[i] for i in range(len(d)))
        cache = {}
        for i, a in enumerate(seq):
            if a in AROM and i in vals[key]:
                data[key[0]].append((key[1], features(seq, i, cache), vals[key][i]))
    return data


def load_exp():
    if not os.path.exists(EXP):
        return []
    by = defaultdict(dict)
    meta = {}
    for r in csv.DictReader(open(EXP)):
        by[r['entry']][int(r['pos'])] = (r['aa'], float(r['z_exp_fixed']) - float(r['z_null_fixed']),
                                         float(r['z_exp_size']) - float(r['z_null_size']))
        meta[r['entry']] = r['protein']
    rows = []
    for e, d in by.items():
        seq = ''.join(d[i][0] for i in range(len(d)))
        cache = {}
        for i, (a, ef, es) in d.items():
            if a in AROM:
                rows.append((e, meta[e], features(seq, i, cache), ef, es))
    return rows


def main():
    data = load_sim()
    exp = load_exp()
    keys = sorted(next(iter(data.values()))[0][1].keys())
    p('Context of aromatic hubs: excess centrality over a residue-neutral null (position-adjusted z; simulation - null)')
    models = {}
    for rule_i, rule in enumerate(('fixed', 'size')):
        p(f'\n=== Rule: {rule} ===')
        uni = {}
        for ff, rows in sorted(data.items()):
            rr = [(g, f, v[rule_i]) for g, f, v in rows if v[rule_i] is not None]
            y = [v for g, f, v in rr]
            feats = [f for g, f, v in rr]
            groups = [g for g, f, v in rr]
            p(f'\n  {ff}: {len(rr)} aromatic residues in {len(set(groups))} proteins; mean excess {sum(y) / len(y):+.2f}, '
              f'share > 0 {sum(v > 0 for v in y) / len(y):.0%}')
            uni[ff] = {k: spearman([f[k] for f in feats], y) for k in keys}
            best = (None, -1e9)
            for alpha in (1.0, 10.0, 100.0, 1000.0):
                r2 = cv_r2(feats, y, groups, keys, alpha)
                if r2 > best[1]:
                    best = (alpha, r2)
            X, st = standardise(feats, keys)
            beta, my = ridge_fit(X, y, best[0])
            models[(ff, rule)] = (beta, my, st)
            p(f'    ridge (alpha {best[0]:g}) cross-validated R^2 by protein = {best[1]:.3f}')
        p('\n  Univariate Spearman of each feature with the aromatic excess (* = same sign and |rho| >= 0.10 in all three force fields)')
        p('  ' + f'{"feature":14s}' + ''.join(f'{ff:>18s}' for ff in sorted(data)))
        order = sorted(keys, key=lambda k: -abs(sum(uni[ff][k] for ff in data) / len(data)))
        for k in order:
            v = [uni[ff][k] for ff in sorted(data)]
            star = '*' if (all(x >= 0.10 for x in v) or all(x <= -0.10 for x in v)) else ' '
            p(f'  {k:14s}' + ''.join(f'{x:+18.2f}' for x in v) + f'  {star}')
        # experimental test
        p(f'\n  Experiment ({rule} rule): aromatic residues of the PRE ensembles (one ensemble per protein), observed excess over '
          f'their own null vs excess predicted from sequence context')
        er = [(e, prot, f, (ef if rule == 'fixed' else es)) for e, prot, f, ef, es in exp if e in PRIMARY]
        if not er:
            continue
        obs = [o for e, prot, f, o in er]
        line = f'    {len(er)} aromatic residues in {len(set(e for e, *_ in er))} proteins. Spearman(observed, predicted):'
        preds = {}
        for ff in sorted(data):
            beta, my, st = models[(ff, rule)]
            Xe, _ = standardise([f for e, prot, f, o in er], keys, st)
            preds[ff] = [my + sum(b * x for b, x in zip(beta, row)) for row in Xe]
            line += f'  {ff} {spearman(preds[ff], obs):+.2f}'
        cons = [sum(preds[ff][i] for ff in preds) / len(preds) for i in range(len(obs))]
        line += f'  | consensus {spearman(cons, obs):+.2f}'
        p(line)
        # permutation p for the consensus correlation, shuffling observed values within the whole set
        rng = random.Random(0)
        r0 = spearman(cons, obs)
        ge = sum(1 for _ in range(2000) if spearman(cons, rng.sample(obs, len(obs))) >= r0)
        p(f'    consensus permutation p (one-sided) = {(ge + 1) / 2001:.3f}')
        byp = defaultdict(lambda: [[], []])
        for (e, prot, f, o), c in zip(er, cons):
            byp[prot][0].append(o)
            byp[prot][1].append(c)
        pm = [(prot, sum(v[0]) / len(v[0]), sum(v[1]) / len(v[1]), len(v[0])) for prot, v in byp.items()]
        p('    per protein (observed mean excess / predicted mean excess, n aromatics): ' + '  '.join(
            f'{prot} {o:+.2f}/{c:+.2f} ({n})' for prot, o, c, n in sorted(pm, key=lambda t: -t[1])))
        p(f'    protein-level Spearman(observed, predicted) = {spearman([t[1] for t in pm], [t[2] for t in pm]):+.2f} (n = {len(pm)})')
    open(os.path.join(HERE, 'context_stickers_summary.txt'), 'w').write('\n'.join(out) + '\n')


if __name__ == '__main__':
    main()
