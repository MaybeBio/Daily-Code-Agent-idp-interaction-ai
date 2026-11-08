"""
ped_null.py — residue-identity-neutral baseline for the PED ensembles: are aromatic residues more central (and acidic
residues less central) in PRE-restrained ensembles than geometry alone predicts?

Null: for each ensemble's exact sequence, self-avoiding Calpha chains (3.8 Å bonds, virtual bond angle uniform in
85-150°, 100-150° or 115-155° (a third of the pool each, to cover compact and expanded ensembles), random dihedrals) whose only residue-specific property is excluded volume (hard core 0.8 * sigma_ij from
CALVADOS-2 bead sizes; phosphoresidues use S/T sizes). No attraction of any kind. N_POOL chains are grown (with local
back-tracking); then, for every experimental model, the pool chain with the nearest radius of gyration is taken
(Rg quantile matching), so the null has the experimental ensemble's compaction.
Both ensembles are analysed identically (contact graphs with the fixed 0.8 nm rule and the bead-size rule, exact
betweenness, mean over models, position-adjusted z). Per class: excess = experimental mean z - null mean z.
Validation: the same comparison against PED's IDPConformerGenerator initial pools (α-syn, Sic1, drkN).
Standard library only.  python3 ped_null.py  ->  ped_null_summary.txt, ped_null.csv
"""
import csv
import json
import math
import os
import random
import sys
from collections import defaultdict
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import exp_anchor as ea  # noqa: E402
import hub_frames_tests as ht  # noqa: E402

PED = ea.PED
N_POOL = 4000
BOND = 3.8
ENSEMBLES = [  # entry, file, protein, restraints
    ('PED00024', 'PED00024/PED00024e001.pdb', 'α-synuclein', 'PRE'),
    ('PED00003', 'PED00003/PED00003e001.pdb', 'β-synuclein', 'PRE'),
    ('PED00006', 'PED00006/PED00006e001.pdb', 'α/β-synuclein hybrid', 'PRE'),
    ('PED00022', 'PED00022/PED00022e001.pdb', 'drkN', 'NOE/PRE'),
    ('PED00156', 'PED00156/PED00156e002.pdb', 'drkN', 'NOE/PRE/smFRET'),
    ('PED00157', 'PED00157/PED00157e002.pdb', 'drkN', 'NOE/PRE/RDC/SAXS (X-EISD)'),
    ('PED00023', 'PED00023/PED00023e001.pdb', 'Sic1', 'PRE (16 models)'),
    ('PED00159', 'PED00159/PED00159e002.pdb', 'Sic1', 'PRE/SAXS'),
    ('PED00160', 'PED00160/PED00160e002.pdb', 'Sic1', 'PRE/SAXS/CS'),
    ('PED00161', 'PED00161/PED00161e002.pdb', 'pSic1 (7 phosphosites)', 'PRE/SAXS'),
    ('PED00225', 'PED00225/PED00225e000.pdb', 'eIF4G1 1-249', 'NOE/PRE/RDC/SAXS'),
    ('PED00498', 'PED00498/pdbfile.pdb', 'V2R C-tail', 'NOE/PRE'),
    ('PED00499', 'PED00499/pdbfile.pdb', 'GHSR C-tail', 'NOE/PRE'),
    ('PED00500', 'PED00500/pdbfile.pdb', 'β2AR C-tail', 'NOE/PRE'),
    ('PED00511', 'PED00511/pdbfile.pdb', 'PDEγ', 'NOE/PRE'),
    ('PED00013', 'PED00013/PED00013e001.pdb', 'PaaA2', 'NOE (no PRE)'),
    ('PED00004', 'PED00004/PED00004e001.pdb', 'p27 KID', 'NOE (no PRE)'),
    ('PED00494', 'PED00494/pdbfile.pdb', 'FUS NTD', 'DEER (no PRE)'),
]
POOLS = [('PED00425', 'PED00425/pdbfile.pdb', 'α-synuclein'), ('PED00423', 'PED00423/pdbfile.pdb', 'Sic1'),
         ('PED00427', 'PED00427/pdbfile.pdb', 'drkN')]
# one ensemble per protein for the per-protein tally (the most restrained / largest)
PRIMARY = ['PED00024', 'PED00003', 'PED00006', 'PED00157', 'PED00160', 'PED00161', 'PED00225', 'PED00498',
           'PED00499', 'PED00500', 'PED00511']
PHOS = {'SPO': 'x', 'SEP': 'x', 'TPO': 'y', 'PTR': 'z'}           # phosphoresidues get their own codes
CLASS = {**dict.fromkeys('FWY', 'FWY'), 'H': 'His', **dict.fromkeys('DE', 'DE'), **dict.fromkeys('KR', 'KR'),
         **dict.fromkeys('GSA', 'GSA'), 'x': 'phospho', 'y': 'phospho', 'z': 'phospho'}
CLASSES = ['FWY', 'DE', 'KR', 'GSA', 'His', 'phospho', 'other']
SIZE_AS = {'x': 'S', 'y': 'T', 'z': 'Y'}
out = []


def p(s=''):
    out.append(s)
    print(s, flush=True)


def load(path):
    """Like exp_anchor.ped_models, but phosphoresidues keep their own one-letter codes."""
    seq, models, cur, first = [], [], None, True
    for line in open(os.path.join(PED, path)):
        if line.startswith('MODEL'):
            if cur is not None and cur[0]:
                models.append(cur)
                first = False
            cur = ([], [], [])
        elif line.startswith(('ATOM', 'HETATM')) and line[12:16].strip() == 'CA' and line[16] in ' A':
            if cur is None:
                cur = ([], [], [])
            cur[0].append(float(line[30:38]))
            cur[1].append(float(line[38:46]))
            cur[2].append(float(line[46:54]))
            if first:
                rn = line[17:20]
                seq.append(PHOS.get(rn) or ea.AA3_VARIANTS.get(rn, 'X'))
        elif line.startswith('ENDMDL') and cur is not None:
            models.append(cur)
            cur, first = None, False
    if cur and cur[0]:
        models.append(cur)
    seq = ''.join(seq)
    return seq, [m for m in models if len(m[0]) == len(seq)]


def rg(m):
    n = len(m[0])
    cx, cy, cz = sum(m[0]) / n, sum(m[1]) / n, sum(m[2]) / n
    return math.sqrt(sum((x - cx) ** 2 + (y - cy) ** 2 + (z - cz) ** 2 for x, y, z in zip(*m)) / n)


def sigma_vec(seq):
    s = {r['one']: float(r['sigmas']) * 10.0 for r in csv.DictReader(open(ht.CAL_PARAMS))}
    return [s[SIZE_AS.get(a, a)] for a in seq]


def grow_chain(sig, rng, max_tries=60, back=5, angles=(85.0, 150.0)):
    """One self-avoiding Calpha chain with residue-specific hard cores and no attraction."""
    n = len(sig)
    X, Y, Z = [0.0], [0.0], [0.0]
    X.append(BOND)
    Y.append(0.0)
    Z.append(0.0)
    fails = 0
    while len(X) < n:
        i = len(X)
        ok = False
        for _ in range(max_tries):
            # previous two bonds define a local frame; sample bond angle and dihedral
            if i >= 2:
                ax, ay, az = X[i - 1] - X[i - 2], Y[i - 1] - Y[i - 2], Z[i - 1] - Z[i - 2]
            else:
                ax, ay, az = 1.0, 0.0, 0.0
            na = math.sqrt(ax * ax + ay * ay + az * az)
            ax, ay, az = ax / na, ay / na, az / na
            # any vector not parallel to a
            tx, ty, tz = (0.0, 1.0, 0.0) if abs(ax) > 0.9 else (1.0, 0.0, 0.0)
            ux, uy, uz = ay * tz - az * ty, az * tx - ax * tz, ax * ty - ay * tx
            nu = math.sqrt(ux * ux + uy * uy + uz * uz)
            ux, uy, uz = ux / nu, uy / nu, uz / nu
            vx, vy, vz = ay * uz - az * uy, az * ux - ax * uz, ax * uy - ay * ux
            theta = math.radians(rng.uniform(*angles))      # Calpha virtual bond angle
            phi = rng.uniform(-math.pi, math.pi)
            # new bond direction makes angle (180 - theta) with the previous bond
            ca, sa = math.cos(math.pi - theta), math.sin(math.pi - theta)
            dx = ca * ax + sa * (math.cos(phi) * ux + math.sin(phi) * vx)
            dy = ca * ay + sa * (math.cos(phi) * uy + math.sin(phi) * vy)
            dz = ca * az + sa * (math.cos(phi) * uz + math.sin(phi) * vz)
            x, y, z = X[i - 1] + BOND * dx, Y[i - 1] + BOND * dy, Z[i - 1] + BOND * dz
            si = sig[i]
            clash = False
            for j in range(i - 2):          # |i - j| >= 3 (i-2 is set by the bond angle)
                dmin = 0.4 * (si + sig[j])   # 0.8 * (sigma_i + sigma_j) / 2
                ex = x - X[j]
                if ex > dmin or ex < -dmin:
                    continue
                ey = y - Y[j]
                if ey > dmin or ey < -dmin:
                    continue
                ez = z - Z[j]
                if ex * ex + ey * ey + ez * ez < dmin * dmin:
                    clash = True
                    break
            if not clash:
                X.append(x)
                Y.append(y)
                Z.append(z)
                ok = True
                break
        if not ok:
            fails += 1
            if fails > 200:
                return None
            k = max(2, len(X) - back)
            del X[k:], Y[k:], Z[k:]
    return X, Y, Z


def null_models(seq, exp_models, seed):
    rng = random.Random(seed)
    sig = sigma_vec(seq)
    pool = []
    # three chain stiffnesses (residue-neutral) so that the Rg range covers compact and expanded experimental ensembles
    for angles in ((85.0, 150.0), (100.0, 150.0), (115.0, 155.0)):
        target = len(pool) + N_POOL // 3
        while len(pool) < target:
            c = grow_chain(sig, rng, angles=angles)
            if c is not None:
                pool.append(c)
    prg = sorted((rg(c), k) for k, c in enumerate(pool))
    vals = [v for v, _ in prg]
    chosen = []
    import bisect
    for m in exp_models:
        r = rg(m)
        k = bisect.bisect_left(vals, r)
        k = min(range(max(0, k - 1), min(len(vals), k + 1)), key=lambda t: abs(vals[t] - r))
        chosen.append(pool[prg[k][1]])
    return chosen, sum(rg(m) for m in exp_models) / len(exp_models), sum(rg(c) for c in chosen) / len(chosen), \
        sum(vals) / len(vals)


def class_means(seq, z):
    by = defaultdict(list)
    for a, x in zip(seq, z):
        by[CLASS.get(a, 'other')].append(x)
    return {k: (sum(v) / len(v), len(v)) for k, v in by.items()}


def centrality(seq, models, rule, c):
    seq_size = ''.join(SIZE_AS.get(a, a) for a in seq)        # bead-size rule: phosphoresidues use S/T/Y sizes
    z, used = ea.ensemble_centrality(seq_size, models, rule, c)
    return z


def one(job):
    entry, path, kind, c, seed = job
    seq, models = load(path)
    if len(models) > 1000:
        models = models[:: math.ceil(len(models) / 1000)]
    res = {'entry': entry, 'kind': kind, 'N': len(seq), 'models': len(models)}
    if kind == 'pool':
        for rule in ('fixed', 'size'):
            res[rule] = class_means(seq, centrality(seq, models, rule, c))
        return res
    nm, rg_exp, rg_null, rg_pool = null_models(seq, models, seed)
    res.update(rg_exp=rg_exp, rg_null=rg_null, rg_pool=rg_pool)
    for rule in ('fixed', 'size'):
        ze = centrality(seq, models, rule, c)
        zn = centrality(seq, nm, rule, c)
        res[rule] = (class_means(seq, ze), class_means(seq, zn))
        res[f'per_res_{rule}'] = (seq, ze, zn)
    return res


def main():
    c = json.load(open(ht.C_FILE))['calvados']
    jobs = [(e, f, 'exp', c, 1000 + k) for k, (e, f, prot, r) in enumerate(ENSEMBLES)]
    jobs += [(e, f, 'pool', c, 0) for e, f, prot in POOLS]
    with Pool(9) as pool:
        res = {r['entry']: r for r in pool.imap_unordered(one, jobs)}
    meta = {e: (prot, r) for e, f, prot, r in ENSEMBLES}
    rows = []
    p('Experimental PED ensembles vs a residue-identity-neutral null (self-avoiding chains, bead-size excluded volume, '
      'no attraction, Rg-matched). Excess = experimental - null mean position-adjusted z.')
    for rule in ('fixed', 'size'):
        p(f'\nRule: {rule}')
        p(f'  {"entry":9s} {"protein":24s} {"restraints":24s} {"N":>4s} {"Rg exp/null/pool":>18s}   '
          + ''.join(f'{k:>16s}' for k in ('FWY', 'DE', 'KR', 'GSA', 'phospho')))
        for e, f, prot, r in ENSEMBLES:
            x = res[e]
            ee, nn = x[rule]
            line = f'  {e:9s} {prot[:24]:24s} {r[:24]:24s} {x["N"]:4d} {x["rg_exp"]:6.1f}/{x["rg_null"]:5.1f}/{x["rg_pool"]:5.1f}   '
            for k in ('FWY', 'DE', 'KR', 'GSA', 'phospho'):
                if k in ee:
                    line += f'{ee[k][0] - nn[k][0]:+7.2f} ({ee[k][0]:+.2f})'
                    rows.append([rule, e, prot, r, k, ee[k][1], f'{ee[k][0]:.4f}', f'{nn[k][0]:.4f}', f'{ee[k][0] - nn[k][0]:.4f}'])
                else:
                    line += f'{"-":>16s}'
            p(line)
        for label, ids in (('PRE ensembles, one per protein', PRIMARY),
                           ('all PRE ensembles', [e for e, f, prot, r in ENSEMBLES if 'PRE' in r and 'no PRE' not in r]),
                           ('non-PRE ensembles', [e for e, f, prot, r in ENSEMBLES if 'no PRE' in r])):
            fwy = [res[e][rule][0]['FWY'][0] - res[e][rule][1]['FWY'][0] for e in ids if 'FWY' in res[e][rule][0]]
            de = [res[e][rule][0]['DE'][0] - res[e][rule][1]['DE'][0] for e in ids if 'DE' in res[e][rule][0]]
            sp = lambda k, n: sum(math.comb(n, i) for i in range(k, n + 1)) / 2 ** n if n else float('nan')  # noqa: E731
            kf, kd = sum(v > 0 for v in fwy), sum(v < 0 for v in de)
            p(f'  {label}: aromatic excess > 0 in {kf}/{len(fwy)} (sign p = {sp(kf, len(fwy)):.3f}, mean {sum(fwy) / max(1, len(fwy)):+.2f}); '
              f'acidic excess < 0 in {kd}/{len(de)} (sign p = {sp(kd, len(de)):.3f}, mean {sum(de) / max(1, len(de)):+.2f})')
    p('\nValidation of the null against PED\'s restraint-free IDPConformerGenerator pools (class mean z: pool vs our null '
      'built for the experimental ensemble of the same protein)')
    pairs = {'PED00425': 'PED00024', 'PED00423': 'PED00160', 'PED00427': 'PED00157'}
    for rule in ('fixed', 'size'):
        for pe, ee_ in pairs.items():
            pm = res[pe][rule]
            nm_ = res[ee_][rule][1]
            p(f'  {rule:5s} {pe} ({meta[ee_][0]}): ' + '  '.join(
                f'{k} {pm[k][0]:+.2f}/{nm_[k][0]:+.2f}' for k in ('FWY', 'DE', 'KR', 'GSA') if k in pm and k in nm_))
    with open(os.path.join(HERE, 'ped_null_per_residue.csv'), 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['entry', 'protein', 'restraints', 'pos', 'aa', 'N', 'z_exp_fixed', 'z_null_fixed', 'z_exp_size', 'z_null_size'])
        for e, f, prot, r in ENSEMBLES:
            sf, ef, nf = res[e]['per_res_fixed']
            ss, es, ns = res[e]['per_res_size']
            for i, a in enumerate(sf):
                w.writerow([e, prot, r, i, a, len(sf), f'{ef[i]:.4f}', f'{nf[i]:.4f}', f'{es[i]:.4f}', f'{ns[i]:.4f}'])
    with open(os.path.join(HERE, 'ped_null.csv'), 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['rule', 'entry', 'protein', 'restraints', 'class', 'n', 'z_exp', 'z_null', 'excess'])
        w.writerows(rows)
    open(os.path.join(HERE, 'ped_null_summary.txt'), 'w').write('\n'.join(out) + '\n')


if __name__ == '__main__':
    main()
