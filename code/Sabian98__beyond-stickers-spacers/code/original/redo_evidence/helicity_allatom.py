"""
helicity_allatom.py — per-residue secondary structure from the Cα trace of the all-atom (set C) and expert-built PED
ensembles of Aβ40, ACTR, drkN, PaaA2 and α-synuclein, and whether the residues the frozen model misses are structured.

Cα-only assignment (no DSSP on this Mac): for each quadruple i..i+3, the Cα pseudo-dihedral and the i→i+3 distance.
  helical quadruple   pseudo-dihedral 30 to 80 deg and d(i, i+3) < 6.0 Å
  extended quadruple  pseudo-dihedral |angle| > 150 deg or -150 to -90 deg, and d(i, i+3) > 9.0 Å
A residue is helical in a model when it lies in two consecutive helical quadruples (a turn), extended likewise.
Same models as the labels (ped_v2_labels.load_pooled / ped_null.load, at most 1,000).
  ~/.venvs/bender/bin/python helicity_allatom.py  ->  helicity_allatom.csv, helicity_allatom_summary.txt
"""
import csv
import math
import os
import sys
from collections import defaultdict
from multiprocessing import Pool

import numpy as np
from scipy.stats import rankdata, spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ped_null as pn  # noqa: E402

ENTRIES = {  # entry: (protein, ensemble kind, file)
    'PED00531': ('Aβ40', 'C36m'), 'PED00532': ('Aβ40', 'C22*'), 'PED00533': ('Aβ40', 'a99SB-disp'),
    'PED00534': ('ACTR', 'C36m'), 'PED00535': ('ACTR', 'C22*'), 'PED00536': ('ACTR', 'a99SB-disp'),
    'PED00537': ('drkN', 'C36m'), 'PED00538': ('drkN', 'C22*'), 'PED00539': ('drkN', 'a99SB-disp'),
    'PED00540': ('PaaA2', 'C36m'), 'PED00541': ('PaaA2', 'C22*'), 'PED00542': ('PaaA2', 'a99SB-disp'),
    'PED00543': ('α-synuclein', 'C36m'), 'PED00544': ('α-synuclein', 'C22*'), 'PED00545': ('α-synuclein', 'a99SB-disp'),
    'PED00022': ('drkN', 'expert TraDES/ENSEMBLE'), 'PED00156': ('drkN', 'expert X-EISD'),
    'PED00157': ('drkN', 'expert X-EISD'), 'PED00013': ('PaaA2', 'expert'), 'PED00024': ('α-synuclein', 'expert'),
}
out = []


def p(s=''):
    out.append(s)
    print(s, flush=True)


def files_for(e):
    for r in csv.DictReader(open(os.path.join(HERE, 'ped_v2_set.csv'))):
        if r['id'] == e:
            return r['files']
    for ent, f, _, _ in pn.ENSEMBLES:
        if ent == e:
            return f
    raise KeyError(e)


def pseudo_dihedral(a, b, c, d):
    b0, b1, b2 = a - b, c - b, d - c
    b1n = b1 / np.linalg.norm(b1, axis=-1, keepdims=True)
    v = b0 - (b0 * b1n).sum(-1, keepdims=True) * b1n
    w = b2 - (b2 * b1n).sum(-1, keepdims=True) * b1n
    x = (v * w).sum(-1)
    y = (np.cross(b1n, v) * w).sum(-1)
    return np.degrees(np.arctan2(y, x))


def one(e):
    seq, models = None, []
    for f in files_for(e).split(';'):
        s, m = pn.load(f)
        seq = seq or s
        models += [x for x in m if len(x[0]) == len(seq)]
    if len(models) > 1000:
        models = models[::math.ceil(len(models) / 1000)]
    n = len(seq)
    X = np.array([np.stack(m, axis=1) for m in models])            # models x n x 3
    a, b, c, d = X[:, :-3], X[:, 1:-2], X[:, 2:-1], X[:, 3:]
    phi = pseudo_dihedral(a, b, c, d)
    d14 = np.linalg.norm(d - a, axis=-1)
    hq = (phi > 30) & (phi < 80) & (d14 < 6.0)
    eq = ((np.abs(phi) > 150) | ((phi > -150) & (phi < -90))) & (d14 > 9.0)
    hel = np.zeros((len(models), n), bool)
    ext = np.zeros((len(models), n), bool)
    for q, res in ((hq, hel), (eq, ext)):
        turn = q[:, :-1] & q[:, 1:]                                   # quadruples k and k+1 both set
        for k in range(turn.shape[1]):
            res[:, k:k + 5] |= turn[:, [k]]
    rg = np.sqrt(((X - X.mean(1, keepdims=True)) ** 2).sum(-1).mean(1)) / 10.0
    return e, seq, hel.mean(0), ext.mean(0), len(models), rg


def partial(x, y, z):
    rx, ry, rz = (rankdata(v) for v in (x, y, z))
    A = np.vstack([rz, np.ones_like(rz)]).T
    ex = rx - A @ np.linalg.lstsq(A, rx, rcond=None)[0]
    ey = ry - A @ np.linalg.lstsq(A, ry, rcond=None)[0]
    return float(np.corrcoef(ex, ey)[0, 1])


def main():
    with Pool(8) as pool:
        res = {r[0]: r for r in pool.imap_unordered(one, list(ENTRIES))}
    occ = defaultdict(list)
    for r in csv.DictReader(open(os.path.join(HERE, 'hub_predictor', 'ped_v2_occupancy.csv'))):
        occ[r['entry']].append(float(r['occ_exp']))
    for r in csv.DictReader(open(os.path.join(HERE, 'hub_predictor', 'experimental_occupancy.csv'))):
        occ.setdefault('X' + r['entry'], []).append(float(r['occ_exp']))
    pred = defaultdict(list)
    for r in csv.DictReader(open(os.path.join(HERE, 'gaivi', 'predictions_frozen_2026-10-02_ped_v2.csv'))):
        pred[r['id']].append(float(r['consensus']))
    with open(os.path.join(HERE, 'helicity_allatom.csv'), 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['entry', 'protein', 'kind', 'pos', 'aa', 'helix', 'extended'])
        for e, (prot, kind) in ENTRIES.items():
            _, seq, h, x, _, _ = res[e]
            w.writerows([e, prot, kind, i, a, f'{h[i]:.3f}', f'{x[i]:.3f}'] for i, a in enumerate(seq))

    p('Cα-trace secondary structure of all-atom and expert-built PED ensembles, and the frozen model where it applies')
    p(f'  {"entry":9s} {"protein":12s} {"ensemble":24s} {"models":>6s} {"helix":>6s} {"ext":>5s} {"Rg nm":>6s} '
      f'{"max hel":>7s}  {"model partial":>13s} {"non-helical":>11s} {"err~helix":>9s}')
    for e, (prot, kind) in ENTRIES.items():
        _, seq, h, x, nm, rg = res[e]
        n = len(seq)
        line = (f'  {e:9s} {prot:12s} {kind:24s} {nm:6d} {h.mean():6.2f} {x.mean():5.2f} {np.median(rg):6.2f} '
                f'{h.max():7.2f}')
        if e in pred and e in occ:
            pr, oe = np.array(pred[e]), np.array(occ[e])
            end = np.minimum(np.arange(n), n - 1 - np.arange(n)).astype(float)
            pa = partial(pr, oe, end)
            keep = h < 0.2
            pa_nh = partial(pr[keep], oe[keep], end[keep]) if keep.sum() >= 15 else float('nan')
            err = rankdata(oe) - rankdata(pr)                     # positive where the ensemble ranks higher than predicted
            rho_err = spearmanr(err, h).correlation if h.std() > 0 else float('nan')
            line += f'  {pa:+13.2f} {pa_nh:+11.2f} ({keep.sum():2d}) {rho_err:+9.2f}'
        p(line)
    p('\n  helix/ext = mean fraction of models in which a residue is helical/extended; max hel = most helical residue;')
    p('  non-helical = model partial correlation over residues with helix < 0.2 (count); err~helix = Spearman between the')
    p('  rank error (ensemble rank minus predicted rank) and helicity.')

    p('\nHub occupancy against helicity within each ensemble (Spearman over residues)')
    for e, (prot, kind) in ENTRIES.items():
        _, seq, h, x, _, _ = res[e]
        o = occ.get(e) or occ.get('X' + e)
        if o and len(o) == len(seq) and h.std() > 0:
            p(f'  {e:9s} {prot:12s} {kind:24s} occupancy~helix {spearmanr(o, h).correlation:+.2f}  '
              f'occupancy~extended {spearmanr(o, x).correlation:+.2f}')

    p('\nHelical stretches (residues with helix >= 0.3), 1-based')
    for e, (prot, kind) in ENTRIES.items():
        _, seq, h, _, _, _ = res[e]
        segs, start = [], None
        for i, v in enumerate(list(h) + [0]):
            if v >= 0.3 and start is None:
                start = i
            elif v < 0.3 and start is not None:
                segs.append(f'{seq[start]}{start + 1}-{seq[i - 1]}{i}')
                start = None
        p(f'  {e:9s} {prot:12s} {kind:24s} {", ".join(segs) or "none"}')
    open(os.path.join(HERE, 'helicity_allatom_summary.txt'), 'w').write('\n'.join(out) + '\n')


if __name__ == '__main__':
    main()
