"""
context_neighbourhood_typed.py — the neighbourhood of hubs compared with residues of the same type (2026-10-04 review fix).

context_model_by_ff.py and context_neighbourhood_null.py compared the +-5 window around hubs with the window around all
residues. Hubs are often glycines themselves (about 40% in CALVADOS-2, against 10.5% of residues) and glycines cluster,
so part of that glycine enrichment restated the type of the hub. Here the expected composition of a hub's window is the
mean window composition of all residues of the hub's own type in the same label set, so the enrichment measures the
surroundings beyond the residue's own type. Hubs, windows (residue itself excluded) and classes are as before; 95%
intervals from 1,000 bootstrap resamples of proteins recompute the type baselines in every resample.
  part 'full'  the three force fields (data_v3, all proteins) and the 27 PED ensembles, hubs by position-adjusted occupancy
  part 'null'  the 507 sequences with all three force fields, hubs by occupancy, simulation z, null z and excess z
               (../sim_null_v3_per_residue.csv)
  ~/.venvs/bender/bin/python context_neighbourhood_typed.py [full|null|both]
      ->  context_neighbourhood_typed.csv, context_neighbourhood_null_typed.csv, context_neighbourhood_typed.txt
"""
import csv
import os
import sys
from collections import defaultdict

import numpy as np

import context_model_by_ff as C
import context_neighbourhood_null as NN
import sticker_spacer_by_ff as SB
import train_hub_model as T

AAS = C.AAS
NBOOT = 1000
out = []


def p(s=''):
    out.append(s)
    print(s, flush=True)


def typed(seqs, vals, rng):
    """Window composition around hubs divided by the mean window composition of residues of the hub's own type."""
    P, K = len(seqs), len(C.CLASSES)
    A = np.zeros((P, 20, K))          # summed window composition of all residues of each type
    Nt = np.zeros((P, 20))            # residues of each type
    Hs = np.zeros((P, K))             # summed window composition of the hubs
    Ht = np.zeros((P, 20))            # hubs of each type
    for j, (s, v) in enumerate(zip(seqs, vals)):
        n = len(s)
        r = C.AT.posadj(np.asarray(v, float))
        k = max(1, int(round(0.1 * n)))
        hub = np.zeros(n, bool)
        hub[np.argsort(-(r + np.random.default_rng(n).uniform(0, 1e-9, n)))[:k]] = True
        frac = NN.window_fractions(s)
        for i, a in enumerate(s):
            if a not in AAS or np.isnan(frac[i]).any():
                continue
            t = AAS.index(a)
            A[j, t] += frac[i]
            Nt[j, t] += 1
            if hub[i]:
                Hs[j] += frac[i]
                Ht[j, t] += 1

    def ratio(w):
        base = (np.tensordot(w, A, axes=([w.ndim - 1], [0]))) / np.maximum(w @ Nt, 1e-12)[..., None]      # (B?, 20, K)
        ht = w @ Ht                                                                         # (B?, 20)
        exp_ = (ht[..., None] * base).sum(-2) / ht.sum(-1)[..., None]
        obs = (w @ Hs) / ht.sum(-1)[..., None]
        return obs / exp_
    est = ratio(np.ones(P))
    W = np.stack([np.bincount(rng.integers(0, P, P), minlength=P) for _ in range(NBOOT)]).astype(float)
    bs = ratio(W)
    lo, hi = np.percentile(bs, [2.5, 97.5], axis=0)
    hub_g = Ht[:, AAS.index('G')].sum() / Ht.sum()
    return est, lo, hi, hub_g


def write(rows, name):
    with open(os.path.join(C.HERE, name), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def table(res, cols):
    lab = max(len(c[0]) for c in C.CLASSES)
    p(f'  {"neighbour class":{lab}s}' + ''.join(f'{c:>22s}' for c in cols))
    for ci, (cname, _) in enumerate(C.CLASSES):
        p(f'  {cname:{lab}s}' + ''.join(f'{res[c][0][ci]:>8.2f} ({res[c][1][ci]:.2f}-{res[c][2][ci]:.2f})' for c in cols))


def full(rng):
    prots, y = T.load(os.path.join(C.HERE, 'data_v3'), 'occupancy')
    sets = {C.NICE[ff]: C.label_set(prots, y, ff)[1:] for ff in C.FFS}
    sets['PED'] = SB.ped_sets()
    res, rows = {}, []
    for name, (seqs, Y) in sets.items():
        res[name] = typed(seqs, Y, rng)
        est, lo, hi, hub_g = res[name]
        p(f'{name}: {len(seqs)} proteins, {100 * hub_g:.0f}% of hubs are glycine')
        for ci, (cname, letters) in enumerate(C.CLASSES):
            rows.append({'label_set': name, 'neighbour_class': cname, 'residues': letters, 'enrichment': est[ci],
                         'lo': lo[ci], 'hi': hi[ci], 'proteins': len(seqs)})
    p('\nWindow composition around hubs relative to residues of the same type (all proteins)')
    table(res, list(sets))
    write(rows, 'context_neighbourhood_typed.csv')


def null(rng):
    rows_in = defaultdict(dict)
    for r in csv.DictReader(open(os.path.join(C.EV, 'sim_null_v3_per_residue.csv'))):
        rows_in[(r['force_field'], r['UniProt_ID'])][int(r['pos'])] = (r['aa'], float(r['z_sim_fixed']), float(r['z_null_fixed']))
    prots, y = T.load(os.path.join(C.HERE, 'data_v3'), 'occupancy')
    rows = []
    for ff in C.FFS:
        seqs, occ, zsim, znull = [], [], [], []
        for u in sorted(u for f, u in rows_in if f == ff):
            d = rows_in[(ff, u)]
            s = ''.join(d[i][0] for i in range(len(d)))
            if len(y.get((ff, u), {})) != len(s) or prots[(ff, u)]['sequence'] != s:
                continue
            seqs.append(s)
            occ.append(np.array([y[(ff, u)][i] for i in range(len(s))]))
            zsim.append(np.array([d[i][1] for i in range(len(s))]))
            znull.append(np.array([d[i][2] for i in range(len(s))]))
        variants = {'occupancy': occ, 'simulation': zsim, 'null': znull, 'excess': [a - b for a, b in zip(zsim, znull)]}
        res = {v: typed(seqs, vals, rng) for v, vals in variants.items()}
        p(f'\n{C.NICE[ff]}: {len(seqs)} proteins (same-type baseline)')
        table(res, list(variants))
        for v, (est, lo, hi, _) in res.items():
            for ci, (cname, _) in enumerate(C.CLASSES):
                rows.append({'force_field': C.NICE[ff], 'hubs_by': v, 'neighbour_class': cname, 'enrichment': est[ci],
                             'lo': lo[ci], 'hi': hi[ci], 'proteins': len(seqs)})
    write(rows, 'context_neighbourhood_null_typed.csv')


def main():
    part = sys.argv[1] if len(sys.argv) > 1 else 'both'
    rng = np.random.default_rng(5)
    if part in ('full', 'both'):
        full(rng)
    if part in ('null', 'both'):
        null(rng)
    open(os.path.join(C.HERE, f'context_neighbourhood_typed_{part}.txt'), 'w').write('\n'.join(out) + '\n')


if __name__ == '__main__':
    main()
