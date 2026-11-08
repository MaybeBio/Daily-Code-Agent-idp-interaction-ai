"""
context_neighbourhood_null.py — are the glycine-rich, proline-poor windows around hubs a bead-size effect?

On the 507 sequences simulated with all three force fields, hubs are the top 10% of residues per protein by
  occupancy   hub occupancy (data_v3), position-adjusted as in context_model_by_ff.py, for reference
  simulation  position-adjusted z of time-averaged betweenness in the simulation (sim_null_v3_per_residue.csv)
  null        the same in the residue-neutral null (self-avoiding chains with the force field's own bead sizes and no
              attraction, Rg-matched to the simulation), i.e. geometry alone
  excess      simulation minus null, i.e. what attraction adds beyond geometry
For each, the composition of the +-5 window around hubs (the residue itself excluded) is divided by the composition of
the window around all residues, with 95% bootstrap intervals over proteins (1,000 resamples).
  ~/.venvs/bender/bin/python context_neighbourhood_null.py  ->  context_neighbourhood_null.txt, context_neighbourhood_null.csv
"""
import csv
import os
from collections import defaultdict

import numpy as np

import context_model_by_ff as C
import train_hub_model as T

AT = C.AT
FFS = C.FFS
NBOOT = 1000
out = []


def p(s=''):
    out.append(s)
    print(s, flush=True)


def window_fractions(s):
    n = len(s)
    frac = np.zeros((n, len(C.CLASSES)))
    lo = np.clip(np.arange(n) - 5, 0, n)
    hi = np.clip(np.arange(n) + 6, 0, n)
    cnt = hi - lo - 1
    for ci, (_, letters) in enumerate(C.CLASSES):
        x = np.array([a in letters for a in s], float)
        c = np.concatenate([[0.0], np.cumsum(x)])
        frac[:, ci] = np.where(cnt > 0, (c[hi] - c[lo] - x) / np.maximum(cnt, 1), np.nan)
    return frac


def neighbourhood(seqs, vals, rng):
    Hs, Hn, As, An = [], [], [], []
    for s, v in zip(seqs, vals):
        n = len(s)
        r = AT.posadj(np.asarray(v, float))
        k = max(1, int(round(0.1 * n)))
        hub = np.zeros(n, bool)
        hub[np.argsort(-(r + np.random.default_rng(n).uniform(0, 1e-9, n)))[:k]] = True
        frac = window_fractions(s)
        ok = ~np.isnan(frac).any(1)
        Hs.append(frac[hub & ok].sum(0)), Hn.append((hub & ok).sum())
        As.append(frac[ok].sum(0)), An.append(ok.sum())
    Hs, Hn, As, An = map(np.array, (Hs, Hn, As, An))
    est = (Hs.sum(0) / Hn.sum()) / (As.sum(0) / An.sum())
    W = np.stack([np.bincount(rng.integers(0, len(Hn), len(Hn)), minlength=len(Hn)) for _ in range(NBOOT)]).astype(float)
    bs = ((W @ Hs) / (W @ Hn)[:, None]) / ((W @ As) / (W @ An)[:, None])
    lo, hi = np.percentile(bs, [2.5, 97.5], axis=0)
    return est, lo, hi


def main():
    rows = defaultdict(dict)
    for r in csv.DictReader(open(os.path.join(C.EV, 'sim_null_v3_per_residue.csv'))):
        rows[(r['force_field'], r['UniProt_ID'])][int(r['pos'])] = (r['aa'], float(r['z_sim_fixed']),
                                                                    float(r['z_null_fixed']))
    prots, y = T.load(os.path.join(C.HERE, 'data_v3'), 'occupancy')
    rng = np.random.default_rng(3)
    res, csv_rows = {}, []
    for ff in FFS:
        ids = sorted(u for f, u in rows if f == ff)
        seqs, occ, zs, zn = [], [], [], []
        for u in ids:
            d = rows[(ff, u)]
            if sorted(d) != list(range(len(d))):
                continue
            s = ''.join(d[i][0] for i in range(len(d)))
            k = (ff, u)
            if len(y.get(k, {})) != len(s) or prots[k]['sequence'] != s:
                continue
            seqs.append(s)
            occ.append(np.array([y[k][i] for i in range(len(s))]))
            zs.append(np.array([d[i][1] for i in range(len(s))]))
            zn.append(np.array([d[i][2] for i in range(len(s))]))
        variants = {'occupancy': occ, 'simulation': zs, 'null': zn, 'excess': [a - b for a, b in zip(zs, zn)]}
        for name, vals in variants.items():
            res[(ff, name)] = neighbourhood(seqs, vals, rng)
            est, lo, hi = res[(ff, name)]
            for ci, (cname, letters) in enumerate(C.CLASSES):
                csv_rows.append({'force_field': C.NICE[ff], 'hubs_by': name, 'neighbour_class': cname,
                                 'enrichment': est[ci], 'lo': lo[ci], 'hi': hi[ci], 'proteins': len(seqs)})
        p(f'{C.NICE[ff]}: {len(seqs)} proteins')
    p('\nComposition of the +-5 window around hubs relative to all residues (95% bootstrap interval over proteins)')
    for ff in FFS:
        p(f'\n{C.NICE[ff]}')
        p(f'  {"neighbour class":22s}' + ''.join(f'{v:>22s}' for v in ('occupancy', 'simulation', 'null', 'excess')))
        for ci, (cname, _) in enumerate(C.CLASSES):
            p(f'  {cname:22s}' + ''.join(f'{res[(ff, v)][0][ci]:>8.2f} ({res[(ff, v)][1][ci]:.2f}-{res[(ff, v)][2][ci]:.2f})'
                                       for v in ('occupancy', 'simulation', 'null', 'excess')))
    with open(os.path.join(C.HERE, 'context_neighbourhood_null.csv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(csv_rows[0]))
        w.writeheader()
        w.writerows(csv_rows)
    open(os.path.join(C.HERE, 'context_neighbourhood_null.txt'), 'w').write('\n'.join(out) + '\n')


if __name__ == '__main__':
    main()
