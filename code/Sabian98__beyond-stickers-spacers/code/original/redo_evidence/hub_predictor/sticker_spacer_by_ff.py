"""
sticker_spacer_by_ff.py — parts A and B of sticker_spacer_test.py for every label set, with identical definitions:
  A  hubs = the top 10% of each protein by hub occupancy (ties broken by a tiny uniform jitter, numpy rng seeded 0,
     proteins in data order); enrichment of aromatic F/W/Y, stickers F/W/Y/R and spacers G/S among hubs relative to
     their share of all residues, pooled over proteins.
  B  share of the within-protein variance of position-adjusted occupancy (quadratic trend in relative position removed,
     residual scaled to unit variance per protein) that residue identity explains (between-type sum of squares / total).
Label sets: the three force fields in a dataset folder (default data_v3) and the 27 experiment-built PED ensembles
(18 development + 9 blind, i.e. without the unweighted A-beta40 and ACTR simulations).
  ~/.venvs/bender/bin/python sticker_spacer_by_ff.py [data folder]  ->  sticker_spacer_by_ff_summary.txt, sticker_spacer_by_ff.csv
"""
import csv
import os
import sys
from collections import defaultdict

import numpy as np

import train_hub_model as T

HERE = os.path.dirname(os.path.abspath(__file__))
EV = os.path.dirname(HERE)
AAS = 'ACDEFGHIKLMNPQRSTVWY'
SETS = {'aromatic FWY': set('FWY'), 'stickers FWYR': set('FWYR'), 'spacers GS': set('GS')}
RAW_MD = {'PED00533', 'PED00536'}
out = []


def p(s=''):
    out.append(s)
    print(s, flush=True)


def ped_sets():
    occ = defaultdict(dict)
    for f, col in (('experimental_occupancy.csv', 'occ_exp'), ('fresh_ped_occupancy.csv', 'occ_exp'),
                   ('ped_extra_occupancy.csv', 'occ_exp')):
        for r in csv.DictReader(open(os.path.join(HERE, f))):
            if f == 'ped_extra_occupancy.csv' and r['set'] != 'blind':
                continue                                   # the drkN / Sic1 development entries are already in
            occ[r['entry']][int(r['pos'])] = (r['aa'], float(r[col]))
    seqs, Y = [], []
    for e in sorted(occ):
        if e in RAW_MD:
            continue
        d = occ[e]
        seqs.append(''.join(d[i][0] for i in range(len(d))))
        Y.append(np.array([d[i][1] for i in range(len(d))]))
    return seqs, Y


def analyse(name, seqs, Y):
    rng = np.random.default_rng(0)
    tot = {'res': 0, 'hub': 0}
    cnt = {k: {'res': 0, 'hub': 0} for k in SETS}
    resid, aa_all = [], []
    for s, occ in zip(seqs, Y):
        n = len(s)
        k = max(1, int(round(0.1 * n)))
        hubs = set(np.argsort(-(occ + rng.uniform(0, 1e-9, n)))[:k])
        tot['res'] += n
        tot['hub'] += k
        for key, st in SETS.items():
            idx = [i for i, a in enumerate(s) if a in st]
            cnt[key]['res'] += len(idx)
            cnt[key]['hub'] += sum(i in hubs for i in idx)
        x = np.arange(n) / max(n - 1, 1)
        A = np.vstack([np.ones(n), x, x ** 2]).T
        r = occ - A @ np.linalg.lstsq(A, occ, rcond=None)[0]
        if r.std() > 0:
            resid.append(r / r.std())
            aa_all.append(np.array(list(s)))
    r = np.concatenate(resid)
    a = np.concatenate(aa_all)
    between = sum((a == t).sum() * r[a == t].mean() ** 2 for t in AAS if (a == t).any())
    row = {'label_set': name, 'proteins': len(seqs), 'residues': tot['res'],
           'identity_share': between / (r ** 2).sum()}
    for key, c in cnt.items():
        row[f'enrich_{key.split()[0]}'] = (c['hub'] / tot['hub']) / (c['res'] / tot['res'])
    p(f'  {name:28s} {len(seqs):6d} proteins  identity {100 * row["identity_share"]:5.1f}%   enrichment among hubs: '
      + '  '.join(f'{key} {row["enrich_" + key.split()[0]]:.2f}' for key in SETS))
    return row


def main():
    data = sys.argv[1] if len(sys.argv) > 1 else 'data_v3'
    prots, y = T.load(os.path.join(HERE, data), 'occupancy')
    p(f'Hubs versus stickers and spacers, and the identity share of position-adjusted occupancy ({data})')
    rows = []
    for ff in ('CALVADOS-2', 'MPIPI-GG', 'Mpipi-Recharged'):
        keys = [k for k in prots if k[0] == ff and len(y.get(k, {})) == len(prots[k]['sequence'])]
        seqs = [prots[k]['sequence'] for k in keys]
        Y = [np.array([y[k][i] for i in range(len(s))]) for k, s in zip(keys, seqs)]
        rows.append(analyse(ff, seqs, Y))
    rows.append(analyse('PED, built from experiment', *ped_sets()))
    with open(os.path.join(HERE, 'sticker_spacer_by_ff.csv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    open(os.path.join(HERE, 'sticker_spacer_by_ff_summary.txt'), 'w').write('\n'.join(out) + '\n')


if __name__ == '__main__':
    main()
