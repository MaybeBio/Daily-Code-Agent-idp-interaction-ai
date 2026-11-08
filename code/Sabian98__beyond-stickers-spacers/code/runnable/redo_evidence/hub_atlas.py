"""
hub_atlas.py — predicted hubs in designed sticker-spacer sequences and across the human IDRome (hub_atlas_plan.md,
hashed 2026-10-11 before any prediction). Exploratory and descriptive: frozen consensus model only, no simulations.
  OMP_NUM_THREADS=1 ~/.venvs/bender/bin/python hub_atlas.py  ->  hub_atlas/*.csv, hub_atlas/summary.txt
"""
import csv
import os
import pickle
import sys
from collections import defaultdict
from multiprocessing import Pool

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, 'hub_predictor'))
sys.path.insert(0, HERE)
import train_hub_model as T  # noqa: E402
import clinvar_hub as CH  # noqa: E402

OUT = os.path.join(HERE, 'hub_atlas')
EDGE = 15
PLDS = {'FUS': 'P35637_1_285', 'hnRNPA1': 'P09651_187_372', 'TDP-43': 'Q13148_252_414', 'TAF15': 'Q92804_1_234',
        'EWSR1': 'Q01844_1_360', 'DDX4': 'Q9NQI0_1_252'}
AAS = 'ACDEFGHIKLMNPQRSTVWY'
ARO = set('FWY')
_M = None


def _init():
    global _M
    _M = pickle.load(open(CH.MODELS, 'rb'))['consensus']


def _pred(job):
    name, seq = job
    return name, _M.predict(T.features(seq, 10))


def designed():
    """(name, L, spacer pair, sticker, sequence) for the three designed series of the plan."""
    out = []

    def build(L, pair, sticker):
        k = max(1, round(120 / (L + 1)))
        spacer = (pair * L)[:L]
        return (spacer + sticker) * k
    for L in range(2, 17):
        out.append((f'len_{L}', L, 'GS', 'Y', build(L, 'GS', 'Y')))
    for pair in ('GG', 'SS', 'GS', 'NQ', 'PG', 'AA', 'EK'):
        out.append((f'spacer_{pair}', 6, pair, 'Y', build(6, pair, 'Y')))
    for st in 'YFWRKL':
        out.append((f'sticker_{st}', 6, 'GS', st, build(6, 'GS', st)))
    return out


def offsets(seq, sticker):
    pos = [i for i, a in enumerate(seq) if a == sticker]
    return np.array([min(abs(i - p) for p in pos) for i in range(len(seq))])


def main():
    os.makedirs(OUT, exist_ok=True)
    posm = CH.positional_baseline()
    D = CH.idrome()
    have = defaultdict(list)
    for r in csv.DictReader(open(os.path.join(HERE, 'clinvar_hub', 'wt_profiles.csv'))):
        have[r['idr']].append((float(r['pred_consensus']), float(r['pred_position'])))
    todo = [(n, D[n]['seq']) for n in D if n not in have]
    des = designed()
    print(f'IDRome {len(D)} IDRs, {len(have)} already predicted, {len(todo)} to predict; {len(des)} designed', flush=True)
    jobs = [(d[0], d[4]) for d in des] + todo
    with Pool(8, initializer=_init) as pool:
        new = dict(pool.imap_unordered(_pred, jobs, chunksize=16))
    lines = [f'hub atlas (hub_atlas_plan.md), frozen consensus model; excess = prediction - positional baseline; '
             f'residues within {EDGE} of the ends excluded in the designed summaries', '']

    # part 1: designed sequences
    with open(os.path.join(OUT, 'designed.csv'), 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['name', 'L', 'spacer', 'sticker', 'pos', 'aa', 'offset', 'pred', 'excess'])
        summ = {}
        for name, L, pair, st, seq in des:
            p = new[name]
            b = posm.predict(CH.pos_features(len(seq)))
            off = offsets(seq, st)
            for i, a in enumerate(seq):
                w.writerow([name, L, pair, st, i, a, off[i], f'{p[i]:.5f}', f'{p[i] - b[i]:.5f}'])
            keep = (np.arange(len(seq)) >= EDGE) & (np.arange(len(seq)) < len(seq) - EDGE)
            summ[name] = {o: float(np.mean((p - b)[keep & (off == o)])) for o in range(0, 9) if (keep & (off == o)).any()}
    lines.append('(a) GS spacer of length L, Y sticker: mean excess by offset from the nearest sticker (0 = sticker)')
    lines.append('   L  ' + ' '.join(f'{o:>7d}' for o in range(9)))
    for L in range(2, 17):
        s = summ[f'len_{L}']
        lines.append(f'  {L:2d}  ' + ' '.join(f'{s[o]:+.4f}' if o in s else '       ' for o in range(9)))
    lines.append('')
    lines.append('(b) L = 6, Y sticker, spacer pair: excess at offsets 0..3')
    for pair in ('GG', 'SS', 'GS', 'NQ', 'PG', 'AA', 'EK'):
        s = summ[f'spacer_{pair}']
        lines.append(f'  {pair}  ' + ' '.join(f'{s[o]:+.4f}' for o in range(4)))
    lines.append('')
    lines.append("(b') L = 6, GS spacer, sticker: excess at the sticker, offset 1 and offset 3")
    for st in 'YFWRKL':
        s = summ[f'sticker_{st}']
        lines.append(f'  {st}  {s[0]:+.4f} {s[1]:+.4f} {s[3]:+.4f}')
    with open(os.path.join(OUT, 'designed_offsets.csv'), 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['name', 'offset', 'mean_excess'])
        for name, s in summ.items():
            w.writerows([name, o, f'{v:.5f}'] for o, v in s.items())

    # part 2: human IDRome
    hub_n, all_n = defaultdict(int), defaultdict(int)
    grid = defaultdict(list)
    n_res = 0
    with open(os.path.join(OUT, 'idrome_new_profiles.csv'), 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['idr', 'pos', 'aa', 'pred_consensus', 'pred_position'])
        for n in sorted(D):
            seq = D[n]['seq']
            if n in have:
                p = np.array([x[0] for x in have[n]])
                b = np.array([x[1] for x in have[n]])
            else:
                p = new[n]
                b = posm.predict(CH.pos_features(len(seq)))
                w.writerows([n, i, a, f'{p[i]:.5f}', f'{b[i]:.5f}'] for i, a in enumerate(seq))
            assert len(p) == len(seq), n
            k = max(1, round(0.1 * len(seq)))
            top = set(np.lexsort((np.arange(len(seq)), -p))[:k])
            for i, a in enumerate(seq):
                all_n[a] += 1
                if i in top:
                    hub_n[a] += 1
                win = seq[max(0, i - 4):i] + seq[i + 1:i + 5]
                grid[(win.count('G'), min(3, sum(c in ARO for c in win)))].append(p[i] - b[i])
            n_res += len(seq)
            if n in PLDS.values():
                with open(os.path.join(OUT, f'pld_{n}.csv'), 'w', newline='') as f2:
                    w2 = csv.writer(f2)
                    w2.writerow(['pos', 'aa', 'pred', 'baseline'])
                    w2.writerows([i + 1, a, f'{p[i]:.5f}', f'{b[i]:.5f}'] for i, a in enumerate(seq))
    tot_h, tot = sum(hub_n.values()), sum(all_n.values())
    enr = {a: np.log2((hub_n[a] / tot_h) / (all_n[a] / tot)) for a in AAS if all_n[a]}
    with open(os.path.join(OUT, 'enrichment.csv'), 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['aa', 'n_all', 'n_hub', 'log2_enrichment'])
        w.writerows([a, all_n[a], hub_n[a], f'{enr[a]:.4f}'] for a in AAS)
    with open(os.path.join(OUT, 'context_grid.csv'), 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['n_gly_pm4', 'n_aro_pm4', 'n', 'mean_excess'])
        for (g, ar), v in sorted(grid.items()):
            w.writerow([g, ar, len(v), f'{np.mean(v):.5f}'])
    lines += ['', f'(c) IDRome: {len(D)} IDRs, {n_res} residues; log2 enrichment among predicted hubs (top 10% per IDR)',
              '  ' + '  '.join(f'{a} {enr[a]:+.2f}' for a in sorted(enr, key=lambda a: -enr[a])), '',
              '(d) mean excess by glycine (rows, 0-8) and aromatic (columns, 0-3+) residues within +-4; n >= 500']
    for g in range(9):
        row = []
        for ar in range(4):
            v = grid.get((g, ar), [])
            row.append(f'{np.mean(v):+.4f}' if len(v) >= 500 else '   --  ')
        lines.append(f'  G{g}  ' + ' '.join(row))
    open(os.path.join(OUT, 'summary.txt'), 'w').write('\n'.join(lines) + '\n')
    print('\n'.join(lines))


if __name__ == '__main__':
    main()
