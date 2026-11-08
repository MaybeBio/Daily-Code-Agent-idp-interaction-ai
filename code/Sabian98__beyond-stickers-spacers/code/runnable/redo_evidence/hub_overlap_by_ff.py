"""
hub_overlap_by_ff.py — how often the force fields pick the same residues of a sequence as hubs.
  For every protein labelled in two force fields, take the 10% of residues with the highest hub occupancy in each
  (ties broken at random, seed 0) and report the share of one force field's hubs that are also hubs in the other
  (chance 0.10), the same share for hubs in the middle 70% of the chain, and the share with a hub of the other force
  field within two residues, against random hub sets of the same size (5 draws per protein).
  ~/.venvs/bender/bin/python hub_overlap_by_ff.py  ->  hub_overlap_by_ff_summary.txt
Planned use (Results 3.3, at the Recharged update): "For a given sequence, a third of the hubs of one force field are
also hubs in another, against a tenth expected by chance, and two thirds lie within two residues of one."
"""
import csv
import os
from collections import defaultdict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SETS = {'CALVADOS-2': 'hub_frames_v3_cal', 'Mpipi-GG': 'hub_frames_v3_mpipi', 'Mpipi-Recharged': 'hub_frames_v3_recharged'}


def load(d):
    done = {r['UniProt_ID']: int(r['N']) for r in csv.DictReader(open(os.path.join(HERE, d, 'per_protein.csv')))}
    occ = defaultdict(dict)
    with open(os.path.join(HERE, d, 'per_residue.csv')) as fh:
        for r in csv.DictReader(fh):
            if r['UniProt_ID'] in done and r.get('occupancy', '') != '':
                occ[r['UniProt_ID']][int(r['pos'])] = float(r['occupancy'])
    return {u: np.array([v[i] for i in sorted(v)]) for u, v in occ.items() if len(v) == done[u]}


def main():
    rng = np.random.default_rng(0)
    L = {ff: load(d) for ff, d in SETS.items()}

    def top(x):
        k = max(1, int(round(0.1 * len(x))))
        return np.lexsort((rng.random(len(x)), -x))[:k]

    def near(ha, hb, w=2):
        hb = np.asarray(hb)
        return np.mean([np.any(np.abs(hb - i) <= w) for i in ha])

    out = ['Hub overlap between force fields (top 10% of residues per protein; chance for the exact overlap is 0.10)']
    for a, b in (('CALVADOS-2', 'Mpipi-Recharged'), ('CALVADOS-2', 'Mpipi-GG'), ('Mpipi-GG', 'Mpipi-Recharged')):
        same, mid, win, win_rand = [], [], [], []
        for u in sorted(set(L[a]) & set(L[b])):
            x, y = L[a][u], L[b][u]
            if len(x) != len(y) or len(x) < 20:
                continue
            ha, hb = top(x), top(y)
            same.append(len(set(ha) & set(hb)) / len(ha))
            n = len(x)
            hm = [i for i in ha if 0.15 * n <= i < 0.85 * n]
            if hm:
                mid.append(len(set(hm) & set(hb)) / len(hm))
            win.append(near(ha, hb))
            win_rand.append(np.mean([near(ha, rng.choice(n, len(hb), replace=False)) for _ in range(5)]))
        out.append(f'{a} vs {b}: n = {len(same)}; same residue {np.median(same):.2f} (mean {np.mean(same):.2f}); '
                   f'middle 70% {np.median(mid):.2f}; within 2 residues {np.median(win):.2f} (random {np.median(win_rand):.2f})')
    open(os.path.join(HERE, 'hub_overlap_by_ff_summary.txt'), 'w').write('\n'.join(out) + '\n')
    print('\n'.join(out))


if __name__ == '__main__':
    main()
