"""
arrangement_ff_table.py — force-field rows of Supplementary Table tab:arrangement, with the functions of
arrangement_test.py (matched-position Spearman: raw, position-partial, and after removing the position-adjusted mean
of each residue type; median over proteins).
Upper block: every protein labelled in both force fields, type means from all proteins of each force field.
Lower block: the sequences labelled in all three force fields, type means from that subset, and the correlation
between the 20 residue-type means of each pair.
  ~/.venvs/bender/bin/python arrangement_ff_table.py  ->  arrangement_ff_table.txt
"""
import os

import numpy as np

import arrangement_test as A

HERE = os.path.dirname(os.path.abspath(__file__))
AAS = 'ACDEFGHIKLMNPQRSTVWY'


def main():
    D = {'CALVADOS-2': A.read('hub_frames_v3_cal'), 'Mpipi-GG': A.read('hub_frames_v3_mpipi'),
         'Mpipi-Recharged': A.read('hub_frames_v3_recharged')}
    tm = {k: A.type_means(v) for k, v in D.items()}
    lines = ['Force-field rows of tab:arrangement (median per-protein Spearman: raw / position removed / identity also '
             'removed)', '', 'Upper block: all proteins labelled in both force fields']
    pairs = [('CALVADOS-2', 'Mpipi-GG'), ('CALVADOS-2', 'Mpipi-Recharged'), ('Mpipi-GG', 'Mpipi-Recharged')]
    for a, b in pairs[:2]:
        keys = sorted(set(D[a]) & set(D[b]))
        r, pa, re = A.agreement(D[a], D[b], [(u, u) for u in keys], tm[a], tm[b])
        lines.append(f'  {a} versus {b}: n = {len(r)}  {np.nanmedian(r):.2f}  {np.nanmedian(pa):.2f}  {np.nanmedian(re):.2f}')
    three = sorted(set(D['CALVADOS-2']) & set(D['Mpipi-GG']) & set(D['Mpipi-Recharged']))
    sub = {k: {u: v[u] for u in three} for k, v in D.items()}
    tms = {k: A.type_means(v) for k, v in sub.items()}
    lines += ['', f'Lower block: the {len(three)} sequences labelled in all three force fields']
    for a, b in pairs:
        r, pa, re = A.agreement(sub[a], sub[b], [(u, u) for u in three], tms[a], tms[b])
        x = [tms[a][t] for t in AAS if t in tms[a] and t in tms[b]]
        y = [tms[b][t] for t in AAS if t in tms[a] and t in tms[b]]
        lines.append(f'  {a} versus {b}: n = {len(r)}  {np.nanmedian(r):.2f}  {np.nanmedian(pa):.2f}  '
                     f'{np.nanmedian(re):.2f}  type means r = {np.corrcoef(x, y)[0, 1]:.2f}')
    print('\n'.join(lines))
    open(os.path.join(HERE, 'arrangement_ff_table.txt'), 'w').write('\n'.join(lines) + '\n')


if __name__ == '__main__':
    main()
