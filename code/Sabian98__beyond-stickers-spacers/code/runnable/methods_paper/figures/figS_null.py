"""
figS_null.py — Supplementary figure (2026-10-14): the two null-model panels formerly in Fig. 2 (d, e of
fig2_flicker_null.py), drawn with the same code and inputs; no new statistics.
  a  mean position-adjusted z per residue type in residue-neutral nulls against CALVADOS-2 bead diameter, size-scaled
     against fixed 0.8 nm cutoff (sim_null_by_type_all.csv)
  b  excess of each force field over the matched null per residue type, fixed cutoff, 907 sequences simulated with all
     three force fields (sim_null_v3_by_type.csv)
  ~/.venvs/finches/bin/python figS_null.py  ->  figS_null.pdf / .png
"""
import csv
import os

import matplotlib.pyplot as plt
import numpy as np

import figstyle as fs

CAL_TABLE = os.path.expanduser((__import__('os').environ['BSS_DATA'] + '/Scripts/CALVADOS IDP/input/residues_CALVADOS2.csv'))
FFS = ['CALVADOS-2', 'MPIPI-GG', 'Mpipi-Recharged']


def main():
    fig = plt.figure(figsize=(7.2, 1.75))
    gs = fig.add_gridspec(1, 2, width_ratios=[1, 2], wspace=0.25, left=0.07, right=0.98, top=0.86, bottom=0.2)
    # C: null bias vs bead size, fixed vs size rule (CALVADOS-2 bead sizes)
    axc = fig.add_subplot(gs[0, 0])
    sig = {r['one']: float(r['sigmas']) for r in csv.DictReader(open(CAL_TABLE))}
    tab = list(csv.DictReader(open(os.path.join(fs.EV, 'sim_null_by_type_all.csv'))))
    for rule, lab in (('size', 'cutoff scaled to bead size'), ('fixed', 'fixed 0.8 nm cutoff')):
        pts = [(sig[r['aa']], float(r['z_null']), r['aa']) for r in tab if r['rule'] == rule and r['force_field'] == 'CALVADOS-2']
        x, y = np.array([p[0] for p in pts]), np.array([p[1] for p in pts])
        col = fs.RULE_COLORS[rule]
        axc.scatter(x, y, s=10, color=col, label=f'{lab} (r = {np.corrcoef(x, y)[0, 1]:+.2f})', zorder=3)
        k, b = np.polyfit(x, y, 1)
        xx = np.linspace(x.min(), x.max(), 10)
        axc.plot(xx, k * xx + b, color=col, lw=0.8)
        OFFS = {'G': (3, 1), 'F': (-8, 2), 'K': (-8, -6), 'Y': (-3, 5), 'R': (4, -3)}   # keep the clustered labels apart
        if rule == 'size':
            for xi, yi, a in pts:
                if a in 'GASWYFRK':
                        axc.annotate(a, (xi, yi), xytext=OFFS.get(a, (-7, 1)), textcoords='offset points',
                                 fontsize=5.5, color=col)
    axc.axhline(0, color='0.6', lw=0.5)
    axc.set_ylim(-1.45, 2.25)
    axc.set_yticks([-1, -0.5, 0, 0.5, 1, 1.5])
    axc.set_xlabel('bead diameter σ (nm)')
    axc.set_ylabel('null mean z')
    fs.panel(axc, 'a', x=-0.26, title='size-scaled cutoffs create hubs in residue-neutral chains')
    axc.legend(frameon=False, loc='upper left', handletextpad=0.1, borderaxespad=0.1)

    # D: excess over the null, fixed rule, three force fields
    axd = fig.add_subplot(gs[0, 1])
    # the three force fields on the same sequences (the paired set also simulated with Mpipi-Recharged), v3 labels
    tab_d = list(csv.DictReader(open(os.path.join(fs.EV, 'sim_null_v3_by_type.csv'))))
    ex = {(r['force_field'], r['aa']): float(r['excess']) for r in tab_d if r['rule'] == 'fixed'}
    aas = sorted('ACDEFGHIKLMNPQRSTVWY', key=lambda a: np.mean([ex[(f, a)] for f in FFS]))
    x = np.arange(20)
    for j, ff in enumerate(FFS):
        axd.scatter(x + (j - 1) * 0.22, [ex[(ff, a)] for a in aas], s=11, color=fs.FF_COLORS[ff],
                    label={'MPIPI-GG': 'Mpipi-GG'}.get(ff, ff), zorder=3)
    axd.axhline(0, color='0.6', lw=0.5)
    axd.set_xticks(x)
    axd.set_xticklabels(aas)
    fs.color_ticks(axd)
    axd.set_ylabel('simulation − null')
    fs.panel(axd, 'b', x=-0.1, title='what each force field adds over the null (fixed cutoff)')
    axd.legend(frameon=False, loc='upper left', handletextpad=0.1)

    fs.save(fig, 'figS_null')


if __name__ == '__main__':
    main()
