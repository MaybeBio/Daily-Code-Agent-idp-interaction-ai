"""
figS_nmr_maps.py — Supplementary figure: Fig. 5f-style panels for the nine other pairs of the NMR case series (HMGA1a is Fig. 5f),
the five of the original plan and the four of Amendment 2 (2026-10-15, figS_nmr_maps_prep2.py). For each
protein, the change in contact probability on phosphorylation (mean of three phosphorylated minus mean of three
unmodified CALVADOS-2 replicates) above the diagonal and the difference between the two most different unmodified
replicates (noise) below it, with the measured amide CSP (orange) and the predicted contact rewiring E (blue), each
scaled to its maximum, along the top, and the phosphosites in red. Maps from cache/nmr_maps_<protein>.npz
(figS_nmr_maps_prep.py); CSP, E and the pre-registered partial correlation from redo_evidence/phospho_nmr
(per_residue.csv, per_protein.csv). Descriptive; no new statistic.
  ~/.venvs/finches/bin/python figS_nmr_maps.py  ->  figS_nmr_maps.pdf / .png
"""
import os
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import figstyle as fs

CACHE = os.path.join(fs.HERE, 'cache')
NAMES = {'PAGE4': 'PAGE4', 'p53TAD': 'p53 TAD', 'Sic1': 'Sic1', 'NHE1tail': 'NHE1 tail', 'HMGA1a': 'HMGA1a',
         'H1-0CTD': 'H1.0 CTD', 'aSynS129': 'αSyn pS129†', 'aSynY125': 'αSyn pY125†', 'CD79b': 'CD79b†', 'Tir': 'Tir†'}
# original plan by partial rho (HMGA1a is Fig. 5f), then the four pairs of Amendment 2 (2026-10-15)
ORDER = ['NHE1tail', 'Sic1', 'H1-0CTD', 'PAGE4', 'p53TAD', 'aSynS129', 'aSynY125', 'CD79b', 'Tir']
LIM = 0.12


def main():
    A2 = os.path.join(fs.EV, 'phospho_nmr_amend2', 'results_cal')
    per = pd.concat([pd.read_csv(os.path.join(fs.EV, 'phospho_nmr', 'per_residue.csv')),
                     pd.read_csv(os.path.join(A2, 'per_residue.csv'))])
    pp = pd.concat([pd.read_csv(os.path.join(fs.EV, 'phospho_nmr', 'per_protein.csv')).set_index('protein'),
                    pd.read_csv(os.path.join(A2, 'per_protein.csv')).set_index('protein')])
    fig = plt.figure(figsize=(7.2, 7.3))
    outer = fig.add_gridspec(3, 3, hspace=0.32, wspace=0.3, left=0.06, right=0.98, top=0.96, bottom=0.11)
    im = None
    for k, prot in enumerate(ORDER):
        mp = np.load(os.path.join(CACHE, f'nmr_maps_{prot}.npz'))
        P, W = mp['pP_reps'], mp['pW_reps']
        n = P.shape[1]
        dp = P.mean(0) - W.mean(0)
        a_, b_ = max(((0, 1), (0, 2), (1, 2)), key=lambda q: np.abs(W[q[0]] - W[q[1]]).max())
        half = np.where(np.triu(np.ones((n, n), bool), 1).T, dp, W[a_] - W[b_])
        sites = [int(re.sub(r'\D', '', s)) for s in pp.loc[prot, 'sites'].split()]
        g = outer[k // 3, k % 3].subgridspec(2, 1, height_ratios=[0.28, 1], hspace=0.1)
        ax, tr = fig.add_subplot(g[1]), fig.add_subplot(g[0])
        im = ax.imshow(half, origin='lower', cmap='PRGn', vmin=-LIM, vmax=LIM, extent=(0.5, n + 0.5, 0.5, n + 0.5),
                       interpolation='nearest', aspect='auto')
        ax.plot([0.5, n + 0.5], [0.5, n + 0.5], color='0.75', lw=0.4)
        for s in sites:
            ax.plot([s], [0.5], marker='^', ms=2.4, color='#b2182b', clip_on=False, zorder=5)
        ax.set_xlim(0.5, n + 0.5)
        ax.set_ylim(0.5, n + 0.5)
        ax.tick_params(labelsize=5.5, length=1.5, pad=1)
        if k % 3 == 0:
            ax.set_ylabel('residue', fontsize=6, labelpad=1)
        if k // 3 == 2:
            ax.set_xlabel('residue', fontsize=6, labelpad=1)
        d = per[per.protein == prot].sort_values('pos')
        tr.bar(d.pos, d.CSP / d.CSP.max(), width=0.9, color='#d95f02', lw=0)
        pos, e = d.pos.to_numpy().astype(float), d.E.to_numpy().astype(float)
        brk = np.where(np.diff(pos) > 6)[0] + 1
        tr.plot(np.insert(pos, brk, np.nan), np.insert(e / np.nanmax(e), brk, np.nan), color='#1b6ca8', lw=0.8)
        for s in sites:
            tr.axvline(s, color='#b2182b', lw=0.5)
        tr.set_xlim(0.5, n + 0.5)
        tr.set_ylim(0, 1.15)
        tr.set_axis_off()
        r = pp.loc[prot]
        fs.panel(tr, 'abcdefghi'[k], x=-0.08, y=1.12,
                 title=f'{NAMES[prot]}, {len(sites)} site{"s" if len(sites) > 1 else ""}, partial ρ = {r.rho_p:+.2f}')
        print(f'{prot}: N {n}, sites {sites}, max|dp| {np.abs(dp).max():.3f}, max|noise| {np.abs(W[a_] - W[b_]).max():.3f}')
    cax = fig.add_axes([0.08, 0.035, 0.28, 0.012])
    cb = fig.colorbar(im, cax=cax, orientation='horizontal', ticks=[-0.1, 0, 0.1])
    cb.set_label('change in contact probability', fontsize=6, labelpad=1)
    cb.ax.tick_params(labelsize=5.5, length=1.5)
    cb.outline.set_linewidth(0.4)
    fig.text(0.42, 0.05, 'above the diagonal: phosphorylated − unmodified\nbelow the diagonal: unmodified replicate − replicate (noise)',
             fontsize=5.6, va='top', color='#333333')
    fig.text(0.78, 0.05, 'measured CSP', fontsize=5.8, color='#C0682F', va='top')
    fig.text(0.78, 0.035, 'predicted rewiring, E', fontsize=5.8, color='#2F7DB5', va='top')
    fig.text(0.78, 0.02, 'phosphosite', fontsize=5.8, color='#b2182b', va='top')
    fig.text(0.42, 0.018, '† added later, chosen before the runs', fontsize=5.4, color='#555555', va='top')
    fs.save(fig, 'figS_nmr_maps')


if __name__ == '__main__':
    main()
