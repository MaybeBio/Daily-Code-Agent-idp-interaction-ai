"""
figS_predictor.py — Supplementary figure with the two predictor panels that left the main Fig. 3 in the restructure of
2026-10-04 (formerly Fig. 3B and 3D of fig3_predictor.py, unchanged in content):
  A  One held-out protein near the median score: simulated occupancy, prediction, position-only baseline.
  B  Change in predicted occupancy when the centre residue of four 100-residue hosts is replaced (blind predictions
     frozen on 2026-10-02, consensus model), to be tested by the design batch.
  ~/.venvs/bender/bin/python figS_predictor.py  ->  figS_predictor.pdf / .png
"""
import csv
import os
from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import spearmanr

import figstyle as fs

HP = os.path.join(fs.EV, 'hub_predictor')
GV = os.path.join(fs.EV, 'gaivi')


def main():
    fig = plt.figure(figsize=(7.2, 2.5))
    gs = fig.add_gridspec(1, 2, wspace=0.3, left=0.07, right=0.91, top=0.86, bottom=0.2, width_ratios=[1, 1.15])
    # B
    ax = fig.add_subplot(gs[0, 0])
    ex = list(csv.DictReader(open(os.path.join(HP, 'cv_full_example.csv'))))
    pos = np.array([int(r['pos']) for r in ex]) + 1
    t, p1, p2 = (np.array([float(r[k]) for r in ex]) for k in ('occupancy', 'pred_hgb', 'pred_position'))
    ax.plot(pos, t, color='black', lw=1.0, label='simulation (1,000 snapshots)')
    ax.plot(pos, p1, color=fs.FF_COLORS['CALVADOS-2'], lw=1.0,
            label=f'predicted from sequence (ρ = {spearmanr(p1, t).correlation:.2f})')
    ax.plot(pos, p2, color='0.6', lw=0.9, ls='--', label=f'position only (ρ = {spearmanr(p2, t).correlation:.2f})')
    ax.axhline(0.10, color='0.8', lw=0.5)
    ax.set_xlabel('residue')
    ax.set_ylabel('hub occupancy')
    ax.set_xlim(1, pos.max())
    ax.set_ylim(0, max(t.max(), p1.max()) * 1.45)
    fs.panel(ax, 'a', x=-0.13, y=1.05, title=f'Held-out protein near the median ({ex[0]["protein_id"]}, {len(ex)} residues)')
    ax.legend(frameon=False, loc='upper left', ncol=1)

    # D
    ax = fig.add_subplot(gs[0, 1])
    pr = defaultdict(dict)
    for r in csv.DictReader(open(os.path.join(GV, 'predictions_frozen_2026-10-02.csv'))):
        if r['id'].startswith(('HOST_', 'GUEST_')) and r['pos'] == '50':
            pr[r['id']] = float(r['consensus'])
    hosts = [('GS', '(GS)$_{50}$'), ('AVG', 'mean IDR'), ('CHG', 'E/K-rich'), ('PLD', 'prion-like')]
    aas = 'ACDEFGHIKLMNPQRSTVWY'
    M = np.full((20, 4), np.nan)
    for j, (h, _) in enumerate(hosts):
        base = pr[f'HOST_{h}']
        for i, a in enumerate(aas):
            M[i, j] = pr.get(f'GUEST_{h}_{a}50', base) - base
    order = np.argsort(-np.nanmean(M, axis=1))
    lim = np.nanmax(np.abs(M))
    im = ax.imshow(M[order].T, cmap='RdBu_r', vmin=-lim, vmax=lim, aspect='auto')
    ax.set_xticks(range(20))
    ax.set_xticklabels([aas[i] for i in order])
    fs.color_ticks(ax)
    ax.set_yticks(range(4))
    ax.set_yticklabels([h[1] for h in hosts])
    ax.tick_params(length=0)
    for sp in ax.spines.values():
        sp.set_visible(False)
    cb = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
    cb.set_label('Δ predicted occupancy', fontsize=6.5)
    cb.ax.tick_params(labelsize=6)
    cb.ax.yaxis.set_major_formatter(plt.FormatStrFormatter('%+.2f'))
    ax.set_xlabel('residue placed at the center of the host')
    fs.panel(ax, 'b', x=-0.2, y=1.05, title='Predicted effect of one residue depends on context')

    fs.save(fig, 'figS_predictor')


if __name__ == '__main__':
    main()
