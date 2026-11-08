"""
fig_atlas.py — main-text figure (2026-10-11, ALBATROSS Fig. 3/4 and IDRome Figs. 1-3 style): predicted hubs in designed
sticker-spacer sequences and across the human IDRome. Frozen consensus model only (redo_evidence/hub_atlas.py,
hub_atlas_plan.md hashed before any prediction); exploratory and descriptive. Excess = prediction minus the positional
baseline.
  a  GS spacer of length L between tyrosine stickers: mean excess at each offset from the nearest sticker
  b  six-residue spacers of different composition between tyrosines: excess against offset
  c  six-residue GS spacers with different stickers: excess at the sticker and one residue away
  d  log2 enrichment of each residue type among the predicted hubs (top 10% per IDR) of all human IDRs
  e  mean excess by the number of glycine and aromatic residues within four positions, all human IDR residues
  figS_atlas_pld: predicted profiles of six prion-like domains, aromatic residues as ticks and glycine runs shaded
  ~/.venvs/bender/bin/python fig_atlas.py  ->  fig_atlas.pdf / .png, figS_atlas_pld.pdf / .png
"""
import csv
import os
from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np

import figstyle as fs
from figS_case_fus import STICK, g_runs

A = os.path.join(fs.EV, 'hub_atlas')
# ALBATROSS-style blue-yellow-red heat scale (their Fig. 3a), saturated at both ends
HEAT = plt.get_cmap('Spectral_r')
PLDS = [('FUS', 'P35637_1_285'), ('hnRNPA1', 'P09651_187_372'), ('TDP-43', 'Q13148_252_414'),
        ('TAF15', 'Q92804_1_234'), ('EWSR1', 'Q01844_1_360'), ('DDX4', 'Q9NQI0_1_252')]
SPACERS = [('GG', '#00843d', 'o'), ('GS', '#2bb673', 's'), ('SS', '#0096c7', 'D'), ('NQ', '#6a3d9a', '^'),
           ('AA', '#5f5f5f', 'v'), ('PG', '#e7298a', 'P'), ('EK', '#ff7f00', 'X')]


def smooth(x, w=3):
    return np.array([x[max(0, i - w):i + w + 1].mean() for i in range(len(x))])


def main():
    off = defaultdict(dict)
    for r in csv.DictReader(open(os.path.join(A, 'designed_offsets.csv'))):
        off[r['name']][int(r['offset'])] = float(r['mean_excess'])
    fig = plt.figure(figsize=(7.2, 2.05))
    row = fig.add_gridspec(1, 5, width_ratios=[1.12, 1.0, 0.85, 1.3, 0.85], wspace=0.85, left=0.065,
                           right=0.955, top=0.85, bottom=0.25)
    G = np.full((9, 4), np.nan)
    for r in csv.DictReader(open(os.path.join(A, 'context_grid.csv'))):
        if int(r['n']) >= 500:
            G[int(r['n_gly_pm4']), int(r['n_aro_pm4'])] = float(r['mean_excess'])

    # a  spacer length heatmap
    ax = fig.add_subplot(row[0, 0])
    Ls = list(range(2, 17))
    M = np.full((len(Ls), 9), np.nan)
    for i, L in enumerate(Ls):
        for o, v in off[f'len_{L}'].items():
            M[i, o] = v
    im = ax.imshow(M, aspect='auto', cmap=HEAT, vmin=0, vmax=0.025, origin='lower',
                   extent=(-0.5, 8.5, Ls[0] - 0.5, Ls[-1] + 0.5))
    cb = fig.colorbar(im, ax=ax, fraction=0.07, pad=0.03, ticks=[0, 0.01, 0.02], extend='min')
    cb.ax.set_yticklabels(['0', '.01', '.02'])
    cb.ax.tick_params(labelsize=5, length=1.5, pad=1)
    ax.set_xlabel('distance from Y')
    ax.set_ylabel('GS spacer length')
    ax.set_xticks([0, 2, 4, 6, 8])
    ax.set_yticks([2, 6, 10, 14])
    fs.panel(ax, 'a', x=-0.36, title='sticker spacing')

    # b  spacer composition
    ax = fig.add_subplot(row[0, 1])
    for pair, col, mk in SPACERS:
        s = off[f'spacer_{pair}']
        o = sorted(s)
        ax.plot(o, [s[k] for k in o], '-', marker=mk, color=col, ms=3, lw=1.4, mec='white', mew=0.4, label=pair)
    ax.axhline(0, color='0.7', lw=0.5)
    ax.set_xlabel('distance from Y')
    ax.set_ylabel('excess occupancy')
    ax.set_xticks(range(0, 4))
    ax.set_ylim(-0.085, 0.055)
    ax.set_yticks([-0.04, 0, 0.04])
    ax.legend(frameon=False, fontsize=5, ncol=4, handlelength=0.8, columnspacing=0.4, handletextpad=0.3,
              loc='lower center', borderaxespad=0.0)
    fs.panel(ax, 'b', x=-0.4, title='spacer type')

    # c  sticker type
    ax = fig.add_subplot(row[0, 2])
    sts = 'YFWRKL'
    x = np.arange(len(sts))
    ax.bar(x - 0.2, [off[f'sticker_{s}'][0] for s in sts], width=0.38, color='#d6006f', label='sticker')
    ax.bar(x + 0.2, [off[f'sticker_{s}'][1] for s in sts], width=0.38, color='#00a651', label='next residue')
    ax.axhline(0, color='0.6', lw=0.5)
    ax.set_xticks(x, list(sts))
    fs.color_ticks(ax)
    ax.set_ylim(-0.03, 0.022)
    ax.set_yticks([-0.02, 0, 0.02])
    ax.legend(frameon=False, fontsize=5, handlelength=0.8, handletextpad=0.3, loc='lower left', borderaxespad=0.0)
    fs.panel(ax, 'c', x=-0.3, title='sticker type')

    # d  proteome residue enrichment
    ax = fig.add_subplot(row[0, 3])
    enr = {r['aa']: float(r['log2_enrichment']) for r in csv.DictReader(open(os.path.join(A, 'enrichment.csv')))}
    order = sorted(enr, key=lambda a: enr[a])
    cols = ['#d6006f' if a in 'FWY' else '#00a651' if a in 'GS' else '#0072b2' if a in 'DEKR' else '#b5b5b5' for a in order]
    ax.bar(range(20), [enr[a] for a in order], color=cols, width=0.75)
    ax.axhline(0, color='0.5', lw=0.5)
    ax.set_xticks(range(20), order, fontsize=5)
    ax.tick_params(axis='x', length=1.5, pad=1.5)
    ax.set_xlim(-0.7, 19.7)
    ax.set_ylabel('log$_2$ enrichment\namong predicted hubs')
    fs.panel(ax, 'd', x=-0.27, title='human IDRome')

    # e  context grid
    ax = fig.add_subplot(row[0, 4])
    im = ax.imshow(G, aspect='auto', cmap=HEAT, vmin=-0.01, vmax=0.08, origin='lower', extent=(-0.5, 3.5, -0.5, 8.5))
    ax.set_xticks(range(4), ['0', '1', '2', '3+'])
    ax.set_xlabel('F, W, Y within ±4')
    ax.set_ylabel('G within ±4')
    fs.panel(ax, 'e', x=-0.5, title='neighborhood')
    cb = fig.colorbar(im, ax=ax, fraction=0.1, pad=0.05, ticks=[0, 0.04, 0.08])
    cb.ax.set_yticklabels(['0', '.04', '.08'])
    cb.ax.tick_params(labelsize=5, length=1.5, pad=1)
    fs.save(fig, 'fig_atlas')


def plds():
    """Supplementary figure: predicted profiles of six prion-like domains (panel e of hub_atlas_plan.md)."""
    fig = plt.figure(figsize=(7.2, 2.2))
    sub = fig.add_gridspec(2, 3, hspace=0.6, wspace=0.1, left=0.02, right=0.99, top=0.9, bottom=0.1)
    for k, (name, idr) in enumerate(PLDS):
        ax = fig.add_subplot(sub[k // 3, k % 3])
        rows = list(csv.DictReader(open(os.path.join(A, f'pld_{idr}.csv'))))
        seq = ''.join(r['aa'] for r in rows)
        p = np.array([float(r['pred']) for r in rows])
        xx = np.arange(1, len(seq) + 1)
        for a_, b_ in g_runs(seq):
            ax.axvspan(a_ - 0.5, b_ + 0.5, color='#7f7f7f', alpha=0.15, lw=0, zorder=0)
        ax.plot(xx, smooth(p), color='#E69F00', lw=0.9)
        tr = ax.get_xaxis_transform()
        for i, c in enumerate(seq):
            if c in 'FWY':
                ax.plot([i + 1, i + 1], [0, 0.15], color=STICK, lw=0.4, transform=tr)
        first = int(idr.split('_')[1])
        ax.set_xlim(1, len(seq))
        ax.set_yticks([])
        ax.tick_params(axis='x', labelsize=5.5, length=1.5, pad=1)
        ax.set_xticks([1, len(seq)], [str(first), str(first + len(seq) - 1)])
        ax.set_title(name, fontsize=7, loc='left', pad=2)
    fs.save(fig, 'figS_atlas_pld')


if __name__ == '__main__':
    main()
    plds()
