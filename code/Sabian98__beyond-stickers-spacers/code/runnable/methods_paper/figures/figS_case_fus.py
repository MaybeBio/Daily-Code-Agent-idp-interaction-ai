"""
figS_case_fus.py — supplementary panels of the FUS N-terminal domain case study (PED00494, 270 residues). The model
with its close-up and the full-length profiles are in main-text Fig. 5 (fig5_case_fus.py), which imports
the helpers here.
a  gallery of 10 models of the experimental ensemble, coloured by hub occupancy (render_ensemble.png)
b  zoom on residues 150-240: hub occupancy in experiment (bars), simulations and prediction (lines)
c  hubs of 150 of the 551 models ordered by radius of gyration (case_fus/raster_hubs.csv), glycine runs shaded
  ~/.venvs/bender/bin/python figS_case_fus.py  ->  figS_case_fus.pdf / .png
"""
import csv
import os

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import ConnectionPatch
from scipy.stats import spearmanr

import figstyle as fs

D = os.path.join(fs.HERE, 'case_fus')
CMAP = LinearSegmentedColormap.from_list('occ', ['#f2f2f2', '#c6dbef', '#6baed6', '#2171b5', '#08306b'])
HUB = '#E69F00'
STICK = '#c51b7d'   # tyrosine / aromatic stickers (magenta, distinct from the blue occupancy scale)
TRACKS = [('occ_experiment', 'experiment (PED00494)', '#222222'), ('occ_calvados2', 'CALVADOS-2', fs.FF_COLORS['CALVADOS-2']),
          ('occ_recharged', 'Mpipi-Recharged', fs.FF_COLORS['Mpipi-Recharged']), ('pred_consensus', 'sequence model', HUB)]
ROLE = {**dict.fromkeys('FWY', STICK), 'G': '#6f6f6f', 'S': '#a9a9a9', **dict.fromkeys('DEKR', '#b2182b')}


def smooth(x, w=3):
    return np.array([x[max(0, i - w):i + w + 1].mean() for i in range(len(x))])


def show_image(ax, path, text, alpha_min=0.02):
    ax.set_axis_off()
    if not os.path.exists(path):
        ax.add_patch(plt.Rectangle((0, 0), 1, 1, transform=ax.transAxes, fc='#f2f2f2', ec='#bbbbbb', lw=0.6))
        ax.text(0.5, 0.5, text, ha='center', va='center', transform=ax.transAxes, fontsize=6.5, color='#777777')
        return
    img = plt.imread(path)
    if img.shape[-1] == 4:
        mask = img[..., 3] > alpha_min
    else:
        mask = img[..., :3].min(-1) < 0.97
    ys, xs = np.nonzero(mask)
    pad = 12
    img = img[max(0, ys.min() - pad):ys.max() + pad, max(0, xs.min() - pad):xs.max() + pad]
    ax.imshow(img)


def g_runs(seq, k=4):
    runs, i = [], 0
    while i < len(seq):
        j = i
        while j < len(seq) and seq[j] == 'G':
            j += 1
        if j - i >= k:
            runs.append((i + 1, j))
        i = max(j, i + 1)
    return runs


def main():
    r = list(csv.DictReader(open(os.path.join(D, 'case_fus_residues.csv'))))
    seq = ''.join(x['aa'] for x in r)
    n = len(seq)
    V = {k: np.array([float(x[k]) for x in r]) for k in r[0] if k not in ('pos', 'aa')}
    x = np.arange(1, n + 1)
    fig = plt.figure(figsize=(7.2, 4.5))
    gs = fig.add_gridspec(2, 3, width_ratios=[1, 1, 0.82], height_ratios=[1.35, 1], wspace=0.38, hspace=0.45,
                          left=0.06, right=0.98, top=0.95, bottom=0.08)

    # B: gallery + colour bar
    ax = fig.add_subplot(gs[0, :2])
    show_image(ax, os.path.join(D, 'render_ensemble.png'), 'ensemble gallery (ChimeraX)')
    fs.panel(ax, 'a', x=0.0, y=1.0, title='ten models of the experimental ensemble')
    cax = ax.inset_axes([0.0, -0.02, 0.24, 0.04])
    cb = fig.colorbar(plt.cm.ScalarMappable(cmap=CMAP, norm=plt.Normalize(0, 0.25)), cax=cax, orientation='horizontal')
    cb.set_ticks([0, 0.125, 0.25])
    cb.ax.tick_params(labelsize=5.8, length=2)
    cb.set_label('hub occupancy (tubes, bars)', fontsize=6, labelpad=1)

    # D: zoom on residues 150-240 (right of the gallery)
    lo, hi = 150, 240
    ax = fig.add_subplot(gs[0, 2])
    xs = np.arange(lo, hi + 1)
    sl = slice(lo - 1, hi)
    ax.bar(xs, V['occ_experiment'][sl], width=0.85, color=[CMAP(min(v / 0.25, 1)) for v in V['occ_experiment'][sl]], lw=0)
    for k, lab, col in TRACKS[1:]:
        ax.plot(xs, smooth(V[k], 1)[sl], color=col, lw=0.8)
    for a_, b_ in g_runs(seq):
        if b_ >= lo and a_ <= hi:
            ax.axvspan(max(a_, lo) - 0.5, min(b_, hi) + 0.5, color='#7f7f7f', alpha=0.12, lw=0, zorder=0)
    tr = ax.get_xaxis_transform()
    for i in xs:
        if seq[i - 1] in 'FWY':
            ax.plot([i, i], [0, 0.07], color=STICK, lw=0.9, transform=tr)
    ax.set_xlim(lo - 0.6, hi + 0.6)
    ax.set_xticks(range(lo, hi + 1, 30))
    ax.set_ylim(0, ax.get_ylim()[1] * 1.32)                  # room for the key above the bars
    ax.set_ylabel('hub occupancy')
    ax.set_xlabel('residue')
    fs.panel(ax, 'b', x=-0.32, title='residues 150-240')
    ax.text(0.02, 0.98, 'bars, experiment\nlines, simulations and model\nticks, Y   gray, G runs', transform=ax.transAxes,
            va='top', fontsize=5.4, color='#444444')

    # c: hubs of 150 of the 551 models ordered by radius of gyration (moved from Fig. 5 on 2026-10-11)
    ras = list(csv.DictReader(open(os.path.join(D, 'raster_hubs.csv'))))
    axr = fig.add_subplot(gs[1, :])
    order = sorted(range(len(ras)), key=lambda k: float(ras[k]['rg_A']))
    show = [order[int(q)] for q in np.linspace(0, len(order) - 1, 150)]
    for row, k in enumerate(show):
        hs = [int(h) for h in ras[k]['hubs'].split()]
        axr.scatter(hs, [row] * len(hs), s=0.9, marker='s', color=HUB, lw=0)
    for a_, b_ in g_runs(seq):
        axr.axvspan(a_ - 0.5, b_ + 0.5, color='#7f7f7f', alpha=0.12, lw=0, zorder=0)
    axr.set_xlim(1, n)
    axr.set_ylim(-1, len(show))
    axr.set_yticks([0, len(show) - 1], ['compact', 'expanded'])
    axr.tick_params(axis='y', labelsize=6, length=0)
    axr.set_xlabel('residue')
    fs.panel(axr, 'c', x=-0.05, y=1.04, title='hubs of 150 models ordered by radius of gyration, glycine runs shaded')
    fs.save(fig, 'figS_case_fus')


if __name__ == '__main__':
    main()
