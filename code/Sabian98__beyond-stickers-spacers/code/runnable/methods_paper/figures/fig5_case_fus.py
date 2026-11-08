"""
fig5_case_fus.py — main-text case study of the FUS N-terminal low-complexity domain (PED00494, 270 residues).
a  one model of the experimental ensemble with its contact network (light grey) and the hubs of that snapshot
   (orange), full view with a box around residues 205-240 (render_full.png, render_full_region.png)
b  the hand-framed close-up of residues 217-233 (render_zoom_hand.png, or render_zoom.png)
c  hub occupancy along the whole chain in experiment, both simulations and the frozen model (7-residue mean for display;
   rho on raw profiles), glycine runs shaded and aromatic residues as ticks
The gallery, the hub raster and the residue 150-240 zoom are in figS_case_fus.py.
  ~/.venvs/bender/bin/python fig5_case_fus.py  ->  fig5_case_fus.pdf / .png
"""
import csv
import os

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import ConnectionPatch
from scipy.stats import spearmanr

import figstyle as fs
from figS_case_fus import CMAP, D, HUB, STICK, TRACKS, g_runs, show_image, smooth


def main():
    r = list(csv.DictReader(open(os.path.join(D, 'case_fus_residues.csv'))))
    seq = ''.join(x['aa'] for x in r)
    n = len(seq)
    V = {k: np.array([float(x[k]) for x in r]) for k in r[0] if k not in ('pos', 'aa')}
    x = np.arange(1, n + 1)
    fig = plt.figure(figsize=(7.2, 2.7))
    gs = fig.add_gridspec(2, 1, height_ratios=[1.75, 1], hspace=0.3, left=0.09, right=0.99, top=0.93, bottom=0.13)

    # A: full conformer with a box around the zoom region, close-up alongside
    full_p, mask_p = os.path.join(D, 'render_full.png'), os.path.join(D, 'render_full_region.png')
    zoom_p = os.path.join(D, 'render_zoom_hand.png')
    if not os.path.exists(zoom_p):
        zoom_p = os.path.join(D, 'render_zoom.png')
    gA = gs[0].subgridspec(1, 2, width_ratios=[0.8, 1.2], wspace=0.06)
    gB = gA[0, 1].subgridspec(2, 1, height_ratios=[0.17, 1], hspace=0.28)
    axf, axs, axz = fig.add_subplot(gA[0, 0]), fig.add_subplot(gB[0]), fig.add_subplot(gB[1])
    full, mask = plt.imread(full_p), plt.imread(mask_p)
    ys, xs = np.nonzero(full[..., 3] > 0.02)
    y0, y1, x0, x1 = max(0, ys.min() - 12), ys.max() + 12, max(0, xs.min() - 12), xs.max() + 12
    full, mask = full[y0:y1, x0:x1], mask[y0:y1, x0:x1]
    my, mx = np.nonzero(mask[..., 3] > 0.02)
    bx0, bx1, by0, by1 = mx.min() - 25, mx.max() + 25, my.min() - 25, my.max() + 25
    axf.imshow(full)
    axf.set_axis_off()
    axf.add_patch(plt.Rectangle((bx0, by0), bx1 - bx0, by1 - by0, fill=False, ec='#333333', lw=0.8, ls=(0, (3, 2))))
    show_image(axz, zoom_p, '', alpha_min=0.4)
    for ya, yb in ((by0, 0.0), (by1, 1.0)):
        fig.add_artist(ConnectionPatch(xyA=(bx1, ya), coordsA=axf.transData, xyB=(0.0, 1.0 - yb),
                                       coordsB=axz.transAxes, color='#777777', lw=0.5, ls=(0, (3, 2))))
    fs.panel(axf, 'a', x=0.0, y=1.0, title='one model of the experimental ensemble')
    # sequence track of the close-up, styled like the linear sequence plots of IDR papers: plain letters, a continuous
    # hub-occupancy track below them on the colour scale of the tube, triangles above the hubs of this model and ticks
    # every five residues
    lo_s, hi_s = 217, 233
    hubs_model = {int(h) for h in open(os.path.join(D, 'highlight_hubs.txt')).read().split()}
    xs_s = np.arange(lo_s, hi_s + 1)
    occ_s = V['occ_experiment'][lo_s - 1:hi_s]
    axs.imshow(occ_s[None, :], cmap=CMAP, vmin=0, vmax=0.25, aspect='auto', extent=(lo_s - 0.5, hi_s + 0.5, -0.05, 0.3))
    axs.add_patch(plt.Rectangle((lo_s - 0.5, -0.05), hi_s - lo_s + 1, 0.35, fill=False, ec='0.4', lw=0.4))
    for i in xs_s:
        axs.text(i, 0.62, seq[i - 1], ha='center', va='center', fontsize=6.5,
                 color=STICK if seq[i - 1] in 'FWY' else '#b2182b' if seq[i - 1] in 'RK' else 'black')
        if i in hubs_model:
            axs.plot(i, 1.08, marker='v', ms=4.2, color=HUB, mec='white', mew=0.3, clip_on=False)
    for t in range(220, hi_s + 1, 5):
        axs.plot([t, t], [-0.05, -0.17], color='0.3', lw=0.5, clip_on=False)
        axs.text(t, -0.22, str(t), ha='center', va='top', fontsize=5.3)
    axs.set_xlim(lo_s - 0.6, hi_s + 0.6)
    axs.set_ylim(-0.05, 1.0)
    axs.set_axis_off()
    fs.panel(axs, 'b', x=0.0, y=1.45, title='residues 217-233, hubs of this model marked by triangles')
    cax = axf.inset_axes([-0.13, 0.86, 0.3, 0.04])
    cb = fig.colorbar(plt.cm.ScalarMappable(cmap=CMAP, norm=plt.Normalize(0, 0.25)), cax=cax, orientation='horizontal')
    cb.set_ticks([0, 0.25])
    cb.ax.tick_params(labelsize=5.8, length=2)
    cb.set_label('hub occupancy (tube)', fontsize=6, labelpad=1)

    # d: full-length profiles
    ax = fig.add_subplot(gs[1])
    for k, lab, col in TRACKS:
        rho = '' if k == 'occ_experiment' else f'  (ρ = {spearmanr(V[k], V["occ_experiment"]).correlation:.2f})'
        ax.plot(x, smooth(V[k]), color=col, lw=1.3 if k == 'occ_experiment' else 0.9, label=lab + rho)
    tr = ax.get_xaxis_transform()
    for i, c in enumerate(seq):
        if c in 'FWY':
            ax.plot([i + 1, i + 1], [0.0, 0.07], color=STICK, lw=0.6, transform=tr)
    for a_, b_ in g_runs(seq):
        ax.axvspan(a_ - 0.5, b_ + 0.5, color='#7f7f7f', alpha=0.12, lw=0, zorder=0)
    ax.set_xlim(1, n)
    ax.set_xlabel('residue')
    ax.set_ylabel('hub occupancy')
    ax.legend(frameon=False, fontsize=5.5, ncol=4, loc='lower right', bbox_to_anchor=(1.0, 0.98), handlelength=1.2,
              columnspacing=0.8, borderaxespad=0.0)
    fs.panel(ax, 'c', x=-0.085, y=1.06, title='hub occupancy along the chain')
    fs.save(fig, 'fig5_case_fus')


if __name__ == '__main__':
    main()
