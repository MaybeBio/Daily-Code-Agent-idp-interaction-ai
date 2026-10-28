#!/usr/bin/env python3
"""Figure 2A-D: per-pair concentration decomposition.

A/C = stacked probability densities (combined, IDR-1, IDR-2).
B/D = IDR-1 vs IDR-2 voxel scatter with marginal KDEs.

Layout follows the original panels E-H of Figure 1, from which this figure
was split off.  Data are the trajectory-averaged voxel maps written by
analyze_heterotypic_slabs.py: avg_tot and avg_s1 are already in
residues/voxel, so the two species sum to the total and no unit conversion
or box geometry is needed here.
"""
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.stats import gaussian_kde
from scipy.signal import find_peaks

plt.rcParams.update({
    'font.size': 20, 'axes.linewidth': 2, 'axes.spines.top': False,
    'axes.spines.right': False, 'xtick.major.size': 8, 'xtick.major.width': 2,
    'ytick.major.size': 8, 'ytick.major.width': 2,
    'legend.frameon': False, 'figure.dpi': 300,
})

COLORS = {'comp1': '#a559aa', 'comp2': '#59a89c', 'compalt': '#A89759'}

RESULTS_DIR = Path("../DDM_Heterotypic/results_het")
DATA_DIR = Path("figure2ad_data")

PAIRS = [('Q96K19_200_258_Q3LI60_1_44', [COLORS['compalt'], COLORS['comp2']], 'A', 'B'),
         ('Q8WTT2_1_202_Q3LI60_1_44',   [COLORS['comp1'],   COLORS['comp2']], 'C', 'D')]

PEAK_HEIGHT = 6e-8       # minimum KDE height for an annotated peak
KDE_GRID = 500
MARGINAL_GRID = 200
MARGINAL_FRAC = 0.18     # marginal KDE depth, as a fraction of the scatter box


def save_csv(data_dict, filename):
    """Save data to CSV, padding arrays to equal length."""
    DATA_DIR.mkdir(exist_ok=True)
    max_len = max(len(v) if hasattr(v, '__len__') and not isinstance(v, str) else 1
                  for v in data_dict.values())
    for k, v in data_dict.items():
        if not hasattr(v, '__len__') or isinstance(v, str):
            data_dict[k] = [v] * max_len
        elif len(v) < max_len:
            data_dict[k] = list(v) + [np.nan] * (max_len - len(v))
    pd.DataFrame(data_dict).to_csv(DATA_DIR / filename, index=False)


def closest_5050(pair_id):
    """Path to the composition whose split is nearest 50/50."""
    candidates = []
    for f in RESULTS_DIR.glob(f"{pair_id}_*_*.npz"):
        try:
            c1, c2 = (int(t) for t in f.stem.rsplit("_", 2)[1:])
        except ValueError:
            continue
        candidates.append((abs(c1 / (c1 + c2) - 0.5), f))
    if not candidates:
        raise FileNotFoundError(f"no results_het/*.npz for {pair_id}")
    return min(candidates)[1]


def load_densities(pair_id):
    """Total, species-1 and species-2 voxel concentrations (residues/voxel)."""
    path = closest_5050(pair_id)
    with np.load(path, allow_pickle=True) as npz:
        total = npz['avg_tot'].astype(np.float64).ravel()
        d1 = npz['avg_s1'].astype(np.float64).ravel()
        name1, name2 = str(npz['seq_name1']), str(npz['seq_name2'])
    print(f"{pair_id}: {path.name}, {total.size} voxels")
    return total, d1, total - d1, name1, name2


def plot_pair_composite(fig, ax, d, d1, d2, colors):
    """Panels A, C: stacked probability densities for a pair."""
    x_range = np.linspace(0, d.max(), KDE_GRID)
    kdes = [gaussian_kde(v, bw_method='scott')(x_range) for v in (d, d1, d2)]
    peaks, _ = find_peaks(kdes[0], height=PEAK_HEIGHT)
    peaks = sorted(peaks, key=lambda p: kdes[0][p], reverse=True)

    ax.set_axis_off()
    gs_sub = gridspec.GridSpecFromSubplotSpec(
        3, 1, subplot_spec=ax.get_subplotspec(), hspace=0.4)

    labels = ['Combined', 'IDR-1', 'IDR-2']
    fill_colors = [None, colors[0], colors[1]]
    line_colors = ['k', colors[0], colors[1]]

    for i, (kde_y, fill, line, label) in enumerate(
            zip(kdes, fill_colors, line_colors, labels)):
        ax_kde = fig.add_subplot(gs_sub[i])
        if fill:
            ax_kde.fill_between(x_range, 0, kde_y, color=fill, alpha=0.6)
        ax_kde.plot(x_range, kde_y, color=line, lw=1.5)
        if i == 0:
            for pk in peaks:
                ax_kde.axvline(x_range[pk], color='k', ls='--', lw=1.5, alpha=0.8)
                ax_kde.text(x_range[pk], kde_y.max() * 1.05, f'{x_range[pk]:.1f}',
                            ha='center', va='bottom', fontweight='bold')
        ax_kde.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, p: f'{v:.2f}'))
        ax_kde.set_ylabel(label, fontsize=15, rotation=90, labelpad=20, va='center')
        if i == 2:
            ax_kde.set_xlabel('Concentration (Res./Voxel)')


def plot_pair_scatter(fig, ax_main, d, d1, d2, s1, s2, colors):
    """Panels B, D: scatter with marginal KDEs.

    The marginals are insets anchored to ax_main rather than separate
    gridspec cells: aspect='equal' shrinks ax_main inside its cell, and only
    an inset tracks that shrunk box, so the marginals stay flush against the
    top and right spines.
    """
    maxi = np.max(d)
    ax_main.scatter(d1, d2, alpha=0.4, s=40, color='grey',
                    edgecolors='black', lw=0.5, rasterized=True)
    ax_main.set_aspect('equal', adjustable='box')
    ax_main.set_xlim(0, maxi)
    ax_main.set_ylim(0, maxi)
    ax_main.set_yticks(ax_main.get_xticks())
    for t in ax_main.get_xticks():
        if 0 <= t <= maxi:
            ax_main.plot([t, t], [0, max(0, maxi - t)],
                         color='gray', ls='--', alpha=0.3, lw=0.8)
            ax_main.plot([0, max(0, maxi - t)], [t, t],
                         color='gray', ls='--', alpha=0.3, lw=0.8)
    ax_main.plot([0, maxi], [maxi, 0], 'k--', alpha=0.5, lw=2, zorder=10)
    ax_main.set_xlim(0, maxi)
    ax_main.set_ylim(0, maxi)
    ax_main.set_xlabel(f'{s1} (Res./Voxel)')
    ax_main.set_ylabel(f'{s2} (Res./Voxel)')

    rng = np.linspace(0, maxi, MARGINAL_GRID)

    ax_top = ax_main.inset_axes([0, 1, 1, MARGINAL_FRAC])
    kde1 = gaussian_kde(d1)(rng)
    ax_top.fill_between(rng, kde1, alpha=0.5, color=colors[0])
    ax_top.plot(rng, kde1, color=colors[0], lw=2)
    ax_top.set_xlim(0, maxi)
    ax_top.set_ylim(0, kde1.max() * 1.05)
    ax_top.axis('off')

    ax_right = ax_main.inset_axes([1, 0, MARGINAL_FRAC, 1])
    kde2 = gaussian_kde(d2)(rng)
    ax_right.fill_betweenx(rng, kde2, alpha=0.5, color=colors[1])
    ax_right.plot(kde2, rng, color=colors[1], lw=2)
    ax_right.set_ylim(0, maxi)
    ax_right.set_xlim(0, kde2.max() * 1.05)
    ax_right.axis('off')


def create_figure2ad():
    fig = plt.figure(figsize=(30, 7))
    gs = gridspec.GridSpec(1, 4, figure=fig, wspace=0.4)
    add_label = lambda ax, lbl: ax.text(-0.1, 1.2, lbl, transform=ax.transAxes,
                                        fontsize=24, fontweight='bold',
                                        va='top', ha='right')

    for i, (pair_id, clrs, lbl1, lbl2) in enumerate(PAIRS):
        d, d1, d2, s1, s2 = load_densities(pair_id)

        ax = fig.add_subplot(gs[0, 2 * i])
        add_label(ax, lbl1)
        plot_pair_composite(fig, ax, d, d1, d2, clrs)

        ax = fig.add_subplot(gs[0, 2 * i + 1])
        add_label(ax, lbl2)
        plot_pair_scatter(fig, ax, d, d1, d2, s1, s2, clrs)

        save_csv({'densities_total': d, 'densities_s1': d1, 'densities_s2': d2},
                 f'panel{lbl1}{lbl2}_{pair_id}_data.csv')

    plt.tight_layout()
    plt.savefig('Figure2AD.pdf', bbox_inches='tight')
    print("Figure 2A-D saved as Figure2AD.pdf")


if __name__ == '__main__':
    create_figure2ad()