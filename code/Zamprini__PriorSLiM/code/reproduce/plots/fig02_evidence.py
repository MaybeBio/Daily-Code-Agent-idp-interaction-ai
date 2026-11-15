from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.text import Text
from matplotlib.colors import LinearSegmentedColormap, Normalize, ListedColormap, to_rgb, to_rgba
FEATURES = [('local_conservation_peak', 'Conservation peak'), ('mean_consensus_percentile_IDR', 'Conservation percentile'), ('mean_pssm_information', 'PSSM information'), ('motif_min_abs_offset', 'Motif centre distance'), ('motif_center_enrichment', 'Motif centre enrichment'), ('motif_count_log', 'Motif count'), ('density_positive_percentile_IDR', 'Positive-charge density'), ('center_enrich_polar_percentile_IDR', 'Polar centre enrichment'), ('density_aromatic_percentile_IDR', 'Aromatic density'), ('access_occlusion_cis', 'Within-IDR occlusion')]

def axstyle(ax):
    ax.spines[['top', 'right']].set_visible(False)
    ax.tick_params(length=2, width=0.55, pad=1.5)
    for spine in ax.spines.values():
        spine.set_linewidth(0.55)

def draw_d(axes, d):
    bases = ['Explicit_only', 'MotifAE_only']
    for row, base in enumerate(bases):
        p = d['paired'][d['paired'].baseline == base].pivot(index='fold', columns='model', values='ap')
        delta = (p.full - p.baseline) * 100
        ax = axes[0]
        ax.barh(row, delta.mean(), height=0.44, color=BLUE if row == 0 else GOLD, alpha=0.42)
        ax.scatter(delta, row + np.linspace(-0.12, 0.12, 10), s=5, c=BLUE if row == 0 else GOLD, linewidth=0)
        ax.text(delta.mean() + 0.6, row, f'{delta.mean():.1f}', va='center', fontsize=6.5)
        for col, names in enumerate([['FN_rescued', 'TP_lost'], ['FP_corrected', 'TN_harmed']], 1):
            ax = axes[col]
            for j, name in enumerate(names):
                g = d['transitions'][(d['transitions'].baseline == base) & (d['transitions'].transition == name)]
                v = g.fraction.to_numpy() * 100 * (1 if j == 0 else -1)
                ax.barh(row, v.mean(), height=0.44, color=TEAL if j == 0 else RED, alpha=0.85)
                ax.scatter(v, row + np.linspace(-0.12, 0.12, 10), s=3.5, c='#333', linewidth=0, alpha=0.65)
    for ax in axes:
        ax.axvline(0, c='#777', lw=0.55)
        ax.set_ylim(1.5, -0.5)
        axstyle(ax)
    axes[0].set(yticks=[0, 1], yticklabels=['Explicit', 'MotifAE'], xlabel='AP gain (pp)', xlim=(0, 16), xticks=[0, 5, 10, 15])
    axes[0].set_title('Add the other source', pad=4)
    for ax, title, xlabel in [(axes[1], 'Rescued / lost', '% of positives'), (axes[2], 'Corrected / new FP', '% of negatives')]:
        ax.set(yticks=[], xlabel=xlabel)
        ax.set_title(title, pad=4)

def draw_e(ax, cax, d):
    grid, stable = ([], [])
    for base in ['Explicit_only', 'MotifAE_only']:
        for kind in ['positive', 'negative']:
            g = d['effects'][(d['effects'].baseline == base) & (d['effects'].analysis == kind)]
            st = g.groupby('feature').standardized_difference.agg(['mean', 'min', 'max', 'count']).reindex([x for x, _ in FEATURES])
            grid.append(st['mean'].to_numpy())
            stable.append(((st['min'] > 0) | (st['max'] < 0)) & (st['count'] == 10))
    grid, stable = (np.array(grid).T, np.array(stable).T)
    lim = max(0.6, np.ceil(abs(grid).max() * 10) / 10)
    im = ax.imshow(grid, cmap='RdBu_r', vmin=-lim, vmax=lim, aspect='auto')
    for i, j in zip(*np.where(stable)):
        ax.plot(j, i, '.', c='white' if abs(grid[i, j]) > lim * 0.55 else '#222', ms=2.7)
    ax.set(yticks=range(10), yticklabels=[n for _, n in FEATURES], xticks=range(4), xticklabels=['E→F\nPos.', 'E→F\nNeg.', 'M→F\nPos.', 'M→F\nNeg.'])
    ax.tick_params(length=0, labelsize=6.1, pad=1.5)
    ax.axvline(1.5, c='white', lw=2)
    ax.set_title('Features of corrected windows', pad=4)
    for sp in ax.spines.values():
        sp.set_visible(False)
    cb = ax.figure.colorbar(im, cax=cax, orientation='vertical', ticks=[-lim, 0, lim])
    cb.ax.tick_params(labelsize=6, length=1, pad=1)
    cb.set_label('Standardized difference', fontsize=6, labelpad=2)

def draw_g(ax, d, detailed=False):
    allg, broad = (d['allg'], d['broad'])
    ax.scatter(allg.M_center_bce, allg.B_center_bce, s=2.4, c='#afb6bc', alpha=0.26, linewidth=0, rasterized=True)
    ax.scatter(broad.M_center_bce, broad.B_center_bce, s=3.2, c=BLUE, alpha=0.32, linewidth=0, rasterized=True)
    high = broad[broad.near_zero_stable.astype(bool)]
    ax.scatter(high.M_center_bce, high.B_center_bce, s=11, c=GOLD, edgecolor='white', lw=0.25, zorder=4)
    ax.set_xscale('symlog', linthresh=1e-05)
    ax.set_yscale('symlog', linthresh=1e-05)
    lim = max(abs(allg[['M_center_bce', 'B_center_bce']].to_numpy(float)).max(), 0.001) * 1.1
    ax.plot([-lim, lim], [-lim, lim], c='#666', ls='--', lw=0.55)
    ax.axhline(0, c='#ccc', lw=0.4)
    ax.axvline(0, c='#ccc', lw=0.4)
    ax.set(xlabel='MotifAE-only deletion ΔBCE', ylabel='Joint deletion ΔBCE')
    ax.set_xticks([0, 0.0001, 0.01])
    ax.set_yticks([0, 0.0001, 0.01])
    ax.set_xticklabels(['0', '10⁻⁴', '10⁻²'])
    ax.set_yticklabels(['0', '10⁻⁴', '10⁻²'])
    ax.set_title('Dependence after fusion', pad=4)
    r = broad[broad.feature_id == 27779].iloc[0]
    ax.annotate('27779', (r.M_center_bce, r.B_center_bce), xytext=(8, -11), textcoords='offset points', arrowprops=dict(arrowstyle='-', lw=0.45, color=GOLD), fontsize=6.2, color='#80541E')
    ax.text(0.02, 0.97, f'{len(allg):,} IDs; {len(high)} near-zero', ha='left', va='top', transform=ax.transAxes, fontsize=5.8)
    axstyle(ax)
BLUE = TEAL = '#276D91'
GOLD = RED = '#BE8345'
