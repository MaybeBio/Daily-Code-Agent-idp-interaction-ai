"""
Deck figure for the lineage test.
  left:  leftover by kingdom (organism means, 95% bootstrap CI over organisms), raw and adjusted for pi and D+E+K,
         CALVADOS-2, all chains
  right: residue-age score vs leftover, Spearman by chain (raw) and rank-based partial correlation adjusted for pi and
         D+E+K, with 95% Fisher CIs, in three sets: CALVADOS-2 all chains, CALVADOS-2 common set, Mpipi-GG common set
Inputs: Analysis_FDT/lineage_test/{kingdoms.csv, per_protein.csv, calvados_common/per_protein.csv, mpipi_common/per_protein.csv}
Output: Analysis_FDT/lineage_test/fig_lineage.{svg,png}
"""
import os
import sys

sys.dont_write_bytecode = True

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import matplotlib  # noqa: E402
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
from scipy import stats  # noqa: E402

D = '/Volumes/User Homes/Nu Project/BENDER-MT/Traj Level Analysis/Analysis_FDT/lineage_test'
BG, INK, MUTED, GRID = '#F7F4EE', '#1B2432', '#5B6472', '#E2DCCF'
ACC, ADJ, CAL, MP = '#8F5E14', '#7C8BA0', '#2E6FA3', '#A9502A'
plt.rcParams.update({'font.family': 'Helvetica Neue', 'font.size': 15, 'svg.fonttype': 'path', 'axes.edgecolor': INK,
                     'axes.labelcolor': INK, 'xtick.color': INK, 'ytick.color': INK, 'mathtext.default': 'regular',
                     'mathtext.fontset': 'dejavusans'})


def rank_normal(x):
    return stats.norm.ppf((stats.rankdata(x) - 0.5) / len(x))


def resid(y, X):
    X1 = np.column_stack([np.ones(len(y)), X])
    return y - X1 @ np.linalg.lstsq(X1, y, rcond=None)[0]


def fisher_ci(r, n):
    z, se = np.arctanh(r), 1 / np.sqrt(n - 3)
    return np.tanh(z - 1.96 * se), np.tanh(z + 1.96 * se)


def main():
    K = pd.read_csv(os.path.join(D, 'kingdoms.csv'), index_col=0)
    order = ['Viruses', 'Bacteria', 'Protists', 'Fungi', 'Plants', 'Mammals']
    lab = {'Viruses': 'Viruses', 'Bacteria': 'Bacteria  (uni)', 'Protists': 'Protists  (uni)', 'Fungi': 'Fungi  (mixed)',
           'Plants': 'Plants  (multi)', 'Mammals': 'Mammals  (multi)'}
    fig, axes = plt.subplots(1, 2, figsize=(16.6, 5.4), facecolor=BG, gridspec_kw={'width_ratios': [1.1, 1]})
    ax = axes[0]
    ax.set_facecolor(BG)
    y = np.arange(len(order))[::-1]
    ax.axvline(0, color=INK, lw=1.2)
    for k, (col, c, name, dy) in enumerate((('L', ACC, 'as simulated', 0.14), ('L_adj', ADJ, 'adjusted for π and D+E+K', -0.14))):
        m = K.loc[order, f'{col}_mean'].to_numpy()
        lo, hi = K.loc[order, f'{col}_lo'].to_numpy(), K.loc[order, f'{col}_hi'].to_numpy()
        for yy, a, b in zip(y + dy, lo, hi):
            ax.plot([a, b], [yy, yy], color=c, lw=3, solid_capstyle='round')
        ax.plot(m, y + dy, 'o', color=c, ms=9, label=name)
    ax.set_yticks(y)
    ax.set_yticklabels([f'{lab[k]}   {int(K.at[k, "organisms"])} org.' for k in order], fontsize=14)
    ax.set_xlabel('Leftover, organism mean (SD units; 95% CI)', fontsize=15)
    ax.set_title('Kingdoms, CALVADOS-2 (10,282 chains, 2,156 organisms)', fontsize=16, loc='left', color=INK, pad=10)
    ax.legend(frameon=False, fontsize=13, loc='upper left', bbox_to_anchor=(0.0, 1.0))
    for s_ in ('top', 'right', 'left'):
        ax.spines[s_].set_visible(False)
    ax.tick_params(axis='y', length=0)
    ax.set_xlim(-0.6, 0.45)
    ax = axes[1]
    ax.set_facecolor(BG)
    sets = [('CALVADOS-2, all chains', os.path.join(D, 'per_protein.csv'), CAL),
            ('CALVADOS-2, common set', os.path.join(D, 'calvados_common', 'per_protein.csv'), CAL),
            ('Mpipi-GG, common set', os.path.join(D, 'mpipi_common', 'per_protein.csv'), MP)]
    yy = np.arange(len(sets))[::-1]
    ax.axvline(0, color=INK, lw=1.2)
    res = {}
    for y0, (name, path, c) in zip(yy, sets):
        p = pd.read_csv(path)
        n = len(p)
        r_raw = stats.spearmanr(p.age, p.L)[0]
        ra = resid(rank_normal(p.age.to_numpy()), np.column_stack([rank_normal(p.pi.to_numpy()), rank_normal(p.dek.to_numpy())]))
        r_par = stats.pearsonr(ra, p.L_adj)[0]
        res[name] = (round(r_raw, 3), round(r_par, 3), n)
        for r, dy, mk, fill in ((r_raw, 0.14, 'o', c), (r_par, -0.14, 's', 'white')):
            lo, hi = fisher_ci(r, n)
            ax.plot([lo, hi], [y0 + dy, y0 + dy], color=c, lw=3, solid_capstyle='round')
            ax.plot(r, y0 + dy, mk, color=c, mfc=fill, mew=2.2, ms=9)
    ax.plot([], [], 'o', color=INK, ms=9, label='raw (Spearman)')
    ax.plot([], [], 's', color=INK, mfc='white', mew=2.2, ms=9, label='beyond π and D+E+K (partial)')
    ax.set_yticks(yy)
    ax.set_yticklabels([f'{s[0]}   n = {res[s[0]][2]:,}' for s in sets], fontsize=14)
    ax.set_xlabel('Correlation of residue-age score with leftover (95% CI)', fontsize=15)
    ax.set_title('Later-arriving residues, both models', fontsize=16, loc='left', color=INK, pad=10)
    ax.legend(frameon=False, fontsize=13, loc='center right', bbox_to_anchor=(1.0, 0.33))
    ax.set_xlim(-0.1, 0.5)
    ax.set_xticks([-0.1, 0, 0.1, 0.2, 0.3])
    for s_ in ('top', 'right', 'left'):
        ax.spines[s_].set_visible(False)
    ax.tick_params(axis='y', length=0)
    fig.tight_layout(w_pad=3)
    for ext, kw in (('svg', {}), ('png', {'dpi': 110})):
        fig.savefig(os.path.join(D, f'fig_lineage.{ext}'), facecolor=BG, **kw)
    print(res)


if __name__ == '__main__':
    main()
