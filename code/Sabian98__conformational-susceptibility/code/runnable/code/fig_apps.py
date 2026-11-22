"""
fig_apps.py — manuscript Fig. 5 in one file (2026-10-15; replaces the three-file panel stack fig_heat_regulation.png,
fig_mutation_effects.png and fig_tracks.pdf in results.tex). CALVADOS-2, exact reweighting at 300 K, stored results only.
  a  change of the total response of Rg against the change of chain size, on phosphorylating one Ser/Thr (left; sites
     with effective sample fraction > 0.1) and on doubling the ionic strength (right; chains with salt effective sample
     > 0.1, coloured by |net charge per residue|), with the Spearman correlation (as fig_total_response.heat_regulation)
  b  median change of the total response by the residue introduced, over the reliable single substitutions of
     mutation_scan.py (as fig_total_response.mutation_effects)
  c, d  sequence tracks of HMGA1 (human) and myelin basic protein (pig): fig_extra.sequence_track
Inputs: Analysis_FDT/total_response/heat_crosses/{sites,chains}.csv, mutation_scan/analysis_by_mut.csv, BENDER_BIO.csv
    BENDER_REANALYSIS_DIR=... python fig_apps.py            (openmpipi env)
Output: BENDER-MT/manuscript/figs/fig_apps.{pdf,png}
"""
import sys

sys.dont_write_bytecode = True

import os  # noqa: E402

import matplotlib  # noqa: E402

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402
from scipy import stats  # noqa: E402

HERE = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, HERE)
import fig_extra as FX  # noqa: E402

TR = FX.TR
HCX = f'{TR}/heat_crosses'
BIO = (__import__('os').environ['CS_DATA'] + '/BENDER-BIO/BENDER_BIO.csv')
OUT = f'{FX.ROOT}/manuscript/figs/fig_apps'
C2 = FX.COL['CALVADOS-2']
CLS = {**dict.fromkeys('DE', ('acidic', FX.ACIDIC)), **dict.fromkeys('KR', ('basic', FX.BASIC)),
       **dict.fromkeys('FWY', ('aromatic', '#C9A86A')), **dict.fromkeys('NQSTHG', ('polar, Gly', '#3E8E7E')),
       **dict.fromkeys('AVLIMCP', ('hydrophobic, Pro', '#8A8F98'))}


def crosses():
    """Reliable phosphosites and salt chains, exactly as fig_total_response.heat_regulation selects them."""
    S = pd.read_csv(f'{HCX}/sites.csv')
    C = pd.read_csv(f'{HCX}/chains.csv').drop_duplicates('UniProt_ID')
    bio = pd.read_csv(BIO, usecols=['UniProt_ID', 'NCPR'], low_memory=False).set_index('UniProt_ID')
    C = C.join(bio, on='UniProt_ID')
    return S[S.ess > 0.1], C[C.salt_ess > 0.1]


def phospho_panel(ax, S):
    x, y = S.phospho_dRg_pct, S.phospho_dtotal
    ext = (-12.6, 11.0, -0.62, 0.45)                       # every site inside the axes (capped model, 2026-10-25)
    assert x.min() > ext[0] and x.max() < ext[1] and y.min() > ext[2] and y.max() < ext[3], 'sites outside the axes'
    ax.hexbin(x, y, gridsize=(34, 40), bins='log', cmap=FX.ff_cmap(C2), mincnt=1, linewidths=0, extent=ext)
    ax.axhline(0, color='0.55', lw=0.5, zorder=0)
    ax.axvline(0, color='0.55', lw=0.5, zorder=0)
    rho = stats.spearmanr(x, y)[0]
    ax.text(0.97, 0.04, f'ρ = {rho:.2f}\n{len(S):,} sites', transform=ax.transAxes, ha='right', va='bottom')
    ax.set_xlim(ext[:2])
    ax.set_ylim(ext[2:])
    ax.set_xticks([-10, -5, 0, 5, 10])
    ax.set_title('Phosphorylation of one Ser/Thr', loc='left')
    ax.set_xlabel(f'Change of {FX.RG} (%)')
    ax.set_ylabel('Change of total\nresponse (% per K)')
    return rho


def salt_panel(ax, C, fig, cax):
    x, y = C.salt_dRg_pct, C.salt_dtotal
    lim = (-11.5, 23.5), (-0.56, 0.5)                      # every chain inside the axes (capped model, 2026-10-25)
    assert x.min() > lim[0][0] and x.max() < lim[0][1] and y.min() > lim[1][0] and y.max() < lim[1][1]
    cm = LinearSegmentedColormap.from_list('ncpr', plt.get_cmap('viridis_r')(np.linspace(0.1, 1.0, 256)))
    sc = ax.scatter(x, y, c=C.NCPR.abs(), cmap=cm, s=5, lw=0, vmin=0, vmax=0.4, zorder=3)
    ax.axhline(0, color='0.55', lw=0.5, zorder=0)
    ax.axvline(0, color='0.55', lw=0.5, zorder=0)
    rho = stats.spearmanr(x, y)[0]
    ax.text(0.97, 0.04, f'ρ = {rho:.2f}\n{len(C):,} chains', transform=ax.transAxes, ha='right', va='bottom')
    ax.set_xlim(lim[0])
    ax.set_ylim(lim[1])
    ax.set_xticks([-10, 0, 10, 20])
    ax.set_title('Ionic strength doubled', loc='left')
    ax.set_xlabel(f'Change of {FX.RG} (%)')
    cb = fig.colorbar(sc, cax=cax)
    cb.set_label('|net charge per residue|', fontsize=6.5)
    cb.set_ticks([0, 0.1, 0.2, 0.3, 0.4])
    cb.ax.tick_params(labelsize=6.5, width=0.5, length=2)
    cb.outline.set_visible(False)
    return rho


def mutation_panel(ax):
    t = pd.read_csv(f'{TR}/mutation_scan/analysis_by_mut.csv', index_col=0).sort_values('median_change')
    x = np.arange(len(t))
    ax.bar(x, t.median_change, color=[CLS[a][1] for a in t.index], width=0.75, lw=0)
    ax.axhline(0, color=FX.INK, lw=0.6)
    ax.set_xticks(x)
    ax.set_xticklabels(t.index)
    ax.tick_params(axis='x', length=0, pad=2)
    ax.set_xlim(-0.6, len(t) - 0.4)
    lo, hi = min(t.median_change.min(), 0), max(t.median_change.max(), 0)      # every bar inside the axes
    ax.set_ylim(lo - 0.12 * (hi - lo), hi + 0.12 * (hi - lo))
    ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.005))
    ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: '0' if v == 0 else FX.minus(f'{v:.3f}')))
    ax.set_xlabel('Residue introduced')
    ax.set_ylabel('Median change of total\nresponse (% per K)')
    seen = {}
    for a in t.index:
        seen.setdefault(CLS[a][0], CLS[a][1])
    ax.legend(handles=[matplotlib.patches.Patch(color=c, label=lab) for lab, c in seen.items()], loc='lower right',
              handlelength=0.9, handletextpad=0.35, labelspacing=0.2, borderaxespad=0.3)   # empty corner, below zero
    return int(t.n.sum())


def main():
    FX.house()
    S, C = crosses()
    fig = plt.figure(figsize=(FX.WIDTH, 6.35))
    r1 = fig.add_gridspec(1, 6, width_ratios=[0.86, 0.27, 0.86, 0.04, 0.8, 1.32], left=0.1, right=0.99, top=0.945,
                          bottom=0.705, wspace=0.05)
    aa = fig.add_subplot(r1[0, 0])
    ab = fig.add_subplot(r1[0, 2])
    cax = fig.add_subplot(r1[0, 3])
    am = fig.add_subplot(r1[0, 5])
    r_ph = phospho_panel(aa, S)
    r_salt = salt_panel(ab, C, fig, cax)
    n_mut = mutation_panel(am)
    gs = fig.add_gridspec(2, 1, hspace=0.45, top=0.6, bottom=0.085, left=0.2, right=0.99)
    tr = [FX.sequence_track(fig, gs[k], uid, t) for k, (uid, t) in enumerate(FX.TRACK_EX)]
    FX.track_legend(fig, 'lower center', (0.595, -0.004))
    FX.add_letters(fig, [(aa, 'a', {'row': 1}), (am, 'b', {'row': 1})] +
                   [(t[0], 'cd'[k], {'xaxes': t}) for k, t in enumerate(tr)])
    FX.save(fig, OUT)
    out = {'phospho_rho': round(r_ph, 2), 'salt_rho': round(r_salt, 2), 'sites': len(S), 'chains': len(C),
           'mutants': n_mut}
    print('written', OUT, out)
    return out


if __name__ == '__main__':
    main()
