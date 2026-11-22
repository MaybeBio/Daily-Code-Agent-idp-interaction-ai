"""
fig_manuscript.py — validation panels of manuscript Fig. 1 (d-f; assembled by fig_manuscript2.fig1), also
saved stand-alone.
  a  integrated susceptibility (278-320 K) vs measured change of <Rg>, CALVADOS-2, all chains
  b  measured / predicted change for the eight observables, CALVADOS-2 and Mpipi-GG
  c  median chi(Rg) at 300 K by own nu bin, three force fields
    python fig_manuscript.py
Output: BENDER-MT/manuscript/figs/fig1_validation.{png,pdf}
"""
import sys

sys.dont_write_bytecode = True

import os  # noqa: E402

import matplotlib  # noqa: E402

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import figfont  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

ROOT = '/Volumes/User Homes/Nu Project/BENDER-MT'
TLA = f'{ROOT}/Traj Level Analysis'
OUT = f'{ROOT}/manuscript/figs'
COL = {'CALVADOS-2': '#2E6FA3', 'Mpipi-GG': '#A9502A', 'Mpipi-Recharged': '#6B4FA0'}   # as fig_total_response
OBS = [('rg', 'Rg'), ('ree', 'Ree'), ('delta', 'Asphericity'), ('fragmentation_index', 'Fragmentation'),
       ('global_efficiency', 'Global eff.'), ('avg_clustering', 'Clustering'), ('transitivity', 'Transitivity'),
       ('degree_assortativity', 'Assortativity')]


def panel_scatter(ax):
    """Integrated susceptibility (278-320 K) against the measured change of <Rg>, CALVADOS-2, all chains."""
    v = pd.read_csv(f'{TLA}/Analysis_FDT/fdt_validation.csv', usecols=['UniProt_ID', 'metric', 'integral', 'observed_diff'])
    r = v[v.metric == 'rg']
    ax.scatter(r.integral, r.observed_diff, s=1.5, alpha=0.25, color=COL['CALVADOS-2'], rasterized=True, lw=0)
    lim = np.nanpercentile(np.r_[r.integral, r.observed_diff], [0.2, 99.8])
    ax.plot(lim, lim, color='k', lw=0.7)
    ax.set(xlim=lim, ylim=lim, xlabel=r'Predicted $\Delta\langle R_g\rangle$, 278$\rightarrow$320 K (nm)',
           ylabel=r'Measured $\Delta\langle R_g\rangle$ (nm)')
    ax.text(0.04, 0.96, f'CALVADOS-2\n{r.UniProt_ID.nunique():,} chains', transform=ax.transAxes, va='top')


def panel_ratio(ax):
    """Measured / predicted change for the eight observables, CALVADOS-2 and Mpipi-GG."""
    sc = pd.read_csv(f'{TLA}/Analysis_FDT/fdt_validation_summary.csv')
    sc = sc[sc['T'] == 'integral'].set_index('metric')
    sm = pd.read_csv(f'{TLA}/Analysis_FDT_MPIPI/mpipi_fdt_validation_summary.csv')
    sm = sm[sm['T'] == 'integral'].set_index('metric')
    x = np.arange(len(OBS))
    rc = [sc.loc[m, 'mean_slope'] / sc.loc[m, 'mean_chi'] for m, _ in OBS]
    rm = [sm.loc[m, 'ratio_observed_to_integral'] for m, _ in OBS]
    ax.bar(x - 0.2, rc, 0.38, color=COL['CALVADOS-2'], label='CALVADOS-2')
    ax.bar(x + 0.2, rm, 0.38, color=COL['Mpipi-GG'], label='Mpipi-GG')
    ax.axhline(1, color='k', lw=0.7)
    ax.set(ylim=(0.8, 1.2), ylabel='Measured / predicted change', xticks=x)
    ax.set_xticklabels([n for _, n in OBS], rotation=45, ha='right')
    ax.legend(frameon=False, loc='upper left', fontsize=7, handlelength=1.2)
    return rc, rm


def panel_chi_by_nu(ax):
    """Median chi(Rg) at 300 K by own nu bin, three force fields."""
    a = pd.read_csv(f'{TLA}/Analysis_FDT_RECHARGED/alongside/A_chi_by_nu.csv')
    a = a[a.own_nu_bin != 'all']
    mids = {'(0.3, 0.45]': 0.43, '(0.45, 0.5]': 0.475, '(0.5, 0.54]': 0.52, '(0.54, 0.58]': 0.56, '(0.58, 0.63]': 0.605,
            '(0.63, 1.0]': 0.65}
    for m, g in a.groupby('model'):
        ax.plot(g.own_nu_bin.map(mids), g.chi_median, 'o-', ms=3, lw=1.2, color=COL[m], label=m)
    ax.axhline(0, color='0.6', lw=0.6)
    ax.axvline(0.5, color='0.6', lw=0.6, ls=':')
    ax.set(xlabel=r'Compaction, $\nu$ (bin)', ylabel=r'$\chi_T(R_g)$ at 300 K (% per K)')
    ax.set_xticks(list(mids.values()))
    ax.set_xticklabels(['<0.45', '0.47', '0.52', '0.56', '0.60', '>0.63'], fontsize=6.5)
    ax.legend(frameon=False, loc='upper right', bbox_to_anchor=(1.04, 1.04), fontsize=7, handlelength=1.2)


def main():
    """Stand-alone version of the three validation panels (now panels d-f of Fig. 1, fig_manuscript2.fig1)."""
    os.makedirs(OUT, exist_ok=True)
    plt.rcParams.update({'font.family': 'Arial', 'font.size': 8, 'axes.spines.top': False, 'axes.spines.right': False})
    figfont.apply()                                   # Arial for all figures
    fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.4), gridspec_kw={'width_ratios': [1, 1.25, 1.15]})
    panel_scatter(ax[0])
    rc, rm = panel_ratio(ax[1])
    panel_chi_by_nu(ax[2])
    for i, p in enumerate(ax):
        p.text(-0.3, 1.06, 'abc'[i], transform=p.transAxes, fontweight='bold', fontsize=10)
    fig.tight_layout(w_pad=2.5)
    for ext in ('png', 'pdf'):
        fig.savefig(f'{OUT}/fig1_validation.{ext}', dpi=300)
    print('ratios C2', np.round(rc, 3), 'GG', np.round(rm, 3))


if __name__ == '__main__':
    main()
