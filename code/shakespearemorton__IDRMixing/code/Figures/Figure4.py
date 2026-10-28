#!/usr/bin/env python3
"""Figure 4: FUS-A1 analysis, demixing predictions, and network topology

Also produces Figure4_extra.pdf: a pairwise matrix of the sabari proteins
coloured by demixing index (D_norm), with pairs that do not form a mixed
condensate (dGij > DG_MIX_CUTOFF, or no dGij at all) greyed out and crossed
through.

Panels A/B read the trajectory-averaged voxel maps written by DDM.py
(DDM_Heterotypic/results/known); the tie lines and their errors come from
DDM.py via Figure4_tielines.py (run it first).  avg_tot and avg_s1 are already in
residues/voxel, so the two species sum to the total and no unit conversion
or box geometry is needed here.  Panels C/D and Figure4_extra read the
DDM_Heterotypic CSVs, with the same 5% two-phase gate as Figure3.py.
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.colors as mcolors
from pathlib import Path
from collections import defaultdict
from scipy.stats import spearmanr, linregress
from matplotlib.patches import Patch
from pycirclize import Circos

plt.rcParams.update({
    'font.size': 20, 'axes.linewidth': 2, 'axes.spines.top': False, 'axes.spines.right': False,
    'xtick.major.size': 8, 'xtick.major.width': 2, 'ytick.major.size': 8, 'ytick.major.width': 2,
    'legend.frameon': False, 'figure.dpi': 300
})

DATA_DIR = Path("figure4_data")
RESULTS_DIR = Path("../DDM_Heterotypic/results/known")
TIELINES_CSV = DATA_DIR / "fus_a1_tielines.csv"   # written by Figure4_tielines.py
EXP_CSV = "../DDM_Heterotypic/heterotypic_exp.csv"
MIN_DENSE_FRAC = 0.05   # two phases only if >= 5% of voxels sit at or above c_den
MAX_DENSE_FRAC = 0.11   # ... and <= 11%: above that the condensate fills the box

# Panel C label offsets (points) and alignment for the crowded labels;
# every other protein is labelled below-right at (9, -15).
PANEL_C_LABELS = {'CBX2': (-10, 8, 'right'), 'KDM1A': (10, 6, 'left'),
                  'Ki67_R12_pm9x2': (10, -18, 'left'), 'NELFA': (-10, -18, 'right'),
                  'RCOR3': (-10, 8, 'right'), 'LEO1': (-10, 8, 'right')}

# phase-concentration points in panels A and B (and their insets)
PHASE_POINT = dict(ms=14, mec='black', mew=2, elinewidth=2, capsize=5, capthick=2, zorder=6)

COLORS = {'ratio_2':'#A6736A', 'ratio_7': '#5A708C', 'normal': '#2c3e50',
          'ewsr1_fus': '#27ae60', 'ctr9': '#e67e22'}

# --- Figure4_extra settings -------------------------------------------------
SABARI_CSV = "../DDM_Heterotypic/heterotypic_sabari.csv"
DG_MIX_CUTOFF = -3.0      # dGij above this -> not mixed -> grey + X
GREY = '#c8c8c8'          # fill for masked (not-mixed) cells
EMPTY = '#ffffff'         # fill for pairs that were never simulated

PROTEIN_ID_MAP = {
    'CCNT1_HUMAN': 'Cyclin T1', 'DDX3X_HUMAN': 'Ddx3', 'DDX4_HUMAN': 'Ddx4', 'DYRK3_HUMAN': 'DYRK3',
    'FUS_HUMAN': 'FUS', 'G3BP1_HUMAN': 'G3BP1', 'GATA3_HUMAN': 'GATA3', 'ROA1_HUMAN': 'hnRNPA1',
    'CBX5_HUMAN': 'HP1α', 'NPHN_HUMAN': 'Nephrin', 'NPM_HUMAN': 'NPM1', 'NUP98_HUMAN': 'Nup98',
    'RBM14_HUMAN': 'RBM-14', 'SYN1_HUMAN': 'Syn-1', 'SYN2_HUMAN': 'Syn-2', 'RBP56_HUMAN': 'TAF-15',
    'TADBP_HUMAN': 'TDP-43', 'UBQL2_HUMAN': 'UBQ2', 'HSPB8_HUMAN': 'HSPB8', 'ESR1_HUMAN': 'Erα'
}


def save_csv(data_dict, filename):
    DATA_DIR.mkdir(exist_ok=True)
    max_len = max(len(v) if hasattr(v, '__len__') and not isinstance(v, str) else 1 for v in data_dict.values())
    for k, v in data_dict.items():
        if not hasattr(v, '__len__') or isinstance(v, str): data_dict[k] = [v] * max_len
        elif len(v) < max_len: data_dict[k] = list(v) + [np.nan] * (max_len - len(v))
    pd.DataFrame(data_dict).to_csv(DATA_DIR / filename, index=False)

def apply_two_phase_gate(df):
    """A voxel map resolves two phases only when between MIN_DENSE_FRAC and
    MAX_DENSE_FRAC of its voxels sit at or above c_den: too few and there is
    no dense phase, too many and the condensate fills the box, leaving no
    clear dilute phase. Otherwise c_dil is the mean over all voxels and no
    dense value (and so no dGij) is defined."""
    df = df.copy()
    df['two_phase'] = df['frac_vox_ge_cden'].between(MIN_DENSE_FRAC, MAX_DENSE_FRAC)
    one = ~df['two_phase']
    for dil in ('cdil', 'cdil1', 'cdil2'):
        df.loc[one, dil] = df.loc[one, f'{dil}_onephase']
    df.loc[one, ['cden', 'cden1', 'cden2', 'dGij', 'dGij_err']] = np.nan
    return df

def get_pair_data(pair_id):
    """Load voxel maps for a protein pair from RESULTS_DIR.

    The npz filename is the base_name ("{name1}_{name2}_{frac1}_{frac2}"),
    so the last two underscore-separated tokens give the composition.  The
    trailing underscore in the glob keeps e.g. FUS_A1 from matching
    FUS_A1p12D.
    """
    compositions = []
    s1 = s2 = None
    for f in sorted(RESULTS_DIR.glob(f"{pair_id}_*_*.npz")):
        try:
            c1, c2 = (int(t) for t in f.stem.rsplit("_", 2)[1:])
        except ValueError:
            continue
        if c1 == 0 or c2 == 0:
            continue

        with np.load(f, allow_pickle=True) as npz:
            d = npz['avg_tot'].astype(np.float64).ravel()
            phi1 = npz['avg_s1'].astype(np.float64).ravel()
            s1, s2 = str(npz['seq_name1']), str(npz['seq_name2'])
        phi2 = d - phi1

        compositions.append({'comp1': c1, 'comp2': c2, 'phi1': phi1, 'phi2': phi2})

    if not compositions:
        print(f"  ! no {RESULTS_DIR}/*.npz for {pair_id}")
        return None, None
    return sorted(compositions, key=lambda x: x['comp1'])[:5], (s1, s2)

def get_experimental_data():
    """FUS-A1 experimental data from Farag et al."""
    df_dil = pd.DataFrame([(100, 1.780, 0.194, 0.000, 0.000), (75, 0.385, 0.011, 0.036, 0.000),
                           (50, 0.143, 0.042, 0.117, 0.000), (25, 0.007, 0.000, 0.319, 0.000),
                           (0, 0.000, 0.000, 0.763, 0.000)], columns=['Ratio', 'FUS', 'FUS_err', 'A1', 'A1_err'])
    df_den = pd.DataFrame([(75, 194.059, 66.806, 84.549, 27.480), (50, 150.107, 29.506, 134.733, 23.796),
                           (25, 59.757, 22.436, 147.483, 44.416)], columns=['Ratio', 'FUS', 'FUS_err', 'A1', 'A1_err'])
    return df_dil, df_den

def plot_fus_a1_panels(fig, gs_main):
    """Panels A-B: FUS-A1 simulation vs experiment."""
    from matplotlib.ticker import MaxNLocator

    comps_1, _ = get_pair_data('FUS_A1')
    comps_2, _ = get_pair_data('FUS_A1p12D')
    if comps_1 is None: print("FUS-A1 data not found"); return

    gs_sub = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs_main, wspace=0.3)
    ax_sim, ax_exp = fig.add_subplot(gs_sub[0, 0]), fig.add_subplot(gs_sub[0, 1])
    axins_sim, axins_exp = ax_sim.inset_axes([0.55, 0.55, 0.4, 0.4]), ax_exp.inset_axes([0.55, 0.55, 0.4, 0.4])

    for ax, lbl in [(ax_sim, 'A'), (ax_exp, 'B')]:
        ax.text(-0.1, 1.1, lbl, transform=ax.transAxes, fontsize=24, fontweight='bold', va='top', ha='right',
                gid=f'label_{lbl}')

    # Plot simulations: voxels, plus the DDM tie lines from Figure4_tielines.py
    tie = pd.read_csv(TIELINES_CSV)
    all_phi = []
    tot_phi = []
    for pair, comps, color in [('FUS_A1', comps_1, COLORS['ratio_2']),
                               ('FUS_A1p12D', comps_2, COLORS['ratio_7'])]:
        if comps is None: continue
        for c in comps:
            phi1, phi2 = c['phi1'], c['phi2']
            all_phi.extend([phi1, phi2])
            tot_phi.extend(phi1+phi2)
            ax_sim.scatter(phi1, phi2, alpha=0.2, s=40, color=color, rasterized=True)
            t = tie[(tie['pair'] == pair) & (tie['comp1'] == c['comp1'])].iloc[0]
            dil, den = (t['dil1'], t['dil2']), (t['den1'], t['den2'])
            ax_sim.plot([dil[0], den[0]], [dil[1], den[1]], color=color, lw=2, alpha=0.8)
            ax_sim.errorbar(den[0], den[1], xerr=t['den1_sem'], yerr=t['den2_sem'], fmt='o',
                            color=color, **PHASE_POINT)
            axins_sim.errorbar(dil[0], dil[1], xerr=t['dil1_sem'], yerr=t['dil2_sem'], fmt='o',
                               color=color, **PHASE_POINT)
            axins_sim.plot([dil[0], den[0]], [dil[1], den[1]], color=color, lw=2, alpha=0.8, zorder=4)

    maxi = max(tot_phi) * 1.1
    # inset (FUS, A1) limits: just past the simulated dilute points and their errors
    sim_inset = (1.15 * (tie['dil1'] + tie['dil1_sem']).max(),
                 1.15 * (tie['dil2'] + tie['dil2_sem']).max())

    # Plot experimental
    df_dil, df_den = get_experimental_data()
    color = COLORS['ratio_2']
    for _, r in df_dil.iterrows():
        for ax in [ax_exp, axins_exp]:
            ax.errorbar(r['FUS'], r['A1'], xerr=r['FUS_err'], yerr=r['A1_err'], fmt='o', color=color,
                       **PHASE_POINT)
    for _, r in df_den.iterrows():
        ax_exp.errorbar(r['FUS'], r['A1'], xerr=r['FUS_err'], yerr=r['A1_err'], fmt='o', color=color,
                       **PHASE_POINT)
    for ratio in [75, 50, 25]:
        d, D = df_dil[df_dil['Ratio'] == ratio].iloc[0], df_den[df_den['Ratio'] == ratio].iloc[0]
        for ax in [ax_exp, axins_exp]:
            ax.plot([d['FUS'], D['FUS']], [d['A1'], D['A1']], color=color, ls='--', lw=2, alpha=0.8)

    # Format
    for ax, ins, lim, lim_in, xl, yl, ttl in [
        (ax_sim, axins_sim, maxi, sim_inset, 'FUS (Res./Voxel)', 'A1 (Res./Voxel)', 'Simulated'),
        (ax_exp, axins_exp, 350, (2, 2), 'FUS (mg/mL)', 'A1 (mg/mL)', 'Experimental [Farag et al.]')]:

        # Main Axes
        ax.set_xlim(0, lim); ax.set_ylim(0, lim); ax.set_aspect('equal')
        ax.plot([0, lim], [lim, 0], color='gray', ls='--', alpha=0.5, lw=2)

        # Inset Formatting
        # start just below 0 so points sitting on an axis are drawn in full
        ins.set_xlim(-0.04 * lim_in[0], lim_in[0]); ins.set_ylim(-0.04 * lim_in[1], lim_in[1])
        ins.xaxis.set_major_locator(MaxNLocator(nbins=2))
        ins.yaxis.set_major_locator(MaxNLocator(nbins=2))

        ax.set_xlabel(xl); ax.set_ylabel(yl); ax.set_title(ttl, fontweight='bold')

    # Legend - Positioned 1/3 up (0.35) and far right (1.0)
    ax_sim.legend(handles=[Patch(fc=COLORS['ratio_2'], ec='black', label='A1$_{WT}$'),
                           Patch(fc=COLORS['ratio_7'], ec='black', label='A1$_{p12D}$')],
                  loc='center right', bbox_to_anchor=(1.1, 0.35))

    def add_gridlines(ax, limit):
        # Determine ticks (will be integers for insets due to MaxNLocator)
        x_ticks = [tick for tick in ax.get_xticks() if 0 <= tick <= limit]

        # Plot grid for X axis
        for x in x_ticks:
            ax.plot([x, x], [0, max(0, limit - x)], color='gray',
                    linestyle='--', alpha=0.3, linewidth=0.8)
        # Plot grid for Y axis (using x_ticks to ensure square/even grid)
        for y in x_ticks:
            ax.plot([0, max(0, limit - y)], [y, y], color='gray',
                    linestyle='--', alpha=0.3, linewidth=0.8)

    add_gridlines(ax_sim, maxi)
    add_gridlines(ax_exp, 350)
    add_gridlines(axins_exp, .2)

def plot_demixing_difference(fig, ax):
    """Panel C: ΔDemixing vs MED1 pellet ratio."""
    ax.text(-0.1, 1.1, 'C', transform=ax.transAxes, fontsize=24, fontweight='bold', va='top', ha='right')

    if not os.path.exists(SABARI_CSV):
        ax.text(0.5, 0.5, 'Data not found', ha='center', va='center', transform=ax.transAxes); return
    df_demix = apply_two_phase_gate(pd.read_csv(SABARI_CSV))
    df_pellet = pd.read_csv("../DATASETS/sabari_all_pellet.csv")

    # Calculate ΔDemixing for each protein
    data = []
    for prot in df_demix['seq_name1'].unique():
        if prot in ['MED1', 'NPM1']: continue
        med1 = df_demix[((df_demix['seq_name1'] == prot) & (df_demix['seq_name2'] == 'MED1')) |
                        ((df_demix['seq_name1'] == 'MED1') & (df_demix['seq_name2'] == prot))]['D_norm']
        npm1 = df_demix[((df_demix['seq_name1'] == prot) & (df_demix['seq_name2'] == 'NPM1')) |
                        ((df_demix['seq_name1'] == 'NPM1') & (df_demix['seq_name2'] == prot))]['D_norm']
        if len(med1) > 0 and len(npm1) > 0:
            data.append({'protein': prot, 'ddemix': npm1.values[0] - med1.values[0]})

    df_pellet['protein'] = df_pellet.apply(lambda x: str(x['Gene name']).strip() if pd.notna(x['Gene name']) else str(x['Uniprot ID']), axis=1)
    df_plot = pd.merge(pd.DataFrame(data), df_pellet[['protein', 'Log2(average_P/S)']].rename(columns={'Log2(average_P/S)': 'pellet_ratio'}))

    if df_plot.empty: ax.text(0.5, 0.5, 'No data', ha='center', va='center', transform=ax.transAxes); return

    outliers = {'EWSR1': ('ewsr1_fus', 's'), 'FUS': ('ewsr1_fus', 's'), 'CTR9': ('ctr9', '^')}
    normal = df_plot[~df_plot['protein'].isin(outliers)]

    ax.scatter(normal['ddemix'], normal['pellet_ratio'], s=300, c=COLORS['normal'], edgecolors='k', alpha=0.7)
    for _, r in df_plot[df_plot['protein'].isin(outliers)].iterrows():
        c, m = outliers[r['protein']]
        ax.scatter(r['ddemix'], r['pellet_ratio'], s=350, c=COLORS[c], edgecolors='k', marker=m, zorder=3)

    if len(normal) > 1:
        res = linregress(normal['ddemix'], normal['pellet_ratio'])
        x_rng = np.array([df_plot['ddemix'].min(), df_plot['ddemix'].max()])
        ax.plot(x_rng, res.intercept + res.slope * x_rng, 'k--', alpha=0.5, linewidth=3)

        rho_filt, _ = spearmanr(normal['ddemix'], normal['pellet_ratio'])
        rho_full, _ = spearmanr(df_plot['ddemix'], df_plot['pellet_ratio'])

        # Display both values
        ax.text(0.05, 0.95, f'Filtered ρ = {rho_filt:.2f}\nFull ρ = {rho_full:.2f}',
                transform=ax.transAxes, va='top', bbox=dict(fc='white', alpha=0.8))

    ax.set_xlabel('ΔD (NPM1 - MED1)'); ax.set_ylabel('MED1 Pellet / Supernatant'); ax.set_ylim(-5, 15)

    for _, r in df_plot.iterrows():
        dx, dy, ha = PANEL_C_LABELS.get(r['protein'], (9, -15, 'left'))
        ax.annotate(r['protein'], (r['ddemix'], r['pellet_ratio']), xytext=(dx, dy),
                    textcoords='offset points', fontsize=18, ha=ha)

    save_csv(df_plot.to_dict('list'), 'panelC_demixing_difference.csv')

def plot_network(fig, ax):
    """Panel D: Network topology of mixing interactions."""

    if not os.path.exists(EXP_CSV): ax.text(0.5, 0.5, 'Data not found', ha='center', va='center', transform=ax.transAxes); return

    df = apply_two_phase_gate(pd.read_csv(EXP_CSV))
    df['seq_name1'], df['seq_name2'] = df['seq_name1'].replace(PROTEIN_ID_MAP), df['seq_name2'].replace(PROTEIN_ID_MAP)
    all_proteins = sorted(set(df['seq_name1']) | set(df['seq_name2']))

    adj_df = pd.DataFrame(0.0, index=all_proteins, columns=all_proteins)
    degree, edge_strength = defaultdict(int), {}

    for _, r in df.iterrows():
        if r.get('dGij', 0) < -4 and r.get('D_norm', 1) < 0.13:
            u, v = tuple(sorted([r['seq_name1'], r['seq_name2']]))
            adj_df.at[u, v] = 1 - r['D_norm']
            edge_strength[(u, v)] = abs(r['dGij'])
            degree[r['seq_name1']] += 1; degree[r['seq_name2']] += 1

    sorted_proteins = sorted(all_proteins, key=lambda x: degree.get(x, 0))
    adj_df = adj_df.reindex(index=sorted_proteins, columns=sorted_proteins)
    for p in sorted_proteins:
        adj_df.at[p, p] = 1

    s_min, s_max = (min(edge_strength.values()), max(edge_strength.values())) if edge_strength else (0, 1)
    def link_kws(f, t):
        if f == t: return dict(ec='none', lw=0, alpha=0)
        raw = edge_strength.get(tuple(sorted([f, t])), 0)
        alpha = 0.1 + 0.7 * ((raw - s_min) / (s_max - s_min) if s_max > s_min else 0.5)
        return dict(ec='k', lw=0.2, alpha=alpha)

    active = {p for p in all_proteins if degree.get(p, 0) > 0}
    degs = [degree[p] for p in active]
    norm = mcolors.Normalize(vmin=min(degs), vmax=max(degs)) if degs else mcolors.Normalize(0, 1)
    cmap = plt.get_cmap('magma')
    colors = {n: mcolors.to_hex(cmap(norm(degree[n]))) for n in active}

    pos = ax.get_position(); ax.remove()
    # the label goes on the figure: ax is replaced by the polar axes below
    fig.text(pos.x0, pos.y1, 'D', fontsize=24, fontweight='bold', va='top', ha='right', gid='label_D')
    ax_polar = fig.add_axes(pos, projection='polar')
    circos = Circos.chord_diagram(adj_df, space=2, order=sorted_proteins, cmap=colors,
                                   label_kws=dict(r=110, size=17), link_kws_handler=link_kws)
    circos.plotfig(ax=ax_polar)

# --------------------------------------------------------------- Figure4_extra
def build_sabari_matrices(df, dg_cutoff=DG_MIX_CUTOFF):
    """Symmetric D_norm and dGij matrices over every simulated sabari pair.

    Pairs are unordered, so a row is written into both (i, j) and (j, i).
    Where the same pair appears more than once (e.g. several compositions),
    the values are averaged.
    """
    df = df.copy()
    df['seq_name1'] = df['seq_name1'].replace(PROTEIN_ID_MAP)
    df['seq_name2'] = df['seq_name2'].replace(PROTEIN_ID_MAP)

    proteins = sorted(set(df['seq_name1']) | set(df['seq_name2']))
    d_mat = pd.DataFrame(np.nan, index=proteins, columns=proteins, dtype=float)
    g_mat = pd.DataFrame(np.nan, index=proteins, columns=proteins, dtype=float)

    for (a, b), grp in df.groupby([df[['seq_name1', 'seq_name2']].min(axis=1),
                                   df[['seq_name1', 'seq_name2']].max(axis=1)]):
        d = grp['D_norm'].mean(skipna=True)
        g = grp['dGij'].mean(skipna=True) if 'dGij' in grp else np.nan
        for i, j in [(a, b), (b, a)]:
            d_mat.at[i, j] = d
            g_mat.at[i, j] = g

    # not mixed: dGij above the cutoff, or no dGij at all (single phase)
    masked = ~(g_mat < dg_cutoff)
    return d_mat, g_mat, masked


def plot_demixing_matrix(out_pdf='Figure4_extra.pdf', dg_cutoff=DG_MIX_CUTOFF,
                         annotate=False):
    """Figure4_extra: pairwise demixing matrix for the sabari proteins.

    Cells are coloured by demixing index D_norm.  Pairs that do not form a
    mixed condensate (dGij > dg_cutoff, or dGij missing) are drawn grey with
    an X through them.  Pairs that were never simulated are left blank.
    """
    if not os.path.exists(SABARI_CSV):
        print(f"Sabari data not found ({SABARI_CSV})")
        return

    df = apply_two_phase_gate(pd.read_csv(SABARI_CSV))
    d_mat, g_mat, masked = build_sabari_matrices(df, dg_cutoff)
    proteins = list(d_mat.index)
    n = len(proteins)
    if n == 0:
        print("No sabari pairs found"); return

    simulated = d_mat.notna() | g_mat.notna()
    shown = d_mat.where(~masked)          # only unmasked cells get a colour

    size = max(8., 0.55 * n + 4.)
    fig, ax = plt.subplots(figsize=(size, size * 0.92))

    cmap = plt.get_cmap('viridis').copy()
    cmap.set_bad(EMPTY)
    vals = shown.values.astype(float)
    finite = vals[np.isfinite(vals)]
    vmin, vmax = (float(finite.min()), float(finite.max())) if finite.size else (0., 1.)
    if vmax - vmin < 1e-9:
        vmin, vmax = vmin - 0.05, vmax + 0.05

    im = ax.imshow(np.ma.masked_invalid(vals), cmap=cmap, vmin=vmin, vmax=vmax,
                   origin='upper', interpolation='none')

    # grey fill + X for the not-mixed pairs; blank for never-simulated pairs
    for i in range(n):
        for j in range(n):
            if not simulated.iat[i, j]:
                ax.add_patch(plt.Rectangle((j - .5, i - .5), 1, 1, fc=EMPTY,
                                           ec='none', zorder=2))
                continue
            if masked.iat[i, j]:
                ax.add_patch(plt.Rectangle((j - .5, i - .5), 1, 1, fc=GREY,
                                           ec='none', zorder=2))
                ax.plot([j - .5, j + .5], [i - .5, i + .5], color='0.35',
                        lw=1.8, zorder=3, solid_capstyle='round')
                ax.plot([j - .5, j + .5], [i + .5, i - .5], color='0.35',
                        lw=1.8, zorder=3, solid_capstyle='round')
            elif annotate and np.isfinite(vals[i, j]):
                lum = cmap((vals[i, j] - vmin) / (vmax - vmin))[:3]
                ax.text(j, i, f"{vals[i, j]:.2f}", ha='center', va='center',
                        fontsize=max(6, 150 // max(n, 1)), zorder=4,
                        color='white' if sum(lum) < 1.5 else 'black')

    # cell borders
    ax.set_xticks(np.arange(-.5, n, 1), minor=True)
    ax.set_yticks(np.arange(-.5, n, 1), minor=True)
    ax.grid(which='minor', color='white', lw=1.5, zorder=6)
    ax.set_axisbelow(False)
    ax.tick_params(which='minor', length=0)

    ax.set_xticks(range(n)); ax.set_yticks(range(n))
    ax.set_xticklabels(proteins, rotation=90, fontsize=max(8, 220 // max(n, 1)))
    ax.set_yticklabels(proteins, fontsize=max(8, 220 // max(n, 1)))
    for s in ax.spines.values():
        s.set_visible(True)

    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
    cbar.set_label('Demixing index $D/D_{max}$')

    fig.legend(handles=[Patch(fc=GREY, ec='0.35',
                              label=f'Not mixed ($\\Delta G_{{ij}} > {dg_cutoff:g}$ $k_BT$)')],
               loc='upper center', bbox_to_anchor=(0.5, 0.0), fontsize=16)
    ax.set_title('Pairwise demixing, sabari proteins', fontweight='bold', pad=15)

    fig.tight_layout()
    fig.savefig(out_pdf, bbox_inches='tight')
    plt.close(fig)

    # long-format source data
    recs = []
    for i, p1 in enumerate(proteins):
        for j, p2 in enumerate(proteins):
            if j <= i or not simulated.iat[i, j]:
                continue
            recs.append({'protein1': p1, 'protein2': p2,
                         'D_norm': d_mat.iat[i, j], 'dGij': g_mat.iat[i, j],
                         'mixed': bool(~masked.iat[i, j])})
    if recs:
        save_csv(pd.DataFrame(recs).to_dict('list'), 'figure4_extra_demixing_matrix.csv')

    n_mask = int((masked.values & simulated.values).sum() // 2)
    print(f"{out_pdf} saved: {n} proteins, "
          f"{int(np.triu(simulated.values, 1).sum())} pairs "
          f"({n_mask} greyed out at dGij > {dg_cutoff:g})")


def align_label_heights(fig, ref_gid, gid):
    """Put the figure-level panel label `gid` at the same height as the
    axes-level label `ref_gid`, after the layout has settled."""
    fig.canvas.draw()
    labels = {t.get_gid(): t for t in fig.findobj(lambda o: o.get_gid() in (ref_gid, gid))}
    if len(labels) < 2:
        return
    ref = labels[ref_gid]
    y = fig.transFigure.inverted().transform(ref.get_transform().transform(ref.get_position()))[1]
    labels[gid].set_y(y)


def create_figure4():
    """Create Figure 4."""
    DATA_DIR.mkdir(exist_ok=True)

    fig = plt.figure(figsize=(30, 15))
    gs = gridspec.GridSpec(3, 7, figure=fig, hspace=0, wspace=0.15)

    plot_fus_a1_panels(fig, gs[0:2, 0:4])
    plot_demixing_difference(fig, fig.add_subplot(gs[2, 0:4]))
    plot_network(fig, fig.add_subplot(gs[0:3, 4:7]))

    plt.tight_layout()
    align_label_heights(fig, 'label_A', 'label_D')
    plt.savefig('Figure4.pdf', bbox_inches='tight')
    print(f"Figure 4 saved as Figure4.pdf, data in {DATA_DIR}/")

if __name__ == '__main__':
    create_figure4()