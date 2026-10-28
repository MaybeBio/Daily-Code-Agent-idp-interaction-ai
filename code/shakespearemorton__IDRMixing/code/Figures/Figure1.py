import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import gaussian_kde
from scipy.signal import find_peaks

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import ScalarFormatter
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.colors import LinearSegmentedColormap, to_rgb, to_hex

import MDAnalysis as mda

DATASETS_DIR = Path("../DATASETS")
DDM_HOMOTYPIC_DIR = Path("../DDM_Homotypic")
TRAJECTORY_DIR = Path("../../calv_comp/simulations_training")
OUTPUT_DATA_DIR = Path("figure1_data")

CHOSEN_VOXEL_SIZE = 20
N_BOOTSTRAP = 5000
RNG = np.random.default_rng(0)

RC_BASE = {
    'axes.linewidth': 2, 'axes.spines.top': False, 'axes.spines.right': False,
    'xtick.major.size': 8, 'xtick.major.width': 2,
    'ytick.major.size': 8, 'ytick.major.width': 2,
    'legend.frameon': False, 'figure.dpi': 300,
}
RC_FIGURE1 = {**RC_BASE, 'font.size': 24}
plt.rcParams.update(RC_BASE)

COLORS = {'dense': '#e02b35', 'mid': '#EBECED', 'dilute': '#2066a8'}
PANEL_COLORS = {'c_den': COLORS['dense'], 'c_dil': COLORS['dilute']}

FINITE_SIZE_PROTEIN_COLORS = ['#84a98c', '#dda15e', '#8a817c']
FINITE_SIZE_BOX_LABELS = {170: '17\u00d717\u00d7170 nm', 200: '20\u00d720\u00d7200 nm', 250: '25\u00d725\u00d7250 nm'}

A1_CONDITION_COLORS = ('#2166ac', '#b2182b', '#1a9850', '#762a83', '#e08214', '#4d4d4d')


def discover_voxel_sizes(root=DDM_HOMOTYPIC_DIR):
    """Voxel sizes S with a slab_results_<S>.csv, read from disk."""
    sizes = []
    for f in sorted(Path(root).glob('slab_results_*.csv')):
        m = re.fullmatch(r'slab_results_(\d+(?:\.\d+)?)\.csv', f.name)
        if m:
            sizes.append(int(m.group(1)))
    return sorted(sizes)


def load_slab_csv(voxel_size, root=DDM_HOMOTYPIC_DIR):
    path = Path(root) / f'slab_results_{voxel_size}.csv'
    df = pd.read_csv(path)
    return df.drop_duplicates('seq_name', keep='last')

def voxel_count_column(df):
    """Name of the voxel-count column, or None if absent."""
    for col in ('num_voxels', 'n_voxels'):
        if col in df.columns:
            return col
    return None


def drop_untrusted_boxes(df):
    col = voxel_count_column(df)
    if col is None:
        return df
    return df[df[col] > df[col].min()]


def merge_with_reference(sim, ref):
    """Join simulation (DDM) and reference (GDS) data; flag rows where
    either method failed to resolve a dense phase."""
    df = ref.merge(sim, on='seq_name', suffixes=('', '_sim'))
    df = df.replace([np.inf, -np.inf], np.nan).dropna(subset=['cdil', 'cdil_sim'])
    df['single_phase'] = df['cden'].isna() | df['cden_sim'].isna()
    return df


def compute_agreement(ref_vals, sim_vals):
    """Deviation, spread, and rank/scale agreement between reference and
    simulated concentrations, computed on log10(sim/ref) since
    concentrations span orders of magnitude and errors are multiplicative."""
    ref_vals = np.asarray(ref_vals, float)
    sim_vals = np.asarray(sim_vals, float)
    ok = np.isfinite(ref_vals) & np.isfinite(sim_vals) & (ref_vals > 0) & (sim_vals > 0)
    x, y = ref_vals[ok], sim_vals[ok]

    result = dict.fromkeys(['dev_pct', 'dev_lo_pct', 'dev_hi_pct', 'spread_pct', 'rho', 'slope'], np.nan)
    result['n'] = len(x)
    if len(x) < 3:
        return result

    log_ratio = np.log10(y / x)
    median = np.median(log_ratio)
    bootstrap_medians = np.median(RNG.choice(log_ratio, size=(N_BOOTSTRAP, len(log_ratio)), replace=True), axis=1)
    lo, hi = np.percentile(bootstrap_medians, [2.5, 97.5])
    to_pct = lambda v: 100 * (10 ** v - 1)

    result['dev_pct'], result['dev_lo_pct'], result['dev_hi_pct'] = to_pct(median), to_pct(lo), to_pct(hi)
    result['spread_pct'] = to_pct(1.4826 * np.median(np.abs(log_ratio - median)))
    result['rho'] = stats.spearmanr(x, y).statistic
    result['slope'] = stats.linregress(np.log10(x), np.log10(y)).slope
    return result


def agreement_by_metric(df):
    """Agreement for c_den (two-phase systems only) and c_dil (all resolved systems)."""
    two_phase = df[~df.single_phase]
    return {
        'c_den': compute_agreement(two_phase['cden'], two_phase['cden_sim']),
        'c_dil': compute_agreement(df['cdil'], df['cdil_sim']),
    }


def load_panel_bcd_data(bulow_df):
    voxel_sizes = discover_voxel_sizes()
    results = {S: load_slab_csv(S) for S in voxel_sizes}
    # Exclude the untrusted boxes here, before any statistics are computed,
    # so panels B-G and the reported rho / n all describe the same set.
    merged = {S: drop_untrusted_boxes(merge_with_reference(df, bulow_df))
              for S, df in results.items()}

    agreement_table = pd.DataFrame([
        {'S': S, 'metric': metric, **stats_}
        for S in sorted(merged)
        for metric, stats_ in agreement_by_metric(merged[S]).items()
    ])

    df = merged.get(CHOSEN_VOXEL_SIZE)
    df_stats = agreement_by_metric(df) if df is not None else None
    chosen_slab_df = results.get(CHOSEN_VOXEL_SIZE)
    if chosen_slab_df is not None:
        chosen_slab_df = drop_untrusted_boxes(chosen_slab_df)
    return df, df_stats, agreement_table, chosen_slab_df


def build_panel_f_data(bulow_df, slab_df):
    """dG (and its error) comes from the GDS reference; frac_two_peak is a
    simulation-only quantity from the slab results. Matched on seq_name."""
    ref_cols = [c for c in ('seq_name', 'dG', 'dG_err') if c in bulow_df.columns]
    sim_cols = [c for c in ('seq_name', 'frac_two_peak', 'n_voxels') if c in slab_df.columns]
    return bulow_df[ref_cols].merge(slab_df[sim_cols], on='seq_name')


def load_panel_g_data(chosen_slab_df):
    """dG and its error for three proteins at three box sizes: the two
    explicit sizes in finite_size.csv, plus the standard 20x20x200 nm box
    already simulated for the main dataset, matched on seq_name/base."""
    finite_size_df = pd.read_csv(DATASETS_DIR / "finite_size.csv")
    proteins = sorted(finite_size_df['base'].unique())

    rows = []
    for protein in proteins:
        for box_x in (170, 250):
            match = finite_size_df[(finite_size_df.base == protein) & (finite_size_df.box_x == box_x)]
            if len(match):
                rows.append({'protein': protein, 'box_x': box_x,
                             'dG': match['dG'].iloc[0], 'dG_err': match['dG_err'].iloc[0]})
        match = chosen_slab_df[chosen_slab_df.seq_name == protein]
        if len(match):
            rows.append({'protein': protein, 'box_x': 200,
                         'dG': match['dG'].iloc[0], 'dG_err': match['dG_err'].iloc[0]})
    return pd.DataFrame(rows)


def load_panel_h_data():
    """A1 homotypic phase diagram data, read directly from A1_phase.csv
    (written by A1_Phase/run_ddm.py) -- no merging with other datasets and
    no recomputation."""
    return pd.read_csv(Path("../A1_Phase/A1_phase.csv"))


def color_shades(base_color, n):
    """n evenly spaced tints of base_color, from light to the full color."""
    rgb = np.array(to_rgb(base_color))
    fractions = np.linspace(0.7, 0.0, n)
    return [to_hex((1 - f) * rgb + f * np.ones(3)) for f in fractions]


def save_csv(data_dict, filename):
    """Save panel data to CSV, padding arrays to equal length."""
    OUTPUT_DATA_DIR.mkdir(exist_ok=True)
    max_len = max(len(v) if hasattr(v, '__len__') and not isinstance(v, str) else 1 for v in data_dict.values())
    for k, v in data_dict.items():
        if not hasattr(v, '__len__') or isinstance(v, str):
            data_dict[k] = [v] * max_len
        elif len(v) < max_len:
            data_dict[k] = list(v) + [np.nan] * (max_len - len(v))
    pd.DataFrame(data_dict).to_csv(OUTPUT_DATA_DIR / filename, index=False)


def plot_panel_A(ax_map, ax_kde):
    """Spatial concentration map (Z vs X, averaged over Y) + KDE of voxel
    concentrations for Q3LI60_1_44, both in mM."""
    npz_path = DDM_HOMOTYPIC_DIR / f'results_{CHOSEN_VOXEL_SIZE}/Q3LI60_1_44.npz'
    with np.load(npz_path) as z:
        active_dens = z['active_dens'].astype(np.float64)
        conv = float(z['conv'])

    kde = gaussian_kde(active_dens)
    extra = 0.05 * (active_dens.max() - active_dens.min())
    x_range = np.linspace(active_dens.min() - extra, active_dens.max() + extra, 4096)
    kde_y = kde(x_range)
    peaks, _ = find_peaks(kde_y, height=1e-8)
    c_dilute = x_range[peaks[0]]
    c_dense = x_range[peaks[-1]] if len(peaks) > 1 else x_range[peaks[0]]

    direc = TRAJECTORY_DIR / "Q3LI60_1_44"
    u = mda.Universe(str(direc / "top.pdb"), str(list(direc.glob("*.dcd"))[0]))
    u.trajectory[-1]
    
    # Save the last frame of the trajectory as a PDB file
    OUTPUT_DATA_DIR.mkdir(exist_ok=True)
    pdb_output_path = OUTPUT_DATA_DIR / "Q3LI60_1_44_last_frame.pdb"
    u.atoms.write(str(pdb_output_path))

    box = u.dimensions[:3]
    n_chunks = np.maximum(np.round(box / CHOSEN_VOXEL_SIZE).astype(int), 1)
    chunk_dims = box / n_chunks
    idx = np.clip((u.atoms.positions // chunk_dims).astype(int), 0, n_chunks - 1)

    counts_3d = np.zeros(n_chunks, dtype=np.float64)
    for i, j, k in idx:
        counts_3d[i, j, k] += 1
    concentration_zx = np.mean(counts_3d, axis=1).T * conv

    cmap = LinearSegmentedColormap.from_list(
        'custom', [(0.0, COLORS['dilute']), (0.15, COLORS['mid']), (1.0, COLORS['dense'])])

    ax_map.imshow(concentration_zx.T, origin='lower', extent=[0, box[2], 0, box[0]],
                  aspect='equal', cmap=cmap, vmin=x_range.min(), vmax=x_range.max())
    for i in range(n_chunks[2] + 1):
        ax_map.axvline(i * chunk_dims[2], color='k', ls=':', lw=0.5, alpha=0.5)
    for i in range(n_chunks[0] + 1):
        ax_map.axhline(i * chunk_dims[0], color='k', ls=':', lw=0.5, alpha=0.5)
    ax_map.set_xlabel('Z (\u00c5)')
    ax_map.set_ylabel('X (\u00c5)')

    for i in range(len(x_range) - 1):
        norm = (x_range[i] - x_range.min()) / (x_range.max() - x_range.min())
        ax_kde.fill_between(x_range[i:i + 2], 1e-8, kde_y[i:i + 2], color=cmap(norm), alpha=0.6)
    ax_kde.plot(x_range, kde_y, 'k-', lw=1.5)
    ax_kde.axvline(c_dilute, color=COLORS['dilute'], ls='--', lw=3, alpha=0.8)
    if len(peaks) > 1:
        ax_kde.axvline(c_dense, color=COLORS['dense'], ls='--', lw=3, alpha=0.8)

    #ax_kde.set_yscale('log')
    ax_kde.set_ylim(bottom=1e-8, top=kde_y.max() *1.2)

    y_max = kde_y.max()
    ax_kde.text(c_dilute + 15, y_max * 1.2, f'{c_dilute:.1f}', ha='center', va='bottom',
                color=COLORS['dilute'], fontweight='bold')
    if len(peaks) > 1:
        ax_kde.text(c_dense - 25, y_max * 1.2, f'{c_dense:.1f}', ha='center', va='bottom',
                    color=COLORS['dense'], fontweight='bold')
    ax_kde.set_xlabel('Concentration (mM)')
    ax_kde.set_ylabel('Prob. Density')

    save_csv({'concentration_zx_flat': concentration_zx.flatten(), 'kde_x': x_range, 'kde_y': kde_y,
              'c_dilute': c_dilute, 'c_dense': c_dense}, 'panelA_composite_data.csv')


def error_bar_pair(part, value_col, err_col):
    """Asymmetric (lower, upper) error array, or None if the column is absent.
    The lower bar is capped just short of the point so a log axis never
    sees a non-positive limit."""
    if err_col not in part or not len(part):
        return None
    err = np.nan_to_num(part[err_col].to_numpy(float), nan=0.0)
    if not np.any(err > 0):
        return None
    val = part[value_col].to_numpy(float)
    return np.vstack([np.minimum(err, 0.999 * val), err])


def plot_panel_B(ax, table, chosen_voxel_size):
    """Median deviation from the GDS vs voxel side length."""
    ax.axhline(0, color='black', lw=1.5, alpha=0.6)
    for metric, colour in PANEL_COLORS.items():
        d = table[table.metric == metric].sort_values('S')
        ax.fill_between(d.S, d.dev_lo_pct, d.dev_hi_pct, color=colour, alpha=0.2, lw=0)
        ax.plot(d.S, d.dev_pct, 'o-', color=colour, lw=3, ms=10, mec='black', mew=0.8)
    ax.axvline(chosen_voxel_size, color='grey', ls='--', lw=3)
    ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.set_xlabel('Voxel Side Length (\u00c5)')
    ax.set_ylabel('Deviation from GDS (%)')
    ax.set_ylim(-3, 3)
    ax.legend(handles=[
        Line2D([0], [0], color=PANEL_COLORS['c_den'], lw=3, label='$c_{den}$'),
        Line2D([0], [0], color=PANEL_COLORS['c_dil'], lw=3, label='$c_{dil}$'),
        Line2D([0], [0], color='grey', ls='--', lw=2, label=f'{chosen_voxel_size:g} \u00c5'),
    ], fontsize=20, loc='upper right')


def plot_panel_scatter(ax, df, col, colour, label, stats_, grey=True, errorbars=True):
    """Log-log scatter of GDS vs DDM concentrations (panels C, D)."""
    finite = df[[col, col + '_sim']].dropna()
    finite = finite[(finite > 0).all(axis=1)]
    span = np.log10([finite.min().min(), finite.max().max()])
    lim = 10 ** (span + 0.05 * np.ptp(span) * np.array([-1, 1]))
    ax.plot(lim, lim, 'k--', alpha=0.5, lw=2, zorder=1)

    err_ref = 'e' + col[1:]
    err_sim = col + '_sem'

    parts = [(df[df.single_phase], 'grey', 2)] if grey else []
    parts.append((df[~df.single_phase], colour, 3))
    for part, c, z in parts:
        if not len(part):
            continue
        ax.errorbar(part[col], part[col + '_sim'],
                    xerr=error_bar_pair(part, col, err_ref) if errorbars else None,
                    yerr=error_bar_pair(part, col + '_sim', err_sim) if errorbars else None,
                    fmt='o', color=c, ecolor=c, alpha=0.75, ms=13, mec='black', mew=1,
                    elinewidth=1, capsize=3, ls='none', zorder=z)

    ax.text(0.05, 0.95, f'$\\rho$ = {stats_["rho"]:.3f}   (n={stats_["n"]})',
            transform=ax.transAxes, va='top', fontsize=20,
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    ax.set_xlabel(f'${label}$ using GDS (mM)')
    ax.set_ylabel(f'${label}$ using DDM (mM)')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_aspect('equal', adjustable='box')
    ax.set_xlim(*lim)
    ax.set_ylim(*lim)


def plot_panel_E(ax, df):
    """Fraction of frames in which the per-frame KDE resolves two peaks,
    plotted against dG. Points are uniform grey."""
    df_clean = df.dropna(subset=['dG', 'frac_two_peak', 'n_voxels']).copy()
    df_clean['frac_two_peak'] = 100 * df_clean['frac_two_peak']
    ax.errorbar(df_clean['dG'], df_clean['frac_two_peak'],
                xerr=None,
                fmt='none', ecolor='black', elinewidth=1, capsize=3, zorder=2)
    points = ax.scatter(df_clean['dG'], df_clean['frac_two_peak'], c='grey',  s=170,
                         edgecolor='black', linewidth=1, alpha=0.85, zorder=3)
    ax.axvline(-3, color='black', ls='--', lw=2, alpha=0.6)
    ax.set_xlabel(r'$\Delta G$')
    ax.axvspan(-3, 0, color='#99b9d3', alpha=0x4d / 255, lw=0, zorder=0)
    ax.set_xticks([-10, -8, -6, -4, -2, 0])
    ax.set_ylabel('Frames with 2 peaks (%)')
    ax.set_xlim(-10.5, 0.5)
    ax.set_ylim(40, 105)


def plot_panel_F(ax, df):
    """Bar chart of dG per protein across three simulation box sizes.
    Each protein gets a hue; box size is encoded as a light-to-dark shade
    of that hue. Protein names are colour-matched, bold x-tick labels
    rather than a second legend, to keep the panel compact."""
    proteins = sorted(df['protein'].unique())
    box_sizes = sorted(df['box_x'].unique())
    bar_width = 0.8 / len(box_sizes)
    offsets = (np.arange(len(box_sizes)) - (len(box_sizes) - 1) / 2) * bar_width

    for p_idx, protein in enumerate(proteins):
        colours = color_shades(FINITE_SIZE_PROTEIN_COLORS[p_idx], len(box_sizes))
        for b_idx, box_x in enumerate(box_sizes):
            row = df[(df.protein == protein) & (df.box_x == box_x)]
            if not len(row):
                continue
            ax.bar(p_idx + offsets[b_idx], row['dG'].iloc[0], width=bar_width * 0.9,
                   yerr=row['dG_err'].iloc[0], color=colours[b_idx], edgecolor='black',
                   linewidth=1, capsize=3, ecolor='black')

    ax.set_xticks(range(len(proteins)))
    ax.set_xticklabels(proteins, fontsize=18, fontweight='bold', rotation=20, ha='right')
    for tick, colour in zip(ax.get_xticklabels(), FINITE_SIZE_PROTEIN_COLORS):
        tick.set_color(colour)

    ax.axhline(0, color='black', lw=1.5)
    ax.set_ylim(0, -11.9)
    ax.set_ylabel(r'$\Delta G$')

    legend_shades = color_shades('#808080', len(box_sizes))
    ax.legend(handles=[Patch(facecolor=c, edgecolor='black', label=FINITE_SIZE_BOX_LABELS[box_x])
                        for box_x, c in zip(box_sizes, legend_shades)],
              title='Box size', fontsize=16, title_fontsize=16, loc='upper left')


def plot_panel_G(ax, df, colors=A1_CONDITION_COLORS, min_dg=0.0):
    is_hottest = df['temp_K'] == df.groupby('condition')['temp_K'].transform('max')
    coex = df['two_phase'] & (df['dG_kT'] <= -abs(min_dg)) & ~is_hottest

    volumes = {c: float(g['box_x_nm'].iloc[0] ** 2 * g['box_z_nm'].iloc[0])
               for c, g in df.groupby('condition')}
    conditions = sorted(volumes, key=volumes.get)

    for condition, colour in zip(conditions, colors):
        g = df[df.condition == condition].sort_values('temp_K')
        paired = g[coex.loc[g.index]]
        single = g[~coex.loc[g.index]].copy()
        hottest = single['temp_K'] == g['temp_K'].max()
        single.loc[hottest, 'c_dil_mM'] = single.loc[hottest, 'c_mean_mM']
        single.loc[hottest, 'c_dil_err'] = 0.0

        for _, row in paired.iterrows():
            ax.plot([row.c_dil_mM, row.c_den_mM], [row.temp_K, row.temp_K],
                    '-', color=colour, lw=1, alpha=0.45, zorder=1)
        if len(paired) > 1:
            ax.plot(paired.c_dil_mM, paired.temp_K, '-', color=colour, lw=2.5, zorder=2)
            ax.plot(paired.c_den_mM, paired.temp_K, '-', color=colour, lw=2.5, zorder=2)
        for col, err in (('c_dil_mM', 'c_dil_err'), ('c_den_mM', 'c_den_err')):
            ax.errorbar(paired[col], paired.temp_K, xerr=paired[err].fillna(0),
                        fmt='o', color=colour, ecolor=colour, ms=13, mec='black',
                        mew=1, elinewidth=1, capsize=3, zorder=3)

        if len(single):
            ax.errorbar(single.c_dil_mM, single.temp_K, xerr=single.c_dil_err.fillna(0),
                        fmt='o', color=colour, ecolor=colour, ms=13, mfc='white',
                        mec=colour, mew=2, elinewidth=1, capsize=3, zorder=3)

    handles = [Line2D([0], [0], marker='o', ls='-', color=c, label=cond)
               for cond, c in zip(conditions, colors)]
    handles.append(Line2D([0], [0], marker='o', ls='none', color='0.3',
                          mfc='white', mew=2, label='Single phase'))
    ax.legend(handles=handles, fontsize=16, loc='upper right')
    ax.set_ylim(303,320)
    #ax.set_xscale('log')
    ax.set_xlabel('Concentration (mM)')
    ax.set_ylabel('Temperature (K)')


def add_label(ax, letter):
    ax.text(-0.1, 1.05, letter, transform=ax.transAxes, fontsize=28,
             fontweight='bold', va='bottom', ha='right')


def report_missing_dense(df):
    """Count and median dG (from the GDS reference and from the DDM) for
    systems where a dense phase is missing from the GDS, from the DDM, or
    from both."""
    dG_cols = {'GDS': 'dG', 'DDM': 'dG_sim'}
    groups = {
        'no cden from GDS': df['cden'].isna(),
        'no cden from DDM': df['cden_sim'].isna(),
        'no cden from either': df['cden'].isna() & df['cden_sim'].isna(),
    }
    for label, mask in groups.items():
        parts = []
        for source, col in dG_cols.items():
            if col in df.columns:
                vals = df.loc[mask, col].dropna()
                med = f"{vals.median():.2f}" if len(vals) else "n/a"
                parts.append(f"median dG ({source})={med} [n={len(vals)}]")
            else:
                parts.append(f"median dG ({source})=column missing")
        print(f"{label}: n={int(mask.sum())}, " + ", ".join(parts))


def build_figure(panel_f_df, panel_g_df, panel_h_df, panel_bcd_df, panel_bcd_stats, panel_bcd_table):
    plt.rcParams.update(RC_FIGURE1)

    fig = plt.figure(figsize=(30, 14))
    gs = gridspec.GridSpec(2, 4, figure=fig, hspace=0.5, wspace=0.4)

    ax_map = fig.add_subplot(gs[0, 0])
    ax_map.set_anchor('N')
    add_label(ax_map, 'A')
    ax_kde = fig.add_subplot(gs[1, 0])
    plot_panel_A(ax_map, ax_kde)

    ax_b = fig.add_subplot(gs[0, 1])
    ax_b.set_anchor('N')
    add_label(ax_b, 'B')
    plot_panel_B(ax_b, panel_bcd_table, CHOSEN_VOXEL_SIZE)

    ax_c = fig.add_subplot(gs[0, 2])
    ax_c.set_anchor('N')
    add_label(ax_c, 'C')
    plot_panel_scatter(ax_c, panel_bcd_df, 'cden', PANEL_COLORS['c_den'], 'c_{den}',
                        panel_bcd_stats['c_den'], grey=False)

    ax_d = fig.add_subplot(gs[0, 3])
    ax_d.set_anchor('N')
    add_label(ax_d, 'D')
    plot_panel_scatter(ax_d, panel_bcd_df, 'cdil', PANEL_COLORS['c_dil'], 'c_{dil}',
                        panel_bcd_stats['c_dil'])

    ax_f = fig.add_subplot(gs[1, 1])
    ax_f.set_anchor('N')
    add_label(ax_f, 'E')
    plot_panel_E(ax_f, panel_f_df)

    ax_g = fig.add_subplot(gs[1, 2])
    ax_g.set_anchor('N')
    add_label(ax_g, 'F')
    plot_panel_F(ax_g, panel_g_df)

    ax_h = fig.add_subplot(gs[1, 3])
    ax_h.set_anchor('N')
    add_label(ax_h, 'G')
    plot_panel_G(ax_h, panel_h_df)

    return fig


def main():
    bulow_df = pd.read_csv(DATASETS_DIR / "bulow_df.csv")
    panel_bcd_df, panel_bcd_stats, panel_bcd_table, chosen_slab_df = load_panel_bcd_data(bulow_df)
    print(f"panel C/D: {len(panel_bcd_df)} simulations retained "
          f"(n_cden={panel_bcd_stats['c_den']['n']}, n_cdil={panel_bcd_stats['c_dil']['n']}, "
          f"rho_cden={panel_bcd_stats['c_den']['rho']:.3f}, rho_cdil={panel_bcd_stats['c_dil']['rho']:.3f})")
    report_missing_dense(panel_bcd_df)

    panel_f_df = build_panel_f_data(bulow_df, chosen_slab_df)
    panel_g_df = load_panel_g_data(chosen_slab_df)
    panel_h_df = load_panel_h_data()
    fig = build_figure(panel_f_df, panel_g_df, panel_h_df, panel_bcd_df, panel_bcd_stats, panel_bcd_table)
    fig.savefig("Figure1.pdf", bbox_inches='tight',dpi = 300)


if __name__ == "__main__":
    main()