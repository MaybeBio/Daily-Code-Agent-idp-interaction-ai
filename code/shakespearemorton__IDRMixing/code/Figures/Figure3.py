#!/usr/bin/env python3
"""Figure 3: Chemical grammar and demixing analysis (2x3 layout A-F)"""
import numpy as np
import pandas as pd
import colormaps as cmaps
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
from pathlib import Path
from scipy.stats import linregress
from scipy import stats
from sklearn.mixture import GaussianMixture
from sklearn.preprocessing import StandardScaler

plt.rcParams.update({
    'font.size': 24, 'axes.linewidth': 2, 'axes.spines.top': False, 'axes.spines.right': False,
    'xtick.major.size': 8, 'xtick.major.width': 2, 'ytick.major.size': 8, 'ytick.major.width': 2,
    'legend.frameon': False, 'figure.dpi': 300
})

DATA_DIR = Path("figure3_data")
POSITIVE_AA, NEGATIVE_AA = set('KR'), set('DE')
MIN_DENSE_FRAC = 0.05   # two phases only if >= 5% of voxels sit at or above c_den

GIN_CATEGORY_MAP = {
    7: 'C', 9: 'C', 18: 'C', 19: 'C', 23: 'C', 25: 'C', 26: 'C', 3: 'C', 8: 'C', 17: 'C', 24: 'C', 29: 'C',
    0: 'P', 1: 'P', 2: 'P', 4: 'P', 5: 'P', 6: 'P', 11: 'P', 12: 'P', 13: 'P', 14: 'P', 15: 'P', 16: 'P', 21: 'P', 22: 'P', 27: 'P',
    10: 'H', 20: 'H', 28: 'H'
}
INTERACTION_PALETTE = {'C-C': '#2E86AB', 'C-H': '#A23B72', 'C-P': '#F18F01', 'H-H': '#C73E1D', 'H-P': '#708238', 'P-P': '#1B998B'}


def assign_bin(comp):
    """Assign composition to bins."""
    if 85 <= comp <= 100: return '85-95'
    if 65 <= comp < 85: return '65-85'
    if 35 <= comp < 65: return '50'
    if 15 <= comp < 35: return '15-35'
    if 0 <= comp < 15: return '5-15'
    return None

def get_significance_stars(p_value):
    """Convert p-value to significance stars."""
    if p_value < 0.001: return '***'
    if p_value < 0.01: return '**'
    if p_value < 0.05: return '*'
    return ''

def cohens_d(x, y):
    """Calculate Cohen's d for effect size."""
    nx, ny = len(x), len(y)
    dof = nx + ny - 2
    return (np.mean(x) - np.mean(y)) / np.sqrt(((nx-1)*np.std(x, ddof=1)**2 + (ny-1)*np.std(y, ddof=1)**2) / dof)

def save_csv(data_dict, filename):
    DATA_DIR.mkdir(exist_ok=True)
    max_len = max(len(v) if hasattr(v, '__len__') and not isinstance(v, str) else 1 for v in data_dict.values())
    for k, v in data_dict.items():
        if not hasattr(v, '__len__') or isinstance(v, str): data_dict[k] = [v] * max_len
        elif len(v) < max_len: data_dict[k] = list(v) + [np.nan] * (max_len - len(v))
    pd.DataFrame(data_dict).to_csv(DATA_DIR / filename, index=False)

def get_interaction_type(GIN_group_p1, GIN_group_p2):
    if pd.isna(GIN_group_p1) or pd.isna(GIN_group_p2): return 'Other'
    cat1, cat2 = GIN_CATEGORY_MAP.get(int(GIN_group_p1), 'Other'), GIN_CATEGORY_MAP.get(int(GIN_group_p2), 'Other')
    return 'Other' if 'Other' in [cat1, cat2] else '-'.join(sorted([cat1, cat2]))

def calc_net_charge(seq):
    return sum(1 if aa in POSITIVE_AA else -1 for aa in str(seq).upper() if aa in POSITIVE_AA | NEGATIVE_AA) if pd.notna(seq) else 0

def calc_ocn(seq1,seq2):
    seq1 = str(seq1.upper())
    seq2 = str(seq2.upper())
    seq1_pos = sum(1 for aa in seq1 if aa in POSITIVE_AA)
    seq1_neg = sum(1 for aa in seq1 if aa in NEGATIVE_AA)
    seq2_pos = sum(1 for aa in seq2 if aa in POSITIVE_AA)
    seq2_neg = sum(1 for aa in seq2 if aa in NEGATIVE_AA)
    return((seq1_pos*seq2_neg + seq1_neg*seq2_pos)-(seq1_pos*seq2_pos + seq1_neg*seq2_neg))


def calc_total_charged(seq):
    return sum(1 for aa in str(seq).upper() if aa in POSITIVE_AA | NEGATIVE_AA) if pd.notna(seq) else 0

def apply_two_phase_gate(df):
    """A voxel map resolves two phases only when at least MIN_DENSE_FRAC of
    its voxels sit at or above c_den. Otherwise c_dil is the mean over all
    voxels and no dense value (and so no dGij) is defined."""
    df = df.copy()
    df['two_phase'] = df['frac_vox_ge_cden'] >= MIN_DENSE_FRAC
    one = ~df['two_phase']
    for dil in ('cdil', 'cdil1', 'cdil2'):
        df.loc[one, dil] = df.loc[one, f'{dil}_onephase']
    df.loc[one, ['cden', 'cden1', 'cden2', 'dGij', 'dGij_err']] = np.nan
    return df

def load_data():
    df_dem = apply_two_phase_gate(pd.read_csv("../DDM_Heterotypic/heterotypic_main.csv"))
    print(f"{df_dem['two_phase'].sum()} of {len(df_dem)} systems two-phase "
          f"(>= {MIN_DENSE_FRAC:.0%} of voxels at or above c_den)")
    df_dem = df_dem[df_dem['two_phase']]  # gated-out (single-phase) maps are not analysed
    df_prot = pd.read_csv("../DATASETS/gin_samples.csv")
    df_prot['n_charged'] = df_prot['fasta'].apply(lambda s: sum(1 if aa in 'KR' else -1 for aa in s if aa in 'KRED') if pd.notna(s) else np.nan)
    df_dem['composition_bin'] = df_dem['comp1'].apply(assign_bin)
    features = ['mean_lambda', 'scd', 'shd', 'faro']
    feature_maps = {f: dict(zip(df_prot['seq_name'], df_prot[f])) for f in features}

    df = df_dem.copy()
    for p, suf in [('seq_name1', '_p1'), ('seq_name2', '_p2')]:
        df = pd.merge(df, df_prot, left_on=p, right_on='seq_name', how='left', suffixes=('', suf))
        df = df.drop(columns=['seq_name'])
        df = df.rename(columns={c: c + suf for c in df_prot.columns if c != 'seq_name' and c in df.columns and not c.endswith(('_p1', '_p2'))})
    df = df[(df['comp1'] < 95) & (df['comp2'] < 95)].copy()
    return df, feature_maps

DF = None
FEATURE_MAPS = None

def init_data():
    """Load the merged demixing/sequence dataset into the module-level globals."""
    global DF, FEATURE_MAPS
    DF, FEATURE_MAPS = load_data()


def print_pp_samples():
    """Print entry names, dGij, dGi, dGj, and Dnorm for all P-P interaction category samples."""
    if DF is None: init_data()
    df = DF.copy()
    df['interaction_type'] = df.apply(lambda r: get_interaction_type(r['GIN_group_p1'], r['GIN_group_p2']), axis=1)
    pp_df = df[df['interaction_type'] == 'P-P'].copy()

    print("\n" + "="*85)
    print("P-P Interaction Category Samples")
    print("="*85)
    print(f"{'seq_name1':<22} {'seq_name2':<22} {'dGij':<10} {'dGi':<10} {'dGj':<10} {'Dnorm':<10}")
    print("-"*85)

    for _, row in pp_df.iterrows():
        s1 = str(row['seq_name1'])
        s2 = str(row['seq_name2'])
        dGij = row['dGij']
        dGi = row['dG_p1']
        dGj = row['dG_p2']
        Dnorm = row['D_norm']
        print(f"{s1:<22} {s2:<22} {dGij:<10.4f} {dGi:<10.4f} {dGj:<10.4f} {Dnorm:<10.4f}")

    print("="*85 + "\n")


def analyze_gin_demixing():
    """Print top 10 individual GIN groups with lowest median demixing index across all their pairs."""
    if DF is None: init_data()

    # comp1 < 95 and comp2 < 95 are already enforced by load_data(); only the
    # dGij cutoff (matching Panel A/C's selection) needs to be applied here.
    df_filt = DF[DF['dGij'] < -3].copy()
    df_filt = df_filt.dropna(subset=['GIN_group_p1', 'GIN_group_p2', 'D_norm'])
    df_filt['interaction_type'] = df_filt.apply(lambda r: get_interaction_type(r['GIN_group_p1'], r['GIN_group_p2']), axis=1)
    df_filt = df_filt[df_filt['interaction_type'] != 'Other']
    df_filt = df_filt.drop_duplicates(subset=['seq_name1', 'seq_name2', 'comp1'])

    df1 = df_filt[['GIN_group_p1', 'D_norm']].rename(columns={'GIN_group_p1': 'gin'})
    df2 = df_filt[['GIN_group_p2', 'D_norm']].rename(columns={'GIN_group_p2': 'gin'})
    df_long = pd.concat([df1, df2], ignore_index=True)
    df_long['gin'] = df_long['gin'].astype(int)

    gin_stats = df_long.groupby('gin')['D_norm'].agg(['median', 'count'])
    gin_stats = gin_stats.sort_values('median').head(10)

    print("\n" + "="*60)
    print("Top 10 GIN Groups with LOWEST Median Demixing Index")
    print("(across all pairs involving each GIN group)")
    print("="*60)
    print(f"{'GIN Group':<12} {'Median':<12} {'N Samples':<10}")
    print("-"*60)
    for gin, row in gin_stats.iterrows():
        print(f"{gin:<12} {row['median']:<12.4f} {int(row['count']):<10}")
    print("="*60 + "\n")

    # Second table: Interaction type groups (C-C, C-H, C-P, H-H, H-P, P-P)
    stats_interaction = df_filt.groupby('interaction_type')['D_norm'].agg(['median', 'count'])
    stats_interaction = stats_interaction.sort_values('median')

    print("="*60)
    print("Median Demixing Index by Interaction Type")
    print("(C=Charged, H=Hydrophobic, P=Polar)")
    print("="*60)
    print(f"{'Interaction':<15} {'Median':<12} {'N Samples':<10}")
    print("-"*60)
    for interaction, row in stats_interaction.iterrows():
        print(f"{interaction:<15} {row['median']:<12.4f} {int(row['count']):<10}")
    print("="*60 + "\n")


def panel_A(ax):
    """Panel A: Demixing index distributions by interaction type (dGij < -3)."""
    df_A = DF[DF['dGij'] < -3].copy()
    df_A['interaction_type'] = df_A.apply(lambda r: get_interaction_type(r['GIN_group_p1'], r['GIN_group_p2']), axis=1)
    df_A = df_A[df_A['interaction_type'] != 'Other'].drop_duplicates(subset=['seq_name1', 'seq_name2', 'comp1'])

    order = df_A.groupby('interaction_type')['D_norm'].median().sort_values().index.tolist()

    parts = ax.violinplot([df_A[df_A['interaction_type'] == it]['D_norm'].values for it in order],
                        positions=range(len(order)), widths=0.6, showmeans=False, showmedians=True)

    for pc, it in zip(parts['bodies'], order):
        pc.set_facecolor(INTERACTION_PALETTE.get(it, 'lightgrey'))
        pc.set_alpha(0.8)
        pc.set_edgecolor('black')
        pc.set_linewidth(2)

    parts['cmedians'].set_edgecolor('white')
    parts['cmedians'].set_linewidth(5)

    for partname in ('cbars', 'cmaxes', 'cmins'):
        parts[partname].set_edgecolor('black')
        parts[partname].set_linewidth(2)

    ax.set_ylabel('Demixing Index'); ax.set_xlabel('Interaction Type')
    ax.set_xticks(range(len(order))); ax.set_xticklabels(order); ax.set_ylim(-0.05, 1.05); ax.grid(False)

    for i, it in enumerate(order):
        ax.text(i, df_A[df_A['interaction_type'] == it]['D_norm'].max() + 0.05,
                f'n={len(df_A[df_A["interaction_type"] == it])}', ha='center', va='bottom', fontsize=22)

    save_csv({'interaction_type': df_A['interaction_type'].values, 'demixing': df_A['D_norm'].values}, 'panelC_interaction.csv')


def panel_B(ax):
    """Panel B: Equimolar vs non-equimolar heterotypic dGij for matched pairs, colored by GMM cluster."""
    equimolar = DF[(DF['composition_bin'] == '50') & DF['dGij'].notna()]
    alt = DF[(DF['composition_bin'] != '50') & DF['dGij'].notna()]

    plot_data = []
    for _, line in equimolar.iterrows():
        matches = alt[(alt['seq_name1'] == line['seq_name1']) & (alt['seq_name2'] == line['seq_name2'])]
        for _, match in matches.iterrows():
            ncpr_sum = abs(line['ncpr_p1'] + line['ncpr_p2'])
            ocn_pair = calc_ocn(line['fasta_p1'], line['fasta_p2'])
            if not np.isfinite(ncpr_sum):
                continue
            plot_data.append([line['dGij'], match['dGij'], ncpr_sum, line['D_norm'],
                               line['dG_p1'], line['dG_p2'], line['seq_name1'], line['seq_name2'],ncpr_sum])

    df_B = np.array(plot_data, dtype=object)

    components = 2
    colors = ['#84a98c', '#354f52']
    gmix = GaussianMixture(n_components=components).fit_predict(df_B[:, :2].astype(float))

    ax.set_xlabel(r'$\Delta G_{ij}$ (50:50) ($k_B T$)')
    ax.set_ylabel(r'$\Delta G_{ij}$ (Unequal) ($k_B T$)')
    ax.set_xlim(0.1, -10)
    ax.set_ylim(0.1, -10)
    ax.plot([-10, 0], [-10, 0], 'k--', linewidth=1.0, alpha=0.5)

    point_colors = [colors[i] for i in gmix]
    ax.scatter(df_B[:, 0].astype(float), df_B[:, 1].astype(float),
            c=point_colors, s=100, edgecolor='k', linewidth=0.3, alpha=0.8)
    ax.grid(False)
    # cluster labels as text in each cluster's colour, in the empty
    # strong-at-50:50, weak-elsewhere corner
    for comp in range(components):
        cluster_ncpr = np.mean(df_B[gmix == comp][:, -1].astype(float))
        cluster_std = np.std(df_B[gmix == comp][:, -1].astype(float))
        print(f'Cluster {comp}: mean NCPR sum = {cluster_ncpr:.2f} ± {cluster_std:.2f}')
        ax.annotate(f'$|\Sigma NCPR|$ = {cluster_ncpr:.2f} ± {cluster_std:.2f}',
                    xy=(0.97, 0.04 + comp * 0.08), xycoords='axes fraction', ha='right', va='bottom',
                    color=colors[comp], fontweight='bold')


def panel_C(ax):
    """Panel C: Feature differences (z-scored) between mixing vs demixing pairs at equimolar composition."""
    df_C = DF[DF['comp1'] == 50].copy()

    features = ['mean_lambda', 'scd', 'shd', 'faro']
    feature_data = {f: [] for f in features}
    group1_mask = []
    group2_mask = []

    for _, row in df_C.iterrows():
        p1, p2 = row['seq_name1'], row['seq_name2']
        dG_p1, dG_p2 = row['dG_p1'], row['dG_p2']

        for feature in features:
            f1 = FEATURE_MAPS[feature][p1]
            f2 = FEATURE_MAPS[feature][p2]
            feature_data[feature].append(abs(f1 - f2))

        demix_idx = row['D_norm']
        dGij = row['dGij']

        if demix_idx < 0.13 and dGij < -3:
            group1_mask.append(True)
            group2_mask.append(False)
        elif demix_idx >= 0.13 and (dG_p1 < -3 or dG_p2 < -3):
            group1_mask.append(False)
            group2_mask.append(True)
        else:
            group1_mask.append(False)
            group2_mask.append(False)

    group1_mask = np.array(group1_mask)
    group2_mask = np.array(group2_mask)

    print(f"Group 1 (will mix): {np.sum(group1_mask)} pairs")
    print(f"Group 2 (will demix): {np.sum(group2_mask)} pairs")
    scaler = StandardScaler()
    plot_data_list = []

    for feature in features:
        data = np.array(feature_data[feature])
        g1_data = data[group1_mask]
        g2_data = data[group2_mask]

        combined = np.concatenate([g1_data, g2_data]).reshape(-1, 1)
        scaler.fit(combined)

        g1_scaled = scaler.transform(g1_data.reshape(-1, 1)).flatten()
        g2_scaled = scaler.transform(g2_data.reshape(-1, 1)).flatten()

        for val in g1_scaled:
            plot_data_list.append([feature.upper(), 'Mix', val])
        for val in g2_scaled:
            plot_data_list.append([feature.upper(), 'Demix', val])

    plot_df = pd.DataFrame(plot_data_list, columns=['Feature', 'Group', 'Value'])

    mix_color = '#EF7071'
    demix_color = '#904D4E'
    palette = {'Mix': mix_color, 'Demix': demix_color}

    sns.boxplot(data=plot_df, x='Feature', y='Value', hue='Group',
                palette=palette, ax=ax, width=0.6, linewidth=1,
                fliersize=10, showcaps=True, boxprops=dict(alpha=1))

    for i, feature in enumerate([f.upper() for f in features]):
        for j, group in enumerate(['Mix', 'Demix']):
            data_subset = plot_df[(plot_df['Feature'] == feature) & (plot_df['Group'] == group)]['Value']
            median = data_subset.median()
            x_pos = i + (j - 0.5) * 0.3
            ax.plot([x_pos - 0.12, x_pos + 0.12], [median, median],
                    color='white', linewidth=2, zorder=10)

    p_values = {}
    effect_sizes = {}

    for feature in features:
        data = np.array(feature_data[feature])
        g1_data = data[group1_mask]
        g2_data = data[group2_mask]

        t_stat, p = stats.ttest_ind(g1_data, g2_data)
        p_values[feature] = p

        d = cohens_d(g1_data, g2_data)
        effect_sizes[feature] = d

    for i, feature in enumerate(features):
        sig = get_significance_stars(p_values[feature])
        d = effect_sizes[feature]

        feat_data = plot_df[plot_df['Feature'] == feature.upper()]['Value']
        y_max = feat_data.max()

        if sig:
            x1, x2 = i - 0.15, i + 0.15
            bracket_y = y_max + 0.2
            ax.plot([x1, x1, x2, x2], [bracket_y - 0.05, bracket_y, bracket_y, bracket_y - 0.05],
                    'k-', linewidth=1)
            ax.text(i, bracket_y + 0.1, sig, ha='center', va='bottom', fontweight='bold')

        ax.text(i, plot_df['Value'].min() - 0.4, f'd={d:.2f}', ha='center', va='top',
                 style='italic', color='#555555')

    features_labels = ['$\hat\lambda$', 'SCD', 'SHD', '$F_{aro}$']
    ax.set_ylabel('Feature Difference (z-score)')
    ax.set_xlabel('')
    ax.set_xticklabels(features_labels, ha='center')
    ax.grid(False)
    ax.axhline(y=0, color='k', linewidth=0.8, alpha=0.4, linestyle=':')
    ax.set_ylim(plot_df['Value'].min() - 1, plot_df['Value'].max() + 2)
    handles, labels = ax.get_legend_handles_labels()
    count_labels = {'Mix': f'Mix (n={group1_mask.sum()})', 'Demix': f'Demix (n={group2_mask.sum()})'}
    ax.legend(handles, [count_labels.get(l, l) for l in labels], loc='upper right', framealpha=0.95, edgecolor='black')


def panel_D(ax):
    """Panel D: Distribution of |dG_ii - dG_jj| for mixing vs demixing pairs across all compositions,
    using the same mix/demix criteria as Panel C: mix = D_norm <= 0.13 and dGij < -3;
    demix = D_norm > 0.13 and (dG_p1 < -3 or dG_p2 < -3). Pairs matching neither are excluded."""
    dfx_D = DF.dropna(subset=['dG_p1', 'dG_p2', 'D_norm', 'dGij']).copy()
    dfx_D['dG_deviation'] = abs(dfx_D['dG_p1'] - dfx_D['dG_p2'])

    mix_color = '#EF7071'
    demix_color = '#904D4E'

    mix_mask = (dfx_D['D_norm'] <= 0.13) & (dfx_D['dGij'] < -3)
    demix_mask = (dfx_D['D_norm'] > 0.13) & ((dfx_D['dG_p1'] < -3) | (dfx_D['dG_p2'] < -3))

    mix_df = dfx_D[mix_mask]
    demix_df = dfx_D[demix_mask]

    # Empirical CDF of |dG_ii - dG_jj| over mixing pairs: smallest x with CDF(x) >= 0.5
    mix_dev = np.sort(mix_df['dG_deviation'].values)
    mix_cdf = np.arange(1, len(mix_dev) + 1) / len(mix_dev)
    x50 = mix_dev[np.searchsorted(mix_cdf, 0.5)]
    demix_frac_below = (demix_df['dG_deviation'] < x50).mean()
    print(f"Panel D: 50% of mixing pairs (n={len(mix_df)}) have |dG_ii - dG_jj| < {x50:.2f} kBT "
          f"({demix_frac_below:.0%} of demixing pairs, n={len(demix_df)}, fall below this)")
    print(f"  -> well-mixed condensates predominately form when IDRs have similar "
          f"phase-separation propensities (difference < {x50:.1f} kBT)")

    bins = np.arange(0, 10.5, 0.5)
    demix_counts, _ = np.histogram(demix_df['dG_deviation'], bins=bins)
    mix_counts, _ = np.histogram(mix_df['dG_deviation'], bins=bins)
    centers = (bins[:-1] + bins[1:]) / 2

    ax.bar(centers, demix_counts, width=0.45, color=demix_color, edgecolor='black',
        label=f'Demix, n={len(demix_df)}')
    ax.bar(centers, -mix_counts, width=0.45, color=mix_color, edgecolor='black',
        label=f'Mix, n={len(mix_df)}')
    ax.axhline(0, color='black', linewidth=1)
    ax.set_yticklabels([f'{abs(int(t))}' for t in ax.get_yticks()])
    ax.set_xlabel('|ΔG$_{ii}$ - ΔG$_{jj}$| ($k_BT$)'); ax.set_ylabel('Count')
    ax.set_xlim(-0.1, 10)
    ax.legend(loc='lower right', framealpha=0.95, edgecolor='black')

    # Panels E and F both reserve space for a colorbar, which shrinks their plot area.
    # Reserve the same space here (then hide it) so this panel renders at the same
    # physical size as its neighbors, and only then fix the box to a square.
    dummy_mappable = plt.cm.ScalarMappable(cmap='Greys')
    dummy_mappable.set_array([])
    cbar_d = plt.colorbar(dummy_mappable, ax=ax, pad=0.02, shrink=0.5)
    cbar_d.ax.set_visible(False)
    ax.set_box_aspect(1)

    save_csv({
        'dG_deviation': dfx_D['dG_deviation'].values,
        'D_norm': dfx_D['D_norm'].values,
        'demix_group': np.select([mix_mask.values, demix_mask.values], ['Mix', 'Demix'], default='Excluded'),
        'seq_name1': dfx_D['seq_name1'].values,
        'seq_name2': dfx_D['seq_name2'].values
    }, 'panelD_dG_diff.csv')


def panel_E(ax):
    """Panel E: Heterotypic dGij vs average homotypic dG for mixing pairs (D_norm < 0.13, dGij < 0)."""
    df_E = DF[DF['D_norm'] < 0.13].copy()
    #df_E = df_E[df_E['comp1'] == 50]
    df_E['average_dG'] = (df_E['dG_p1'] + df_E['dG_p2']) / 2
    df_E = df_E[df_E['dGij'] < 0]

    sns.histplot(data=df_E, x='average_dG', y='dGij', binwidth=(0.2, 0.2), cmap=cmaps.blues_dark_r, cbar=True, ax=ax,
                 cbar_kws={'shrink': 0.5, 'format': '%d'})
    cbar_ax = ax.figure.axes[-1]  # histplot appends the colorbar as the last axes on the figure
    cbar_ax.yaxis.set_label_position('left')
    cbar_ax.set_ylabel('Count', rotation=90, labelpad=8)

    slope, intercept, r_val, _, _ = linregress(df_E['average_dG'], df_E['dGij'])
    ax.plot([-10, 0], [slope * -10 + intercept, intercept], 'k--', lw=3)
    ax.plot([-10, 0], [-10, 0], color='grey', ls='--', lw=2)
    ax.text(0.95, 0.05, f'Slope = {slope:.2f}\n$R^2$ = {r_val**2:.2f}', transform=ax.transAxes, fontsize=22, va='bottom', ha='right',
            bbox=dict(boxstyle='round,pad=0.5', fc='white', alpha=0.7))
    ax.set_xlabel('(ΔG$_{ii}$ + ΔG$_{jj}$)/2 ($k_BT$)'); ax.set_ylabel('Heterotypic ΔG$_{ij}$ ($k_BT$)')
    ax.set_xlim(0.5, -10); ax.set_ylim(0.5, -10); ax.set_aspect('equal'); ax.grid(False)
    save_csv({'average_dG': df_E['average_dG'].values, 'dGij': df_E['dGij'].values, 'slope': slope, 'r_squared': r_val**2}, 'panelE_dG_avg.csv')

def panel_F(ax):
    """Panel F: Heterotypic dGij vs total net charge in the simulation box, for repeated pairs (D_norm < 0.13)."""
    df_F = DF[DF['D_norm'] < 0.13].copy()
    df_F = df_F.dropna(subset=['fasta_p1', 'fasta_p2', 'n1', 'n2', 'dGij'])
    df_F['net_charge1'] = df_F['fasta_p1'].apply(calc_net_charge) * df_F['n1']
    df_F['net_charge2'] = df_F['fasta_p2'].apply(calc_net_charge) * df_F['n2']
    df_F['pair_id'] = df_F['seq_name1'] + '_' + df_F['seq_name2']
    df_F['abs_total_net_charge'] = (df_F['net_charge1'] + df_F['net_charge2']).abs()
    df_F['total_charged'] = (df_F['fasta_p1'].apply(calc_total_charged) * df_F['n1'] +
                             df_F['fasta_p2'].apply(calc_total_charged) * df_F['n2'])
    pair_counts = df_F['pair_id'].value_counts()
    df_F = df_F[df_F['pair_id'].isin(pair_counts[pair_counts > 1].index)]

    # every tested (two-phase) system, as a faint background
    tested = DF.dropna(subset=['fasta_p1', 'fasta_p2', 'n1', 'n2', 'dGij'])
    tested_net = (tested['fasta_p1'].apply(calc_net_charge) * tested['n1'] +
                  tested['fasta_p2'].apply(calc_net_charge) * tested['n2']).abs()
    ax.scatter(tested_net, tested['dGij'], color='grey', s=60, alpha=0.2, lw=0, zorder=0)

    for pid in df_F['pair_id'].unique():
        pdata = df_F[df_F['pair_id'] == pid].sort_values('abs_total_net_charge')
        ax.plot(pdata['abs_total_net_charge'], pdata['dGij'], color='gray', alpha=0.3, lw=1.5, zorder=1)
    # only phase-separating pairs (dGij < -3) are coloured; the rest are grey like the tested systems
    ps = df_F['dGij'] < -3
    ax.scatter(df_F.loc[~ps, 'abs_total_net_charge'], df_F.loc[~ps, 'dGij'], color='grey', s=60,
               alpha=0.2, lw=0, zorder=0)
    scatter = ax.scatter(df_F.loc[ps, 'abs_total_net_charge'], df_F.loc[ps, 'dGij'],
                         c=df_F.loc[ps, 'total_charged'], cmap='YlOrBr',
                         s=100, alpha=0.7, edgecolors='black', lw=1, zorder=2)
    ax.set_xlabel('|Net Charge in Simulation|'); ax.set_ylabel('Heterotypic ΔG$_{ij}$ ($k_BT$)')
    ax.set_ylim(0.5, -10); ax.set_box_aspect(1); ax.grid(False)
    cbar = plt.colorbar(scatter, ax=ax, pad=0.02, shrink=0.5)
    cbar.ax.yaxis.set_label_position('left')
    cbar.set_label('Num. of Charged Residues', rotation=90, labelpad=8)
    save_csv({'abs_net_charge': df_F['abs_total_net_charge'].values, 'dGij': df_F['dGij'].values,
              'total_charged': df_F['total_charged'].values, 'pair_id': df_F['pair_id'].values}, 'panelF_charge.csv')


def create_figure3():
    if DF is None: init_data()

    fig = plt.figure(figsize=(33, 22))
    gs = gridspec.GridSpec(2, 3, figure=fig, hspace=-0.1, wspace=0.15, height_ratios=[0.6, 1])
    add_label = lambda ax, lbl: ax.text(-0.1, 1.12, lbl, transform=ax.transAxes, fontsize=28, fontweight='bold', va='top', ha='right')

    ax = fig.add_subplot(gs[0, 0]); add_label(ax, 'A'); panel_A(ax)
    ax = fig.add_subplot(gs[0, 1]); add_label(ax, 'B'); panel_B(ax)
    ax = fig.add_subplot(gs[0, 2]); add_label(ax, 'C'); panel_C(ax)
    ax = fig.add_subplot(gs[1, 0]); add_label(ax, 'D'); panel_D(ax)
    ax = fig.add_subplot(gs[1, 1]); add_label(ax, 'E'); panel_E(ax)
    ax = fig.add_subplot(gs[1, 2]); add_label(ax, 'F'); panel_F(ax)

    plt.tight_layout()
    plt.savefig('Figure3.pdf', bbox_inches='tight')
    print("Figure 3 saved as Figure3.pdf")

    df = load_data()[0]
    df = df[df['comp1'] == 50].copy()
    group1_mask = []
    group2_mask = []

    for _, row in df.iterrows():
        dG_p1, dG_p2 = row['dG_p1'], row['dG_p2']

        demix_idx = row['D_norm']
        dGij = row['dGij']

        if demix_idx < 0.13 and dGij < -3:
            group1_mask.append(True)
            group2_mask.append(False)
        elif demix_idx >= 0.13 and (dG_p1 < -3 or dG_p2 < -3):
            group1_mask.append(False)
            group2_mask.append(True)
        else:
            group1_mask.append(False)
            group2_mask.append(False)

    group1_mask = np.array(group1_mask)
    group2_mask = np.array(group2_mask)
    df['ocn'] = df.apply(lambda r: calc_ocn(r['fasta_p1'], r['fasta_p2']), axis=1)
    plt.figure()
    plt.scatter(df[group1_mask]['ocn'], df[group1_mask]['dGij'], color='blue', label='Mixing', alpha=0.5)
    plt.savefig('test.pdf')


if __name__ == '__main__':
    init_data()
    print_pp_samples()
    analyze_gin_demixing()
    create_figure3()