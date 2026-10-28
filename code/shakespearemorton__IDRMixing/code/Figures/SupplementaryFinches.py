import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import linregress
from finches import CALVADOS_frontend

plt.rcParams.update({
    "axes.prop_cycle": (plt.cycler('color', ['k', 'r', 'b', 'g']) + 
                        plt.cycler('ls', ['-', '--', ':', '-.'])),
    "figure.figsize": (3.3, 2.5),
    "figure.dpi": 600,
    "font.size": 8,
    "axes.labelsize": 8,
    "legend.fontsize": 7,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "lines.linewidth": 1  # Standard for IEEE
})

# Data paths - adjust BULOW_CSV if bulow_df is named/located differently.
# bulow_df is the new version of df_training.csv, and must contain a 'seq_name'
# column plus the raw sequence ('fasta') so it can be used as a lookup table.
MIX_CSV = '../DDM_Heterotypic/heterotypic_main.csv'
BULOW_CSV = '../DATASETS/bulow_df.csv'
MIN_DENSE_FRAC = 0.05  # two phases only if >= 5% of voxels sit at or above c_den


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


def calculate_regression_line(x_data, y_data, x_range):
    """Calculates linregress and returns the predicted y values."""
    res = linregress(x_data, y_data)
    return x_range, res.intercept + res.slope * x_range, res.rvalue


def attach_sequences(mixdata, bulow_df):
    """heterotypic_main.csv does not carry raw sequences directly (seq1/seq2),
    only seq_name1/seq_name2. Look the raw fasta strings up in bulow_df keyed on
    seq_name - the same pattern Supplementary3.py uses to attach GIN_group/dG
    from gin_samples.csv onto the heterotypic table."""
    seq_lookup = bulow_df[['seq_name', 'fasta']].drop_duplicates('seq_name')

    for name_col, seq_col in [('seq_name1', 'seq1'), ('seq_name2', 'seq2')]:
        mixdata = pd.merge(mixdata, seq_lookup, left_on=name_col, right_on='seq_name', how='left')
        mixdata = mixdata.drop(columns=['seq_name']).rename(columns={'fasta': seq_col})

    missing = mixdata[mixdata['seq1'].isna() | mixdata['seq2'].isna()]
    if len(missing):
        print(f"Warning: dropping {len(missing)} rows with no sequence match in bulow_df "
              f"(unmatched names: {sorted(set(missing['seq_name1']) | set(missing['seq_name2']))[:10]}...)")

    return mixdata.dropna(subset=['seq1', 'seq2']).copy()


# 1. Data Loading & Filtering
cf = CALVADOS_frontend()

bulow_df = pd.read_csv(BULOW_CSV)

# Process Mix Data (heterotypic) - pull seq1/seq2 in from bulow_df by name
mixdata = apply_two_phase_gate(pd.read_csv(MIX_CSV))
mixdata = mixdata[mixdata['two_phase']]  # gated-out (single-phase) maps are not analysed
mixdata = mixdata[mixdata['comp1'] == 50].copy()
mixdata['dGij'] = mixdata['dGij'].fillna(0)
mixdata = attach_sequences(mixdata, bulow_df)
mixdata['epsilon'] = [cf.epsilon(s1, s2) for s1, s2 in zip(mixdata['seq1'], mixdata['seq2'])]

# Process PS Data (homotypic) - bulow_df already has fasta + dG directly, no merge needed
ps_data = bulow_df.copy()
ps_data['epsilon'] = ps_data['fasta'].apply(lambda x: cf.epsilon(x, x))

# 2. Plotting
fig, axes = plt.subplots(2, 2, figsize=(5.2, 5.0), constrained_layout=True)

# Top-left: Plot dG Data
x_dg = np.linspace(-10, 0, 100)
_, y_dg, r_dg = calculate_regression_line(ps_data['dG'], ps_data['epsilon'], x_dg)
axes[0, 0].scatter(ps_data['dG'], ps_data['epsilon'], s=25, alpha=0.6, edgecolors='none')
axes[0, 0].plot(x_dg, y_dg, c='#ca0020', linestyle='--', linewidth=3, label=f'$R^2$={r_dg**2:.2f}')
axes[0, 0].set_xlabel(r'Homotypic $\Delta$G ($k_BT$)')
axes[0, 0].set_ylabel('Finches $\epsilon_{ii}$')

# Top-right: Plot Demixing Data
x_mix = np.linspace(0, 1, 100)
_, y_mix, r_mix = calculate_regression_line(mixdata['D_norm'], mixdata['epsilon'], x_mix)
axes[0, 1].scatter(mixdata['D_norm'], mixdata['epsilon'], s=25, alpha=0.6)
axes[0, 1].plot(x_mix, y_mix, c='#ca0020', linestyle='--', linewidth=3, label=f'$R^2$={r_mix**2:.2f}')
axes[0, 1].set_xlabel('Demixing Index')
axes[0, 1].set_ylabel('Finches $\epsilon_{ij}$')

# Calculate chi for mixdata
mixdata['epsilon_ii'] = [cf.epsilon(s1, s1) for s1 in mixdata['seq1']]
mixdata['epsilon_jj'] = [cf.epsilon(s2, s2) for s2 in mixdata['seq2']]
mixdata['chi'] = 2 * mixdata['epsilon'] - mixdata['epsilon_ii'] - mixdata['epsilon_jj']

# Bottom-left: chi versus demixing
x_chi_demix = np.linspace(mixdata['D_norm'].min(), mixdata['D_norm'].max(), 100)
_, y_chi_demix, r_chi_demix = calculate_regression_line(mixdata['D_norm'], mixdata['chi'], x_chi_demix)
axes[1, 0].scatter(mixdata['D_norm'], mixdata['chi'], s=25, alpha=0.6)
axes[1, 0].plot(x_chi_demix, y_chi_demix, c='#ca0020', linestyle='--', linewidth=3, label=f'$R^2$={r_chi_demix**2:.2f}')
axes[1, 0].set_xlabel('Demixing Index')
axes[1, 0].set_ylabel(r'2$\epsilon_{ij}$ - $\epsilon_{ii}$ - $\epsilon_{jj}$')

# Bottom-right: chi versus dGij
x_chi_dg = np.linspace(mixdata['dGij'].min(), mixdata['dGij'].max(), 100)
_, y_chi_dg, r_chi_dg = calculate_regression_line(mixdata['dGij'], mixdata['chi'], x_chi_dg)
axes[1, 1].scatter(mixdata['dGij'], mixdata['chi'], s=25, alpha=0.6)
axes[1, 1].plot(x_chi_dg, y_chi_dg, c='#ca0020', linestyle='--', linewidth=3, label=f'$R^2$={r_chi_dg**2:.2f}')
axes[1, 1].set_xlabel(r'Heterotypic $\Delta$G$_{ij}$ ($k_BT$)')
axes[1, 1].set_ylabel(r'2$\epsilon_{ij}$ - $\epsilon_{ii}$ - $\epsilon_{jj}$')

for ax in axes.flat:
    ax.legend()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

# Add panel labels
for i, label in enumerate(['A', 'B', 'C', 'D']):
    row, col = divmod(i, 2)
    axes[row, col].text(-0.15, 1.05, label, transform=axes[row, col].transAxes,
                        fontsize=10, fontweight='bold', va='top', ha='right')

plt.savefig('SupplementaryFinches.pdf', dpi=800, bbox_inches='tight', transparent=True, format='pdf')