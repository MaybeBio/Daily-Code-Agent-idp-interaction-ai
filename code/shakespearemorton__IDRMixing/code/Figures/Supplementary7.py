#!/usr/bin/env python3
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

plt.rcParams.update({
    'font.size': 16, 'axes.linewidth': 2, 'axes.spines.top': False, 'axes.spines.right': False,
    'xtick.major.size': 6, 'xtick.major.width': 2, 'ytick.major.size': 6, 'ytick.major.width': 2,
    'legend.frameon': False, 'figure.dpi': 300, 'figure.facecolor': 'white', 'axes.facecolor': 'white'
})

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

def load_and_prepare_data():
    """Load datasets and reshape from pair-wise to per-protein format."""
    df_dem = apply_two_phase_gate(pd.read_csv("../DDM_Heterotypic/heterotypic_main.csv"))
    df_dem = df_dem[df_dem['two_phase']]  # gated-out (single-phase) maps are not analysed
    df_prot = pd.read_csv("../DATASETS/gin_samples.csv")

    # Merge GIN group info for protein 1 and protein 2
    df = pd.merge(df_dem, df_prot[['seq_name', 'GIN_group']], 
                  left_on='seq_name1', right_on='seq_name', how='left')
    df = df.rename(columns={'GIN_group': 'GIN_group_p1'}).drop(columns=['seq_name'])
    
    df = pd.merge(df, df_prot[['seq_name', 'GIN_group']], 
                  left_on='seq_name2', right_on='seq_name', how='left')
    df = df.rename(columns={'GIN_group': 'GIN_group_p2'}).drop(columns=['seq_name'])

    # Filter for simulations that underwent phase separation (using dGij < -3 threshold)
    # and ensure Demixing Index (D_norm) is recorded
    df_phase_sep = df[(df['dGij'] < -3) & (df['D_norm'].notna())].copy()

    # Melt the pair data so each protein gets a row for every phase-separated simulation it is in
    df1 = df_phase_sep[['seq_name1', 'GIN_group_p1', 'D_norm']].rename(
        columns={'seq_name1': 'seq_name', 'GIN_group_p1': 'GIN_group'})
    df2 = df_phase_sep[['seq_name2', 'GIN_group_p2', 'D_norm']].rename(
        columns={'seq_name2': 'seq_name', 'GIN_group_p2': 'GIN_group'})

    df_long = pd.concat([df1, df2], ignore_index=True)
    
    # Create the composite x-axis label combining Sequence Name and GIN Group
    df_long['label'] = df_long.apply(
        lambda x: f"{x['seq_name']}\n(GIN {int(x['GIN_group']) if pd.notna(x['GIN_group']) else '?'})", 
        axis=1
    )
    
    return df_long

def plot_demixing_by_protein():
    """Generate a clean, readable boxplot of demixing indices per unique protein."""
    df = load_and_prepare_data()

    # Sort the x-axis order by median demixing index for a cleaner overall shape
    order = df.groupby('label')['D_norm'].median().sort_values().index

    # Use a wide figure size to comfortably fit ~32 unique proteins
    fig, ax = plt.subplots(figsize=(22, 8))

    # Box plot displaying quartiles, medians, and outliers
    sns.boxplot(
        data=df, 
        x='label', 
        y='D_norm', 
        order=order, 
        ax=ax,
        color='#2E86AB', # Reusing a blue from your INTERACTION_PALETTE
        width=0.5,
        linewidth=1.5,
        fliersize=6,
        flierprops={'marker': 'o', 'markerfacecolor': '#A23B72', 'markeredgecolor': 'none', 'alpha': 0.6}
    )

    ax.set_ylabel('Demixing Index', fontsize=18)
    ax.set_xlabel('Protein Sequence', fontsize=18, labelpad=15)
    ax.set_ylim(-0.05, 1.05)
    
    # Rotate x-axis labels to prevent text collision
    plt.xticks(rotation=90, fontsize=12)
    
    # Add a subtle grid to help track horizontal levels across a wide plot
    ax.grid(True, alpha=0.3, axis='y', linestyle='--')

    plt.tight_layout()
    plt.savefig('Demixing_by_Protein.pdf', bbox_inches='tight')
    print("Plot saved as Demixing_by_Protein.pdf")

if __name__ == '__main__':
    plot_demixing_by_protein()
