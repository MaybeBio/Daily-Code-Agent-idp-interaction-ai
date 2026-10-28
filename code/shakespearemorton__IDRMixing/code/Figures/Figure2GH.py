#!/usr/bin/env python3
"""Figure 2 Panels G-H: G: KDE distribution of demixing index, H: Bar graph by composition"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.stats import gaussian_kde

plt.rcParams.update({
    'font.size': 20, 'axes.linewidth': 2, 'axes.spines.top': False, 'axes.spines.right': False,
    'xtick.major.size': 8, 'xtick.major.width': 2, 'ytick.major.size': 8, 'ytick.major.width': 2,
    'legend.frameon': False, 'figure.dpi': 300
})

BAR_COLORS = ['#cad2c5', '#84a98c', '#52796f', '#354f52', '#2f3e46']
BIN_ORDER = ['5-15', '15-35', '50', '65-85', '85-95']
DATA_DIR = Path("figure2gh_data")

def assign_bin(comp):
    """Assign composition to bins."""
    if 85 <= comp <= 100: return '85-95'
    if 65 <= comp < 85: return '65-85'
    if 35 <= comp < 65: return '50'
    if 15 <= comp < 35: return '15-35'
    if 0 <= comp < 15: return '5-15'
    return None

def load_data():
    """Load and merge datasets."""
    df_dem = pd.read_csv("../DATASETS/heterotypic_slabs.csv")
    df_prot = pd.read_csv("../DATASETS/gin_samples.csv")
    print(f'Total number of Samples: {len(df_dem)}')
    df = df_dem.copy()
    
    cols = ['seq_name', 'faro', 'ah_ij', 'N', 'shd', 'mean_lambda', 'SPR_svr']
    # Fixed column lookup names from 'seq_name_1' -> 'seq_name1'
    for p, sfx in [('seq_name1', '_p1'), ('seq_name2', '_p2')]:
        df = pd.merge(df, df_prot[cols], left_on=p, right_on='seq_name', how='left', suffixes=('', '_temp'))
        df = df.drop(columns=['seq_name']).rename(columns={c: c + sfx for c in cols[1:]})
    
    return df[(df['dGij'] < -3) & (df['comp1'] != 100) & (df['comp2'] != 100)].copy()

def create_panels_GH():
    """Create Panels G and H."""
    DATA_DIR.mkdir(exist_ok=True)
    df = load_data()
    print(f'Filtered Samples: {len(df)}')
    
    # Updated column names: composition1 -> comp1, demixing_index -> D_norm
    df['composition_bin'] = df['comp1'].apply(assign_bin)
    df_h = df[df['composition_bin'].notna()]
    comp_stats = df_h.groupby('composition_bin')['D_norm'].agg(['mean', 'std', 'count']).reindex(BIN_ORDER)
    comp_stats['sem'] = comp_stats['std'] / np.sqrt(comp_stats['count'])

    # Panel G data filtered for comp1 == 50
    df_g = df#[df['comp1'] == 50]
    demix_vals = df_g['D_norm'].dropna()
    
    # Save raw data
    demix_vals.to_csv(DATA_DIR / "panel_G_demixing_values.csv", index=False)
    comp_stats.to_csv(DATA_DIR / "panel_H_composition_stats.csv")
    df_h[['composition_bin', 'D_norm']].to_csv(DATA_DIR / "panel_H_raw_data.csv", index=False)
    
    # Create figure
    fig, (ax_g, ax_h) = plt.subplots(1, 2, figsize=(12, 4))
    add_label = lambda ax, lbl: ax.text(-0.1, 1.12, lbl, transform=ax.transAxes, fontsize=24, fontweight='bold', va='top', ha='right')
    
    # Panel G: KDE distribution plot
    add_label(ax_g, 'G')
    
    kde = gaussian_kde(demix_vals)
    x_range = np.linspace(demix_vals.min() - 0.05, demix_vals.max() + 0.05, 500)
    kde_values = kde(x_range)
    
    ax_g.fill_between(x_range, kde_values, alpha=0.5, color='lightgrey', edgecolor='black', linewidth=2)
    ax_g.plot(x_range, kde_values, color='black', linewidth=2)
    ax_g.set_yscale('log')
    ax_g.set_xlabel('Demixing Index', fontsize=18)
    ax_g.set_ylabel('Density', fontsize=18)
    ax_g.set_xlim(-0.05, 1.05)
    ax_g.grid(True, alpha=0.3, axis='x')
    
    stats_text = f'n = {len(demix_vals)}\nMean = {demix_vals.mean():.3f}\nMedian = {demix_vals.median():.3f}\nStd = {demix_vals.std():.3f}'
    ax_g.text(0.95, 0.95, stats_text, transform=ax_g.transAxes, fontsize=14, va='top', ha='right',
              bbox=dict(boxstyle='round,pad=0.5', fc='white', alpha=0.7))
    
    # Panel H: Bar graph
    add_label(ax_h, 'H')
    x_pos = np.arange(len(BIN_ORDER))
    for i, (pos, color) in enumerate(zip(x_pos, BAR_COLORS)):
        max_val = df_h[df_h["composition_bin"] == BIN_ORDER[i]]["D_norm"].max()
        print(f'Maximum demixing in bin {BIN_ORDER[i]}: {max_val}')
        ax_h.bar(pos, comp_stats['mean'].iloc[i], yerr=comp_stats['sem'].iloc[i],
                 capsize=8, alpha=0.8, color=color, edgecolor='black', linewidth=2)
        ax_h.text(pos, comp_stats['mean'].iloc[i] + comp_stats['sem'].iloc[i] + 0.003,
                  f'n={int(comp_stats["count"].iloc[i])}', ha='center', va='bottom', fontsize=14)
    
    ax_h.set_xlabel('IDR-1 Composition Range (%)', fontsize=18)
    ax_h.set_ylabel('Demixing Index', fontsize=18)
    ax_h.set_xticks(x_pos)
    ax_h.set_xticklabels(BIN_ORDER)
    ax_h.set_xlim(-0.6, len(BIN_ORDER) - 0.4)
    ax_h.set_ylim(0, 0.35)
    ax_h.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    plt.savefig('Figure2GH.pdf', bbox_inches='tight')
    print(f"Figure saved as Figure2GH.pdf, data saved in {DATA_DIR}/")

if __name__ == '__main__':
    create_panels_GH()