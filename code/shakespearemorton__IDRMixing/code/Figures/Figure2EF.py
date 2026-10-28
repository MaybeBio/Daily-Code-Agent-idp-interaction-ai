#!/usr/bin/env python3
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

plt.rcParams.update({
    'font.size': 16, 'axes.linewidth': 1.5, 'axes.spines.top': False, 'axes.spines.right': False,
    'xtick.major.size': 6, 'xtick.major.width': 1.5, 'ytick.major.size': 6, 'ytick.major.width': 1.5,
    'legend.frameon': False, 'figure.dpi': 300, 'figure.facecolor': 'white', 'axes.facecolor': 'white'
})

COLORS = ['#2c0703', '#890620', '#b6465f', '#da9f93', '#ebd4cb',
          '#090F1A', '#2A3A59', '#5A708C', '#D9D4BF', '#E5E3D9']

# ---------------------------------------------------------------- paths
# Updated to match Supplementary5.py methodology
RESULTS_DIR = Path("../DDM_Heterotypic/results_het")
HET_CSV = Path("../DATASETS/heterotypic_slabs.csv")
DATA_DIR = Path("figure2ac_data")

def save_csv(data_dict, filename):
    """Save data to CSV, padding arrays to equal length."""
    DATA_DIR.mkdir(exist_ok=True)
    max_len = max(len(v) if hasattr(v, '__len__') and not isinstance(v, str) else 1 for v in data_dict.values())
    for k, v in data_dict.items():
        if not hasattr(v, '__len__') or isinstance(v, str):
            data_dict[k] = [v] * max_len
        elif len(v) < max_len:
            data_dict[k] = list(v) + [np.nan] * (max_len - len(v))
    pd.DataFrame(data_dict).to_csv(DATA_DIR / filename, index=False)

def analyze_single_pair(pair_id, ax, color_offset=0):
    """Analyze a single protein pair showing all compositions on one plot."""
    if not HET_CSV.exists():
        print(f"Dataset {HET_CSV} not found"); return
        
    df = pd.read_csv(HET_CSV)
    
    # Locate all compositions for this specific pair
    pair_rows = df[df['base_name'].str.startswith(pair_id)]
    if pair_rows.empty:
        print(f"Pair {pair_id} not found in {HET_CSV}"); return
    
    s1 = pair_rows.iloc[0]['seq_name1']
    s2 = pair_rows.iloc[0]['seq_name2']
    
    compositions = []
    
    for _, row in pair_rows.iterrows():
        base_name = row['base_name']
        
        # Extract comps from the base_name format {s1}_{len1}_{s2}_{len2}_{comp1}_{comp2}
        parts = base_name.split('_')
        try:
            c1, c2 = int(parts[-2]), int(parts[-1])
        except ValueError:
            continue
            
        if c1 == 0 or c2 == 0: 
            continue
            
        # Load the voxel maps from the updated results directory
        npz_path = RESULTS_DIR / f"{base_name}.npz"
        if not npz_path.exists():
            print(f"Missing {npz_path}"); continue
            
        with np.load(npz_path, allow_pickle=True) as p:
            tot = p['avg_tot'].astype(np.float64).ravel()
            phi1 = p['avg_s1'].astype(np.float64).ravel()
        
        phi2 = tot - phi1
        demix = row['D_norm']
        
        compositions.append({'comp1': c1, 'comp2': c2, 'phi1': phi1, 'phi2': phi2, 'demixing': demix})
    
    compositions = sorted(compositions, key=lambda x: x['comp1'])[:5]
    if not compositions:
        print(f"No valid compositions for {pair_id}"); return
    
    # Save raw data
    for comp in compositions:
        save_csv({'phi1': comp['phi1'], 'phi2': comp['phi2'], 'demixing': comp['demixing']},
                 f"{pair_id}_{comp['comp1']}_{comp['comp2']}_data.csv")
    
    # Plotting logic
    all_phi = np.concatenate([np.concatenate([c['phi1'], c['phi2']]) for c in compositions])
    maxi = all_phi.max() * 1.1
    
    for idx, comp in enumerate(compositions):
        label = f"{comp['comp1']}:{comp['comp2']} (D = {comp['demixing']:.3f})"
        ax.scatter(comp['phi1'], comp['phi2'], alpha=0.7, s=40, color=COLORS[idx + color_offset],
                   edgecolors='black', linewidth=0.2, rasterized=True, label=label)
    
    ax.set_xlim(0, maxi); ax.set_ylim(0, maxi)
    ax.set_aspect('equal', adjustable='box')
    ax.yaxis.set_major_locator(ax.xaxis.get_major_locator())
    
    # Dashed line across diagonal
    ax.plot([0, maxi], [maxi, 0], color='gray', ls='--', alpha=0.5, lw=2)
    
    for t in ax.get_xticks():
        if 0 <= t <= maxi:
            ax.plot([t, t], [0, max(0, maxi - t)], color='gray', ls='--', alpha=0.3, lw=0.8)
            ax.plot([0, max(0, maxi - t)], [t, t], color='gray', ls='--', alpha=0.3, lw=0.8)
    
    ax.set_xlabel(f'{s1} (Res. / Voxel)'); ax.set_ylabel(f'{s2} (Res. / Voxel)')
    ax.legend(loc='best', fontsize=15)
    print(f"Analyzed {len(compositions)} compositions for {pair_id}")

if __name__ == "__main__":
    pair_ids = ['O60229_1_34_Q96RR1_1_152', 'Q3LI60_1_44_Q92581_1_59']
    color_offsets = [5, 0]
    
    fig, axes = plt.subplots(2, 1, figsize=(8, 16))
    for pair_id, ax, offset in zip(pair_ids, axes, color_offsets):
        analyze_single_pair(pair_id, ax, color_offset=offset)
    
    plt.savefig('Figure2AC.pdf', facecolor='white', dpi=300)
    print(f"\nFigure saved as Figure2AC.pdf, data saved in {DATA_DIR}/")
