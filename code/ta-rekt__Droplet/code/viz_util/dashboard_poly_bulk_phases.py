"""
Dashboard for poly_bulk_phases experiment.

Visualizes polymer collapse phase transition induced by bulk-polymer interactions.
"""

import json
import numpy as np
import matplotlib.pyplot as plt


def dashboard_poly_bulk_phases(data_dir, param_range):
    """
    Create visualization dashboard for poly_bulk_phases experiment.
    
    Args:
        data_dir: Path to data directory (e.g., "data/20260124/poly_bulk_phases_16_n1-1")
        param_range: List of polymer self-interaction values (e.g., np.arange(-1.0, 1.25, 0.25))
    
    Returns:
        fig: Matplotlib figure with 4 subplots
    """
    # Load data from JSON files
    ncpus = 9
    results = {}
    for cpu in range(1, ncpus + 1):
        filepath = f"{data_dir}/output_cpu{cpu}.json"
        with open(filepath, 'r') as f:
            data = json.load(f)
            results.update(data)
    
    # Extract raw metrics
    sorted_keys = sorted(results.keys(), key=int)
    
    rg_series = []
    hams_poly = []
    hams_bulk = []
    acceptance_rates = []
    
    for key in sorted_keys:
        rg_series.append(results[key]['1']['radius_of_gyration'])
        hams_poly.append(results[key]['1']['hamiltonian_poly'])
        hams_bulk.append(results[key]['1']['hamiltonian_bulk'])
        acceptance_rates.append(results[key]['1']['acceptance_rate'])
    
    # Compute derived metrics (post burn-in)
    burn_in_idx = len(rg_series[0]) // 2  # 50% burn-in
    
    mean_rg = [np.mean(rg[burn_in_idx:]) for rg in rg_series]
    std_rg = [np.std(rg[burn_in_idx:]) for rg in rg_series]
    
    # Create figure with 4 subplots
    fig = plt.figure(figsize=(14, 10))
    
    # Plot 1: Mean Rg vs self-interaction
    ax1 = plt.subplot(2, 2, 1)
    ax1.plot(param_range, mean_rg, 'o-', linewidth=2, markersize=8, color='navy')
    ax1.axhline(y=mean_rg[0], color='red', linestyle='--', alpha=0.5, label='Collapsed state')
    ax1.axvline(x=0, color='gray', linestyle='--', alpha=0.5, label='ε=0')
    ax1.set_xlabel('Polymer Self-Interaction (ε_poly)', fontsize=12)
    ax1.set_ylabel('⟨Radius of Gyration⟩', fontsize=12)
    ax1.set_title('Polymer Size vs Self-Interaction', fontsize=14, fontweight='bold')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Rg fluctuations (phase transition indicator)
    ax2 = plt.subplot(2, 2, 2)
    ax2.plot(param_range, std_rg, 'o-', linewidth=2, markersize=8, color='darkgreen')
    ax2.axvline(x=0, color='gray', linestyle='--', alpha=0.5)
    ax2.set_xlabel('Polymer Self-Interaction (ε_poly)', fontsize=12)
    ax2.set_ylabel('σ(Rg)', fontsize=12)
    ax2.set_title('Rg Fluctuations (Phase Transition Indicator)', fontsize=14, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    
    # Plot 3: Hamiltonian equilibration check
    ax3 = plt.subplot(2, 2, 3)
    # Show key parameter values
    key_indices = [0, 4, 6, 8]  # ε = -1.0, 0.0, 0.5, 1.0
    for idx in key_indices:
        eps_val = param_range[idx]
        ax3.plot(hams_poly[idx], label=f'Poly (ε={eps_val:.2f})', alpha=0.7)
    ax3.set_xlabel('Frame', fontsize=12)
    ax3.set_ylabel('Polymer Hamiltonian', fontsize=12)
    ax3.set_title('Equilibration Check (Polymer)', fontsize=14, fontweight='bold')
    ax3.legend(fontsize=9)
    ax3.grid(True, alpha=0.3)
    
    # Plot 4: Bulk hamiltonian equilibration
    ax4 = plt.subplot(2, 2, 4)
    for idx in key_indices:
        eps_val = param_range[idx]
        ax4.plot(hams_bulk[idx], label=f'Bulk (ε_poly={eps_val:.2f})', alpha=0.7)
    ax4.set_xlabel('Frame', fontsize=12)
    ax4.set_ylabel('Bulk Hamiltonian', fontsize=12)
    ax4.set_title('Equilibration Check (Bulk)', fontsize=14, fontweight='bold')
    ax4.legend(fontsize=9)
    ax4.grid(True, alpha=0.3)
    
    plt.tight_layout()
    return fig
