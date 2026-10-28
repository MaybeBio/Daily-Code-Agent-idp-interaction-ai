"""
Dashboard for prewet_enhancers experiment.

Visualizes TF condensation between parallel DNA rods as function of rod separation.
"""

import json
import numpy as np
import matplotlib.pyplot as plt


def dashboard_prewet_enhancers(data_dir, param_range):
    """
    Create visualization dashboard for prewet_enhancers experiment.
    
    Args:
        data_dir: Path to data directory (e.g., "data/20251218/prewet_enhancers_32_1-16")
        param_range: List of rod separation values (e.g., list(range(1, 17)))
    
    Returns:
        fig: Matplotlib figure with 4 subplots
    """
    # Load data from JSON files
    ncpus = 16
    results = {}
    for cpu in range(1, ncpus + 1):
        filepath = f"{data_dir}/output_cpu{cpu}.json"
        with open(filepath, 'r') as f:
            data = json.load(f)
            results.update(data)
    
    # Extract raw metrics
    sorted_keys = sorted(results.keys(), key=int)
    
    hamiltonians = []
    density_lattices = []
    acceptance_rates = []
    
    for key in sorted_keys:
        hamiltonians.append(results[key]['1']['hamiltonian'])
        density_lattices.append(np.array(results[key]['1']['density_lattice']))
        acceptance_rates.append(results[key]['1']['acceptance_rate'])
    
    # Compute derived metrics
    mean_hams = [np.mean(h) for h in hamiltonians]
    mean_force = np.gradient(mean_hams, param_range)
    
    # Create figure with 4 subplots
    fig = plt.figure(figsize=(14, 10))
    
    # Plot 1: Mean force vs distance
    ax1 = plt.subplot(2, 2, 1)
    ax1.plot(param_range, mean_force, 'o-', linewidth=2, markersize=6)
    ax1.set_xlabel('Rod Distance (lattice units)', fontsize=12)
    ax1.set_ylabel('Mean Force (dH/dr)', fontsize=12)
    ax1.set_title('Force between rods', fontsize=14, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: YZ cross-sections for selected separations
    ax2 = plt.subplot(2, 2, 2)
    # Show 4 representative separations
    selected_alphas = [0, 4, 9, 15]  # Indices for alpha = 1, 5, 10, 16
    x_slice = 16  # Middle of lattice
    
    for i, idx in enumerate(selected_alphas):
        density_yz = density_lattices[idx][x_slice, :, :]
        
        # Create simple visualization (can be enhanced)
        extent = [0, 32, 0, 32]
        alpha_val = param_range[idx]
        
        # Plot as small subpanel (simplified - could use inset axes)
        # For now, just show the last one
        if i == len(selected_alphas) - 1:
            im = ax2.imshow(density_yz.T, origin='lower', cmap='viridis', 
                           extent=extent, aspect='auto')
            ax2.set_xlabel('Y (lattice units)', fontsize=12)
            ax2.set_ylabel('Z (lattice units)', fontsize=12)
            ax2.set_title(f'TF Density YZ slice (α={alpha_val})', fontsize=14, fontweight='bold')
            plt.colorbar(im, ax=ax2, label='TF Density')
    
    # Plot 3: XZ cross-section at midplane
    ax3 = plt.subplot(2, 2, 3)
    # Show separation = 8 as example
    mid_idx = 7  # alpha = 8
    y_slice = 16  # Approximate midplane between rods
    density_xz = density_lattices[mid_idx][:, y_slice, :]
    
    im = ax3.imshow(density_xz.T, origin='lower', cmap='viridis', aspect='auto')
    ax3.set_xlabel('X (lattice units)', fontsize=12)
    ax3.set_ylabel('Z (lattice units)', fontsize=12)
    ax3.set_title(f'TF Density XZ slice (α={param_range[mid_idx]})', fontsize=14, fontweight='bold')
    plt.colorbar(im, ax=ax3, label='TF Density')
    
    # Plot 4: Hamiltonian equilibration check
    ax4 = plt.subplot(2, 2, 4)
    # Show first and last parameter values
    ax4.plot(hamiltonians[0], label=f'α={param_range[0]}', alpha=0.7)
    ax4.plot(hamiltonians[-1], label=f'α={param_range[-1]}', alpha=0.7)
    ax4.set_xlabel('Frame', fontsize=12)
    ax4.set_ylabel('Hamiltonian', fontsize=12)
    ax4.set_title('Equilibration Check', fontsize=14, fontweight='bold')
    ax4.legend()
    ax4.grid(True, alpha=0.3)
    
    plt.tight_layout()
    return fig
