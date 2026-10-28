"""
Plot concentration profile at fixed distance from sphere surface as a function of phi_inf.

This script solves the BVP for spherical geometry at various phi_inf values
and plots phi(z=0.2R) vs phi_inf.

KEY: Uses continuation (previous solution as initial guess) to stay on the same
solution branch as phi_inf is varied.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_bvp
import json
import os
from viz_util.bvp_solver_membrane import (Parameters, solve_auto, make_ode, make_bc, 
                            compute_mu_inf, clip01)


def solve_with_previous(p: Parameters, prev_sol, tol: float = 1e-6):
    """
    Solve BVP using previous solution as initial guess.
    This keeps us on the same solution branch.
    """
    if prev_sol is None:
        return solve_auto(p, verbose=False)
    
    mu_inf = compute_mu_inf(p)
    ode = make_ode(p, mu_inf)
    bc = make_bc(p, mu_inf)
    
    # Interpolate previous solution onto new grid
    z_new = np.linspace(0, p.L, 500)
    
    # Scale z if domain size changed
    scale = prev_sol.x[-1] / p.L if prev_sol.x[-1] > 0 else 1.0
    z_scaled = np.minimum(z_new * scale, prev_sol.x[-1] * 0.999)
    y_new = prev_sol.sol(z_scaled)
    y_new[0] = clip01(y_new[0], p.eps)
    
    # Try with previous solution as initial guess
    for current_tol in [tol, tol*10, tol*100]:
        try:
            sol = solve_bvp(ode, bc, z_new, y_new, p=prev_sol.p.copy(),
                           tol=current_tol, max_nodes=50000, verbose=0)
            if sol.success:
                phi_min = np.min(sol.y[0])
                phi_max = np.max(sol.y[0])
                phim = sol.p[0]
                if phi_min > 0 and phi_max < 1 and 0 < phim < 1:
                    return sol
        except Exception:
            continue
    
    # Fallback to automatic solver
    return solve_auto(p, verbose=False)


def main():
    # Sphere radius
    R = 5.6
    
    # Domain size
    L = 15.0
    
    # Distance from surface (as fraction of R)
    z_frac = 0.2
    z_probe = 6 - R  # Actual distance
    
    # Range of phi_inf values to scan (more points for smoother curve)
    phi_inf_values = np.linspace(0.01, 0.244, 50)
    
    # Storage for results
    phi_at_probe = []
    phi_at_surface = []
    phi_m_values = []
    successful_phi_inf = []
    solutions = []  # Store full solutions for profile plotting
    
    print(f"Scanning phi_inf from {phi_inf_values[0]:.2f} to {phi_inf_values[-1]:.2f}")
    print(f"Sphere radius R = {R}")
    print(f"Probe distance z = {z_frac}R = {z_probe}")
    print("-" * 50)
    
    prev_sol = None  # Track previous solution for continuation
    
    date = "20260218"
    run = "_enhancer2_collapsed"
    
    for phi_inf in phi_inf_values:
        p = Parameters(
            chi_b=2.2,
            chi_m=1.0,
            omega_m=-0.1,
            omega_b=-0.7,
            chi_bb=-0.0,
            chi_bm=-0.0,
            phi_inf=phi_inf,
            kappa=1.2,
            L=L,
            R=R,
        )
        
        # Use continuation: previous solution as initial guess
        sol = solve_with_previous(p, prev_sol)
        
        if sol is not None and sol.success:
            # Update previous solution for next iteration
            prev_sol = sol
            
            # Evaluate solution at probe distance
            phi_probe, _ = sol.sol(z_probe)
            phi_surface, _ = sol.sol(0)
            
            phi_at_probe.append(phi_probe[0] if hasattr(phi_probe, '__len__') else phi_probe)
            phi_at_surface.append(phi_surface[0] if hasattr(phi_surface, '__len__') else phi_surface)
            phi_m_values.append(sol.p[0])
            successful_phi_inf.append(phi_inf)
            solutions.append(sol)  # Store solution for profile plotting
            print(f"  phi_inf={phi_inf:.3f}: phi(0)={phi_at_surface[-1]:.4f}, phi({z_frac}R)={phi_at_probe[-1]:.4f}, phi_m={phi_m_values[-1]:.4f}")
        else:
            print(f"  phi_inf={phi_inf:.3f}: FAILED - stopping sweep")
            # Stop the sweep at first failure
            break
    
    # Convert to arrays
    successful_phi_inf = np.array(successful_phi_inf)
    phi_at_probe = np.array(phi_at_probe)
    phi_at_surface = np.array(phi_at_surface)
    phi_m_values = np.array(phi_m_values)
    
    # Check if we have any solutions
    if len(solutions) == 0:
        print("\nNo solutions found!")
        return
    
    # Plotting
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    # Left panel: phi at various locations vs phi_inf
    ax1 = axes[0]
    # ax1.plot(successful_phi_inf, phi_at_surface, 'b-o', markersize=4, label=r'$\phi(z=0)$ (surface)')
    ax1.plot(successful_phi_inf, phi_at_probe, 'r-s', markersize=4, label=rf'$\phi(z={z_frac}R)$')
    # ax1.plot(successful_phi_inf, phi_m_values, 'g-^', markersize=4, label=r'$\phi_m$ (membrane)')
    ax1.plot(successful_phi_inf, successful_phi_inf, 'k--', alpha=0.5, label=r'$\phi_\infty$ (reference)')
    
    ax1.set_xlabel(r'$\phi_\infty$', fontsize=12)
    ax1.set_ylabel(r'$\phi$', fontsize=12)
    ax1.set_title(f'Concentration vs bulk value (R={R})', fontsize=12)
    ax1.legend(loc='best')
    ax1.set_xlim(phi_inf_values[0]*0.93, phi_inf_values[-1] *1.07)
    # ax1.set_ylim(-0.05, 1.05)
    ax1.grid(True, alpha=0.3)
    
    # Right panel: concentration profiles phi(z) for selected phi_inf values
    ax2 = axes[1]
    z_plot = np.linspace(0, L, 200)
    
    # Plot every 5th solution (so ~10 curves for 50 points)
    cmap = plt.cm.viridis
    indices_to_plot = range(0, len(solutions), 5)
    n_curves = len(list(indices_to_plot))
    
    for i, idx in enumerate(range(0, len(solutions), 5)):
        sol = solutions[idx]
        phi_inf = successful_phi_inf[idx]
        phi_profile, _ = sol.sol(z_plot)
        color = cmap(i / max(n_curves - 1, 1))
        ax2.plot(z_plot, phi_profile, '-', color=color, linewidth=1.5, 
                label=rf'$\phi_\infty={phi_inf:.2f}$')
        # Mark phi_inf with horizontal dashed line
        ax2.axhline(phi_inf, color=color, linestyle=':', alpha=0.3)
    
    ax2.set_xlabel(r'$z$ (distance from surface)', fontsize=12)
    ax2.set_ylabel(r'$\phi(z)$', fontsize=12)
    ax2.set_title('Concentration profiles', fontsize=12)
    ax2.axvline(x=z_probe, linestyle="dashed", color="black", alpha=0.5, lw=1)
    
    # First legend for solution curves
    leg1 = ax2.legend(loc='upper right', fontsize=8, ncol=2)
    
    # Add a second legend for model parameters (text only, no lines)
    params_text = (
        r"$\chi_b={}, \chi_m={}$ ".format(p.chi_b, p.chi_m) +
        r"$\omega_m={}, \omega_b={}$ ".format(p.omega_m, p.omega_b) +
        r"$\chi_bb={}, \chi_bm={}, \kappa={}$".format(p.chi_bb, p.chi_bm, p.kappa)
    )
    plt.suptitle(params_text, fontsize=14)
    
    # Re-add the first legend after text (since adding text removes it)
    ax2.add_artist(leg1)
    
    ax2.set_xlim(0, L)
    ax2.set_ylim(0, 1)
    ax2.grid(True, alpha=0.3)

    # Save results to JSON

    data = {
        "parameters": {
            "chi_b": p.chi_b,
            "chi_m": p.chi_m,
            "omega_m": p.omega_m,
            "omega_b": p.omega_b,
            "chi_bb": p.chi_bb,
            "chi_bm": p.chi_bm,
            "kappa": p.kappa,
            "L": L,
            "R": R,
            "z_frac": z_frac,
        },
        "phi_inf_values": successful_phi_inf.tolist(),
        "phi_at_probe": phi_at_probe.tolist(),
        "concentration_profiles": {
            str(successful_phi_inf[idx]): sol.sol(np.linspace(0, L, 200))[0].tolist()
            for idx, sol in enumerate(solutions)
        }
    }

    output_dir = "../data/{}".format(date)
    os.makedirs(output_dir, exist_ok=True)
    output_file = f'{output_dir}/spherical_phi_vs_phi_inf_run{run}.json'
    with open(output_file, 'w') as f:
        json.dump(data, f, indent=2)
    print(f"Saved data to {output_file}")

    plt.tight_layout()
    plt.savefig('../media/figures/{}/spherical_phi_vs_phi_inf_run{}.png'.format(date,run), dpi=150, bbox_inches='tight')
    print(f"\nSaved plot to spherical_phi_vs_phi_inf.png")
    plt.show()


if __name__ == "__main__":
    main()
