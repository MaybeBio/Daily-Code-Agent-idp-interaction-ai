"""
Demo: run the BVP solver and plot the bulk concentration profile φ(z)
alongside the equilibrium membrane composition φ_m.

Output: media/figures/20260508/bvp_concentration_profile.png
"""

import sys, os
import json
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from viz_util.bvp_solver import (
    Parameters, solve_auto, solve_bvp_robust, compute_spinodal, compute_binodal,
    dfb_dphi, d2fb_dphi2, compute_mu_inf, excess_free_energy
)

# ── Parameters ───────────────────────────────────────────────────────────────
p = Parameters(
    nb=1.0,
    chi_b=4,      # moderate Flory-Huggins bulk interaction (above spinodal threshold)
    omega_b=-1.3,   # surface attraction to bulk species
    kappa=1,
    phi_inf=0.010,    # far-field volume fraction (below spinodal)
    L=5.0,
    R=10.0,
    eps=1e-8,
    conserve_volume=False
)

# ── Thierry's formula ─────────────────────────────────────────────────────────

def phi_linear(z, p):

    alpha = d2fb_dphi2(p.phi_inf, p)

    h0 = p.omega_b
    lam = np.sqrt(alpha/p.kappa)

    return ( - h0 * p.R * np.exp(- lam * z) 
            / ((z + p.R) * p.kappa * (lam + 1 / p.R)))

# ── Solve ─────────────────────────────────────────────────────────────────────
print("Solving BVP …")
phi_infs = np.linspace(0.0095, 0.0105, 10)
z = np.linspace(0, p.L, 500)
phi_profiles = np.zeros((10, 500))
phi_profiles_lin = np.zeros((10, 500))

for k, phi_inf in enumerate(phi_infs):
    p.phi_inf = phi_inf
    sol = solve_bvp_robust(p, verbose=True)
    # sol = solve_auto(p, verbose=True)

    if sol is None or not sol.success:
        print("Solver failed – try adjust parameters.")
        sys.exit(1)

    if p.conserve_volume:
        phi_profiles[k], dphi, iphi = sol.sol(z)
    else:
        phi_profiles[k], dphi = sol.sol(z)

    phi_profiles_lin[k] = phi_profiles[k, -1] + phi_linear(z, p)

print(f"\nSolution summary:")
print(f"  φ(0)  = {phi_profiles[0][0]:.4f}   (surface bulk composition)")
print(f"  φ(L)  = {phi_profiles[0][-1]:.4f}  (should ≈ φ_inf = {p.phi_inf})")
print("mu_infinity = {}".format(-dfb_dphi(p.phi_inf, p)))

# ── Plot ──────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 1, figsize=(5, 3.5))

# Left: concentration profile
ax = axes
ax.plot(z, phi_profiles.T, color="#2563EB", lw=1)  #label=r"$\phi(z)$"
ax.plot(z, phi_profiles[0], color="#2563EB", lw=1, label=r"exact solution")
ax.plot(z, phi_profiles_lin.T, color="orange", lw=1, linestyle="dashed")  #label="analytical"
ax.plot(z, phi_profiles_lin[0], color="orange", lw=1, linestyle="dashed", label="linear approx.")

# ax.axhline(p.phi_inf, color="gray", ls="--", lw=1, label=r"$\phi_\infty$")
phi_spin_low, phi_spin_high = compute_spinodal(p)
if phi_spin_low is not None:
    ax.axhspan(phi_spin_low, phi_spin_high, alpha=0.10, color="red",
               label=f"spinodal ({phi_spin_low:.2f}, {phi_spin_high:.2f})")
phi_bin_low, phi_bin_high = compute_binodal(p)
if phi_bin_low is not None:
    ax.axhline(phi_bin_low, color="purple", ls="-.", lw=1.2,
               label=f"binodal ({phi_bin_low:.2f}, {phi_bin_high:.2f})")
    ax.axhline(phi_bin_high, color="purple", ls="-.", lw=1.2)
ax.axvline(0, color="k", lw=0.8, ls=":")
# ax.scatter(np.zeros(10), np.array([phi_profiles[:, 0]]), color="#2563EB", zorder=5, s=50)
ax.set_xlabel(r"Distance from surface $z$", fontsize=12)
ax.set_ylabel(r"Volume fraction $\phi$", fontsize=12)
ax.set_title("Bulk concentration profile", fontsize=13)
# ax.set_ylim(0.05, 0.2)
ax.set_xlim(0, p.L)
# ax.set_xlim(0, 1.5)
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)
fig.suptitle(rf"$R={p.R}$, $\chi_b={p.chi_b}$, $\omega_b={p.omega_b}$")

# ── Config ───────────────────────────────────────────────────────────────────
EXP_NAME = "fisher_two_spheres_test"
OUT_DATE = "20260722"
OUT_PNG  = f"media/figures/{OUT_DATE}/{EXP_NAME}.png"
OUT_JSON = f"media/figures/{OUT_DATE}/{EXP_NAME}.json"

fig.tight_layout()
os.makedirs(os.path.dirname(OUT_PNG), exist_ok=True)
fig.savefig(OUT_PNG, dpi=150, bbox_inches="tight")
print(f"\nSaved figure → {OUT_PNG}")

with open(OUT_JSON, "w") as f:
    json.dump(phi_profiles, f, indent=2)
print(f"wrote -> {OUT_JSON}")