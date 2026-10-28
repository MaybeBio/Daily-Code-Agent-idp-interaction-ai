"""
Demo: run the BVP solver and plot the bulk concentration profile φ(z)
alongside the equilibrium membrane composition φ_m.

Output: media/figures/20260508/bvp_concentration_profile.png
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from viz_util.bvp_solver_membrane import Parameters, solve_auto, compute_spinodal, d2fb_dphi2

# ── Parameters ───────────────────────────────────────────────────────────────
p = Parameters(
    nb=1.0,
    nm=1.0,
    chi_b=2.5,      # moderate Flory-Huggins bulk interaction (above spinodal threshold)
    chi_m=0.0,      # membrane self-interaction
    omega_m=0.0,   # membrane self-affinity
    omega_b=-0.5,   # surface attraction to bulk species
    chi_bb=0.0,    # bulk-bulk coupling at surface
    chi_bm=0.0,    # bulk-membrane coupling
    kappa=1.0,
    phi_inf=0.1,    # far-field volume fraction (below spinodal)
    L=10.0,
    R=5
)

# ── Thierry's formula ─────────────────────────────────────────────────────────

def phi_profile(z, p):

    alpha = d2fb_dphi2(p.phi_inf, p)

    h0 = p.omega_b
    lam = np.sqrt(alpha/p.kappa)

    return (p.phi_inf - h0 * p.R * np.exp(- lam * z) 
            / ((z + p.R) * p.kappa * (lam + 1 / p.R)))

# ── Solve ─────────────────────────────────────────────────────────────────────
print("Solving BVP …")
sol = solve_auto(p, verbose=True)

if sol is None or not sol.success:
    print("Solver failed – try adjust parameters.")
    sys.exit(1)

phi_m = sol.p[0]
z = np.linspace(0, p.L, 500)
phi, dphi = sol.sol(z)

print(f"\nSolution summary:")
print(f"  φ(0)  = {phi[0]:.4f}   (surface bulk composition)")
print(f"  φ_m   = {phi_m:.4f}   (membrane composition)")
print(f"  φ(L)  = {phi[-1]:.4f}  (should ≈ φ_inf = {p.phi_inf})")

# ── Plot ──────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 1, figsize=(7, 4.5))

# Left: concentration profile
ax = axes
ax.plot(z, phi-p.phi_inf, color="#2563EB", lw=2, label=r"$\phi(z)$")

# ax.axhline(p.phi_inf, color="gray", ls="--", lw=1, label=r"$\phi_\infty$")
phi_spin_low, phi_spin_high = compute_spinodal(p)
if phi_spin_low is not None:
    ax.axhspan(phi_spin_low, phi_spin_high, alpha=0.10, color="red",
               label=f"spinodal ({phi_spin_low:.2f}, {phi_spin_high:.2f})")
ax.axvline(0, color="k", lw=0.8, ls=":")
ax.scatter([0], [phi[0]], color="#2563EB", zorder=5, s=50)
ax.scatter([0], [phi_m], color="#DC2626", zorder=5, s=70,
           marker="*", label=rf"$\phi_m = {phi_m:.3f}$")
ax.plot(z, phi_profile(z, p)-p.phi_inf, color="orange", linestyle="dashed", label="analytical")

ax.set_xlabel(r"Distance from surface $z$", fontsize=12)
ax.set_ylabel(r"Volume fraction $\phi$", fontsize=12)
ax.set_title("Bulk concentration profile", fontsize=13)
ax.set_xlim(0, 4)
ax.set_ylim(1e-9, None)

ax.set_yscale('log')

ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)

# Right: gradient profile
# ax2 = axes[1]
# ax2.plot(z, dphi, color="#059669", lw=2)
# ax2.axhline(0, color="gray", ls="--", lw=1)
# ax2.set_xlabel(r"Distance from surface $z$", fontsize=12)
# ax2.set_ylabel(r"$\partial_z\phi$", fontsize=12)
# ax2.set_title("Composition gradient", fontsize=13)
# ax2.set_xlim(0, p.L)
# ax2.grid(True, alpha=0.3)

# Parameter annotation
param_text = (
    rf"$\chi_b={p.chi_b}$, $\chi_m={p.chi_m}$, $\omega_m={p.omega_m}$"
    "\n"
    rf"$\omega_b={p.omega_b}$, $\chi_{{bb}}={p.chi_bb}$, $\chi_{{bm}}={p.chi_bm}$"
    "\n"
    rf"$\kappa={p.kappa}$, $\phi_\infty={p.phi_inf}$, $L={p.L}$"
)
fig.text(0.5, -0.04, param_text, ha="center", fontsize=9, color="0.4")

fig.tight_layout()
out = "media/figures/20260508/bvp_concentration_profile.png"
os.makedirs(os.path.dirname(out), exist_ok=True)
fig.savefig(out, dpi=150, bbox_inches="tight")
print(f"\nSaved figure → {out}")
