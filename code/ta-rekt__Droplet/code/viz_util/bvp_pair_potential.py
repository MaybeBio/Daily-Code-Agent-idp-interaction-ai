"""
Plot: compute the pair potential using the linear superposition approx.
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from viz_util.bvp_solver import (
    Parameters, solve_auto, d2fb_dphi2, dfb_dphi, fb, compute_spinodal,
    compute_binodal, psi_linear, F, compute_mu_inf
)

# ── Parameters ───────────────────────────────────────────────────────────────
p = Parameters(
    nb=1.0,
    chi_b=2.5,      # moderate Flory-Huggins bulk interaction (above spinodal threshold)
    omega_b=-0.1,   # surface attraction to bulk species
    kappa=1,
    phi_inf=0.15,    # far-field volume fraction (below spinodal)
    L=15.0,
    eps = 1e-6,
    R=5,
    conserve_volume=True
)

def sph_lap(z_arr, dphi_arr, r_arr):
    """Spherical Laplacian d²φ/dz² + (2/r)·dφ/dz, from dphi array and radii r_arr."""
    d2phi = np.gradient(dphi_arr, z_arr) if len(z_arr) >= 2 else np.zeros_like(dphi_arr)
    return d2phi + 2.0 * dphi_arr / r_arr


def sph_lap_lin_left(z_arr, p):
    """Analytical spherical Laplacian of psi_linear for left sphere (r = R+z)."""
    psi, dpsi, lam = psi_linear(z_arr, p)
    r = p.R + z_arr
    d2psi = dpsi * (-lam - 1.0/r) + psi / r**2
    return d2psi + 2.0 * dpsi / r


def sph_lap_lin_right(u_arr, p):
    """Analytical spherical Laplacian of psi_linear for right sphere (r2 = R+u).

    u_arr = d - z : distance from right sphere's surface.
    Returns Laplacian w.r.t. z: d²ψ/du² - 2·(dψ/du)/r₂, since dψ/dz = -dψ/du.
    """
    psi, dpsi_du, lam = psi_linear(u_arr, p)
    r2 = p.R + u_arr
    d2psi = dpsi_du * (-lam - 1.0/r2) + psi / r2**2
    return d2psi - 2.0 * dpsi_du / r2


# ── Solve ─────────────────────────────────────────────────────────────────────
print("Solving BVP …")
sol = solve_auto(p, verbose=True)

if sol is None or not sol.success:
    print("Solver failed – try adjust parameters.")
    sys.exit(1)

nsteps = 1000
z = np.linspace(0, p.L, nsteps)
distances = np.linspace(0, p.L, nsteps)

# full solution, left sphere
phi, dphi, iphi = sol.sol(z)
psi = phi - phi[-1]
dpsi = dphi
psi_lin, dpsi_lin, lam = psi_linear(z, p)

# Laplacians for the full left-sphere profiles
lap_psi     = sph_lap(z, dpsi, p.R + z)
lap_psi_lin = sph_lap_lin_left(z, p)

print("mu_inf: ", -dfb_dphi(p.phi_inf, p))
print("F (linear): ", F(psi_lin, lap_psi_lin, p))
print("surface term: ", 2 * np.pi * p.R**2 * p.omega_b * psi_lin[0])

# ── Compute potential ─────────────────────────────────────────────────────────

Fvals = np.zeros(nsteps)
Fvals_lin = np.zeros(nsteps)

Uvals = np.zeros(nsteps)
Uvals_lin = 4 * np.pi * p.R**2 * psi_lin * p.omega_b


# ── Plot ──────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(3,1, figsize=(5, 9))

# Left: concentration profile
ax = axes[0]
ax.plot(z, phi, color="#2563EB", lw=2, label=r"$\phi(z)$")
ax.plot(z, phi[-1]+psi_lin, color="orange", linestyle="dashed", label="linear approx.")
# ax.set_yscale('log')

ax.axhline(p.phi_inf, color="gray", ls="--", lw=1, label=r"$\phi_\infty$")
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
ax.scatter([0], np.array([phi[0]]), color="#2563EB", zorder=5, s=50)
ax.set_xlabel(r"Distance from surface $z$", fontsize=12)
ax.set_ylabel(r"Volume fraction $\phi$", fontsize=12)
ax.set_title("Bulk concentration profile", fontsize=13)
ax.set_xlim(0, p.L)
# ax.set_ylim(0.0, 0.2)
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)

ax2 = axes[1]

# minmax = np.where(np.diff(np.sign(np.diff(Fvals))))[0]

# dFvals = np.diff(Fvals) * nsteps / p.L
# d2Fvals = np.diff(dFvals) * nsteps / p.L

# ax.axvline(x=distances[minmax[0]], linestyle="dashed", color="black", lw=1)
# ax.axvline(x=distances[minmax[1]], linestyle="dashed", color="black", lw=1)

# ax2.plot(distances, Uvals)
ax2.plot(distances, Uvals_lin)
# ax2.set_ylim(-5, 5)
ax2.axhline(y=0, linestyle="dashed", color="black")

ax3 = axes[2]

z_full = np.linspace(0, p.L, nsteps)
ax3.plot(z_full, psi_lin, color="blue")

for i, d in enumerate([1,3]):
    nd_pts = int((d/p.L)*nsteps)
    z_plot = np.linspace(0, d, nd_pts)
    psi_lin2, _, _ = psi_linear(d - z_plot, p)

    ax3.plot(z_plot, psi_lin2, color="C{}".format(d))
    ax3.plot(z_plot, psi_lin[:nd_pts] + psi_lin2, color="C{}".format(d), linestyle="dashed")

# ax.set_yscale('log')
ax3.set_xlim(0, p.L)

fig.tight_layout()
out = "media/figures/20260512/bvp_pair_potential.png"
os.makedirs(os.path.dirname(out), exist_ok=True)
fig.savefig(out, dpi=150, bbox_inches="tight")
print(f"\nSaved figure → {out}")
