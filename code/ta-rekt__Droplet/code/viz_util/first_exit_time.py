"""
Polymer prewetting as a concentration sensor 

A. linear approximation (analytical)
    1. plotting the equilibrium distance distribution between poly beads
    2. integrating over it for \psi_probe
    3. computing the fisher rate from \psi_probe

B. nonlinear profile from solver
    4. get free energy profile from numerical solution of \psi(z), Derjaguin
    5. repeat 1-3
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import numpy as np
import matplotlib
from dataclasses import dataclass
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from viz_util.bvp_solver import (
    Parameters, solve_auto, d2fb_dphi2, dfb_dphi, fb, compute_spinodal,
    compute_binodal, psi_linear, F, compute_mu_inf
)

@dataclass
class Probe:
    kb: float = 1.0
    ku: float = 1.0

# ── Parameters ────────────────────────────────────────────────────────────
p = Parameters(
    nb=1.0,
    chi_b=2.5,     
    omega_b=-0.3,   
    kappa=1,
    phi_inf=0.01,   
    L=20.0,
    eps=1e-6,
    R=5,
    conserve_volume=False,
    parallel_plates = False
)

# polymer params
s = 150   # inter-bead distance
a = 1   # Kuhn length

# ── Equilibrium distance distribution ─────────────────────────────────────
nsteps = 10000
r = np.linspace(0, p.L, nsteps)     # distance from center
d = np.linspace(0, p.L-p.R, nsteps) # distance from surface

def compute_P_eq(p: Parameters, d: np.ndarray) -> np.ndarray:
    """Compute the equilibrium distribution of bead separation on a polymer"""
    
    # profile
    psi_lin = psi_linear(d, p)[0]

    # potentials
    U_lin = 4 * np.pi * p.R**2 * p.omega_b * psi_lin
    U_ent = 3*(d+p.R)**2/(2*s*a**2)

    # equilibrium distribution
    P_eq_unnorm = (d+p.R)**2 * np.exp(- U_ent - U_lin)
    P_eq = P_eq_unnorm / np.sum(P_eq_unnorm)

    # concentration at the probe as a function of distance
    phi_boundary = p.phi_inf
    Phi_d = phi_boundary + psi_lin[0] + psi_lin 

    return P_eq, Phi_d, U_lin


# ── Plot ──────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(1,1, figsize=(5, 3.5))

phi_infs = np.linspace(0.001, 0.1, 15)

P_eq, Phi_d, U_lin = compute_P_eq(p, d)
phi_probe_ = np.sum(P_eq * Phi_d)
mean_dist = np.sum(P_eq * d)
U_entropy = 3*(d+p.R)**2/(2*s*a**2)
U_geometry = -2*np.log(d+p.R)
U = U_lin+U_entropy+U_geometry

ax.set_ylabel(r"$U(d)$", color="black")
ax.set_xlabel(r"distance $d$", color="black")
# ax.set_xlim(0,8)
ax.plot(d, U)

####### computing the MFPTs out of the potential using Texier eq 4.74

local_extrema = np.where(np.diff(np.sign(np.diff(U))))[0]
x2, x1 = local_extrema[:2]
d2U = np.gradient(np.gradient(U))
delta1 = 1/np.sqrt(d2U[x1])
delta2 = 1/np.sqrt(-d2U[x2])

for l in local_extrema:
    ax.axvline(d[l], color="black", alpha=0.5, linestyle="dotted")

U_at_extrema = U[local_extrema]
U2, U1 = U_at_extrema[:2]

D = 1.0
mean_hit_time = 2*np.pi*delta1*delta2*np.exp((U2-U1)/D)
print(mean_hit_time)


# ── Save ──────────────────────────────────────────────────────────────────
fig.tight_layout()
out = "media/figures/20260522/first_exit_time.png"
os.makedirs(os.path.dirname(out), exist_ok=True)
fig.savefig(out, dpi=150, bbox_inches="tight")
print(f"\nSaved figure → {out}")



###### potential #######

# Left: concentration profile
# ax = ax[0]
# ax.set_title("Linear")
# ax.set_ylabel(r"$V(d)$", color="black")
# ax.plot(d, U, color="black", lw=2, label=r"$U_{ads}$")
# # ax.plot(d, U_entropy, color="blue", lw=1, label=r"$U_{ent}$")
# # ax.plot(d, U_geometry, color="red", lw=1, label=r"$U_{geo}$")
# ax.set_xlabel("d (distance between surfaces)", fontsize=12)
# # ax.set_ylim(-1.10*np.max(np.abs(U)), 1.10*np.max(np.abs(U)))



# ── Save ──────────────────────────────────────────────────────────────────
# fig.tight_layout()
# out = "media/figures/20260522/polymer_sensitivity.png"
# os.makedirs(os.path.dirname(out), exist_ok=True)
# fig.savefig(out, dpi=150, bbox_inches="tight")
# print(f"\nSaved figure → {out}")


