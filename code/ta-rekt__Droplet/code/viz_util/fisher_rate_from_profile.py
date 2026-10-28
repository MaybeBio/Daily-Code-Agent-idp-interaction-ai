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

from viz_util.bvp_solver import (
    Parameters, solve_auto, d2fb_dphi2, dfb_dphi, fb, compute_spinodal,
    compute_binodal, psi_linear, F, compute_mu_inf
)
import matplotlib.pyplot as plt
from dataclasses import dataclass
import matplotlib
import numpy as np
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
matplotlib.use("Agg")


@dataclass
class Probe:
    kb: float = 1.0
    ku: float = 1.0


# ── Parameters ────────────────────────────────────────────────────────────
p = Parameters(
    nb=1.0,
    chi_b=2.5,
    omega_b=-0.6,
    kappa=1,
    phi_inf=0.1,
    L=20.0,
    eps=1e-6,
    R=2.5,
    conserve_volume=False,
    parallel_plates=False
)

# polymer params
s = 150   # inter-bead distance
a = 1   # Kuhn length

# ── Equilibrium distance distribution ─────────────────────────────────────
nsteps = 1000
r = np.linspace(0, p.L, nsteps)     # distance from center
d = np.linspace(0, p.L-p.R, nsteps)  # distance from surface


def compute_P_eq(p: Parameters, d: np.ndarray) -> np.ndarray:
    """Compute the equilibrium distribution of bead separation on a polymer"""

    # sol = solve_auto(p, verbose=False)

    # profile
    # phi = sol.sol(d)[0]
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


# ── Mean probe concentration varying total phi volume ─────────────────────
n_params = 100
phi_probe = np.zeros(n_params)
phi_tf_vals = np.linspace(0.003, 0.15, n_params)
dr = (p.L-p.R) / nsteps

for i, _phi in enumerate(phi_tf_vals):
    p.phi_inf = _phi
    P_eq, Phi_d, U_lin = compute_P_eq(p, d)
    phi_probe[i] = np.sum(P_eq * Phi_d)


# ── Fisher rate ───────────────────────────────────────────────────────────
def fisher_rate(p2: Probe,
                phi_probe: np.ndarray,
                phi_tf: np.ndarray) -> np.ndarray:

    # fisher rate without microcondensation, 2 state probe
    fisher_base = phi_tf * p2.kb * p2.ku / (p2.kb + p2.ku)

    # fisher rate accounting for probe-enhancer microcondensation
    dphi_probe = np.gradient(phi_probe, phi_tf)
    inv_taub = phi_probe * p2.kb * p2.ku / (p2.kb + p2.ku)
    fisher_probe = phi_tf**2 * inv_taub * (dphi_probe / phi_probe)**2

    return fisher_base, fisher_probe


# Biophysical parameters
a = 3e-3        # um
D_tf = 7.4      # um^2/s
D_cplx = 0.1    # um^2/s
dl_tf = 4*a*D_tf
ub = 5

probe = Probe(
    ku=ub,
    kb=dl_tf,
)
fisher_base, fisher_probe = fisher_rate(probe, phi_probe, phi_tf_vals)

# ── Fisher rate UB ────────────────────────────────────────────────────────

alpha = d2fb_dphi2(phi_tf_vals, p)
h0 = p.omega_b
lam = np.sqrt(alpha/p.kappa)

dalpha = (1.0 / p.nb) * (-1.0 / phi_tf_vals **
                         2 + 1.0 / (1.0 - phi_tf_vals) ** 2)
dlam = dalpha / (2.0 * np.sqrt(alpha * p.kappa))
phi_R = phi_tf_vals - h0 / (p.kappa * (lam + 1.0 / p.R))
dphi_R = 1.0 + dlam * (phi_tf_vals - phi_R) / (lam + 1.0 / p.R)
gain = dphi_R * phi_tf_vals / phi_R

fisher_surface = phi_R * probe.kb * probe.ku / (probe.kb + probe.ku)

# ── Plot ──────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(1, 1, figsize=(5, 3.5))

# phi_infs = np.linspace(0.001, 0.1, 15)

# for phi_inf in phi_infs:
#     p.phi_inf = phi_inf
#     P_eq, Phi_d, U_lin = compute_P_eq(p, d)
#     phi_probe_ = np.sum(P_eq * Phi_d)
#     mean_dist = np.sum(P_eq * d)
#     U_entropy = 3*(d+p.R)**2/(2*s*a**2)
#     U_geometry = -2*np.log(d+p.R)
#     U = U_lin+U_entropy+U_geometry

#     # axr = axl.twinx()
#     ax[0].plot(d, P_eq, color="#2563EB", lw=1, label=r"$\psi(d)$")
#     ax[0].axvline(mean_dist, color="#2563EB", lw=1, alpha=0.5)
#     # ax[1].plot(d, Phi_d-p.phi_inf, color="green", lw=1)
#     # ax[1].axhline(y=phi_probe_, color="green", alpha=0.5)

phi_a, phi_b = compute_binodal(p)
ax.axvline(x=phi_a, linestyle="dashed", color="black", lw=1, label="binodal")

# SPRT
eps = 0.32
K = np.log((1-eps)/eps)
c_over_dc = 10

T_dec_base = 2*K*np.tanh(K/2)*c_over_dc**2 / fisher_base
T_dec_surface = 2*K*np.tanh(K/2)*c_over_dc**2 / fisher_surface
T_dec_probe = 2*K*np.tanh(K/2)*c_over_dc**2 / fisher_probe

ax.plot(phi_tf_vals, T_dec_base, color="C0",
        linestyle="dotted", label="background")
ax.plot(phi_tf_vals, T_dec_surface, color="C0",
        linestyle="dashed", label="promoter only")
ax.plot(phi_tf_vals, T_dec_probe, color="C0",
        linestyle="solid", label="with enhancer")
ax.set_yscale('log')
ax.set_ylabel("Decision time (s)")
ax.set_xlabel(r"$\phi_{TF}$ (total volume fraction)", fontsize=12)
ax.legend()

# ax[1].plot(phi_tf_vals, fisher_base, color="red", lw=1, label="background",linestyle="dotted")
# ax[1].plot(phi_tf_vals, fisher_surface, color="red", linestyle="dashed", lw=1, label="promoter only")
# ax[1].plot(phi_tf_vals, fisher_probe, color="red", lw=1, label="with enhancer")

# axr = ax[0]
# axr.set_ylabel(r"$P_{eq}(d)$", color="#2563EB")
# axr.set_xlabel(r"distance $d$", color="black")
# axr.set_xlim(0,8)
# axr.set_yscale('log')

# axl = ax[1]
# axl.set_xlim(0,2)
# axl.set_ylabel(r"Phi profile", color="black")

# fig, axl = plt.subplots(1,1, figsize=(5, 3.5))
# axl = ax[1]
# axl.set_ylabel(r"Fisher information rate", color="black")
# axl.set_xlabel(r"$\phi_{TF}$ (total volume fraction)", fontsize=12)
# axl.set_yscale('log')
# axl.legend(fontsize=8)

# ── Save ──────────────────────────────────────────────────────────────────
fig.tight_layout()
out = "media/figures/20260525/polymer_sensitivity.png"
os.makedirs(os.path.dirname(out), exist_ok=True)
fig.savefig(out, dpi=150, bbox_inches="tight")
print(f"\nSaved figure → {out}")


###### potential #######

# Left: concentration profile
# axl = ax[0]
# axl.set_title("Linear")
# axl.set_ylabel(r"$V(d)$", color="black")
# axl.plot(d, U, color="black", lw=2, label=r"$U_{ads}$")
# # axl.plot(d, U_entropy, color="blue", lw=1, label=r"$U_{ent}$")
# # axl.plot(d, U_geometry, color="red", lw=1, label=r"$U_{geo}$")
# axl.set_xlabel("d (distance between surfaces)", fontsize=12)
# # axl.set_ylim(-1.10*np.max(np.abs(U)), 1.10*np.max(np.abs(U)))


###### gradient #######

# axl.set_ylabel(r"$\phi_{probe}$", color="black")
# axl.plot(phi_tf_vals, phi_probe, color="C2", lw=2)
# axl.plot(phi_tf_vals, phi_tf_vals, color="C2", lw=1.5, linestyle="dashed", label=r"$\phi_{TF}$")
# axl.set_xlabel(r"$\phi_{TF}$ (total volume fraction)", fontsize=12)
# axl.legend()

# axl.set_ylabel(r"$d\log \phi / d\log \phi_{\infty}$", color="black")
# axl.plot(np.log(phi_tf_vals), np.gradient(np.log(phi_probe))/np.gradient(np.log(phi_tf_vals)),
#          color="C2", lw=2,
#          label=r"$d\log \phi_{probe}$")
# axl.plot(np.log(phi_tf_vals),
#         np.gradient(np.log(phi_tf_vals))/np.gradient(np.log(phi_tf_vals)), color="C2", lw=1.5,
#          linestyle="dashed", label=r"$d\log \phi_{\infty}$")
# axl.plot(np.log(phi_tf_vals),
#         np.gradient(np.log(phi_R))/np.gradient(np.log(phi_tf_vals)), color="C2", lw=1.5,
#          linestyle="dotted", label=r"$d\log \phi_{R}$")
# axl.set_xlabel(r"$\log \phi_{\infty}$ (total volume fraction)", fontsize=12)
# axl.legend()

# ── Save ──────────────────────────────────────────────────────────────────
# fig.tight_layout()
# out = "media/figures/20260522/polymer_sensitivity.png"
# os.makedirs(os.path.dirname(out), exist_ok=True)
# fig.savefig(out, dpi=150, bbox_inches="tight")
# print(f"\nSaved figure → {out}")
