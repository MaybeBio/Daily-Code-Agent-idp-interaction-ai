import sys, os, glob, json
from datetime import date

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import numpy as np
from scipy.integrate import cumulative_trapezoid
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from viz_util.bvp_solver import (
    Parameters, solve_auto, d2fb_dphi2, dfb_dphi, fb, compute_spinodal,
    compute_binodal, psi_linear, F, compute_mu_inf, regularized_log
)

EXP_NAME = "prewetting_graphical_construction_above_pw"
SLICE_AT  = None                       # value for the third axis; None => median
# OUT_DATE  = date.today().strftime("%Y%m%d")
OUT_DATE = "20260521"
DATA_DIR  = f"data/{OUT_DATE}/{EXP_NAME}"
OUT_PATH  = f"media/figures/{OUT_DATE}/{EXP_NAME}.png"

# ── Model ───────────────────────────────────────────────────────────────────
p = Parameters(
    nb=1.0,
    chi_b=4.0,      # moderate Flory-Huggins bulk interaction (above spinodal threshold)
    omega_b=-1.25,   # surface attraction to bulk species
    kappa=1,
    phi_inf=0.009,    # far-field volume fraction (below spinodal)
    L=5.0,
    R=np.inf,
    eps=1e-8,
    conserve_volume=False
)

mu_inf = compute_mu_inf(p)
phi_low, phi_high = compute_binodal(p)

# ── Plot ─────────────────────────────────────────────────────────────────────

fig, axes = plt.subplots(2,1,figsize=(4.2,7))

eps = 1e-4
phi_0s = np.linspace(eps, 1-eps, 1000)
W = fb(phi_0s, p) - fb(p.phi_inf, p) - mu_inf*(phi_0s - p.phi_inf)
dJ = -p.omega_b - 1e-3 * regularized_log(phi_0s/(1-phi_0s))
gamma = cumulative_trapezoid(np.sqrt(2*p.kappa*W)-dJ, phi_0s, initial=0.0)

integrand = -np.sqrt(2*p.kappa*W) +dJ
sols = np.where(np.diff(np.sign(integrand)))[0]

ax = axes[0]

i=0
for s in sols:
    if i==0:
        ax.axvline(x=phi_0s[s+1], color="black", alpha=0.8, lw=1, 
        label=rf"$\phi_0'$", linestyle="dotted")
    else:
        ax.axvline(x=phi_0s[s+1], color="black", alpha=0.8, lw=1, linestyle="dotted")
    i+=1

ax.plot(phi_0s, dJ, linestyle="dashed", lw=1, color="black", 
           label=r"$-h_0$")
ax.plot(phi_0s, np.sqrt(2*p.kappa*W), color="black", label=r"$\sqrt{2\kappa W}$")
ax.plot(phi_0s, -np.sqrt(2*p.kappa*W), color="black")
ax.set_title("Surface free energy integrand")
ax.set_ylim(0, None)
ax.set_xlabel(r"$\phi|_0$")
ax.legend()

ax = axes[1]

i=0
for s in sols:
    if i==0:
        ax.axvline(x=phi_0s[s+1], color="black", alpha=0.8, lw=1, 
        label=rf"$\phi_0'$", linestyle="dotted")
    else:
        ax.axvline(x=phi_0s[s+1], color="black", alpha=0.8, lw=1, linestyle="dotted")
    i+=1

ax.plot(phi_0s, gamma, color="red")
ax.set_title("Surface free energy")
ax.set_ylabel(r"$\gamma(\phi)$")
ax.set_xlabel(r"$\phi|_0$")

fig.suptitle(
    f"R = {p.R}, " 
        + rf"$\chi = {p.chi_b}$, "
        + rf"Binodal at $\phi_\infty={round(phi_low, 3)}$"
        + "\n"+ rf"$\phi_\infty = {p.phi_inf}, \quad h_0={p.omega_b}$",
    fontsize=13,
)
fig.tight_layout()

os.makedirs(os.path.dirname(OUT_PATH), exist_ok=True)
fig.savefig(OUT_PATH, dpi=150, bbox_inches="tight")
print(f"saved -> {OUT_PATH}")
