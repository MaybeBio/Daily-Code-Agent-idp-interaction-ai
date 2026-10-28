"""
Why the Rouse chain is ~2x slower than the 1D reduction: the coordinate the
reduction calls "diffusing with D_eff = 2 kBT/gamma" is sub-diffusive over the whole
range of first-passage times.

(a) MSD of the end-to-end distance from viz_util/rouse_end_to_end_msd.jl (a free chain
    at the production parameters) against the Markovian 2 D_eff t and against exact
    discrete-Rouse theory; (b) the apparent diffusivity D_app(t) = MSD/(2t) in units of
    D_eff, with the first-passage time scales of data/20260912 marked.

Exact theory: in the dimensionless integrator of src/rouse.jl the chain is an OU
process dx = -A x dtau + sqrt(dtau) xi, so C_stat = A^{-1}/2 and, for the end-to-end
vector u = e_N - e_1,
    <|dR(t)|^2> = 3 posScale^2 * u' A^{-1} (I - exp(-A t/tau_bond)) u,
which is 6 D_eff t as t -> 0.  The radial (scalar-distance) MSD is 1/3 of that while
the displacement is small compared to |R|.

    python viz_util/rouse_end_to_end_msd.py [msd_json] [out_dir]
"""

import json
import os
import sys

import matplotlib.pyplot as plt
import numpy as np

MSD_JSON = sys.argv[1] if len(sys.argv) > 1 else "media/figures/20260914/rouse_end_to_end_msd.json"
OUT_DIR  = sys.argv[2] if len(sys.argv) > 2 else "media/figures/20260914"

d = json.load(open(MSD_JSON))
t, msd = np.array(d["t"]), np.array(d["msd"])
Deff, tb, tr = d["D_eff"], d["tau_bond"], d["tau_rouse"]
s, a = int(d["s"]), d["a"]
N = s + 1

# ── exact discrete-Rouse end-to-end MSD ──────────────────────────────────────
A = (np.diag(np.full(N, 2.0)) + np.diag(np.full(N - 1, -1.0), 1)
     + np.diag(np.full(N - 1, -1.0), -1))
A[0, 0] = A[-1, -1] = 1.0
w, V = np.linalg.eigh(A)
u = np.zeros(N); u[0], u[-1] = 1.0, -1.0
c = V.T @ u
nz = w > 1e-12
posScale2 = 2.0 / (3.0 / a**2)                      # 2 kBT / kappa, kBT = 1


def msd_vec(tt):
    """<|dR(t)|^2> of the end-to-end VECTOR (3 components)."""
    x = np.outer(np.atleast_1d(tt) / tb, w[nz])
    return 3 * posScale2 * ((c[nz]**2 / w[nz]) * (1 - np.exp(-x))).sum(axis=1)


tt = np.logspace(np.log10(t[0] / 3), np.log10(max(t[-1], tr) * 3), 600)
theory_radial = msd_vec(tt) / 3.0                   # radial component

# first-passage scales of the sweep this was built for
MARKS = [(4.61e-5, "solver MFPT$_{\\rm exit}$ (low $c$)", "#1f77b4"),
         (1.15e-4, "Rouse MFPT$_{\\rm exit}$ (low $c$)",  "#cc3311"),
         (5.31e-4, "Rouse MFPT$_{\\rm entry}$",           "#2e7d5b"),
         (tr,      r"$\tau_{\rm Rouse}$",                 "0.35")]
D_STAR = 1 / 2.25                                   # observed exit slowdown, 1/2.25

fig, ax = plt.subplots(1, 2, figsize=(12.4, 4.5))

ax[0].loglog(t, msd, "o", ms=4.5, color="#cc3311", label="Rouse chain (simulated)")
ax[0].loglog(tt, theory_radial, "-", color="k", lw=1.3, label="exact discrete-Rouse theory")
ax[0].loglog(tt, 2 * Deff * tt, "--", color="#1f77b4", lw=1.6,
             label=r"Markovian $2 D_{\rm eff} t$  ($D_{\rm eff}=2k_BT/\gamma$)")
ax[0].set_xlabel("lag $t$  [s]", fontsize=12)
ax[0].set_ylabel(r"MSD of the end-to-end distance", fontsize=11)
ax[0].set_title("(a)  the reduced coordinate is sub-diffusive", fontsize=12, loc="left")
ax[0].set_ylim(msd.min() / 3, 3 * max(theory_radial.max(), (2 * Deff * tt).max() / 30))
ax[0].legend(fontsize=9, loc="upper left")
ax[0].grid(alpha=0.25, which="both")

ax[1].semilogx(t, msd / (2 * Deff * t), "o", ms=4.5, color="#cc3311", label="simulated")
ax[1].semilogx(tt, theory_radial / (2 * Deff * tt), "-", color="k", lw=1.3, label="exact theory")
ax[1].axhline(1.0, color="#1f77b4", ls="--", lw=1.6, label="what the 1D reduction assumes")
ax[1].axhline(D_STAR, color="#cc3311", ls=":", lw=1.4,
              label=f"$D^*/D_{{\\rm eff}}$ implied by the exit times ({D_STAR:.2f})")
for x0, lab, col in MARKS:
    ax[1].axvline(x0, color=col, lw=0.9, alpha=0.65)
    ax[1].text(x0, 0.975, " " + lab, rotation=90, fontsize=8, color=col,
               va="top", ha="center", transform=ax[1].get_xaxis_transform())
ax[1].set_xlabel("lag $t$  [s]", fontsize=12)
ax[1].set_ylabel(r"$D_{\rm app}(t)\,/\,D_{\rm eff}$", fontsize=12)
ax[1].set_title("(b)  apparent diffusivity vs the reduction's constant", fontsize=12, loc="left")
ax[1].set_ylim(0, 1.45)
ax[1].legend(fontsize=8.5, loc="lower left", framealpha=0.92)
ax[1].grid(alpha=0.25)

fig.tight_layout()
os.makedirs(OUT_DIR, exist_ok=True)
for ext in ("pdf", "png"):
    out = os.path.join(OUT_DIR, f"rouse_end_to_end_msd.{ext}")
    fig.savefig(out, format=ext, bbox_inches="tight", dpi=150 if ext == "png" else None)
    print("saved", out)

short = t < 5 * tb
print(f"\nD_app/D_eff at the shortest lag ({t[0]:.2e} s) = {msd[0]/(2*Deff*t[0]):.3f}"
      "   <- the integrator realises D_eff exactly")
print(f"simulated/theory over t < 5 tau_bond: "
      f"{np.mean(msd[short]/(theory_radial_at := msd_vec(t[short])/3)):.3f}")
for x0, lab, _ in MARKS:
    print(f"  D_app/D_eff at {lab:34s} (t={x0:.2e} s) = {msd_vec(x0)[0]/(6*Deff*x0):.3f}")
i = np.argmin(abs(msd_vec(tt) / (6 * Deff * tt) - D_STAR))
print(f"\n  D_app/D_eff = {D_STAR:.3f} is reached at t = {tt[i]:.2e} s = {tt[i]/tb:.1f} tau_bond, "
      f"where the radial rms displacement is {np.sqrt(msd_vec(tt[i])[0]/3):.0f} "
      f"(well width ~12, exit distance ~95)")
