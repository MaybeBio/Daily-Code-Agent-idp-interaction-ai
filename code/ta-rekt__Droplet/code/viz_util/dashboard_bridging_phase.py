"""Capillary-bridging phase diagrams in the (d, phi_inf) plane.

Reads the JSON produced by scratchpad/compute_bridging.jl (one entry per
(chi, omega_b) panel) and renders a 2x2 grid of phase diagrams.

For every (d, phi_inf) cell the Julia side recorded whether each branch of the
two-sphere BVP converged:
  * unbridged (dilute gap)   -> U_exists, with excess free energy F_U
  * bridged   (dense neck)   -> B_exists, with excess free energy F_B

The background colour is the *existence* category, which exposes the metastable
/ spinodal structure:
    none (neither converges) | unbridged only | bistable (both) | bridged only
The black line is the *thermodynamic* transition: within the bistable wedge it
is where F_B = F_U (equilibrium switches from dilute to bridged), i.e. the
capillary-condensation onset.  The dashed line marks that chi's dilute binodal.

Usage:
    python viz_util/dashboard_bridging_phase.py [bridging.json] [out.png]
or from viz.ipynb:
    from viz_util.dashboard_bridging_phase import plot_bridging_phase
    plot_bridging_phase("path/to/bridging.json")
"""

import json
import sys

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Patch

# category codes
NONE, UNBRIDGED, BISTABLE, BRIDGED = -1, 0, 1, 2
_CMAP = ListedColormap(["#e6e6e6", "#4a90d9", "#f0ad4e", "#c0392b"])
_NORM = BoundaryNorm([-1.5, -0.5, 0.5, 1.5, 2.5], _CMAP.N)
_LABELS = [("#e6e6e6", "none (no solution)"),
           ("#4a90d9", "unbridged (dilute)"),
           ("#f0ad4e", "bistable (both)"),
           ("#c0392b", "bridged (dense neck)")]


def _arr(x):
    """Nested list with None -> np.nan into a float ndarray."""
    return np.array([[np.nan if v is None else v for v in row] for row in x],
                    dtype=float)


def _edges(centers):
    """Cell edges for pcolormesh from monotonic centers."""
    c = np.asarray(centers, dtype=float)
    mid = 0.5 * (c[:-1] + c[1:])
    return np.concatenate([[c[0] - (mid[0] - c[0])], mid,
                           [c[-1] + (c[-1] - mid[-1])]])


def _panel(ax, panel):
    d = np.asarray(panel["d_grid"], dtype=float)          # (Nd,)
    phi = np.asarray(panel["phi_inf"], dtype=float)       # (Nphi,)
    U = _arr(panel["U_exists"]).astype(bool)              # (Nphi, Nd)
    B = _arr(panel["B_exists"]).astype(bool)
    FU = _arr(panel["F_U"])
    FB = _arr(panel["F_B"])

    cat = np.full(U.shape, NONE, dtype=float)
    cat[U & ~B] = UNBRIDGED
    cat[~U & B] = BRIDGED
    cat[U & B] = BISTABLE

    ax.pcolormesh(_edges(d), _edges(phi), cat, cmap=_CMAP, norm=_NORM,
                  shading="flat")

    # Equilibrium (min-F) bridged region -> black onset line.
    eqB = np.where(B & ~U, 1.0,
          np.where(U & ~B, 0.0,
          np.where(U & B, (FB < FU).astype(float), np.nan)))
    if np.any(np.isfinite(eqB)) and np.nanmin(eqB) < 0.5 < np.nanmax(eqB):
        ax.contour(d, phi, np.nan_to_num(eqB, nan=0.0), levels=[0.5],
                   colors="k", linewidths=2.0)

    binodal = panel["binodal"]
    ax.axhline(binodal, ls="--", lw=1.0, color="0.25")
    ax.text(d[-1], binodal, "  binodal", va="center", ha="left",
            fontsize=8, color="0.25")

    ax.set_title(rf"$\chi={panel['chi']:.1f}$, $\omega_b={panel['omega_b']:.1f}$"
                 rf"   ($\phi_b={binodal:.3f}$)", fontsize=10)
    ax.set_xlim(d[0], d[-1])
    ax.set_ylim(phi[0], phi[-1] * 1.02)


def plot_bridging_phase(json_path, save_path=None):
    with open(json_path) as f:
        data = json.load(f)
    panels = data["panels"]
    R = data["R"]

    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    for ax, panel in zip(axes.ravel(), panels):
        _panel(ax, panel)
    for ax in axes[-1, :]:
        ax.set_xlabel(rf"surface separation  $d$   ($R={R:.0f}$)")
    for ax in axes[:, 0]:
        ax.set_ylabel(r"background concentration  $\phi_\infty$")

    handles = [Patch(facecolor=c, edgecolor="0.4", label=l) for c, l in _LABELS]
    handles.append(plt.Line2D([0], [0], color="k", lw=2.0,
                              label=r"equilibrium onset ($F_B=F_U$)"))
    fig.legend(handles=handles, loc="lower center", ncol=5, frameon=False,
               fontsize=9, bbox_to_anchor=(0.5, -0.01))
    fig.suptitle("Capillary-bridging phase diagrams  (two-sphere BVP, "
                 "fisher_two_spheres settings)", fontsize=12)
    fig.tight_layout(rect=(0, 0.03, 1, 0.98))

    if save_path:
        fig.savefig(save_path, dpi=140, bbox_inches="tight")
        print("wrote", save_path)
    return fig


if __name__ == "__main__":
    jp = sys.argv[1] if len(sys.argv) > 1 else "bridging.json"
    sp = sys.argv[2] if len(sys.argv) > 2 else jp.rsplit(".", 1)[0] + ".png"
    plot_bridging_phase(jp, sp)
