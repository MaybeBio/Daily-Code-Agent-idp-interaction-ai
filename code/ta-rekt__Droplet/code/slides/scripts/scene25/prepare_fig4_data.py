"""
Reduce the solver output behind Fig. 4 D-G of the droplet paper to what the
slides plot: for each run, the far-field volume fractions, the surface volume
fraction on the stable branch (and on the metastable one, where bistable), the
branch, and whether the solver converged.

Reads, from the repo root's data/ (gitignored; regenerate it with
figure_production/fig4/run_fig4_solver.sh):
    data/20261004/{small_R,medium_R,large_R}/demo_bvp_profiles.json   (D, F)
    data/20261007/run{1..5}/demo_bvp_profiles.json                    (E, G)
the same files, and the same fields, as figure_production/fig4/fig4.ipynb:
D reads the profiles from "phi_solver", E from "phi_profiles_full".

Writes slides/scripts/scene25/fig4_curves.json (small, kept with the scripts).

Usage (from slides/):
    python3 scripts/scene25/prepare_fig4_data.py
"""

import json
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[3]
OUT = Path(__file__).with_name("fig4_curves.json")

V_WATER = 0.0297         # nm^3: the lattice site
SITE = V_WATER ** (1 / 3)  # nm
RUNS = {
    "D": [(f"data/20261004/{name}/demo_bvp_profiles.json", "phi_solver")
          for name in ("small_R", "medium_R", "large_R")],
    "E": [(f"data/20261007/run{k}/demo_bvp_profiles.json", "phi_profiles_full")
          for k in range(1, 6)],
}


def reduce(path, profiles):
    d = json.load(open(REPO / path))
    return dict(
        source=path,
        R_nm=d["R"] * 0.3,                     # as labelled in fig4.ipynb
        h0=-d["params"]["surface"]["omega_b"] / SITE ** 2,  # k_BT / nm^2
        binodal=d["binodal"],
        phi_inf=d["phi_inf"],
        phi_surf=[profile[0] for profile in d[profiles]],
        phi_other=[p[0] if len(p) else None for p in d["phi_other"]],
        branch=d["branch"],
        success=d["success"],
    )


if __name__ == "__main__":
    curves = {panel: [reduce(path, profiles) for path, profiles in runs]
              for panel, runs in RUNS.items()}
    OUT.write_text(json.dumps(curves))
    for panel, runs in curves.items():
        for run in runs:
            print(f"{panel}: R = {run['R_nm']:g} nm, h0 = {run['h0']:.3f} kT/nm^2, "
                  f"{len(run['phi_inf'])} points, bistable: {'thick' in run['branch']}")
    print(f"wrote {OUT} ({OUT.stat().st_size / 1e3:.0f} kB)")
