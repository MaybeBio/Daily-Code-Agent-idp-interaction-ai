#!/usr/bin/env python3
"""Reduce the inputs of Fig. 5 to small files in processed_data/.

    python figure_production/fig5/prepare_fig5_data.py

The figure uses the June toy model (spec 9f78af2:spec/rouse_mfpt_attractive.jl,
re-run as fig5/rouse_toy.jl): dimensionless units, contact at d = R = 2.5, box
edge L = 10. Everything but the histograms is analytic and computed here, with
the formulas of that spec and of viz-rouse.ipynb cell 23.

Reads (gitignored, under the repo root):
  data/*/fig5_rouse_toy/output_cpu*.json   Rouse distance traces (panel D), the
                                           most recent dated folder; from
                                           fig5/run_fig5_rouse.sh
Writes:
  processed_data/fig5_potentials.csv   d, U_ads_k{i}, U_tot_k{i} for all 16 phi_inf
  processed_data/fig5_panel_e.csv      d, U_tot at high c on [R, D_MAX_E]: panel E,
                                       like the original, runs past L
  processed_data/fig5_scalars.json     model parameters, phi_inf, chosen indices,
                                       ridges d*(c), d^o, d^x, run statistics
  processed_data/fig5_hist.csv         pooled distance histograms -- only when the
                                       Rouse output exists
"""
import argparse, glob, json, os, sys
import numpy as np

# the June model
NB, CHI_B, OMEGA_B, KAPPA, EPS, R, L = 1.0, 2.5, -0.5, 1.0, 1e-6, 2.5, 10.0
N, KUHN, KBT = 100, 1.0, 1.0
S = N - 1                                  # beadSep: tracked beads 1 and N
PHI = np.linspace(0.003, 0.15, 16)         # alphas of the June spec
K_LOW, K_HIGH = 0, 9                       # cell 37 (0-based)
K_UADS = (3, 9, 15)                        # panel A: distinct ranges; k = 0 is a mere tick
D_O_FACTOR = 1.2                           # d^o = 1.2 max_c d*(c) (old caption)
N_BINS = 100                               # cell 37: bins=100, range=(R, L)
BURN_IN = 0.0                              # cell 37 used none
D_MAX_E = 15.0                             # panel E's curve reaches past L, as drawn


def resolve_repo(start):
    d = os.path.abspath(start)
    while d != os.path.dirname(d):
        if os.path.isfile(os.path.join(d, "Project.toml")) and os.path.isdir(os.path.join(d, "spec")):
            return d
        d = os.path.dirname(d)
    sys.exit("repo root not found above " + start)


def potentials(phi, d):
    """U_ads, U_tot and the analytic force -dU_tot/dd (cell 23)."""
    p = Parameters(nb=NB, chi_b=CHI_B, omega_b=OMEGA_B, kappa=KAPPA, phi_inf=phi, eps=EPS, R=R, L=L)
    lam = np.sqrt(d2fb_dphi2(phi, p) / KAPPA)
    pref = 4 * np.pi * R**3 * OMEGA_B**2 / (KAPPA * (lam + 1 / R))
    u_ads = -pref * np.exp(-lam * (d - R)) / d
    u_tot = u_ads - 2 * np.log(d) + 3 * d**2 / (2 * S * KUHN**2)
    force = 2 / d - 3 * d / (S * KUHN**2) - pref * np.exp(-lam * (d - R)) * (lam / d + 1 / d**2)
    return u_ads, u_tot, force


def ridge(d, force):
    """First maximum of U_tot out from contact: F = -dU/dd turns from - to +.
    Linear interpolation between the bracketing samples, as in cell 23."""
    sc = np.where((force[:-1] < 0) & (force[1:] > 0))[0]
    assert sc.size, "no ridge found"
    i = sc[0]
    return float(d[i] - force[i] * (d[i + 1] - d[i]) / (force[i + 1] - force[i]))


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    repo = resolve_repo(here)
    sys.path.insert(0, repo)
    global Parameters, d2fb_dphi2
    from viz_util.bvp_solver import Parameters, d2fb_dphi2

    found = sorted(glob.glob(os.path.join(repo, "data", "*", "fig5_rouse_toy")))
    ap = argparse.ArgumentParser()
    ap.add_argument("--rouse", default=found[-1] if found else None,
                    help="folder of output_cpu*.json (default: the latest data/*/fig5_rouse_toy)")
    ap.add_argument("--out", default=os.path.join(os.path.dirname(here), "processed_data"))
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    # ---- potentials, on the spec's own grid ------------------------------------
    d = np.linspace(R, L, 500)
    U_ads, U_tot, ridges = [], [], []
    for phi in PHI:
        ua, ut, f = potentials(phi, d)
        U_ads.append(ua); U_tot.append(ut); ridges.append(ridge(d, f))
    with open(os.path.join(a.out, "fig5_potentials.csv"), "w") as fh:
        fh.write("d," + ",".join("U_ads_k%d" % k for k in range(len(PHI))) + ","
                 + ",".join("U_tot_k%d" % k for k in range(len(PHI))) + "\n")
        for i in range(d.size):
            fh.write("%.17g," % d[i] + ",".join("%.17g" % u[i] for u in U_ads + U_tot) + "\n")

    dE = np.linspace(R, D_MAX_E, 800)
    uE = potentials(PHI[K_HIGH], dE)[1]
    with open(os.path.join(a.out, "fig5_panel_e.csv"), "w") as fh:
        fh.write("d,U_tot\n")
        for x, u in zip(dE, uE):
            fh.write("%.17g,%.17g\n" % (x, u))

    scal = {
        "model": "June toy model, spec 9f78af2:spec/rouse_mfpt_attractive.jl",
        "nb": NB, "chi_b": CHI_B, "omega_b": OMEGA_B, "kappa": KAPPA, "eps": EPS,
        "R": R, "L": L, "N": N, "kuhn_length": KUHN, "kBT": KBT, "bead_separation": S,
        "phi_inf": PHI.tolist(), "K_LOW": K_LOW, "K_HIGH": K_HIGH, "K_UADS": list(K_UADS),
        "x_ridge": ridges,
        # landmarks of panel E, code coordinates (contact at d = R; the figure
        # labels it 2R, as the original did)
        "d_x": R, "d_o": D_O_FACTOR * max(ridges), "d_o_factor": D_O_FACTOR,
    }

    # ---- Rouse histograms, when the runs exist ---------------------------------
    files = sorted(glob.glob(os.path.join(a.rouse, "output_cpu*.json"))) if a.rouse else []
    hist_path = os.path.join(a.out, "fig5_hist.csv")
    if not files:
        print("no Rouse output (data/*/fig5_rouse_toy) -- panel (d) will show theory only.\n"
              "Run: bash figure_production/fig5/run_fig5_rouse.sh")
        scal["rouse"] = None
    else:
        runs = []
        for fpath in files:
            for dd in json.load(open(fpath)).values():
                runs.extend(dd.values())
        assert all(np.isclose(r["alpha"], PHI[[K_LOW, K_HIGH]], rtol=1e-12).any() for r in runs), \
            "a Rouse run's phi_inf is neither of the two plotted concentrations"
        edges = np.linspace(R, L, N_BINS + 1)
        cols, stats = {}, {}
        for tag, k in (("low", K_LOW), ("high", K_HIGH)):
            mine = [r for r in runs if np.isclose(r["alpha"], PHI[k], rtol=1e-12)]
            assert mine, "no Rouse runs at phi_inf = %.6g (k = %d)" % (PHI[k], k)
            # the potential the chain felt is the one we plot: same formula, so
            # agreement to rounding
            dev = max(float(np.max(np.abs(np.asarray(r["total_potential"]) - U_tot[k]))) for r in mine)
            assert dev < 1e-9, "%s: Rouse potential differs from U_tot by %.3g" % (tag, dev)
            samples = np.concatenate([np.asarray(r["distance_traj"], float)[int(BURN_IN * len(r["distance_traj"])):]
                                      for r in mine])
            cnt, _ = np.histogram(samples, bins=edges)
            inside = cnt.sum()
            cols[tag] = cnt / (inside * np.diff(edges))          # density over [R, L]
            cols[tag + "_counts"] = cnt
            stats[tag] = {"n_runs": len(mine), "seeds": [r.get("seed") for r in mine],
                          "n_samples": int(samples.size),
                          "frac_outside_range": float(1 - inside / samples.size),
                          "P_closed": float(np.sum(samples < ridges[k]) / inside),   # over [R, L], as p_th
                          "runtime_min_mean": float(np.mean([r["runtime_min"] for r in mine])),
                          "potential_max_dev": dev}
        with open(hist_path, "w") as fh:
            fh.write("x_lo,x_hi,p_low,p_high,n_low,n_high\n")
            for i in range(N_BINS):
                fh.write("%.17g,%.17g,%.17g,%.17g,%d,%d\n" % (edges[i], edges[i + 1], cols["low"][i],
                         cols["high"][i], cols["low_counts"][i], cols["high_counts"][i]))
        scal["rouse"] = {"dir": os.path.relpath(a.rouse, repo), "burn_in": BURN_IN,
                         "n_bins": N_BINS, **stats}
        for tag, k in (("low", K_LOW), ("high", K_HIGH)):
            s = stats[tag]
            p_eq = np.exp(-(U_tot[k] - U_tot[k].min()))
            m = d < ridges[k]
            trapz = getattr(np, "trapezoid", None) or np.trapz     # numpy 1.x and 2.x
            p_th = trapz(p_eq[m], d[m]) / trapz(p_eq, d)
            print("%-4s %d runs, %d samples, %.2g outside [R, L], P(d < d*) Rouse %.3f vs theory %.3f"
                  % (tag, s["n_runs"], s["n_samples"], s["frac_outside_range"], s["P_closed"], p_th))

    json.dump(scal, open(os.path.join(a.out, "fig5_scalars.json"), "w"), indent=1)
    print("wrote", ", ".join(sorted(os.path.basename(p) for p in glob.glob(os.path.join(a.out, "fig5_*")))))


if __name__ == "__main__":
    main()
