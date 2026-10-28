#!/usr/bin/env python3
"""
prepare_fig3_data.py -- reduce the raw proximal-promoter runs to the small table
Fig. 3 needs, and write it to Droplet/processed_data/.

Fig. 3 has two panels:

  A  sequential probability ratio test: a realization of the log-likelihood
     ratio R(t) against the +/-K decision thresholds, with the first crossing
     T_dec and the drift asymptote V t.
  B  local Fisher information rate: F_o, F_x, F_local and F_free versus
     concentration -- the same quantity as panel (e) of Fig. 2.

Only panel B needs the simulation. Panel A is generated in the notebook from
panel B's own output, because the paper's SPRT is a drift-diffusion process
whose drift is set by the Fisher information rate (lattice.tex, Eq. `drift`):

    V = (1/2) (dc/c)^2 Fdot(c),   with diffusion constant equal to V,
    <T_dec> = (K/V) tanh(K/2),    K = log((1-eps)/eps).

There is deliberately no empirical log-likelihood trace here. The paper is
explicit that P(O(t)|c) is not estimable from these runs -- "the DNA
configuration is a latent variable of the dynamics: it is not included in
O(t) ... and must be marginalized out, which is computationally intractable."
An attempt to build R(t) from the accumulated transcription instead also fails
on its own terms: the across-run variance of the accumulated signal grows like
t^2 rather than t, and a within-run batch-means estimate of the variance rate
has no plateau, because the open/closed polymer state barely switches within a
run. There is no diffusive regime in this data, so no variance rate to divide
by. The honest construction is the paper's: simulate the process, take its
drift from the measured Fisher rate.

Outputs (Droplet/processed_data/)
---------------------------------
  fig3_fisher_input.csv   16 rows, one per concentration  -> panel B
  fig3_scalars.json       constants, SPRT parameters, provenance

Cell provenance
---------------
Transcribed from viz.ipynb, section "Probe-enhancer transcription"; cell indices
are positions in that notebook's `cells` array (its stored execution counts are
duplicated and unusable): cell 110 (loader), 123 (theta), 126 (conditional
means), 128 (mean RNAP off-time -> kb, ku).

Usage
-----
    python3 prepare_fig3_data.py --data-root ~/Archive/Work/droplet_data

`--repo` defaults to this script's directory and may be the Droplet checkout or
a directory containing it. Requires numpy and pandas.
"""

import argparse
import json
import os
import subprocess
import sys
import time
from datetime import datetime, timezone

import numpy as np
import pandas as pd

CONCS = [0.003, 0.004, 0.005, 0.006, 0.007, 0.008,
         0.009, 0.010, 0.011, 0.012, 0.013, 0.014,
         0.015, 0.016, 0.017, 0.018]
NCONC = len(CONCS)
NRUNS = 16
EXP = "proximal_promoter_transc_d_50_c_{:.3f}"
THRESH_CLOSED = 6.5          # pattern_distance < thresh <=> "closed"  (cell 126)

# SPRT parameters, from lattice.tex: a 10% concentration difference decided with
# error 0.32 ("call the shot correctly at least 68% of the time"), and
# K = log((1-eps)/eps).
EPSILON = 0.32
DELTA_C_OVER_C = 0.10


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Reduce the raw runs to the table Fig. 3 needs.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--data-root", default="data")
    p.add_argument("--repo", default=os.path.dirname(os.path.abspath(__file__)))
    p.add_argument("--date", default="20260423")
    p.add_argument("--out", default=None,
                   help="output directory (default: figure_production/processed_data)")
    return p.parse_args(argv)


def log(msg):
    print("[{}] {}".format(time.strftime("%H:%M:%S"), msg), flush=True)


def resolve_repo(repo):
    """Locate the Droplet checkout (the directory holding viz_util/).

    Accepts the checkout itself, a directory containing it, or any ancestor of
    this script -- it lives in figure_production/fig3/ inside the repo."""
    root = os.path.abspath(os.path.expanduser(repo))
    cands = [root, os.path.join(root, "Droplet")]
    p = root
    for _ in range(4):                      # walk up towards the repo root
        p = os.path.dirname(p)
        cands.append(p)
    for cand in cands:
        if os.path.isdir(os.path.join(cand, "viz_util")):
            return cand
    raise SystemExit("no viz_util/ found at or above {}\n"
                     "pass --repo pointing at the Droplet checkout".format(root))

def git_commit(repo):
    try:
        return subprocess.check_output(["git", "-C", repo, "rev-parse", "HEAD"],
                                       stderr=subprocess.DEVNULL).decode().strip()
    except Exception:
        return None


def load_runs(repo, data_root, date):
    sys.path.insert(0, repo)
    import viz_util.load_transcription as lt

    src = os.path.join(os.path.abspath(os.path.expanduser(data_root)), date)
    if not os.path.isdir(src):
        raise SystemExit("no such data directory: {}".format(src))

    groups = lt.make_groups([
        (c, [("{}/{}".format(src, EXP.format(c)), range(1, NRUNS + 1))])
        for c in CONCS])

    log("reading {} concentrations x {} runs from {} ...".format(NCONC, NRUNS, src))
    t0 = time.time()
    (on_times, off_times, bound_times, unbound_times,
     pattern_distances, transcription_activity, binding_activity,
     local_concentration, run_ids) = lt.load(groups, n1=0, n=0)
    log("read in {:.0f} s".format(time.time() - t0))

    # cell 110: rescale the event times to MC sweeps
    params = np.array(CONCS)
    mult = 50000 / ((params + 0.02) * 62 ** 3)
    for i in range(NCONC):
        for k in range(NRUNS):
            off_times[i][k] = off_times[i][k] * mult[i]

    return dict(off_times=off_times, pattern_distances=pattern_distances,
                local_concentration=local_concentration, src=src)


def reduce_runs(R):
    """cell 126's conditional means and occupancies, plus cell 128's kinetics."""
    dists, conc = R["pattern_distances"], R["local_concentration"]

    def split(cidx, ridx, closed):
        v = np.array(conc[cidx][ridx])
        d = np.array(dists[cidx][ridx])
        idx = np.where(d < THRESH_CLOSED) if closed else np.where(d >= THRESH_CLOSED)
        return idx, v[idx]

    phi_closed = np.zeros(NCONC)
    phi_open = np.zeros(NCONC)
    p_closed = np.zeros(NCONC)
    p_open = np.zeros(NCONC)

    for c in range(NCONC):
        pooled_c, pooled_o = [], []
        n_c = n_o = n_tot = 0
        for r in range(NRUNS):
            ic, vc = split(c, r, True)
            io, vo = split(c, r, False)
            pooled_c = np.concatenate((pooled_c, vc))
            pooled_o = np.concatenate((pooled_o, vo))
            n_c += len(ic[0])
            n_o += len(io[0])
            n_tot += len(conc[c][r])
        phi_closed[c] = np.mean(pooled_c)
        phi_open[c] = np.mean(pooled_o)
        p_closed[c] = n_c / n_tot
        p_open[c] = n_o / n_tot

    mean_toff_rnap = float(np.mean(
        [np.mean(np.concatenate(R["off_times"][k])) for k in range(NCONC)]))

    theta = 40 * 1e3 / (6e23 * (25 * 1e-9) ** 3)      # volume fraction -> uM
    params = np.array(CONCS)
    df = pd.DataFrame({
        "c_phi": params, "c_uM": params * theta,
        "phi_local_closed": phi_closed, "phi_local_open": phi_open,
        "p_closed": p_closed, "p_open": p_open,
    })
    return df, theta, mean_toff_rnap


def main(argv=None):
    args = parse_args(argv)
    repo = resolve_repo(args.repo)
    # processed_data/ is shared by the figures, one level up from this
    # script: figure_production/processed_data/
    out_dir = args.out or os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        "processed_data")
    os.makedirs(out_dir, exist_ok=True)

    R = load_runs(repo, args.data_root, args.date)
    fisher_df, theta, mean_toff = reduce_runs(R)

    # %.17g round-trips float64 exactly. This matters: the LOOCV behind panel B
    # is ill-conditioned at the small end of its lambda grid, and a rounded
    # input selects a different smoothing weight.
    fisher_df.to_csv(os.path.join(out_dir, "fig3_fisher_input.csv"),
                     index=False, float_format="%.17g")

    K = float(np.log((1 - EPSILON) / EPSILON))
    scalars = dict(
        theta=theta, theta_units="uM per unit volume fraction",
        mean_toff_rnap=mean_toff, K_D_bands=2.3e-3, K_D_curves=2.33e-3,
        K_D_note=("cell 128 uses 2.3e-3 for the bootstrap bands and cell 132 "
                  "recomputes kb/ku with 2.33e-3 for the central curves"),
        phi_R=0.02, threshold_closed=THRESH_CLOSED,
        epsilon=EPSILON, K=K, delta_c_over_c=DELTA_C_OVER_C,
        sprt_note=("lattice.tex: R(t)=log P(O(t)|c)/P(O(t)|c+dc) is a "
                   "drift-diffusion process with drift V and EQUAL diffusion "
                   "constant V; V=(1/2)(dc/c)^2*Fdot(c) and "
                   "<T_dec>=(K/V)tanh(K/2) with K=log((1-eps)/eps). "
                   "eps=0.32 and dc/c=0.1 are the values quoted in the text."),
        n_conc=NCONC, n_runs=NRUNS,
        source_data=R["src"], source_date_folder=args.date,
        repo_commit=git_commit(repo), generated_by=os.path.basename(__file__),
        generated_at=datetime.now(timezone.utc).isoformat(timespec="seconds"),
        source_cells="viz.ipynb cells 110, 123, 126, 128")
    with open(os.path.join(out_dir, "fig3_scalars.json"), "w") as f:
        json.dump(scalars, f, indent=2)
        f.write("\n")

    log("K = log((1-{})/{}) = {:.4f}".format(EPSILON, EPSILON, K))
    total = 0
    for name in ("fig3_fisher_input.csv", "fig3_scalars.json"):
        n = os.path.getsize(os.path.join(out_dir, name))
        total += n
        print("  {:28s} {:>8.1f} KB".format(name, n / 1024))
    print("  {:28s} {:>8.1f} KB".format("TOTAL", total / 1024))


if __name__ == "__main__":
    main()
