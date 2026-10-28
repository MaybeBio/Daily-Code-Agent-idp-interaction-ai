#!/usr/bin/env python3
"""
prepare_fig2_data.py -- reduce the raw proximal-promoter runs to the small
tables that main Fig. 2 actually needs, and write them to
Droplet/processed_data/.

Run this once, on a machine that has the 3.3 GB of simulation output. Thereafter
`fig2.ipynb` plots the figure from the ~100 KB written here, with no access to
the raw runs.

What is kept, and why
---------------------
Fig. 2 consumes the raw traces in exactly two ways:

  * as 16-element summaries, one value per TF concentration -- conditional
    means, occupancies, mean transcription; and
  * as histograms on fixed bin edges (the two heatmaps).

Nothing else survives into the figure. So the pooled traces themselves (18.1 M
values per observable, ~435 MB for the three) are reduced here to per-
concentration summaries and binned densities. The histograms are stored raw
(`density=True`, unsmoothed, unnormalised), so the notebook still performs the
gaussian blur and the per-column normalisation and both remain adjustable.

The shape fits behind panels (d,e) are deliberately NOT precomputed: they take
~8 s and are where the modelling lives, so they belong in the notebook.

Outputs (Droplet/processed_data/)
---------------------------------
  fig2_concentration_summary.csv     16 rows, one per concentration
  fig2_hist_local_concentration.csv  27 bins x 16 concentrations -> panel (b)
  fig2_hist_pattern_distance.csv     69 bins x 16 concentrations -> panel (a)
  fig2_transcription_bootstrap.csv   200 rows (100 reps x 2 states) -> panel (c)
  fig2_scalars.json                  constants, bin definitions, provenance

This directory is deliberately not `data/`, which .gitignore excludes wholesale
so that the raw runs never enter the repo. These tables are small enough to
commit, and the `fig2_` prefix leaves room for the other figures.

Cell provenance
---------------
Transcribed from viz.ipynb, section "Probe-enhancer transcription". Cell indices
are positions in the notebook's `cells` array (stored execution counts are
duplicated across the notebook and unusable):

    cell 110  run-major loader
    cell 118  bound/unbound times with hysteresis
    cell 123  theta
    cell 126  conditional means + run-level bootstrap
    cell 128  mean RNAP off-time -> kb, ku
    cell 163  the pooling and the histogram bin edges

Usage
-----
    python3 prepare_fig2_data.py --data-root ~/Archive/Work/droplet_data

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

# thresholds, from cells 118 and 126
D_HYST = 6.5          # hysteresis thresholds d_o = d_x
THRESH_CLOSED = 6.5   # pattern_distance < thresh  <=>  "closed"

# bin edges, from cell 163
YE_LOCAL_CONC = (np.arange(28) - 0.5) / 26.0      # the k/26 lattice grid
N_BINS_DIST = 70                                  # edges; 69 bins
DIST_PERCENTILE = 99.5                            # upper edge of the distance axis

NBOOT = 100           # cell 126


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Reduce the raw runs to the tables Fig. 2 needs.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--data-root", default="data",
                   help="directory holding the dated run folders, i.e. the "
                        "parent of 20260423/")
    p.add_argument("--repo", default=os.path.dirname(os.path.abspath(__file__)),
                   help="the Droplet checkout, or a directory containing it")
    p.add_argument("--date", default="20260423",
                   help="run-date folder under --data-root")
    p.add_argument("--out", default=None,
                   help="output directory (default: figure_production/processed_data)")
    p.add_argument("--seed", type=int, default=0,
                   help="seed for the run-level bootstrap, which is unseeded "
                        "in the notebook; -1 leaves the global RNG untouched")
    return p.parse_args(argv)


def log(msg):
    print("[{}] {}".format(time.strftime("%H:%M:%S"), msg), flush=True)


def resolve_repo(repo):
    """Locate the Droplet checkout (the directory holding viz_util/).

    Accepts the checkout itself, a directory containing it, or any ancestor of
    this script -- it lives in figure_production/fig2/ inside the repo."""
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
        return subprocess.check_output(
            ["git", "-C", repo, "rev-parse", "HEAD"],
            stderr=subprocess.DEVNULL).decode().strip()
    except Exception:
        return None


# ----------------------------------------------------------------------------
# load
# ----------------------------------------------------------------------------

def load_runs(repo, data_root, date):
    """cells 110 and 118."""
    sys.path.insert(0, repo)
    import viz_util.load_transcription as lt
    from viz_util.hysteresis_times import compute_hysteresis_times

    src = os.path.join(os.path.abspath(os.path.expanduser(data_root)), date)
    if not os.path.isdir(src):
        raise SystemExit("no such data directory: {}\n"
                         "pass --data-root (the parent of {}/)".format(src, date))

    groups = lt.make_groups([
        (c, [("{}/{}".format(src, EXP.format(c)), range(1, NRUNS + 1))])
        for c in CONCS
    ])

    log("reading {} concentrations x {} runs from {} ...".format(
        NCONC, NRUNS, src))
    t0 = time.time()
    (on_times, off_times, bound_times, unbound_times,
     pattern_distances, transcription_activity, binding_activity,
     local_concentration, run_ids) = lt.load(groups, n1=0, n=0)
    log("read in {:.0f} s".format(time.time() - t0))

    # cell 110: rescale the event times to MC sweeps
    params = np.array(CONCS)
    total_vol_frac = (params + 0.02) * 62 ** 3
    sampling_multiplier = 50000 / total_vol_frac
    for i in range(NCONC):
        for k in range(NRUNS):
            on_times[i][k] = on_times[i][k] * sampling_multiplier[i]
            off_times[i][k] = off_times[i][k] * sampling_multiplier[i]

    # cell 118. NB the notebook unpacks (unbound, bound) from a function that
    # returns (bound, unbound). The swap is invisible in Fig. 2 -- these feed
    # only tb/tu/tau in cell 163, which the plot never uses -- so they are not
    # carried into the saved tables at all.
    _ = compute_hysteresis_times(pattern_distances, d_o=D_HYST, d_x=D_HYST,
                                 nparams=NCONC, nruns=NRUNS)

    return dict(off_times=off_times,
                pattern_distances=pattern_distances,
                transcription_activity=transcription_activity,
                local_concentration=local_concentration,
                src=src)


# ----------------------------------------------------------------------------
# reduce
# ----------------------------------------------------------------------------

def split_closed_open(values, dists, cidx, ridx, thresh=THRESH_CLOSED, closed=True):
    """cell 126's compute_local_conc: values of one observable in one run,
    restricted to the frames where the probe is closed (or open)."""
    v = np.array(values[cidx][ridx])
    d = np.array(dists[cidx][ridx])
    idx = np.where(d < thresh) if closed else np.where(d >= thresh)
    return idx, v[idx]


def reduce_runs(R, seed):
    """cells 123, 126 and 128, plus the pooling and binning of cell 163."""
    dists = R["pattern_distances"]
    conc = R["local_concentration"]
    transc = R["transcription_activity"]

    # ---- cell 123 ---------------------------------------------------------
    theta = 40 * 1e3 / (6e23 * (25 * 1e-9) ** 3)      # volume fraction -> uM

    # ---- cell 126: conditional means and occupancies ----------------------
    phi_local_closed = np.zeros(NCONC)
    phi_local_open = np.zeros(NCONC)
    p_closed = np.zeros(NCONC)
    p_open = np.zeros(NCONC)

    for c in range(NCONC):
        pooled_closed, pooled_open = [], []
        n_closed = n_open = n_total = 0
        for r in range(NRUNS):
            ic, vc = split_closed_open(conc, dists, c, r, closed=True)
            io, vo = split_closed_open(conc, dists, c, r, closed=False)
            pooled_closed = np.concatenate((pooled_closed, vc))
            pooled_open = np.concatenate((pooled_open, vo))
            n_closed += len(ic[0])
            n_open += len(io[0])
            n_total += len(conc[c][r])
        phi_local_closed[c] = np.mean(pooled_closed)
        phi_local_open[c] = np.mean(pooled_open)
        p_closed[c] = n_closed / n_total
        p_open[c] = n_open / n_total

    transcription_closed = np.zeros(NCONC)
    transcription_open = np.zeros(NCONC)
    for c in range(NCONC):
        for arr, is_closed in ((transcription_closed, True),
                               (transcription_open, False)):
            arr[c] = np.mean(np.concatenate([
                split_closed_open(transc, dists, c, r, closed=is_closed)[1]
                for r in range(NRUNS)]))

    # ---- cell 126: run-level (cluster) bootstrap --------------------------
    if seed >= 0:
        np.random.seed(seed)
        log("bootstrap seeded with {}".format(seed))
    else:
        log("bootstrap left unseeded (as in the notebook)")

    boot_closed = np.zeros((NBOOT, NCONC))
    boot_open = np.zeros((NBOOT, NCONC))
    for i in range(NBOOT):
        runs_ = np.random.randint(0, NRUNS, size=(NCONC, NRUNS))
        for c in range(NCONC):
            for arr, is_closed in ((boot_closed, True), (boot_open, False)):
                pooled = np.concatenate([
                    split_closed_open(transc, dists, c, r, closed=is_closed)[1]
                    for r in runs_[c]])
                arr[i, c] = np.mean(pooled) if len(pooled) else np.nan
    log("bootstrap done ({} replicates)".format(NBOOT))

    # ---- cell 128: the one kinetic input the figure needs -----------------
    mean_toff_rnap = float(np.mean(
        [np.mean(np.concatenate(R["off_times"][k])) for k in range(NCONC)]))

    # ---- cell 163: pool across runs, then bin -----------------------------
    def pool(data):
        out = [np.concatenate([np.ravel(data[c][r]) for r in range(NRUNS)])
               for c in range(NCONC)]
        return [x[np.isfinite(x)] for x in out]

    pooled_transc, pooled_conc, pooled_dist = pool(transc), pool(conc), pool(dists)

    # the distance axis runs to a percentile of the pooled data, so the cutoff
    # is frozen here and recorded in scalars.json
    dist_hi = float(np.percentile(np.concatenate(pooled_dist), DIST_PERCENTILE))
    ye_dist = np.linspace(0.0, dist_hi, N_BINS_DIST)

    def histogram_matrix(pooled, edges):
        return np.column_stack([
            np.histogram(p, bins=edges, density=True)[0] for p in pooled])

    H_conc = histogram_matrix(pooled_conc, YE_LOCAL_CONC)
    H_dist = histogram_matrix(pooled_dist, ye_dist)

    params = np.array(CONCS)
    return dict(
        theta=theta, mean_toff_rnap=mean_toff_rnap, dist_hi=dist_hi,
        summary=pd.DataFrame({
            "c_phi": params,
            "c_uM": params * theta,
            "phi_local_closed": phi_local_closed,
            "phi_local_open": phi_local_open,
            "p_closed": p_closed,
            "p_open": p_open,
            "transcription_closed": transcription_closed,
            "transcription_open": transcription_open,
            "mean_transcription": [p.mean() for p in pooled_transc],
            "mean_local_concentration": [p.mean() for p in pooled_conc],
            "mean_pattern_distance": [p.mean() for p in pooled_dist],
            "n_samples": [p.size for p in pooled_dist],
        }),
        hist_conc=hist_frame(H_conc, YE_LOCAL_CONC, params),
        hist_dist=hist_frame(H_dist, ye_dist, params),
        bootstrap=bootstrap_frame(boot_closed, boot_open, params),
    )


def hist_frame(H, edges, params):
    """bins x concentrations, carrying both bin edges so the notebook can
    rebuild the pcolormesh grid without recomputing the percentile cutoff."""
    df = pd.DataFrame(H, columns=["c_{:.3f}".format(c) for c in params])
    df.insert(0, "bin_right", edges[1:])
    df.insert(0, "bin_left", edges[:-1])
    return df


def bootstrap_frame(boot_closed, boot_open, params):
    cols = ["c_{:.3f}".format(c) for c in params]
    a = pd.DataFrame(boot_closed, columns=cols)
    a.insert(0, "state", "closed")
    b = pd.DataFrame(boot_open, columns=cols)
    b.insert(0, "state", "open")
    out = pd.concat([a, b], ignore_index=True)
    out.insert(0, "replicate", list(range(NBOOT)) * 2)
    return out


# ----------------------------------------------------------------------------
# write
# ----------------------------------------------------------------------------

def write_outputs(D, out_dir, args, repo, src):
    os.makedirs(out_dir, exist_ok=True)

    paths = {
        "fig2_concentration_summary.csv": D["summary"],
        "fig2_hist_local_concentration.csv": D["hist_conc"],
        "fig2_hist_pattern_distance.csv": D["hist_dist"],
        "fig2_transcription_bootstrap.csv": D["bootstrap"],
    }
    # Full round-trip precision, deliberately. The LOOCV that picks the
    # smoothing weight in panel (d) is numerically ill-conditioned at the small
    # end of its lambda grid -- there the roughness penalty nearly vanishes, the
    # system approaches interpolation, and only the 1e-12*trace(A)/n ridge
    # inside fit_shape conditions it. Rounding phi_local_closed to 10
    # significant digits was enough to shift the CV score by 20% and move the
    # selected lambda by an order of magnitude (3.2e-6 -> 3.2e-7, edf 6.1 ->
    # 7.6). %.17g round-trips float64 exactly, so the notebook reproduces the
    # published fit rather than a neighbouring one.
    for name, df in paths.items():
        df.to_csv(os.path.join(out_dir, name), index=False,
                  float_format="%.17g")

    scalars = {
        "theta": D["theta"],
        "theta_units": "uM per unit volume fraction",
        "mean_toff_rnap": D["mean_toff_rnap"],
        "K_D_bands": 2.3e-3,
        "K_D_curves": 2.33e-3,
        "K_D_note": ("cell 128 uses 2.3e-3 and cell 132 recomputes kb/ku with "
                     "2.33e-3, so in panel (e) the bootstrap bands and the "
                     "central curves rest on marginally different K_D. "
                     "Reproduced as published."),
        "phi_R": 0.02,
        "threshold_closed": THRESH_CLOSED,
        "d_hysteresis": D_HYST,
        "dist_bin_hi": D["dist_hi"],
        "dist_bin_percentile": DIST_PERCENTILE,
        "dist_n_edges": N_BINS_DIST,
        "local_conc_bin_note": "edges = (arange(28) - 0.5)/26, the k/26 lattice grid",
        "n_conc": NCONC,
        "n_runs": NRUNS,
        "n_bootstrap": NBOOT,
        "bootstrap_seed": args.seed,
        "source_data": src,
        "source_date_folder": args.date,
        "repo_commit": git_commit(repo),
        "generated_by": os.path.basename(__file__),
        "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "source_cells": "viz.ipynb cells 110, 118, 123, 126, 128, 163",
    }
    with open(os.path.join(out_dir, "fig2_scalars.json"), "w") as f:
        json.dump(scalars, f, indent=2)
        f.write("\n")

    total = 0
    for name in list(paths) + ["fig2_scalars.json"]:
        p = os.path.join(out_dir, name)
        n = os.path.getsize(p)
        total += n
        print("  {:32s} {:>8.1f} KB".format(name, n / 1024))
    print("  {:32s} {:>8.1f} KB".format("TOTAL", total / 1024))


def main(argv=None):
    args = parse_args(argv)
    repo = resolve_repo(args.repo)
    # processed_data/ is shared by the figures, one level up from this
    # script: figure_production/processed_data/
    out_dir = args.out or os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        "processed_data")
    log("repo {}".format(repo))

    R = load_runs(repo, args.data_root, args.date)
    D = reduce_runs(R, args.seed)
    log("writing to {}".format(out_dir))
    write_outputs(D, out_dir, args, repo, R["src"])


if __name__ == "__main__":
    main()
