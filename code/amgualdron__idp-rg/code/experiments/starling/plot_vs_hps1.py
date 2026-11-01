#!/usr/bin/env python
"""STARLING's R_g histograms against HPS1 MD, for the finished proteins of a run.

    mamba run -n starling python experiments/starling/plot_vs_hps1.py runs/starling/check_hps1

Writes to the run directory:
    vs_hps1.png   top row: R_g in nm, with the Mpipi-GG mean marked; bottom row: R_g / <R_g>, which removes the size
                  difference and leaves the shape
    vs_hps1.csv   one row per protein: moments of both, W1, shape W1 and split-half noise floors

W1 is the Wasserstein distance. Noise floors are W1 between two halves of the same set. HPS1 is split into the first and
second half of its trajectory, so slow drift counts as noise. These halves are half the size of the full sets, so a floor
overstates the noise of the full comparison by about sqrt(2).
"""
import argparse
import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

ROOT = pathlib.Path(__file__).resolve().parents[2]
HPS1, STARLING = "#2a78d6", "#eb6834"
SURFACE, INK, MUTED, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e6e5e0"


def compare(name, N, starling, hps1, rng):
    half = len(hps1) // 2
    shuffled = rng.permutation(starling)
    return dict(
        name=name, N=N, n_starling=len(starling), n_hps1=len(hps1),
        mean_starling=starling.mean(), mean_hps1=hps1.mean(), mean_diff_pct=100 * (starling.mean() / hps1.mean() - 1),
        sd_starling=starling.std(ddof=1), sd_hps1=hps1.std(ddof=1),
        sd_ratio=starling.std(ddof=1) / hps1.std(ddof=1),
        cv_starling=starling.std(ddof=1) / starling.mean(), cv_hps1=hps1.std(ddof=1) / hps1.mean(),
        skew_starling=stats.skew(starling), skew_hps1=stats.skew(hps1),
        exkurt_starling=stats.kurtosis(starling), exkurt_hps1=stats.kurtosis(hps1),
        w1_nm=stats.wasserstein_distance(starling, hps1),
        w1_shape=stats.wasserstein_distance(starling / starling.mean(), hps1 / hps1.mean()),
        floor_w1_starling_nm=stats.wasserstein_distance(shuffled[:len(shuffled) // 2], shuffled[len(shuffled) // 2:]),
        floor_w1_hps1_nm=stats.wasserstein_distance(hps1[:half], hps1[half:]))


def outline(ax, values, edges, color, label=None):
    density, _ = np.histogram(values, bins=edges, density=True)
    ax.stairs(density, edges, color=color, lw=2, label=label)
    ax.stairs(density, edges, color=color, alpha=0.10, fill=True)


def style(ax):
    ax.set_facecolor(SURFACE)
    ax.set_yticks([])
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.spines["bottom"].set_color(GRID)
    ax.tick_params(colors=MUTED, labelsize=9, length=0)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("run", type=pathlib.Path)
    p.add_argument("--trajectories", type=pathlib.Path, default=ROOT / "runs" / "adopt" / "traj",
                   help="directory of <name>_pooled_timeseries.csv, as for paper/make_figures.py")
    args = p.parse_args()

    files = pd.read_csv(ROOT / "data" / "benchmark.csv").set_index("name").trajectories
    summary = pd.read_csv(args.run / "summary.csv").set_index("name")
    done = summary[summary.status == "ok"].sort_values("N")
    rng = np.random.default_rng(0)
    rows, panels = [], []
    for name, r in done.iterrows():
        paths = [args.trajectories / f"{f}_pooled_timeseries.csv" for f in files[name].split(";")]
        hps1 = np.concatenate([pd.read_csv(path, usecols=["rg"]).rg.to_numpy() for path in paths])
        starling = np.load(args.run / "samples" / f"{name}.npz")["rg_nm"]
        rows.append(compare(name, int(r.N), starling, hps1, rng))
        panels.append((name, int(r.N), starling, hps1, r.mpipi_gg_rg_nm))
    if not panels:
        raise SystemExit(f"{args.run}: no finished proteins")
    table = pd.DataFrame(rows)
    table.to_csv(args.run / "vs_hps1.csv", index=False)

    fig, axes = plt.subplots(2, len(panels), figsize=(3.1 * len(panels), 5.2), squeeze=False, facecolor=SURFACE)
    for col, (name, N, starling, hps1, mpipi) in enumerate(panels):
        top, bottom = axes[0, col], axes[1, col]
        edges = np.histogram_bin_edges(np.concatenate([starling, hps1]), bins="fd")
        outline(top, hps1, edges, HPS1, "HPS1 MD")
        outline(top, starling, edges, STARLING, "STARLING")
        top.axvline(mpipi, color=MUTED, lw=1, ls=(0, (3, 3)), label="Mpipi-GG mean")
        top.set_title(f"{name}\nN = {N}", fontsize=10, color=INK, loc="left")
        w1 = rows[col]
        top.text(0.98, 0.95, f"W1 {w1['w1_nm']:.2f} nm\nSD ratio {w1['sd_ratio']:.2f}", transform=top.transAxes, ha="right",
                 va="top", fontsize=8.5, color=MUTED)

        a, b = starling / starling.mean(), hps1 / hps1.mean()
        edges = np.histogram_bin_edges(np.concatenate([a, b]), bins="fd")
        outline(bottom, b, edges, HPS1)
        outline(bottom, a, edges, STARLING)
        bottom.text(0.98, 0.95, f"shape W1 {w1['w1_shape']:.3f}", transform=bottom.transAxes, ha="right", va="top",
                    fontsize=8.5, color=MUTED)
        for ax in (top, bottom):
            style(ax)
        top.set_xlabel("$R_g$ (nm)", fontsize=9, color=MUTED)
        bottom.set_xlabel(r"$R_g\,/\,\langle R_g\rangle$", fontsize=9, color=MUTED)
    axes[0, 0].set_ylabel("density, in nm", fontsize=9, color=MUTED)
    axes[1, 0].set_ylabel("density, shape only", fontsize=9, color=MUTED)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper right", ncols=3, frameon=False, fontsize=9.5, labelcolor=INK)
    fig.suptitle("STARLING against HPS1 MD: $R_g$ distributions", x=0.01, ha="left", fontsize=12, color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(args.run / "vs_hps1.png", dpi=200, facecolor=SURFACE)

    pd.set_option("display.width", 220)
    print(table[["name", "N", "mean_starling", "mean_hps1", "mean_diff_pct", "sd_starling", "sd_hps1", "sd_ratio", "cv_starling",
                 "cv_hps1", "skew_starling", "skew_hps1", "w1_nm", "w1_shape", "floor_w1_starling_nm",
                 "floor_w1_hps1_nm"]].round(3).to_string(index=False))


if __name__ == "__main__":
    main()
