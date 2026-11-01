#!/usr/bin/env python
"""Shape of the R_g distribution: the released network (model/network.eqx), STARLING and the HPS1 trajectories.

    mamba run -n idp-proteome python experiments/starling/score_shape.py runs/starling/benchmark

Each distribution is divided by its own mean, so size differences are gone, and the pairs are compared with Wasserstein-1:

    W1 = integral over u in (0, 1) of |Q_a(u) - Q_b(u)|,   Q = quantile function of R_g / <R_g>

W1 is reported x 100, so it reads as a percentage of <R_g>. Writes to the run directory:

    shape_w1.png   four scatter plots, one point per protein: network-HPS1 against network-STARLING and against
                   STARLING-HPS1, then network-HPS1 and STARLING-HPS1 against the shared-shape baseline
    shape_w1.csv   one row per protein: the three W1, their noise floors, the shared-shape baseline and the moments

Noise floor: the W1 that sampling alone gives if the two distributions are identical. It draws from the network's
distribution for the protein, 2000 independent samples for STARLING and n_w1 for HPS1. n_eff = frames / integrated
autocorrelation time of R_g (Geyer's initial positive sequence, per trajectory, summed) overstates W1's noise, so it is
calibrated once over all proteins: n_w1 = n_eff / c^2, with c the median of split-half W1 / (2 x the floor at n_eff). The
two halves of one trajectory hold half the data each, so independent noise makes their W1 twice the full set's floor.

Equal noise: the network is an exact curve and STARLING is 2000 samples, so STARLING's W1 carries more sampling noise.
w1_net2000_hps1 is the mean W1 to HPS1 of 20 sets of 2000 samples drawn from the network, for a like-for-like comparison.

Shared-shape baseline: one shape for every protein, the average quantile function of the HPS1 shapes of the other protein
families. A predictor that does not beat it has captured no sequence-specific shape.
"""
import argparse
import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

import idprg
from idprg import distribution as dist

ROOT = pathlib.Path(__file__).resolve().parents[2]
M = 20000
U = (np.arange(M) + 0.5) / M
N_BOOT = 200
MIN_N = 20    # the network's stated domain
POINT, SURFACE, INK, MUTED, GRID = "#2a78d6", "#fcfcfb", "#0b0b0b", "#52514e", "#e6e5e0"


def tau_geyer(x):
    """Integrated autocorrelation time in frames, Geyer's initial positive sequence."""
    x = x - x.mean()
    f = np.fft.rfft(x, 2 * len(x))
    ac = np.fft.irfft(f * np.conj(f))[:len(x)]
    ac /= ac[0]
    tau = 1.0
    for k in range(1, len(x) // 2 - 1, 2):
        g = ac[k] + ac[k + 1]
        if g <= 0:
            break
        tau += 2 * g
    return tau


def q_samples(x):
    """Quantile function of x / <x> on U."""
    x = np.sort(x) / x.mean()
    return x[(U * len(x)).astype(int)]


def w1(a, b):
    return 100 * float(np.abs(a - b).mean())


def noise_floors(q_net, n_hps1, rng):
    """Median and 95th percentile of W1 from sampling alone, for each pair, drawing from the network's shape."""
    draw = lambda n: q_samples(np.interp(rng.random(n), U, q_net))
    n_hps1 = max(int(round(n_hps1)), 2)
    out = {"net_hps1": [], "net_starling": [], "starling_hps1": []}
    for _ in range(N_BOOT):
        s, h = draw(2000), draw(n_hps1)
        out["net_hps1"].append(w1(q_net, h))
        out["net_starling"].append(w1(q_net, s))
        out["starling_hps1"].append(w1(s, h))
    return {f"floor_{k}{tag}": float(f(v)) for k, v in out.items()
            for tag, f in (("", np.median), ("_p95", lambda v: np.percentile(v, 95)))}


def score_protein(r, p, frames_files, starling, rng):
    parts = [pd.read_csv(f, usecols=["rg"]).rg.to_numpy(float) for f in frames_files]
    hps1 = np.concatenate(parts)
    n_eff = sum(len(x) / tau_geyer(x) for x in parts)
    lnA, lnSm, nu = float(p.lnA), float(np.log(p.S_m)), float(p.nu)
    q_net = dist.quantiles(U, lnA, lnSm, nu) / float(p.mean_rg)
    q_hps1 = q_samples(hps1)
    half = len(hps1) // 2
    row = dict(name=r.name, family=r.family, N=int(r.N), n_frames=len(hps1), n_eff_hps1=n_eff,
               w1_net_hps1=w1(q_net, q_hps1),
               w1_nm_net_hps1=w1(dist.quantiles(U, lnA, lnSm, nu), np.sort(hps1)[(U * len(hps1)).astype(int)]) / 100,
               split_half_hps1=w1(q_samples(hps1[:half]), q_samples(hps1[half:])),
               cv_net=float(p.std_rg / p.mean_rg), cv_hps1=hps1.std(ddof=1) / hps1.mean(),
               skew_net=float(p["skew"]), skew_hps1=stats.skew(hps1))
    if starling is not None:
        q_star = q_samples(starling)
        row |= dict(w1_net_starling=w1(q_net, q_star), w1_starling_hps1=w1(q_star, q_hps1),
                    split_half_starling=w1(q_samples(starling[:1000]), q_samples(starling[1000:])),
                    cv_starling=starling.std(ddof=1) / starling.mean(), skew_starling=stats.skew(starling))
    return row, q_net, q_hps1, q_star if starling is not None else None


def family_bootstrap(t, fn, rng, n=5000):
    """95% interval of fn(table), resampling whole protein families, because variants of one protein are not independent."""
    groups = [g for _, g in t.groupby("family")]
    out = [fn(pd.concat([groups[i] for i in rng.integers(len(groups), size=len(groups))])) for _ in range(n)]
    return np.percentile(out, [2.5, 97.5])


def report(t, rng):
    """t: the proteins scored by all three pairs."""
    print(f"\n{len(t)} proteins, {t.family.nunique()} families. Shape W1, % of <Rg>")
    rows = []
    for label, c, floor in (("network - HPS1", "w1_net_hps1", "floor_net_hps1"),
                            ("network - STARLING", "w1_net_starling", "floor_net_starling"),
                            ("STARLING - HPS1", "w1_starling_hps1", "floor_starling_hps1"),
                            ("shared shape - HPS1", "w1_baseline_hps1", "floor_net_hps1"),
                            ("shared shape - STARLING", "w1_baseline_starling", "floor_net_starling")):
        lo, hi = family_bootstrap(t, lambda d: d[c].median(), rng)
        p95 = floor + "_p95"
        rows.append(dict(pair=label, median=t[c].median(), ci_lo=lo, ci_hi=hi, p90=t[c].quantile(.9), max=t[c].max(),
                         median_floor=t[floor].median(), above_floor_p95=(t[c] > t[p95]).mean()))
    print(pd.DataFrame(rows).round(3).to_string(index=False))

    d = t.w1_starling_hps1 - t.w1_net_hps1
    lo, hi = family_bootstrap(t, lambda g: (g.w1_starling_hps1 - g.w1_net_hps1).median(), rng)
    print(f"\nSTARLING-HPS1 minus network-HPS1: median {d.median():+.3f} [{lo:+.3f}, {hi:+.3f}], network closer on "
          f"{(d > 0).sum()}/{len(d)}; Wilcoxon p = {stats.wilcoxon(d).pvalue:.2g} (treats proteins as independent)")
    d = t.w1_starling_hps1 - t.w1_net2000_hps1
    lo, hi = family_bootstrap(t, lambda g: (g.w1_starling_hps1 - g.w1_net2000_hps1).median(), rng)
    print(f"Equal noise, network as 2000 samples: median W1 {t.w1_net2000_hps1.median():.3f}, p90 {t.w1_net2000_hps1.quantile(.9):.3f}; STARLING minus network "
          f"median {d.median():+.3f} [{lo:+.3f}, {hi:+.3f}], network closer on {(d > 0).sum()}/{len(d)}")
    for label, a, b in (("network-HPS1 vs network-STARLING", "w1_net_hps1", "w1_net_starling"),
                        ("network-HPS1 vs STARLING-HPS1", "w1_net_hps1", "w1_starling_hps1")):
        print(f"Spearman {label}: {stats.spearmanr(t[a], t[b]).statistic:.2f}")
    for who, c in (("network", "w1_net_hps1"), ("STARLING", "w1_starling_hps1")):
        print(f"beats shared shape against HPS1, {who}: {(t[c] < t.w1_baseline_hps1).sum()}/{len(t)}")
    print("\nWidth and asymmetry against HPS1:")
    for who, cv, sk in (("network", "cv_net", "skew_net"), ("STARLING", "cv_starling", "skew_starling")):
        ratio = t[cv] / t.cv_hps1
        print(f"  {who}: CV ratio median {ratio.median():.3f}, median |ratio-1| {(ratio - 1).abs().median():.3f}, r {np.corrcoef(t[cv], t.cv_hps1)[0, 1]:.2f}; "
              f"skew diff median {(t[sk] - t.skew_hps1).median():+.3f}, MAE {(t[sk] - t.skew_hps1).abs().mean():.3f}, r {np.corrcoef(t[sk], t.skew_hps1)[0, 1]:.2f}")
    print(f"\nSplit-half W1 / twice the floor, median, after calibration: HPS1 {(t.split_half_hps1 / (2 * t.floor_net_hps1)).median():.2f}, "
          f"STARLING {(t.split_half_starling / (2 * t.floor_net_starling)).median():.2f} (should be about 1)")
    clear = t[t.w1_baseline_hps1 > t.floor_net_hps1_p95]
    print(f"Proteins whose HPS1 shape is measurably away from the shared shape: {len(clear)}/{len(t)}")
    for who, c in (("network", "w1_net_hps1"), ("network as 2000 samples", "w1_net2000_hps1"), ("STARLING", "w1_starling_hps1")):
        print(f"  closer to HPS1 than the shared shape, {who}: {(clear[c] < clear.w1_baseline_hps1).sum()}/{len(clear)}, median {clear[c].median():.3f} vs {clear.w1_baseline_hps1.median():.3f}")
    cols = ["w1_net_hps1", "w1_net_starling", "w1_starling_hps1", "w1_baseline_hps1"]
    print("Median over families:", t.groupby("family")[cols].median().median().round(3).to_dict())
    print("Proteins with N >= 20:", t[t.N >= MIN_N][cols].median().round(3).to_dict())
    print("Largest family sizes:", t.family.value_counts().head(4).to_dict())
    for c, label in (("cv", "CV"), ("skew", "skew")):
        for who in ("net", "starling"):
            print(f"  Spearman {label} {who} vs HPS1, N >= 30: {stats.spearmanr(t[t.N >= 30][f'{c}_{who}'], t[t.N >= 30][f'{c}_hps1']).statistic:.2f}")
    show = ["name", "N", "w1_net_hps1", "w1_net_starling", "w1_starling_hps1", "w1_baseline_hps1", "floor_net_hps1"]
    print("\nLargest network-HPS1:\n" + t.nlargest(6, "w1_net_hps1")[show].round(2).to_string(index=False))
    print("\nLargest STARLING-HPS1:\n" + t.nlargest(6, "w1_starling_hps1")[show].round(2).to_string(index=False))
    print("\nLargest network-STARLING:\n" + t.nlargest(6, "w1_net_starling")[show].round(2).to_string(index=False))


LABEL = r"shape $W_1$ (% of $\langle R_g\rangle$)"


def panel(ax, t, x, y, xlabel, ylabel, floor_x, floor_y, labels, title, above, below, lim=(0.2, 8)):
    ax.plot(lim, lim, color=MUTED, lw=1, zorder=1)
    ax.axvline(floor_x, color=MUTED, lw=1, ls=(0, (3, 3)), zorder=1)
    ax.axhline(floor_y, color=MUTED, lw=1, ls=(0, (3, 3)), zorder=1)
    inside = t.N >= MIN_N
    ax.scatter(t[x][inside], t[y][inside], s=30, color=POINT, alpha=0.75, edgecolor=SURFACE, linewidth=0.6, zorder=3)
    ax.scatter(t[x][~inside], t[y][~inside], s=30, facecolor=SURFACE, edgecolor=POINT, linewidth=1.4, zorder=3)
    for name, (dx, dy, ha) in labels.items():
        r = t[t.name == name].iloc[0]
        ax.annotate(name, (r[x], r[y]), xytext=(dx, dy), textcoords="offset points", ha=ha, va="center", fontsize=8.5,
                    color=INK)
    ax.text(0.03, 0.97, above, transform=ax.transAxes, ha="left", va="top", fontsize=8.5, color=MUTED)
    ax.text(0.97, 0.03, below, transform=ax.transAxes, ha="right", va="bottom", fontsize=8.5, color=MUTED)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_aspect("equal")
    ticks = [0.2, 0.5, 1, 2, 5]
    ax.set_xticks(ticks, [str(v) for v in ticks])
    ax.set_yticks(ticks, [str(v) for v in ticks])
    ax.minorticks_off()
    ax.grid(color=GRID, lw=0.8, zorder=0)
    ax.set_axisbelow(True)
    ax.set_facecolor(SURFACE)
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_color(GRID)
    ax.tick_params(colors=MUTED, labelsize=9.5, length=0)
    ax.set_xlabel(f"{xlabel}, {LABEL}", fontsize=10, color=INK)
    ax.set_ylabel(f"{ylabel}, {LABEL}", fontsize=10, color=INK)
    ax.set_title(title, fontsize=10.5, color=INK, loc="left", pad=8)


def main():
    a = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    a.add_argument("run", type=pathlib.Path)
    a.add_argument("--trajectories", type=pathlib.Path, default=ROOT / "runs" / "adopt" / "traj")
    args = a.parse_args()

    bench = pd.read_csv(ROOT / "data" / "benchmark.csv")
    family = pd.read_csv(ROOT / "results" / "network_benchmark.csv").set_index("name").family
    bench["family"], bench["N"] = family.reindex(bench.name).to_numpy(), bench.sequence.str.len()
    pred = idprg.predict(list(bench.sequence))
    samples = np.load(args.run / "rg_samples.npz")
    rng = np.random.default_rng(2026)

    rows, q = [], {}
    for r, (_, p) in zip(bench.itertuples(), pred.iterrows()):
        files = [args.trajectories / f"{f}_pooled_timeseries.csv" for f in r.trajectories.split(";")]
        star = samples[r.name] if r.name in samples.files else None
        row, q_net, q_hps1, q_star = score_protein(r, p, files, star, rng)
        rows.append(row)
        q[r.name] = dict(net=q_net, hps1=q_hps1, star=q_star)
    table = pd.DataFrame(rows)

    probe = [np.median([w1(q[n]["net"], q_samples(np.interp(rng.random(max(int(k), 2)), U, q[n]["net"]))) for _ in range(40)])
             for n, k in zip(table.name, table.n_eff_hps1)]
    c = float(np.median(table.split_half_hps1 / (2 * np.array(probe))))
    table["n_w1_hps1"] = table.n_eff_hps1 / c ** 2
    print(f"Calibration: c = {c:.3f}, n_w1 median {table.n_w1_hps1.median():.0f}, n_eff median {table.n_eff_hps1.median():.0f}")
    floors = pd.DataFrame([noise_floors(q[n]["net"], k, rng) for n, k in zip(table.name, table.n_w1_hps1)])
    table = pd.concat([table, floors], axis=1)

    table["w1_net2000_hps1"] = [np.mean([w1(q_samples(np.interp(rng.random(2000), U, q[n]["net"])), q[n]["hps1"]) for _ in range(20)])
                                for n in table.name]

    done = table[table.name.isin(samples.files)]
    for i in done.index:
        name, fam = table.name[i], table.family[i]
        others = done[done.family != fam].name
        shared = np.mean([q[n]["hps1"] for n in others], axis=0)
        table.loc[i, "w1_baseline_hps1"] = w1(shared, q[name]["hps1"])
        table.loc[i, "w1_baseline_starling"] = w1(shared, q[name]["star"])
    table.to_csv(args.run / "shape_w1.csv", index=False)

    t = table[table.name.isin(samples.files)].copy()
    report(t, rng)
    print("\nN below the network's domain:", ", ".join(f"{n} ({k})" for n, k in zip(t.name[t.N < MIN_N], t.N[t.N < MIN_N])))
    make_figure(t, args.run / "shape_w1.png")


def make_figure(t, path):
    NH, NS, SH, BH = "w1_net_hps1", "w1_net_starling", "w1_starling_hps1", "w1_baseline_hps1"
    fN, fS, fB = (t[f"floor_{k}"].median() for k in ("net_hps1", "net_starling", "starling_hps1"))
    fig, axes = plt.subplots(2, 2, figsize=(10.2, 10.6), facecolor=SURFACE)
    panel(axes[0, 0], t, NH, NS, "network vs HPS1", "network vs STARLING", fN, fS,
          {"p27_v78": (-7, 0, "right"), "p27_v56": (0, 10, "center"), "stath_NTD": (7, 3, "left")},
          "a  The network against each reference", "network closer to HPS1", "network closer to STARLING")
    panel(axes[0, 1], t, NH, SH, "network vs HPS1", "STARLING vs HPS1", fN, fB,
          {"p27_v78": (-7, 0, "right"), "p27_v56": (-7, 0, "right")},
          "b  Each predictor against HPS1", "network closer to HPS1", "STARLING closer to HPS1")
    panel(axes[1, 0], t, BH, NH, "shared shape vs HPS1", "network vs HPS1", fN, fN,
          {"p27_v78": (-7, 0, "right"), "p27_v56": (-7, 0, "right"), "stath_NTD": (7, -3, "left"), "RS-peptide": (0, -10, "center")},
          "c  The network against one shape for all", "shared shape closer to HPS1", "network closer to HPS1")
    panel(axes[1, 1], t, BH, SH, "shared shape vs HPS1", "STARLING vs HPS1", fN, fB,
          {"p27_v78": (-7, 0, "right"), "p27_v56": (-7, 0, "right"), "stath_NTD": (-7, -2, "right"), "RS-peptide": (0, 10, "center")},
          "d  STARLING against one shape for all", "shared shape closer to HPS1", "STARLING closer to HPS1")
    fig.text(0.01, 0.005, "One point per protein (121). Solid line: equal distance. Dashed lines: median noise floor for each axis. "
             "Open points: N below 20, outside the network's stated domain.\n"
             "The network and the shared shape are exact curves, while STARLING is 2000 samples, so its W1 carries more sampling noise.",
             fontsize=8.5, color=MUTED)
    fig.tight_layout(rect=(0, 0.035, 1, 1))
    fig.savefig(path, dpi=200, facecolor=SURFACE)


if __name__ == "__main__":
    main()
