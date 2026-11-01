"""The paper's figures of the network and its closed form, written to paper/figures/ (PNG at 300 dpi + PDF), with every
number they quote in paper/figure_metrics.json.

    python paper/make_figures.py --frames FRAMES.npy --trajectories DIR      # every figure
    python paper/make_figures.py ... --only fig21                            # just the named figure(s)

--frames: the MD frame store (data/README.md); --trajectories: the benchmark proteins' HPS1 trajectories, one
<name>_pooled_timeseries.csv per entry of data/benchmark.csv's `trajectories` column (column `rg`, nm), needed only
while the predictions are computed; they are cached in paper/cache/ (--refresh recomputes them).

    fig0_intro_example        one validation protein's predicted density, its <Rg> marked
    fig0b_intro_pair          fig0 for two validation proteins with equal predicted <Rg> and different exponent nu
    fig1_val_moments          network vs. MD <Rg>, sigma, skewness, excess kurtosis (validation)
    fig2_val_distributions    single column: the validation W1 distribution + the proteins at its quartiles and maximum
    fig4_real_proteins        W1 against the HPS1 trajectories of the benchmark + four examples
    fig13_experiment_parity   network and HPS1 <Rg> against SAXS, coloured by sequence family
    fig15_surrogate_examples  closed form vs. network vs. MD at the median / 95th percentile / largest KL(network || closed form)
    fig20_surrogate_summary   fig15's examples + the KL distribution + the exponent nu (MD fit, network, closed form)
    fig21_surrogate_response  the closed form's response to each input, one row each: % change of the mean, exponent nu and
                              amplitude A, + that input's example P(Rg) over the reference chain's
    fig23_shape_from_sequence same <Rg>, different distributions: validation pairs with equal MD <Rg>, back to back
    fig24_val_nu              the network's exponent nu vs. the one fitted to each protein's MD frames (validation)
"""
import argparse
import json
import os
import pathlib
import pickle
import re
import sys

import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import idprg  # noqa: E402,F401  (float64 in JAX)
from idprg import closed_form as CF, data, descriptors, distribution as D, metrics, network  # noqa: E402

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.ticker  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
from matplotlib.text import Text  # noqa: E402
from scipy.stats import wasserstein_distance  # noqa: E402

OUT = HERE / "figures"
CACHE = HERE / "cache" / "collected.pkl"

# ── style ─────────────────────────────────────────────────────────────────────
MODEL, SERIES2 = "#2a78d6", "#eb6834"
INK, INK2 = "#0b0b0b", "#52514e"
MOSS, BURGUNDY = "#3ca066", "#D61431"            # MD frames / network
COL1, COL2 = 3.4, 7.0                            # column widths, inches
NET_STYLE = dict(color=BURGUNDY, lw=1.4)
ORANGE, DARK_ORANGE = "#f0913f", "#de5a26"       # fig0b's first protein / fig1's scatter (green is kept for MD)
CF_STYLE = dict(color=INK, lw=1.1, ls=(0, (4, 2)))
FIG20_FONT_SCALE = 1.25
INTRO_ID = 5844                                  # fig0: a strongly skewed validation protein (MD skewness 1.03)
FAMILY_COLORS = (BURGUNDY, MODEL, "#e8a317")
FAMILY_NAMES = {"A1_LCD": "A1-LCD variants", "Tau": "Tau constructs", "PNt": "pertactin variants", "p27": "p27 permutants"}
WGRID = np.linspace(-3.0, 2.0, 10001)            # ln(Rg / MD <Rg>) grid of the exact W1 (as idprg.metrics)


def set_style():
    plt.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 8, "xtick.labelsize": 7, "ytick.labelsize": 7,
        "legend.fontsize": 7, "axes.linewidth": 0.6, "axes.edgecolor": INK2, "axes.labelcolor": INK,
        "xtick.color": INK2, "ytick.color": INK2, "text.color": INK, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "xtick.major.size": 2.5, "ytick.major.size": 2.5, "axes.spines.top": False, "axes.spines.right": False,
        "legend.frameon": False, "mathtext.default": "regular", "pdf.fonttype": 42, "ps.fonttype": 42, "savefig.dpi": 300,
    })


def panel_label(ax, letter, x=-0.02, y=1.04):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=9, fontweight="bold", ha="right", va="bottom", color=INK)


def save(fig, name, tight=True):
    OUT.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"{name}.{ext}", **(dict(bbox_inches="tight", pad_inches=0.03) if tight else {}))
    plt.close(fig)
    print(f"[fig] {OUT / name}.png / .pdf")


def fit_stats(ref, pred):
    ref, pred = np.asarray(ref, float), np.asarray(pred, float)
    ss_res, ss_tot = np.sum((ref - pred) ** 2), np.sum((ref - ref.mean()) ** 2)
    return dict(n=int(len(ref)), r2=float(1 - ss_res / ss_tot), r=float(np.corrcoef(ref, pred)[0, 1]),
                mae=float(np.mean(np.abs(ref - pred))), bias=float(np.mean(pred - ref)))


def display_name(name):
    """Benchmark names are filesystem-safe (`A1_LCD_plus_7R`, `A1_LCD_aro__`); undo that for print."""
    if name.endswith("_plus_"):
        name = name[: -len("_plus_")] + "_plus"
    tail = re.search(r"_+$", name)
    name = re.sub(r"_+$", "", name) + (("−" * len(tail.group(0))) if tail else "")
    name = re.sub(r"_plus$", "+", name)
    name = re.sub(r"_plus_", " +", name)
    return name.replace("_", " ").replace("  ", " ").strip()


# ── two-exponent numerics, one protein p = (ln A, ln S_m, nu) ─────────────────
def pdf(r, p):
    return D.pdf(np.asarray(r, float), *p)


def quantiles(p, n=20_000):
    return D.quantiles((np.arange(n) + 0.5) / n, *p)


def support(p, q=5e-4):
    s = quantiles(p, 4000)
    return s[int(q * 4000)], s[-1 - int(q * 4000)]


def w1_exact(p, frames):
    fr = np.sort(np.asarray(frames, float))
    r = np.exp(np.log(fr.mean()) + WGRID)
    return float(np.trapezoid(np.abs(D.cdf(r, *p) - np.searchsorted(fr, r, "right") / len(fr)), r))


def kl(p, q):
    """KL(p || q) in nats, trapezoid rule in ln Rg."""
    t = p[1] + D.T
    r = np.exp(t)
    lp, lq = D.logpdf(r, *p), D.logpdf(r, *q)
    return float(np.trapezoid(np.exp(lp + t) * (lp - lq), t))


def moments(P):
    m = np.array([D.moments(*p) for p in P])
    return dict(mean=m[:, 0], std=m[:, 1], skew=m[:, 2], exkurt=m[:, 3])


def triples(lnA, lnSm, u):
    return np.column_stack([lnA, lnSm, D.nu_of(u)])


# ── data ──────────────────────────────────────────────────────────────────────
def collect(frames_path, traj_dir, refresh=False):
    """Validation, training and benchmark tables (cached)."""
    if CACHE.exists() and not refresh:
        return pickle.loads(CACHE.read_bytes())
    P, F = data.proteins(), data.frames(frames_path)
    tr, ev = data.split(P, "standard")
    V, T = P.iloc[ev], P.iloc[tr]
    net, spec = network.load()
    coef = CF.load()
    s = metrics.score(*network.parameters(net, spec, V), V, F)
    actual = dict(mean=V.rg_mean.to_numpy(), std=V.rg_std.to_numpy(), skew=V.rg_skew.to_numpy(), exkurt=V.rg_exkurt.to_numpy())
    val = dict(id=V.id.to_numpy(int), length=V.N.to_numpy(int), frame_row=V.frame_row.to_numpy(int),
               net=triples(*network.parameters(net, spec, V)), cf=triples(*CF.parameters(coef, V)), actual=actual)
    val["pred"], val["pred_cf"] = moments(val["net"]), moments(val["cf"])
    val["w1"] = s.w1_pct.to_numpy() / 100 * actual["mean"]
    val["kl_net_cf"] = np.array([kl(p, q) for p, q in zip(val["net"], val["cf"])])
    val["nu_fitted"] = D.nu_of(data.md_fits().set_index("id").u.reindex(val["id"]).to_numpy(float))
    train = dict(id=T.id.to_numpy(int), length=T.N.to_numpy(float), scd_pos=T.scd_pos.to_numpy(float),
                 scd_neg=T.scd_neg.to_numpy(float), mean_hydro=T.mean_hydro.to_numpy(float),
                 contact_energy=T.contact_energy.to_numpy(float))

    B = data.benchmark()
    pb = triples(*network.parameters(net, spec, descriptors.table(B.sequence)))
    real = []
    for b, p in zip(B.itertuples(), pb):
        fr = np.concatenate([pd.read_csv(pathlib.Path(traj_dir) / f"{f}_pooled_timeseries.csv")["rg"].to_numpy(float)
                             for f in b.trajectories.split(";")])
        p = tuple(p)
        real.append(dict(name=b.name, family=b.family, length=int(b.N), p=p, frames=fr, mean=float(D.moments(*p)[0]),
                         sim_mean=float(fr.mean()), sim_std=float(fr.std()), w1=w1_exact(p, fr),
                         nll=-float(np.mean(D.logpdf(fr, *p))), saxs=float(b.saxs_rg_nm)))
    blob = dict(val=val, train=train, real=real)
    CACHE.parent.mkdir(exist_ok=True)
    CACHE.write_bytes(pickle.dumps(blob))
    return blob


class Frames:
    """frames[row]: one protein's stored MD frames (nm), from the memory-mapped store."""
    def __init__(self, path):
        self.F = data.frames(path)

    def __getitem__(self, row):
        return np.asarray(self.F[int(row)], float)


# ── figures ───────────────────────────────────────────────────────────────────
def density_panel(ax, fr, curves, title, note="", legend=False):
    """MD histogram + one or more (p, label, style) curves."""
    lo = min([np.percentile(fr, 0.1)] + [support(p)[0] for p, _, _ in curves])
    hi = max([np.percentile(fr, 99.9)] + [support(p)[1] for p, _, _ in curves])
    grid = np.linspace(max(lo, 1e-3), hi, 400)
    ax.hist(fr, bins=np.linspace(grid[0], grid[-1], 56), density=True, color=MOSS, alpha=0.55, histtype="stepfilled",
            linewidth=0, label="MD frames", zorder=1)
    for z, (p, label, style) in enumerate(curves):
        ax.plot(grid, pdf(grid, p), label=label, zorder=2 + z, **style)
    ax.set_xlim(grid[0], grid[-1]); ax.set_ylim(bottom=0)
    ax.set_title(title, loc="left", pad=3, fontsize=7.5, color=INK)
    if note:
        ax.text(0.97, 0.95, note, transform=ax.transAxes, ha="right", va="top", fontsize=6.5, color=INK2, linespacing=1.3)
    if legend:
        ax.legend(loc="center right", handlelength=1.3, borderaxespad=0.2)


def fig0(val):
    i = int(np.flatnonzero(val["id"] == INTRO_ID)[0])
    p = val["net"][i]
    lo, hi = support(p, 2e-4)
    grid = np.linspace(lo, hi, 500); y = pdf(grid, p)
    mean = float(val["pred"]["mean"][i])
    fig, ax = plt.subplots(figsize=(COL1, COL1 / 1.5))
    ax.fill_between(grid, y, color=BURGUNDY, alpha=0.12, linewidth=0, zorder=1)
    ax.plot(grid, y, color=BURGUNDY, lw=2.4, zorder=3)
    ax.axvline(mean, color=INK, lw=2.0, ls="--", zorder=4)
    ax.annotate(r"$\langle R_g \rangle$", (mean, 0.97), xycoords=("data", "axes fraction"), xytext=(6, 0),
                textcoords="offset points", ha="left", va="top", fontsize=20, color=INK)
    ax.set_xlim(grid[0], grid[-1]); ax.set_ylim(0, 1.15 * y.max())
    ax.set_xlabel(r"$R_g$ [nm]", fontsize=18); ax.set_ylabel(r"$P(R_g)$", fontsize=18)
    ax.set_yticks([]); ax.tick_params(axis="x", labelsize=15, width=1.2, length=4)
    ax.xaxis.set_major_locator(plt.MaxNLocator(4))
    for side in ("left", "bottom"):
        ax.spines[side].set_linewidth(1.2)
    save(fig, "fig0_intro_example")
    a = val["actual"]
    return dict(id=int(val["id"][i]), length=int(val["length"][i]), rel_w1_pct=float(100 * val["w1"][i] / a["mean"][i]),
                **{f"{k}_md": float(a[k][i]) for k in a}, **{f"{k}_pred": float(val["pred"][k][i]) for k in a})


# fig0b's pair, hand-picked: both validation proteins better than the median W1 with |<Rg> error| < 1.5%, network <Rg>
# within 1% and MD <Rg> within 1%, and among such pairs the largest difference in the network's nu (0.73 vs 0.59)
INTRO_PAIR = (1247, 7299)
INTRO_PAIR_COLORS = (ORANGE, MODEL)
INTRO_PAIR_LABEL_AT = (0.75, 0.3)         # each label sits beside its curve's right flank, at this fraction of its peak


def fig0b(val):
    """fig0 for two proteins: equal predicted <Rg>, different exponent nu. Network densities only, each labelled with N
    and nu in its own colour; dotted line = the mean of the two predicted <Rg>."""
    idx = [int(np.flatnonzero(val["id"] == f)[0]) for f in INTRO_PAIR]
    P = [val["net"][i] for i in idx]
    lo = min(support(p, 2e-4)[0] for p in P); hi = max(support(p, 2e-4)[1] for p in P)
    grid = np.linspace(lo, hi, 600)
    fig, ax = plt.subplots(figsize=(3.8, 1.9))   # saved < 2 in tall
    top = 0.0
    for i, p, col, frac in zip(idx, P, INTRO_PAIR_COLORS, INTRO_PAIR_LABEL_AT):
        y = pdf(grid, p); top = max(top, y.max())
        ax.fill_between(grid, y, color=col, alpha=0.15, linewidth=0, zorder=1)
        ax.plot(grid, y, color=col, lw=1.4, zorder=3)
        k = int(np.argmax(y)); k += int(np.flatnonzero(y[k:] <= frac * y[k])[0])
        ax.annotate(rf"$N$ = {int(val['length'][i])},  $\nu$ = {p[2]:.2f}", (grid[k], y[k]), xytext=(5, 0),
                    textcoords="offset points", ha="left", va="center", fontsize=8.5, color=col)
    mean = float(np.mean(val["pred"]["mean"][idx]))
    ax.axvline(mean, color=INK, lw=0.8, ls=(0, (1, 2)), zorder=4)
    ax.annotate(r"$\langle R_g \rangle$", (mean, 0.98), xycoords=("data", "axes fraction"), xytext=(3, 0),
                textcoords="offset points", ha="left", va="top", fontsize=8.5, color=INK)
    ax.set_xlim(grid[0], grid[-1]); ax.set_ylim(0, 1.2 * top)
    ax.set_xlabel(r"$R_g$ [nm]", fontsize=9); ax.set_ylabel(r"$P(R_g)$ [nm$^{-1}$]", fontsize=9)
    ax.tick_params(labelsize=8, width=0.6, length=2.5)
    ax.xaxis.set_major_locator(plt.MaxNLocator(6)); ax.yaxis.set_major_locator(plt.MaxNLocator(4))
    for side in ("left", "bottom"):
        ax.spines[side].set_linewidth(0.6)
    save(fig, "fig0b_intro_pair")
    a, pr = val["actual"], val["pred"]
    rel = lambda v: float(100 * abs(v[idx[0]] - v[idx[1]]) / min(v[idx[0]], v[idx[1]]))
    return dict(ids=list(INTRO_PAIR), length=[int(val["length"][i]) for i in idx], nu_network=[float(p[2]) for p in P],
                mean_line=mean, mean_diff_pct_pred=rel(pr["mean"]), mean_diff_pct_md=rel(a["mean"]),
                rel_w1_pct=[float(100 * val["w1"][i] / a["mean"][i]) for i in idx],
                **{f"{k}_md": [float(a[k][i]) for i in idx] for k in a}, **{f"{k}_pred": [float(pr[k][i]) for i in idx] for k in a})


def parity_panel(ax, ref, pred, label, unit=""):
    ref, pred = np.asarray(ref, float), np.asarray(pred, float)
    st = fit_stats(ref, pred)
    lo, hi = min(ref.min(), pred.min()), max(ref.max(), pred.max())
    pad = 0.05 * (hi - lo)
    lo, hi = lo - pad, hi + pad
    ax.plot([lo, hi], [lo, hi], color=INK2, lw=0.7, ls=(0, (4, 3)), zorder=1)
    ax.scatter(ref, pred, s=5, color=DARK_ORANGE, alpha=0.45, linewidths=0, rasterized=True, zorder=2)
    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
    u = f" [{unit}]" if unit else ""
    ax.set_xlabel(f"MD{u}"); ax.set_ylabel(f"network{u}")
    ax.set_title(label, loc="left", pad=4, fontsize=9.5, color=INK)
    fx, fy = (ref - lo) / (hi - lo), (pred - lo) / (hi - lo)          # stats in the emptier off-diagonal corner
    n_ul, n_lr = np.sum((fx < 0.62) & (fy > 0.62)), np.sum((fx > 0.42) & (fy < 0.36))
    tx, ty, ha, va = (0.05, 0.95, "left", "top") if n_ul <= n_lr else (0.96, 0.05, "right", "bottom")
    ax.text(tx, ty, f"$R^2$ = {st['r2']:.3f}\nMAE = {st['mae']:.3f}{' ' + unit if unit else ''}\nn = {st['n']}",
            transform=ax.transAxes, ha=ha, va=va, fontsize=8, color=INK, linespacing=1.35)
    return st


def fig1(val):
    specs = (("mean", r"$\langle R_g \rangle$", "nm"), ("std", r"$\sigma_{R_g}$", "nm"), ("skew", "skewness", ""),
             ("exkurt", "excess kurtosis", ""))
    out = {}
    with plt.rc_context({"axes.labelsize": 9, "xtick.labelsize": 8, "ytick.labelsize": 8}):
        fig, axes = plt.subplots(2, 2, figsize=(4.4, 3.3), layout="constrained")
        fig.get_layout_engine().set(w_pad=0.04, h_pad=0.04, wspace=0.03, hspace=0.05)
        for ax, letter, (k, label, unit) in zip(axes.ravel(), "abcd", specs):
            out[k] = parity_panel(ax, val["actual"][k], val["pred"][k], label, unit)
            ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(4)); ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(4))
            ax.text(-0.17, 1.05, letter, transform=ax.transAxes, fontsize=11, fontweight="bold", ha="right", va="bottom", color=INK)
        save(fig, "fig1_val_moments", tight=False)
    return out


def fig24(val):
    with plt.rc_context({"axes.labelsize": 9, "xtick.labelsize": 8, "ytick.labelsize": 8}):
        fig, ax = plt.subplots(figsize=(2.6, 2.6), layout="constrained")
        out = parity_panel(ax, val["nu_fitted"], val["net"][:, 2], r"exponent $\nu$")
        ax.set_xlabel("MD (fitted to the frames)")
        ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(5)); ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(5))
        ax.set_aspect("equal")
        save(fig, "fig24_val_nu", tight=False)
    return out


def fig2(val, frames):
    """Single column, as fig4 for the validation split: (a) the W1 distribution, (b-e) the proteins at its 25th, 50th and
    75th percentiles and the largest W1 (MD histogram + network density). The 50th percentile protein is the median."""
    rel = val["w1"] / val["actual"]["mean"]
    order = np.argsort(rel)
    fig = plt.figure(figsize=(COL1, 4.6))
    gs = fig.add_gridspec(3, 2, height_ratios=[0.8, 1.0, 1.0], wspace=0.32, hspace=0.72, top=0.9, bottom=0.08)
    ax_w = fig.add_subplot(gs[0, :])
    ax_w.hist(100 * rel, bins=np.linspace(0, 100 * rel.max() * 1.02, 61), color=MOSS, alpha=0.55, edgecolor="white", linewidth=0.3)
    ax_w.set_xlim(0, 100 * rel.max() * 1.05)
    ax_w.set_xlabel(r"$W_1$ between network and MD [% of MD $\langle R_g \rangle$]")
    ax_w.set_ylabel("proteins")
    med = 100 * float(np.median(rel))
    ax_w.text(0.88, 0.88, f"validation proteins, n = {len(rel)}\nmedian {med:.1f}% (c)", transform=ax_w.transAxes, ha="right",
              va="top", fontsize=6.5, color=INK2, linespacing=1.35)
    panel_label(ax_w, "a", x=-0.1, y=1.1)
    top, picked = ax_w.get_ylim()[1], []
    axes = [fig.add_subplot(gs[r, c]) for r in (1, 2) for c in (0, 1)]
    for ax, pc, letter in zip(axes, (25, 50, 75, 100), "bcde"):
        i = int(order[int(round(pc / 100 * (len(order) - 1)))])
        ax_w.axvline(100 * rel[i], color=BURGUNDY, lw=0.9, zorder=3)
        ax_w.text(100 * rel[i], top * 1.02, letter, ha="center", va="bottom", fontsize=7.5, fontweight="bold", color=INK, clip_on=False)
        rank = f"largest of {len(rel)}" if pc == 100 else f"{pc}th percentile"
        density_panel(ax, frames[val["frame_row"][i]], [(val["net"][i], "network prediction", NET_STYLE)], f"N = {val['length'][i]}",
                      note=f"{rank}\n$W_1$ = {100 * rel[i]:.1f}%")
        ax.set_ylim(0, 1.3 * ax.get_ylim()[1])
        ax.xaxis.set_major_locator(plt.MaxNLocator(4)); ax.yaxis.set_major_locator(plt.MaxNLocator(4))
        if letter in "de":
            ax.set_xlabel(r"$R_g$ [nm]")
        if letter in "bd":
            ax.set_ylabel(r"$P(R_g)$ [nm$^{-1}$]")
        panel_label(ax, letter, x=-0.2)
        picked.append(dict(percentile=pc, id=int(val["id"][i]), length=int(val["length"][i]), w1_nm=float(val["w1"][i]),
                           rel_w1=float(rel[i])))
    fig.legend(*axes[0].get_legend_handles_labels(), loc="lower center", bbox_to_anchor=(0.5, 0.955), ncol=2)
    save(fig, "fig2_val_distributions")
    return dict(examples=picked, n=len(rel), rel_w1_pct={str(q): float(np.percentile(100 * rel, q)) for q in (25, 50, 75, 90, 95, 99, 100)})


def fig4(real):
    fig = plt.figure(figsize=(COL2, 5.3))
    gs = fig.add_gridspec(3, 2, height_ratios=[0.62, 1.0, 1.0], wspace=0.24, hspace=0.62)
    ax_w = fig.add_subplot(gs[0, :])
    rel = np.array([r["w1"] / r["sim_mean"] for r in real])
    order = np.argsort(rel)
    ax_w.hist(100 * rel, bins=np.linspace(0, 100 * rel.max() * 1.02, 41), color=MOSS, alpha=0.55, edgecolor="white", linewidth=0.4)
    med = 100 * float(np.median(rel))
    ax_w.axvline(med, color=INK, lw=0.9, ls=(0, (4, 3)), zorder=4)
    ax_w.set_xlim(0, 100 * rel.max() * 1.05)
    ax_w.set_xlabel(r"$W_1$ between network and MD distributions [% of MD $\langle R_g \rangle$]")
    ax_w.set_ylabel("proteins")
    ax_w.text(0.92, 0.9, f"benchmark proteins, n = {len(real)}\nMD: their HPS1 trajectories\nsolid lines: the proteins in b–e\n"
              f"dashed line: median, {med:.1f}%", transform=ax_w.transAxes, ha="right", va="top", fontsize=7, color=INK2,
              linespacing=1.35, bbox=dict(facecolor="white", edgecolor="none", pad=1.5))
    panel_label(ax_w, "a", x=-0.035, y=1.08)
    top, picked = ax_w.get_ylim()[1], []
    for ax, pc, letter in zip([fig.add_subplot(gs[i, j]) for i in (1, 2) for j in (0, 1)], (10, 50, 90, 100), "bcde"):
        i = int(order[int(round(pc / 100 * (len(order) - 1)))])
        r = real[i]
        rank = f"largest $W_1$ of {len(real)}" if pc == 100 else f"{pc}th percentile of $W_1$"
        ax_w.axvline(100 * rel[i], color=BURGUNDY, lw=1.0, zorder=3)
        ax_w.text(100 * rel[i], top * 1.02, letter, ha="center", va="bottom", fontsize=8, fontweight="bold", color=INK, clip_on=False)
        density_panel(ax, r["frames"], [(r["p"], "network prediction", NET_STYLE)], display_name(r["name"]),
                      note=f"N = {r['length']}\n{rank}\n$W_1$ = {r['w1']:.3f} nm ({100 * rel[i]:.1f}%)", legend=(letter == "b"))
        ax.set_xlabel(r"$R_g$ [nm]")
        if letter in "bd":
            ax.set_ylabel("probability density [nm$^{-1}$]")
        panel_label(ax, letter, x=-0.08)
        picked.append(dict(percentile=pc, name=r["name"], w1_nm=r["w1"], rel_w1_pct=float(100 * rel[i])))
    save(fig, "fig4_real_proteins")
    sim = lambda k: np.array([r[k] for r in real])
    return dict(n=len(real), mean=fit_stats(sim("sim_mean"), sim("mean")), examples=picked,
                rel_w1_pct=dict(median=med, mean=float(100 * rel.mean()), p90=float(100 * np.percentile(rel, 90)),
                                max=float(100 * rel.max())), nll_mean=float(sim("nll").mean()))


def _darken(color, factor=0.62):
    return tuple(factor * c for c in matplotlib.colors.to_rgb(color))


def fig13(real):
    """<Rg> of the network and of the HPS1 simulations against SAXS, coloured by sequence family."""
    exp, net, hps = (np.array([r[k] for r in real]) for k in ("saxs", "mean", "sim_mean"))
    fam = np.array([FAMILY_NAMES.get(r["family"], r["family"]) for r in real])
    pd.DataFrame(dict(protein=[r["name"] for r in real], family=fam, length=[r["length"] for r in real], rg_saxs_nm=exp,
                      rg_network_nm=net, rg_hps1_nm=hps)).to_csv(OUT / "fig13_per_protein.csv", index=False)
    top = list(pd.Series(fam).value_counts().index[:len(FAMILY_COLORS)])
    colour = np.array([FAMILY_COLORS[top.index(x)] if x in top else MOSS for x in fam])
    dark = np.array([_darken(c) for c in colour])
    order = np.argsort([x in top for x in fam], kind="stable")            # named families drawn on top
    fig, ax = plt.subplots(figsize=(0.6 * COL2, 0.6 * COL2))
    hi = 1.06 * np.max(np.concatenate([exp, net, hps]))
    ax.plot([0, hi], [0, hi], color=INK2, lw=0.7, ls=(0, (4, 3)), zorder=1)
    ax.scatter(exp[order], hps[order], s=62, marker="o", c=colour[order], alpha=0.38, linewidths=0, zorder=2)
    ax.scatter(exp[order], net[order], s=19, marker="^", c=dark[order], edgecolors="white", linewidths=0.3, zorder=3)
    st_t, st_s = fit_stats(exp, hps), fit_stats(exp, net)
    shapes = [Line2D([], [], marker="o", ls="", ms=7, color=MOSS, alpha=0.45, mew=0),
              Line2D([], [], marker="^", ls="", ms=4.5, color=_darken(MOSS), mew=0)]
    first = ax.legend(shapes, [f"HPS1 (n = {st_t['n']}):  $R^2$ = {st_t['r2']:.2f},  MAE = {st_t['mae']:.2f} nm",
                               f"network (n = {st_s['n']}):  $R^2$ = {st_s['r2']:.2f},  MAE = {st_s['mae']:.2f} nm"],
                      loc="lower right", handletextpad=0.2, borderaxespad=0.4, fontsize=7)
    ax.add_artist(first)
    swatches = [Line2D([], [], marker="s", ls="", ms=5, color=c, mew=0) for c in FAMILY_COLORS[:len(top)] + (MOSS,)]
    ax.legend(swatches, [f"{x} ({int((fam == x).sum())})" for x in top] + [f"other families ({int(np.sum(~np.isin(fam, top)))})"],
              loc="upper left", handletextpad=0.2, borderaxespad=0.4, fontsize=7)
    ax.set_xlim(0, hi); ax.set_ylim(0, hi); ax.set_aspect("equal", "box")
    ax.set_xlabel(r"SAXS $\langle R_g \rangle$ [nm]"); ax.set_ylabel(r"computed $\langle R_g \rangle$ [nm]")
    save(fig, "fig13_experiment_parity")
    return dict(n=len(exp), network_vs_saxs=st_s, hps1_vs_saxs=st_t, network_vs_hps1=fit_stats(hps, net),
                network_closer_to_saxs=int(np.sum(np.abs(net - exp) < np.abs(hps - exp))))


def fig15(val, frames, axes=None):
    """Median / 95th percentile / largest KL(network || closed form) on validation; `axes` (three) for fig20."""
    k = val["kl_net_cf"]
    order = np.argsort(k)
    picks = (("median KL", order[len(k) // 2]), ("95th-percentile KL", order[int(round(0.95 * (len(k) - 1)))]),
             ("largest KL", order[-1]))
    own = axes is None
    if own:
        fig, axes = plt.subplots(1, 3, figsize=(COL2, 2.35), gridspec_kw=dict(wspace=0.3))
    out = []
    for ax, letter, (rank, i) in zip(axes, "abc", picks):
        fr = frames[val["frame_row"][i]]; pn, pc = val["net"][i], val["cf"][i]
        lo = min(support(pn)[0], support(pc)[0], np.percentile(fr, 0.2))
        hi = max(support(pn)[1], support(pc)[1], np.percentile(fr, 99.8))
        grid = np.linspace(max(lo, 1e-3), hi, 500)
        ax.hist(fr, bins=np.linspace(grid[0], grid[-1], 45), density=True, color=MOSS, alpha=0.35, linewidth=0, zorder=1, label="MD frames")
        ax.plot(grid, pdf(grid, pn), color=BURGUNDY, lw=1.5, zorder=2, label="network")
        ax.plot(grid, pdf(grid, pc), zorder=3, label="closed form", **CF_STYLE)
        ax.set_xlim(grid[0], grid[-1]); ax.set_ylim(0, None); ax.set_ylim(0, 1.8 * ax.get_ylim()[1])
        ax.set_xlabel(r"$R_g$ [nm]")
        ax.set_title(rank, loc="left", pad=3, fontsize=7.5, color=INK)
        ax.legend(loc="upper right", fontsize=6.5, handlelength=1.8, borderaxespad=0.2,
                  title=f"N = {int(val['length'][i])}\nKL = {k[i]:.2g} nats", title_fontsize=6.5, alignment="right")
        panel_label(ax, letter, x=-0.1)
        out.append(dict(rank=rank, id=int(val["id"][i]), length=int(val["length"][i]), kl=float(k[i]),
                        nu_network=float(pn[2]), nu_closed_form=float(pc[2]), mean_rg_network=float(val["pred"]["mean"][i]),
                        mean_rg_closed_form=float(val["pred_cf"]["mean"][i]), mean_rg_md=float(val["actual"]["mean"][i])))
    axes[0].set_ylabel("probability density [nm$^{-1}$]")
    if own:
        save(fig, "fig15_surrogate_examples")
    return dict(n=len(k), examples=out)


def fig20(val, frames):
    fig = plt.figure(figsize=(COL2, 5.6))
    gs = fig.add_gridspec(2, 6, hspace=0.55, wspace=1.6)
    top = [fig.add_subplot(gs[0, 2 * j:2 * j + 2]) for j in range(3)]
    dx, ex = fig.add_subplot(gs[1, 0:3]), fig.add_subplot(gs[1, 3:6])
    fig15(val, frames, axes=top)
    # (d) KL(network || closed form), linear axis to a little past the 95th percentile, tail collected in the last bin
    k = val["kl_net_cf"]
    lim = 0.02 * np.ceil(1.25 * np.percentile(k, 95) / 0.02)
    dx.hist(np.clip(k, None, lim - 1e-12), bins=np.linspace(0.0, lim, 41), color=MOSS, alpha=0.55, linewidth=0)
    dx.set_xlim(0.0, lim)
    step = 0.01 if lim <= 0.06 + 1e-12 else 0.02 * np.ceil(lim / 0.12)
    ticks = np.arange(0.0, lim + 1e-12, step)
    dx.set_xticks(ticks)
    dx.set_xticklabels([("≥" if np.isclose(t, lim) else "") + f"{t:.2f}".replace("0.00", "0") for t in ticks])
    for q, name, ls in ((50, "median", (0, (4, 3))), (95, "95th percentile", (0, (1, 2)))):
        v = np.percentile(k, q)
        dx.axvline(v, ymax=0.7, color=INK, lw=0.9, ls=ls, label=f"{name}: {v:.2g}")
    dx.set_xlabel("KL(network $\\|$ closed form) [nats]"); dx.set_ylabel("proteins")
    dx.set_ylim(0, 1.22 * dx.get_ylim()[1])
    dx.legend(loc="upper right", fontsize=6.5, handlelength=1.6, borderaxespad=0.2)
    # (e) exponent nu: fitted to each protein's MD frames, the network's, the closed form's
    fitted = val["nu_fitted"]
    ok = np.isfinite(fitted)
    series = (fitted[ok], val["net"][ok, 2], val["cf"][ok, 2])
    bins = np.linspace(min(min(s.min() for s in series), 0.40) - 0.005, max(max(s.max() for s in series), 0.62) + 0.005, 41)
    ex.hist(series[0], bins=bins, color=MOSS, alpha=0.55, linewidth=0, zorder=2, label="fitted to each protein's MD frames")
    ex.hist(series[1], bins=bins, histtype="step", color=BURGUNDY, lw=1.1, zorder=3, label="network")
    ex.hist(series[2], bins=bins, histtype="step", zorder=4, label="closed form", **CF_STYLE)
    for v, txt in ((0.5, r"$\theta$ chain "), (0.588, " good\n solvent")):
        ex.axvline(v, color=INK2, lw=0.8, ls=(0, (4, 3)), zorder=1)
        ex.text(v, 0.98, txt, transform=ex.get_xaxis_transform(), ha="right" if v < 0.55 else "left", va="top", fontsize=6.5,
                color=INK2, linespacing=1.2)
    ex.set_xlabel(r"exponent $\nu$"); ex.set_ylabel("proteins")
    ex.set_ylim(0, max(np.histogram(s, bins)[0].max() for s in series) / 0.38)
    ex.legend(loc="upper left", fontsize=6.5, borderaxespad=0.2, bbox_to_anchor=(0.0, 0.78), title=f"validation split, n = {ok.sum()}",
              title_fontsize=6.5, alignment="left", frameon=True, facecolor="white", edgecolor="none", framealpha=1.0)
    ex.yaxis.set_major_locator(plt.MaxNLocator(5))
    panel_label(dx, "d", x=-0.1); panel_label(ex, "e", x=-0.1)
    for t in fig.findobj(Text):
        t.set_fontsize(t.get_fontsize() * FIG20_FONT_SCALE)
    save(fig, "fig20_surrogate_summary")
    pct = lambda v: {str(q): float(np.percentile(v, q)) for q in (5, 25, 50, 75, 95, 99)}
    st = lambda v: dict(min=float(v.min()), median=float(np.median(v)), max=float(v.max()))
    return dict(kl=dict(mean=float(k.mean()), max=float(k.max()), percentiles=pct(k)),
                nu=dict(fitted=st(series[0]), network=st(series[1]), closed_form=st(series[2]),
                        r2_network_vs_fitted=fit_stats(series[0], series[1])["r2"],
                        r2_closed_form_vs_fitted=fit_stats(series[0], series[2])["r2"]))


RESPONSE = (("length", r"$N$", "chain length"), ("scd_neg", r"$s_-$", "opposite-charge blocks"),
            ("scd_pos", r"$s_+$", "net-charge repulsion"), ("contact_energy", "CE", "contact energy"),
            ("mean_hydro", r"$\langle\lambda\rangle$", "mean hydropathy"))
RESPONSE_EXAMPLES = {"scd_neg": 1, "scd_pos": 99, "contact_energy": 99, "mean_hydro": 99}   # right: input at this training pct
# one colour per input; green is left to the MD frames. The reference chain is grey.
RESPONSE_COLORS = {"length": "#9446f0", "scd_neg": "#2a78d6", "scd_pos": "#f5333f", "contact_energy": "#e8772e",
                   "mean_hydro": "#f9bf0f", None: "#8a8984"}
RESPONSE_QUANTITIES = ((r"$\langle R_g\rangle$", "-"), (r"exponent $\nu$", (0, (4, 2))), (r"amplitude $A$", (0, (1, 1.2))))
RESPONSE_RG_MIN = 2.0                    # lower end of the Rg axis on the right, below the lower wall of the most compact chains
MEAN_LINE = dict(color=INK, lw=0.8)
REF_LINE = dict(color=RESPONSE_COLORS[None], lw=0.7, ls=(0, (1, 1.5)))
REF_PDF = dict(color=RESPONSE_COLORS[None], lw=1.2)
RESPONSE_FONT_SCALE = 1.3


def fig21(train):
    """The closed form's response to each input over its training 1st-99th percentile range, in the input's own units, the
    others at the reference chain (training medians of N and <lambda>; CE at the CE threshold c(N), so h = 0 -- along the
    length sweep CE follows c(N); no charge decoration). One row per input (a-e):
    left: change of the mean <Rg>, the exponent nu and the amplitude A relative to the reference chain (circle), in %;
    b-e share their y axis. Right: P(Rg) of the chain at the diamond (the input at its RESPONSE_EXAMPLES percentile), with
    the reference chain's as a grey outline, on one density scale; row a's right holds the legend."""
    coef = CF.load()
    c = coef["ce_threshold"]
    threshold = lambda N: c["c0"] + c["c1"] * np.log(N)

    def at(rows):
        """(ln A, ln S_m, nu) of the closed form for a list of {length, scd_pos, scd_neg, contact_energy, mean_hydro}."""
        t = pd.DataFrame(rows)
        t = t.assign(N=t.length, log_seq_length=np.log(t.length))
        return [tuple(p) for p in triples(*CF.parameters(coef, t))]

    N_ref = float(np.median(train["length"])); cut = float(threshold(N_ref))
    ref = dict(length=N_ref, scd_neg=0.0, scd_pos=0.0, contact_energy=cut, mean_hydro=float(np.median(train["mean_hydro"])))
    vary = lambda key, x: {**ref, key: x, **({"contact_energy": float(threshold(x))} if key == "length" else {})}
    p_ref = at([ref])[0]; mom_ref = D.moments(*p_ref); m_ref = float(mom_ref[0])
    q_ref = np.array([m_ref, p_ref[2], np.exp(p_ref[0])])
    change = lambda p, mo: 100 * (np.array([mo[0], p[2], np.exp(p[0])]) / q_ref - 1)     # % change of <Rg>, nu, A

    fig = plt.figure(figsize=(COL2, 6.2))
    gs = fig.add_gridspec(5, 2, width_ratios=[1.5, 1.0], hspace=1.05, wspace=0.1, left=0.11, right=0.98, top=0.97, bottom=0.08)
    axes = [fig.add_subplot(gs[r, 0]) for r in range(5)]
    rights = [None] + [fig.add_subplot(gs[r, 1]) for r in range(1, 5)]
    out = dict(reference=dict(ref, h=0.0), ce_threshold=c, reference_mean_rg=m_ref, reference_nu=float(p_ref[2]),
               reference_A=float(q_ref[2]), reference_skew=float(mom_ref[2]), inputs={}, examples=[])
    examples = {}
    for k_, (ax, (key, sym, name)) in enumerate(zip(axes, RESPONSE)):
        v, col = train[key], RESPONSE_COLORS[key]
        lo, hi = np.percentile(v, [1, 99])
        xs = np.geomspace(lo, hi, 160) if key == "length" else np.linspace(lo, hi, 160)
        P = at([vary(key, float(x)) for x in xs])
        y = np.array([change(p, D.moments(*p)) for p in P])
        ax.axhline(0, zorder=1, **REF_LINE)
        for k, (_, ls) in enumerate(RESPONSE_QUANTITIES):
            ax.plot(xs, y[:, k], color=col, lw=1.5, ls=ls, zorder=3)
        ax.plot(ref[key], 0, "o", ms=4.5, color=INK, mec="white", mew=0.7, zorder=5, clip_on=False)
        if key in RESPONSE_EXAMPLES:
            pc = RESPONSE_EXAMPLES[key]
            xv = float(np.percentile(v, pc)); pe = at([vary(key, xv)])[0]; me_ = D.moments(*pe)
            ax.plot([xv] * 3, change(pe, me_), "D", ls="", ms=4.5, color=col, mec=INK, mew=0.6, zorder=5, clip_on=False)
            examples[key] = (xv, pe, float(me_[0]))
            out["examples"].append(dict(input=key, percentile=pc, value=xv, mean_rg=float(me_[0]), std_rg=float(me_[1]),
                                        skew=float(me_[2]), nu=float(pe[2]), A=float(np.exp(pe[0])),
                                        mean_change_pct=float(100 * (me_[0] / m_ref - 1))))
        if key == "contact_energy":
            ax.axvline(cut, color=INK2, lw=0.6, ls=(0, (1, 1.5)), zorder=2)
            ax.text(cut, 0.96, " CE threshold", transform=ax.get_xaxis_transform(), ha="left", va="top", fontsize=6.5, color=INK2)
        ax.set_xlim(lo, hi)
        ax.set_xlabel(f"{sym}, {name}", labelpad=1.5)
        panel_label(ax, "abcde"[k_], x=-0.075, y=0.97)
        rows = []
        for pc in (1, 5, 50, 95, 99):
            x = float(np.percentile(v, pc)); p = at([vary(key, x)])[0]; mo = D.moments(*p); q = D.quantiles([0.05, 0.95], *p)
            rows.append(dict(percentile=pc, value=x, mean_rg=float(mo[0]), std_rg=float(mo[1]), skew=float(mo[2]), q05=float(q[0]),
                             q95=float(q[1]), nu=float(p[2]), A=float(np.exp(p[0])), mean_change_pct=float(100 * (mo[0] / m_ref - 1))))
        out["inputs"][key] = rows
    a = axes[0]
    a.set_xscale("log")
    a.xaxis.set_major_locator(matplotlib.ticker.FixedLocator([25, 50, 100, 200, 400]))
    a.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter()); a.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    lim = np.array([ax.get_ylim() for ax in axes[1:]])
    for ax in axes[1:]:
        ax.set_ylim(lim[:, 0].min(), lim[:, 1].max())
    for ax in axes:
        ax.yaxis.set_major_locator(plt.MaxNLocator(3, integer=True))
        ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda t, _: f"{t:+.0f}".replace("+0", "0").replace("-", "−")))
    fig.text(0.015, 0.52, "change relative to the reference chain [%]", rotation=90, ha="left", va="center", fontsize=8)

    # right: each example chain's P(Rg) over the reference chain's, one density scale for all rows
    hi = max(support(p, 1e-3)[1] for p in [p_ref] + [e[1] for e in examples.values()])
    x = np.linspace(RESPONSE_RG_MIN, hi, 600)
    d_ref = pdf(x, p_ref)
    top = 1.1 * max([d_ref.max()] + [pdf(x, p).max() for _, p, _ in examples.values()])
    for rx, (key, sym, _) in zip(rights[1:], RESPONSE[1:]):
        xv, pe, me = examples[key]; d = pdf(x, pe)
        rx.fill_between(x, d, color=RESPONSE_COLORS[key], alpha=0.55, lw=0, zorder=1)
        rx.plot(x, d, color=INK, lw=0.7, zorder=2)
        rx.plot(x, d_ref, zorder=4, **REF_PDF)          # the reference as an outline on top: no mixed fills
        rx.plot([me, me], [0, pdf(np.array([me]), pe)[0]], zorder=5, **MEAN_LINE)
        rx.axvline(m_ref, color=INK, lw=0.7, ls=(0, (4, 3)), zorder=6)
        rx.text(0.98, 0.92, f"{sym} = {xv:.2f}".replace("= -", "= −") + "\n" +
                rf"$\langle R_g\rangle$ {100 * (me / m_ref - 1):+.0f}%".replace("-", "−"), transform=rx.transAxes, ha="right",
                va="top", fontsize=7, color=INK, linespacing=1.3)
        rx.set_xlim(x[0], x[-1]); rx.set_ylim(0, top); rx.set_yticks([])
        rx.spines["left"].set_visible(False)
        if rx is not rights[-1]:
            rx.tick_params(labelbottom=False)
    rights[-1].set_xlabel(rf"$R_g$ at $N = {ref['length']:.0f}$ [nm]", labelpad=1.5)
    rights[1].set_title(r"$P(R_g)$ of the chain at the diamond", fontsize=7.5, pad=3, loc="left")

    from matplotlib.legend_handler import HandlerTuple
    handles = [Line2D([], [], color=INK2, lw=1.5, ls=ls) for _, ls in RESPONSE_QUANTITIES] + [
        (Line2D([], [], **REF_PDF), Line2D([], [], ls="", marker="o", ms=4.5, color=INK, mec="white")),
        Line2D([], [], ls="", marker="D", ms=4.5, color="white", mec=INK, mew=0.6),
        Line2D([], [], color=INK, lw=0.7, ls=(0, (4, 3)))]
    labels = [lab for lab, _ in RESPONSE_QUANTITIES] + ["reference chain", "example chain", "reference mean"]
    lx = fig.add_subplot(gs[0, 1]); lx.axis("off")
    lx.legend(handles, labels, handler_map={tuple: HandlerTuple(ndivide=None)}, loc="center right", bbox_to_anchor=(1.0, 0.5),
              ncol=2, fontsize=6.5, handlelength=1.6, columnspacing=0.8, borderaxespad=0)
    for t in fig.findobj(Text):                  # every label, tick, legend entry and note, scaled up together
        t.set_fontsize(t.get_fontsize() * RESPONSE_FONT_SCALE)
    save(fig, "fig21_surrogate_response")
    return out


PAIR_TOL_MEAN = 0.01      # pairs: MD <Rg> equal within 1%
PAIR_MIN_DSKEW = 0.2      # fig23's examples: MD skewness differs by at least this
PAIR_COLORS = (MODEL, SERIES2)


def mean_matched_pairs(val, frames):
    """Validation pairs whose MD <Rg> agree within PAIR_TOL_MEAN (N free), W1 between their MD distributions (% of <Rg>),
    and whether the network / closed form order their sigma and skewness like MD."""
    a = val["actual"]; m = a["mean"]
    I, J = np.triu_indices(len(m), 1)
    ok = np.abs(m[I] - m[J]) / np.minimum(m[I], m[J]) <= PAIR_TOL_MEAN
    I, J = I[ok], J[ok]
    fr = lambda i: frames[val["frame_row"][i]]
    md_w1 = np.array([100 * wasserstein_distance(fr(i), fr(j)) / (0.5 * (m[i] + m[j])) for i, j in zip(I, J)])
    sgn = lambda v: np.sign(v[I] - v[J])
    right = {f"{who}_{k}": sgn(a[k]) == sgn(pred[k]) for who, pred in (("network", val["pred"]), ("closed_form", val["pred_cf"]))
             for k in ("std", "skew")}
    return dict(i=I, j=J, md_w1_rel=md_w1, d_skew=np.abs(a["skew"][I] - a["skew"][J]),
                std_ratio=np.exp(np.abs(np.log(a["std"][I] / a["std"][J]))), **right)


def fig23(val, frames):
    """Same <Rg>, different shape: among validation pairs with MD <Rg> within 1% (N free) and MD skewness differing by
    >= 0.2, the four whose MD distributions differ most (W1), no protein twice, chosen without looking at the predictions."""
    a, p = val["actual"], val["pred"]
    pr = mean_matched_pairs(val, frames)
    shape = pr["d_skew"] >= PAIR_MIN_DSKEW
    picked, seen = [], set()
    for t in np.flatnonzero(shape)[np.argsort(-pr["md_w1_rel"][shape])]:
        i, j = int(pr["i"][t]), int(pr["j"][t])
        if i in seen or j in seen:
            continue
        if a["skew"][i] > a["skew"][j]:
            i, j = j, i                                      # less skewed protein on top
        picked.append((t, i, j)); seen |= {i, j}
        if len(picked) == 4:
            break
    picked.sort(key=lambda tij: a["mean"][tij[1]])
    fig, axes = plt.subplots(2, 2, figsize=(COL2, 4.6), gridspec_kw=dict(hspace=0.42, wspace=0.12))
    out = dict(tol_mean=PAIR_TOL_MEAN, min_d_skew=PAIR_MIN_DSKEW, n_pairs=int(len(pr["i"])), n_shape_pairs=int(shape.sum()),
               network_skew_order_right=float(pr["network_skew"][shape].mean()), examples=[])
    num = lambda v: f"{v:.2f}".replace("-0.00", "0.00").replace("-", "−")
    for ax, letter, (t, i, j) in zip(axes.ravel(), "abcd", picked):
        fi, fj = frames[val["frame_row"][i]], frames[val["frame_row"][j]]
        lo = min(np.percentile(fi, 0.1), np.percentile(fj, 0.1)); hi = max(np.percentile(fi, 99.9), np.percentile(fj, 99.9))
        x, bins = np.linspace(lo, hi, 400), np.linspace(lo, hi, 50)
        top = 0.0
        for sign, k, f_, col in ((1, i, fi, PAIR_COLORS[0]), (-1, j, fj, PAIR_COLORS[1])):
            h, e = np.histogram(f_, bins=bins, density=True)
            ax.bar(e[:-1], sign * h, width=np.diff(e), align="edge", color=col, alpha=0.3, lw=0)
            y = pdf(x, val["net"][k])
            ax.plot(x, sign * y, color=col, lw=1.5)
            top = max(top, h.max(), y.max())
            ax.text(0.98, 0.97 if sign > 0 else 0.03,
                    f"N = {int(val['length'][k])},  $\\nu$ = {val['net'][k, 2]:.2f} (network)\n"
                    f"$\\sigma$ = {num(a['std'][k])} nm (network {num(p['std'][k])})\n"
                    f"skewness {num(a['skew'][k])} (network {num(p['skew'][k])})",
                    transform=ax.transAxes, ha="right", va="top" if sign > 0 else "bottom", fontsize=6.3, color=col,
                    linespacing=1.35, zorder=6, bbox=dict(facecolor="white", edgecolor="none", alpha=0.75, pad=0.6))
        ax.axhline(0, color=INK2, lw=0.6)
        ax.axvline(0.5 * (a["mean"][i] + a["mean"][j]), color=INK, lw=0.9, ls=(0, (4, 3)), zorder=4)
        ax.set_xlim(lo, hi); ax.set_ylim(-1.8 * top, 1.8 * top)
        ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{abs(v):g}"))
        ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(6, symmetric=True))
        ax.set_title(rf"MD $\langle R_g\rangle$ = {a['mean'][i]:.2f} / {a['mean'][j]:.2f} nm", loc="left", fontsize=7.5, pad=3)
        panel_label(ax, letter, x=-0.04)
        out["examples"].append(dict(ids=[int(val["id"][i]), int(val["id"][j])], length=[int(val["length"][i]), int(val["length"][j])],
                                    md_mean=[float(a["mean"][i]), float(a["mean"][j])], md_w1_between_pct=float(pr["md_w1_rel"][t]),
                                    nu_network=[float(val["net"][i, 2]), float(val["net"][j, 2])],
                                    md_std=[float(a["std"][i]), float(a["std"][j])], pred_std=[float(p["std"][i]), float(p["std"][j])],
                                    md_skew=[float(a["skew"][i]), float(a["skew"][j])], pred_skew=[float(p["skew"][i]), float(p["skew"][j])]))
    for ax in axes[1]:
        ax.set_xlabel(r"$R_g$ [nm]")
    for ax in axes[:, 0]:
        ax.set_ylabel("probability density [nm$^{-1}$]")
    fig.legend(handles=[Patch(color=INK2, alpha=0.3, lw=0, label="MD frames (histogram)"),
                        Line2D([], [], color=INK2, lw=1.5, label="network prediction"),
                        Line2D([], [], color=INK, lw=0.9, ls=(0, (4, 3)), label=r"MD $\langle R_g\rangle$ of the pair")],
               loc="lower center", bbox_to_anchor=(0.5, 0.96), ncol=3, fontsize=6.5)
    save(fig, "fig23_shape_from_sequence")
    for key, thr, fmt in (("d_skew", (0.0, 0.2, 0.4), "skew_diff_ge_{:g}"), ("std_ratio", (1.0, 1.1, 1.2, 1.3), "std_ratio_ge_{:g}")):
        k = "skew" if key == "d_skew" else "std"
        for th in thr:
            s = pr[key] >= th
            out[fmt.format(th)] = dict(n=int(s.sum()), network=float(pr[f"network_{k}"][s].mean()),
                                       closed_form=float(pr[f"closed_form_{k}"][s].mean()))
    out["md_w1_between_pct"] = dict(median=float(np.median(pr["md_w1_rel"])), max=float(pr["md_w1_rel"].max()))
    return out


FIGS = ("fig0", "fig0b", "fig1", "fig2", "fig4", "fig13", "fig15", "fig20", "fig21", "fig23", "fig24")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--frames", default=os.environ.get(data.FRAMES_ENV), help="MD frame store (.npy)")
    ap.add_argument("--trajectories", help="directory of the benchmark trajectories (<name>_pooled_timeseries.csv)")
    ap.add_argument("--only", nargs="+", choices=FIGS, default=None)
    ap.add_argument("--refresh", action="store_true", help="recompute the cached predictions")
    args = ap.parse_args()
    if (not CACHE.exists() or args.refresh) and not args.trajectories:
        sys.exit("computing the predictions needs --trajectories (the benchmark's HPS1 trajectories)")
    want = lambda f: args.only is None or f in args.only
    set_style()
    d = collect(args.frames, args.trajectories, args.refresh)
    val, train, real, frames = d["val"], d["train"], d["real"], Frames(args.frames)
    mpath = HERE / "figure_metrics.json"
    m = json.loads(mpath.read_text()) if mpath.exists() else {}
    m.update(n_val=len(val["id"]), n_train=len(train["id"]), n_benchmark=len(real))
    if want("fig0"): m["fig0"] = fig0(val)
    if want("fig0b"): m["fig0b"] = fig0b(val)
    if want("fig1"): m["fig1"] = fig1(val)
    if want("fig2"): m["fig2"] = fig2(val, frames)
    if want("fig4"): m["fig4"] = fig4(real)
    if want("fig13"): m["fig13"] = fig13(real)
    if want("fig15"): m["fig15"] = fig15(val, frames)
    if want("fig20"): m["fig20"] = fig20(val, frames)
    if want("fig21"): m["fig21"] = fig21(train)
    if want("fig23"): m["fig23"] = fig23(val, frames)
    if want("fig24"): m["fig24"] = fig24(val)
    mpath.write_text(json.dumps(m, indent=1, default=float))
    print(f"[fig] metrics -> {mpath}")


if __name__ == "__main__":
    main()
