"""How much the shape of P(R_g) varies among human disordered chains with the same mean R_g.

    python proteome/rg_groups.py [--centers 1.5 2 2.5 3 3.5 4 5 6] [--min-chains 50] [--keep-extrapolated]
                                    [--model cf|net]

Chains from proteome/results/predictions.parquet, one per distinct sequence. For each centre R0 the group is every chain
with |<R_g>/R0 - 1| <= 1%, widened to 2% then 3% if it has fewer than --min-chains. Writes to proteome/results/:

--model picks the predictor (closed form, cf, or network, net); groups use that predictor's own <R_g>.

    rg_groups[_net].csv   one row per group: window, size, N range and quantiles of nu, skewness and std/mean
    rg_groups[_net].png   nu within each group, and P(R_g) of the chains at the 5th / 50th / 95th percentile of nu
"""
import argparse
import pathlib
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import idprg  # noqa: E402,F401
from idprg import distribution as dist  # noqa: E402

RES = HERE / "results"
WINDOWS = (0.01, 0.02, 0.03)
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
BLUE = "#2a78d6"
NU_SHADES = {"p5": "#86b6ef", "p50": "#2a78d6", "p95": "#0d366b"}   # one hue, light -> dark with nu


MODEL = {"cf": "closed form", "net": "network"}


def chains(keep_extrapolated, model):
    """One row per distinct sequence, with the chosen predictor's columns unprefixed (mean_rg, nu, ...)."""
    D = pd.read_parquet(RES / "predictions.parquet")
    D = D.rename(columns={c: c[len(model) + 1:] for c in D if c.startswith(model + "_")})
    D = D[~(D.kind.eq("idr") & D.whole_protein)].drop_duplicates("sequence")
    return D if keep_extrapolated else D[D.tier != "extrapolated"]


def groups(D, centers, min_chains):
    out = []
    for r0 in centers:
        rel = (D.mean_rg / r0 - 1).abs()
        for w in WINDOWS:
            g = D[rel <= w]
            if len(g) >= min_chains:
                break
        out.append((r0, w, g))
    return out


def summary(gs):
    rows = []
    for r0, w, g in gs:
        q = lambda c, p: g[c].quantile(p)
        cv = g.std_rg / g.mean_rg
        rows.append(dict(rg_center=r0, window_pct=100 * w, n=len(g), N_min=g.N.min(), N_median=g.N.median(),
                         N_max=g.N.max(), nu_p5=q("nu", .05), nu_p25=q("nu", .25), nu_median=q("nu", .5),
                         nu_p75=q("nu", .75), nu_p95=q("nu", .95), nu_sd=g.nu.std(),
                         skew_p5=q("skew", .05), skew_p95=q("skew", .95),
                         cv_p5=cv.quantile(.05), cv_p95=cv.quantile(.95),
                         corr_nu_lnN=np.corrcoef(g.nu, np.log(g.N))[0, 1]))
    return pd.DataFrame(rows)


def representative(g, p):
    """The chain whose nu is closest to the p-quantile of the group."""
    return g.iloc[(g.nu - g.nu.quantile(p)).abs().argmin()]


def style(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(INK2)
    ax.tick_params(colors=INK2, labelsize=8)
    ax.grid(axis="y", color=GRID, lw=0.6)
    ax.set_axisbelow(True)


def figure(gs, S, path, detail, model):
    fig = plt.figure(figsize=(11, 7.6), facecolor="#fcfcfb")
    top = fig.add_gridspec(2, 1, height_ratios=[1.15, 1], hspace=0.42, left=0.07, right=0.98, top=0.9, bottom=0.08)
    ax = fig.add_subplot(top[0])
    style(ax)
    pos = np.arange(len(gs))
    v = ax.violinplot([g.nu for _, _, g in gs], pos, widths=0.8, showextrema=False)
    for b in v["bodies"]:
        b.set_facecolor(BLUE), b.set_edgecolor("none"), b.set_alpha(0.35)
    for i, (_, _, g) in enumerate(gs):
        lo, q1, md, q3, hi = g.nu.quantile([.05, .25, .5, .75, .95])
        ax.plot([i, i], [lo, hi], color=BLUE, lw=1.2, solid_capstyle="round")
        ax.plot([i, i], [q1, q3], color=BLUE, lw=5, solid_capstyle="round")
        ax.plot(i, md, "o", ms=5, mfc="#fcfcfb", mec=BLUE, mew=1.5)
    ax.set_xticks(pos, [f"{r0:g} nm ±{100 * w:.0f}%\nn={len(g)}, N {g.N.min()}–{g.N.max()}"
                        for r0, w, g in gs], fontsize=8, color=INK2)
    ax.set_ylabel(f"ν ({MODEL[model]})", color=INK, fontsize=9)
    ax.set_title("ν among chains with the same ⟨R_g⟩   (bar: 5–95%, thick: IQR, dot: median)",
                 loc="left", fontsize=10, color=INK)

    bot = top[1].subgridspec(1, len(detail), wspace=0.18)
    byc = {r0: g for r0, _, g in gs}
    for k, r0 in enumerate(detail):
        a = fig.add_subplot(bot[k])
        style(a)
        g = byc[r0]
        rg = np.linspace(0.4 * r0, 2.3 * r0, 600)
        for tag, p in (("p5", .05), ("p50", .5), ("p95", .95)):
            c = representative(g, p)
            y = dist.pdf(rg, c.lnA, np.log(c.S_m), c.nu)
            a.plot(rg, y, color=NU_SHADES[tag], lw=2,
                   label=f"ν {c.nu:.3f} · N {c.N} · skew {c['skew']:.2f}\n{c.gene or c.accession}")
        a.axvline(r0, color=INK2, lw=0.8, ls=":")
        a.set_xlabel("R_g (nm)", fontsize=9, color=INK)
        a.set_title(f"⟨R_g⟩ ≈ {r0:g} nm: chains at ν p5 / p50 / p95", loc="left", fontsize=9.5, color=INK)
        a.legend(fontsize=7, frameon=False, loc="upper right", labelcolor=INK2, labelspacing=0.9)
        a.set_xlim(rg[0], rg[-1])
        if k == 0:
            a.set_ylabel("P(R_g) (nm⁻¹)", fontsize=9, color=INK)
    fig.suptitle(f"Shape of P(R_g) at fixed mean size — human disordered chains (idprg {MODEL[model]})",
                 x=0.07, ha="left", fontsize=12, color=INK)
    fig.savefig(path, dpi=160)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--centers", type=float, nargs="+", default=[1.5, 2, 2.5, 3, 3.5, 4, 5, 6])
    ap.add_argument("--detail", type=float, nargs="+", default=[2, 3, 5], help="centres with P(R_g) panels")
    ap.add_argument("--min-chains", type=int, default=50)
    ap.add_argument("--keep-extrapolated", action="store_true")
    ap.add_argument("--model", choices=MODEL, default="cf")
    args = ap.parse_args(argv)

    D = chains(args.keep_extrapolated, args.model)
    tag = "" if args.model == "cf" else "_" + args.model
    gs = groups(D, args.centers, args.min_chains)
    S = summary(gs)
    S.to_csv(RES / f"rg_groups{tag}.csv", index=False, float_format="%.4g")
    figure(gs, S, RES / f"rg_groups{tag}.png", args.detail, args.model)
    print(f"{len(D)} distinct chains", file=sys.stderr)
    print(S.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
