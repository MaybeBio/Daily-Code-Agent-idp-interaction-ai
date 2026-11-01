"""The two chains of equal size that differ most in nu (network), and their P(R_g) overlaid.

    python proteome/nu_pair.py [--tol 0.03] [--min-n 60]

Candidates: one chain per sequence, extrapolated tier and domain outliers (proteome/domain_map.py) excluded. A pair has
equal size when both the network's <R_g> and the normalized size <R_g> / N^nu_F agree within --tol (which also makes N
nearly equal). Picks the pair with the largest |nu_i - nu_j| (always short chains, where the network's nu spreads most),
and the largest among longer chains (N >= --min-n) that are both in the `confident` tier. Writes proteome/results/nu_pair.png and nu_pair.csv (the two pairs, with sequences).
"""
import argparse
import pathlib
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))
from idprg import data, distribution as dist  # noqa: E402
from rg_groups import GRID, INK, INK2, chains, style  # noqa: E402
from size_vs_nu import coloring  # noqa: E402

RES = HERE / "results"
PAIR = ("#c7362f", "#1c5cab")          # low nu, high nu: the red / blue poles of the s colour scale, darker


def widest_pair(D, tol):
    """(i, j) with equal <R_g> and normalized size within tol and the largest nu difference."""
    lr = np.log(D.net_mean_rg.to_numpy())
    ly = lr - data.NU_F * np.log(D.N.to_numpy())
    nu = D.net_nu.to_numpy()
    o, t = np.argsort(lr), np.log1p(tol)
    hi = np.searchsorted(lr[o], lr[o] + t, side="right")
    best = (-1.0, None, None)
    for a, i in enumerate(o):
        j = o[a + 1:hi[a]]
        j = j[np.abs(ly[j] - ly[i]) <= t]
        if len(j):
            k = j[np.argmax(np.abs(nu[j] - nu[i]))]
            if abs(nu[k] - nu[i]) > best[0]:
                best = (abs(nu[k] - nu[i]), i, k)
    _, i, k = best
    return (i, k) if nu[i] < nu[k] else (k, i)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tol", type=float, default=0.03)
    ap.add_argument("--min-n", type=int, default=60, help="second pair: chains at least this long")
    args = ap.parse_args(argv)

    D = chains(False, "none")
    D = D[~D.chain.isin(pd.read_csv(RES / "domain_outliers.csv").chain)].reset_index(drop=True)
    D["s"] = D.scd_pos + D.scd_neg
    C = D[(D.tier == "confident") & (D.N >= args.min_n)].reset_index(drop=True)
    pairs = [("largest Δν", D.loc[list(widest_pair(D, args.tol))]),
             (f"N ≥ {args.min_n}, both confident", C.loc[list(widest_pair(C, args.tol))])]

    fig = plt.figure(figsize=(13, 5.4), facecolor="#fcfcfb")
    gs = fig.add_gridspec(1, 3, width_ratios=[1.15, 1, 1], wspace=0.25, left=0.05, right=0.98, top=0.85, bottom=0.27)
    ax = fig.add_subplot(gs[0])
    style(ax)
    ax.grid(axis="x", color=GRID, lw=0.6)
    norm, cmap, order = coloring(D, "s")
    ax.scatter(D.loc[order, "net_nu"], D.loc[order, "net_mean_rg"] / D.loc[order, "N"] ** data.NU_F,
               c=D.loc[order, "s"], cmap=cmap, norm=norm, s=2, lw=0, alpha=0.5, rasterized=True)
    for k, (title, P) in enumerate(pairs):
        x, y = P.net_nu, P.net_mean_rg / P.N ** data.NU_F
        ax.plot(x, y, color=INK, lw=0.8, zorder=3)
        ax.scatter(x, y, s=40, c=PAIR, edgecolor=INK, lw=0.8, zorder=4, marker="os"[k])
        ax.annotate(str(k + 1), (x.max(), y.iloc[x.argmax()]), xytext=(7, 0), textcoords="offset points",
                    va="center", fontsize=8.5, color=INK, fontweight="bold")
    ax.set_xlabel("predicted ν (network)", fontsize=9, color=INK)
    ax.set_ylabel(f"⟨R_g⟩ / N^ν_F  (nm)", fontsize=9, color=INK)
    ax.set_title("network, coloured by s (red s−, blue s+)", loc="left", fontsize=9.5, color=INK)

    rows = []
    for k, (title, P) in enumerate(pairs):
        a = fig.add_subplot(gs[k + 1])
        style(a)
        r0 = P.net_mean_rg.mean()
        rg = np.linspace(0.45 * r0, 2.0 * r0, 600)
        for (_, c), col in zip(P.iterrows(), PAIR):
            a.plot(rg, dist.pdf(rg, c.net_lnA, np.log(c.net_S_m), c.net_nu), color=col, lw=2,
                   label=f"{c.gene or c.accession} ({c.chain.split('_')[-1]}), N {c.N}, ⟨R_g⟩ {c.net_mean_rg:.3f} nm\n"
                         f"ν {c.net_nu:.3f} · skew {c['net_skew']:.2f} · σ {c.net_std_rg:.3f} nm · s {c.s:+.2f}")
            a.axvline(c.net_mean_rg, color=col, lw=0.8, ls=":")
            rows.append(dict(pair=k + 1, selection=title, **c[["chain", "gene", "description", "N", "tier", "net_mean_rg",
                                                               "net_S_m", "net_std_rg", "net_skew", "net_exkurt",
                                                               "net_nu", "net_lnA", "cf_nu", "s", "h", "fcr",
                                                               "sequence"]].to_dict()))
        a.set_xlim(rg[0], rg[-1])
        a.set_xlabel("R_g (nm)", fontsize=9, color=INK)
        if k == 0:
            a.set_ylabel("P(R_g) (nm⁻¹)", fontsize=9, color=INK)
        dnu = abs(P.net_nu.diff().iloc[-1])
        a.set_title(f"{k + 1}. {title}: Δν = {dnu:.3f}", loc="left", fontsize=9.5, color=INK)
        a.legend(fontsize=7.5, frameon=False, loc="upper left", bbox_to_anchor=(0, -0.16), labelcolor=INK2,
                 labelspacing=0.8)
    fig.suptitle(f"Same size, most different ν: network P(R_g) of chains whose ⟨R_g⟩ and ⟨R_g⟩/N^ν_F agree within "
                 f"{100 * args.tol:g}%", x=0.05, ha="left", fontsize=12, color=INK)
    fig.savefig(RES / "nu_pair.png", dpi=160)
    pd.DataFrame(rows).to_csv(RES / "nu_pair.csv", index=False, float_format="%.5g")
    print(pd.DataFrame(rows).drop(columns=["description", "sequence"]).round(3).to_string(index=False))
    for r in rows:
        print(f"{r['pair']} {r['gene']:>7} {r['sequence']}")


if __name__ == "__main__":
    main()
