"""Figures for proteome/shuffle_search.py (run it first).

    python proteome/shuffle_plot.py

    shuffle_search.png        per candidate: random shuffles in the network's (<R_g>, skewness) plane, native and the
                              four optimized sequences marked
    shuffle_search_pdfs.png   per candidate: network P(R_g) of the native and optimized sequences
"""
import pathlib
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))
from idprg import distribution as dist  # noqa: E402
from rg_groups import GRID, INK, INK2, style  # noqa: E402

RES = HERE / "results"
OBJ = {"rg_max": ("max ⟨R_g⟩", "#2a78d6", "^"), "rg_min": ("min ⟨R_g⟩", "#eb6834", "v"),
       "skew_max": ("max skew", "#1baf7a", "D"), "skew_min": ("min skew", "#eda100", "s")}
GREY = "#a3a29d"


def main():
    S = pd.read_csv(RES / "shuffle_search.csv")
    R = pd.read_parquet(RES / "shuffle_random.parquet")
    chains = list(S[S.kind == "native"].chain)
    nc = len(chains)
    cols = min(3, nc)
    rows = int(np.ceil(nc / cols))

    fig, axs = plt.subplots(rows, cols, figsize=(4.2 * cols, 3.7 * rows + 0.6), facecolor="#fcfcfb", squeeze=False,
                            gridspec_kw=dict(hspace=0.45, wspace=0.28, left=0.07, right=0.98, top=0.88, bottom=0.08))
    for ax, ch in zip(axs.flat, chains):
        style(ax)
        ax.grid(axis="x", color=GRID, lw=0.6)
        g, r = S[S.chain == ch], R[(R.chain == ch) & R.in_domain]
        nat = g[g.kind == "native"].iloc[0]
        ax.scatter(r.net_mean_rg, r.net_skew, s=4, color=GREY, lw=0, alpha=0.6,
                   label=f"random shuffles ({len(r)})")
        for k, (lab, col, mk) in OBJ.items():
            o = g[g.kind == k].iloc[0]
            ax.scatter(o.net_mean_rg, o.net_skew, s=46, marker=mk, color=col, edgecolor="#fcfcfb", lw=0.8, zorder=4,
                       label=lab)
        ax.scatter(nat.net_mean_rg, nat.net_skew, s=110, marker="*", color=INK, edgecolor="#fcfcfb", lw=0.6,
                   zorder=5, label="native")
        ax.set_title(f"{nat.gene} ({ch.split('_')[-1] if '_' in ch else 'protein'}), N {nat.N}, "
                     f"fcr {nat.fcr:.2f}, ⟨λ⟩ {nat.mean_hydro:.2f}", loc="left", fontsize=9, color=INK)
        ax.set_xlabel("⟨R_g⟩ (nm), network", fontsize=8.5, color=INK)
        ax.set_ylabel("skewness, network", fontsize=8.5, color=INK)
    for ax in axs.flat[nc:]:
        ax.axis("off")
    h, l = axs.flat[0].get_legend_handles_labels()
    fig.legend(h, l, loc="upper left", bbox_to_anchor=(0.07, 0.945), ncol=6, frameon=False, fontsize=8.5,
               labelcolor=INK2)
    fig.suptitle("Residue order alone: composition-preserving shuffles of charge- and hydropathy-rich human IDRs",
                 x=0.07, ha="left", fontsize=12, color=INK)
    fig.savefig(RES / "shuffle_search.png", dpi=160)

    fig, axs = plt.subplots(rows, cols, figsize=(4.2 * cols, 3.4 * rows + 0.6), facecolor="#fcfcfb", squeeze=False,
                            gridspec_kw=dict(hspace=0.45, wspace=0.25, left=0.07, right=0.98, top=0.87, bottom=0.08))
    for ax, ch in zip(axs.flat, chains):
        style(ax)
        g = S[S.chain == ch]
        nat = g[g.kind == "native"].iloc[0]
        hi = g.net_mean_rg.max()
        rg = np.linspace(0.3 * g.net_mean_rg.min(), 2.4 * hi, 600)
        ax.plot(rg, dist.pdf(rg, nat.net_lnA, np.log(nat.net_S_m), nat.net_nu), color=INK, lw=1.6, ls="--",
                label="native")
        for k, (lab, col, _) in OBJ.items():
            o = g[g.kind == k].iloc[0]
            ax.plot(rg, dist.pdf(rg, o.net_lnA, np.log(o.net_S_m), o.net_nu), color=col, lw=1.8, label=lab)
        ax.set_xlim(rg[0], rg[-1])
        ax.set_title(f"{nat.gene}, N {nat.N}", loc="left", fontsize=9.5, color=INK)
        ax.set_xlabel("R_g (nm)", fontsize=8.5, color=INK)
        ax.set_ylabel("P(R_g) (nm⁻¹)", fontsize=8.5, color=INK)
    for ax in axs.flat[nc:]:
        ax.axis("off")
    h, l = axs.flat[0].get_legend_handles_labels()
    fig.legend(h, l, loc="upper left", bbox_to_anchor=(0.07, 0.94), ncol=5, frameon=False, fontsize=8.5,
               labelcolor=INK2)
    fig.suptitle("Network P(R_g): native sequence vs composition-preserving optima", x=0.07, ha="left", fontsize=12,
                 color=INK)
    fig.savefig(RES / "shuffle_search_pdfs.png", dpi=160)


if __name__ == "__main__":
    main()
