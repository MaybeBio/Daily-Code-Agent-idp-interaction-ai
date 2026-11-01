"""Normalized size <R_g> / N^nu_F against the predicted nu, for every disordered chain of the human proteome.

    python proteome/size_vs_nu.py [--color N s scd_neg scd_pos h] [--keep-extrapolated]

nu_F = 0.57874 is the Flory fit ln <R_g> = ln A_F + nu_F ln N of the training proteins (idprg.data), so a chain on the
Flory baseline sits at A_F = 0.209 nm. One panel per predictor (closed form, network), each with its own <R_g> and nu;
one point per distinct sequence. One figure per colouring (--color): N, or one of the closed form's descriptors s-, s+,
the signed s = s+ + s- (diverging red/blue, each arm scaled to its own 1st/99th percentile), h (colour scale clipped at the 1st/99th percentile; larger magnitudes darker and drawn on top). Writes
proteome/results/size_vs_nu.png (N) and size_vs_nu_<descriptor>.png.
"""
import argparse
import pathlib
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, LogNorm, Normalize, TwoSlopeNorm

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from idprg import data  # noqa: E402
from rg_groups import GRID, INK, INK2, chains, style  # noqa: E402

RES = HERE / "results"
MODEL = {"cf": "closed form", "net": "network"}
BLUES = LinearSegmentedColormap.from_list("blues", ["#b7d3f6", "#5598e7", "#256abf", "#0d366b"])
# diverging, for the signed s: red (s-, opposite charges segregated) <- grey -> blue (s+, like charges)
CHARGE = LinearSegmentedColormap.from_list("charge", ["#7a1a1a", "#e34948", "#f2a6a0", "#d6d5d0",
                                                      "#86b6ef", "#2a78d6", "#0d366b"])
COLOR = {"N": "N (residues)", "s": "s = s+ + s−\n← segregated ± blocks · like charges →", "scd_neg": "s− = min(asinh(5 SCD / N), 0)", "scd_pos": "s+ = max(asinh(5 SCD / N), 0)",
         "h": "h (contact energy above threshold)"}


def coloring(D, c):
    """(norm, colormap, draw order) for colouring by column c: larger magnitudes darker and on top."""
    v = D[c]
    if c == "s":                                   # each arm on its own range: s- reaches -0.2, s+ 0.8
        return (TwoSlopeNorm(0, v.quantile(0.01), v.quantile(0.99)), CHARGE,
                (v / np.where(v < 0, -v.quantile(0.01), v.quantile(0.99))).abs().sort_values().index)
    if c == "N":
        return LogNorm(v.min(), v.max()), BLUES, v.sort_values().index
    if c == "scd_neg":
        return Normalize(v.quantile(0.01), 0), BLUES.reversed(), v.sort_values(ascending=False).index
    return Normalize(0, v.quantile(0.99)), BLUES, v.sort_values().index


def figure(D, c):
    a_f = np.exp(data.LN_A_F)
    fig, axs = plt.subplots(1, 2, figsize=(11, 5), sharex=True, sharey=True, facecolor="#fcfcfb",
                            gridspec_kw=dict(wspace=0.08, left=0.07, right=0.9, top=0.84, bottom=0.11))
    norm, cmap, order = coloring(D, c)
    for ax, m in zip(axs, MODEL):
        style(ax)
        ax.grid(axis="x", color=GRID, lw=0.6)
        x, y = D.loc[order, f"{m}_nu"], D.loc[order, f"{m}_mean_rg"] / D.loc[order, "N"] ** data.NU_F
        sc = ax.scatter(x, y, c=D.loc[order, c], cmap=cmap, norm=norm, s=3, lw=0, alpha=0.6, rasterized=True)
        ax.axhline(a_f, color=INK2, lw=0.9, ls="--")
        r = np.corrcoef(x, np.log(y))[0, 1]
        ax.set_title(f"{MODEL[m]}   (r = {r:.2f} with ln size, n = {len(D):,})", loc="left", fontsize=10, color=INK)
        ax.set_xlabel("predicted ν", fontsize=9, color=INK)
        if c == "N":
            print(f"{m}: corr(nu, ln size) {r:.3f}; nu 1-99% {x.quantile(.01):.3f}-{x.quantile(.99):.3f}; "
                  f"size/A_F 1-99% {y.quantile(.01) / a_f:.3f}-{y.quantile(.99) / a_f:.3f}", file=sys.stderr)
    axs[0].set_ylabel(f"⟨R_g⟩ / N^ν_F  (nm),  ν_F = {data.NU_F:.3f}", fontsize=9, color=INK)
    axs[0].annotate("Flory baseline (training fit)", (axs[0].get_xlim()[0], a_f), xytext=(4, 4),
                    textcoords="offset points", fontsize=7.5, color=INK2)
    cb = fig.colorbar(sc, cax=fig.add_axes([0.915, 0.11, 0.012, 0.73]))
    if c == "N":
        cb.set_ticks([20, 50, 100, 200, 500], labels=["20", "50", "100", "200", "500"])
        cb.ax.minorticks_off()
    if c == "s":
        cb.set_ticks([t for t in (-0.2, -0.1, 0, 0.2, 0.4, 0.6, 0.8) if norm.vmin <= t <= norm.vmax])
    cb.set_label(COLOR[c], fontsize=9, color=INK)
    cb.ax.tick_params(labelsize=8, colors=INK2)
    cb.outline.set_visible(False)
    cb.solids.set_alpha(1)
    fig.suptitle(f"Normalized size vs predicted ν, coloured by {COLOR[c].split(' =')[0].split(' (')[0]} — "
                 "human disordered chains (idprg)", x=0.07, ha="left",
                 fontsize=12, color=INK)
    out = RES / ("size_vs_nu.png" if c == "N" else f"size_vs_nu_{c}.png")
    fig.savefig(out, dpi=160)
    plt.close(fig)
    print(f"-> {out}", file=sys.stderr)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--color", nargs="+", choices=COLOR, default=list(COLOR))
    ap.add_argument("--keep-extrapolated", action="store_true")
    args = ap.parse_args(argv)
    D = chains(args.keep_extrapolated, "none")
    D["s"] = D.scd_pos + D.scd_neg
    for c in args.color:
        figure(D, c)


if __name__ == "__main__":
    main()
