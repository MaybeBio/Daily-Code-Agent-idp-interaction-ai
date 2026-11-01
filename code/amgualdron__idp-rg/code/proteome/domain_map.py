"""Does the human disordered proteome fall inside the model's domain? The proteome as a density in descriptor space,
with the validation proteins on top and the outlying chains marked.

    python proteome/domain_map.py [--q 0.99]

Descriptor space: ln N and the nine descriptors of MODEL.md's domain table (scd_pos ... mean_hydro), standardized on the
training proteins; every distinct chain is shown, the extrapolated tier included. A chain's distance to the training set is the Euclidean distance to its nearest training protein
there. The validation proteins (never fitted) calibrate it: proteome chains farther than the validation 99th percentile
are outliers. Panels: the first two principal components of the training proteins, and three descriptor planes (the
closed form's s = s+ + s- and h, <lambda>, N). Writes to proteome/results/:

    domain_map.png
    domain_outliers.csv   the outlying chains, farthest first, with their distance and the descriptors that sit
                          outside the training range
"""
import argparse
import pathlib
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, LogNorm
from scipy.spatial import cKDTree

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))
from idprg import closed_form as cf, data, descriptors  # noqa: E402
from rg_groups import GRID, INK, INK2, chains, style  # noqa: E402

RES = HERE / "results"
DIMS = ("log_seq_length",) + descriptors.COLUMNS
GREYS = LinearSegmentedColormap.from_list("greys", ["#e4e3df", "#a3a29d", "#52514e", "#0b0b0b"])
VAL, OUT = "#2a78d6", "#eb6834"


def prepare(T):
    T = T.assign(s=T.scd_pos + T.scd_neg)
    return T.assign(h=cf.regressors(T, **cf.load()["ce_threshold"]).h)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--q", type=float, default=0.99, help="validation quantile of the distance that marks outliers")
    args = ap.parse_args(argv)

    P = data.proteins()
    tr, va = prepare(P[P.split == "train"]), prepare(P[P.split == "val"])
    D = prepare(chains(True, "none"))           # every chain, the extrapolated tier included
    mu, sd = tr[list(DIMS)].mean(), tr[list(DIMS)].std()
    Z = lambda T: ((T[list(DIMS)] - mu) / sd).to_numpy()
    tree = cKDTree(Z(tr))
    d_val, d_pro = tree.query(Z(va))[0], tree.query(Z(D))[0]
    cut = np.quantile(d_val, args.q)
    D["nn_dist"], D["outlier"] = d_pro, d_pro > cut

    # principal components of the training proteins
    _, _, Vt = np.linalg.svd(Z(tr), full_matrices=False)
    for T in (tr, va, D):
        T[["pc1", "pc2"]] = Z(T) @ Vt[:2].T
    ev = np.linalg.svd(Z(tr), compute_uv=False) ** 2
    ev = ev / ev.sum()

    panels = [("pc1", "pc2", f"PC1 ({ev[0]:.0%})", f"PC2 ({ev[1]:.0%})"),
              ("N", "s", "N (residues)", "s = asinh(5 SCD / N)"),
              ("mean_hydro", "h", "⟨λ⟩", "h (contact energy above threshold)"),
              ("fcr", "shd", "fraction of charged residues", "SHD")]
    fig, axs = plt.subplots(2, 2, figsize=(11, 9.4), facecolor="#fcfcfb",
                            gridspec_kw=dict(hspace=0.3, wspace=0.2, left=0.07, right=0.9, top=0.9, bottom=0.06))
    for ax, (x, y, xl, yl) in zip(axs.flat, panels):
        style(ax)
        ax.grid(axis="x", color=GRID, lw=0.6)
        hb = ax.hexbin(D[x], D[y], gridsize=70, cmap=GREYS, norm=LogNorm(), mincnt=1, lw=0,
                       xscale="log" if x == "N" else "linear")
        o = D[D.outlier]
        ax.scatter(o[x], o[y], s=16, marker="x", color=OUT, lw=0.9,
                   label=f"proteome outliers ({len(o)}, > val. {100 * args.q:g}th pct. distance)")
        ax.scatter(va[x], va[y], s=5, color=VAL, lw=0.4, edgecolor="#fcfcfb", alpha=0.9, label=f"validation proteins ({len(va)})")
        if x == "N":
            ax.set_xticks([20, 50, 100, 200, 500], labels=["20", "50", "100", "200", "500"])
            ax.minorticks_off()
        ax.set_xlabel(xl, fontsize=9, color=INK)
        ax.set_ylabel(yl, fontsize=9, color=INK)
    axs[0, 0].legend(fontsize=7.5, frameon=False, loc="upper left", labelcolor=INK2, markerscale=1.5)
    cb = fig.colorbar(hb, cax=fig.add_axes([0.915, 0.06, 0.012, 0.84]))
    cb.set_label("proteome chains per hexagon", fontsize=9, color=INK)
    cb.ax.tick_params(labelsize=8, colors=INK2)
    cb.outline.set_visible(False)
    fig.suptitle(f"Human disordered proteome ({len(D):,} chains) in the model's descriptor space",
                 x=0.07, ha="left", fontsize=12, color=INK)
    fig.savefig(RES / "domain_map.png", dpi=160)

    cols = ["chain", "kind", "gene", "description", "N", "nn_dist", "tier", "out_of_range", "fcr", "scd_pos",
            "scd_neg", "mean_hydro", "h", "contact_energy", "net_mean_rg", "cf_mean_rg", "sequence"]
    D[D.outlier].sort_values("nn_dist", ascending=False)[cols].to_csv(RES / "domain_outliers.csv", index=False,
                                                                     float_format="%.4g")
    print(f"distance to nearest training protein: validation median {np.median(d_val):.2f}, "
          f"{100 * args.q:g}th pct {cut:.2f}; proteome median {np.median(d_pro):.2f}", file=sys.stderr)
    print(f"proteome within the validation {100 * args.q:g}th pct: {1 - D.outlier.mean():.1%} "
          f"({D.outlier.sum()} outliers); by tier:", file=sys.stderr)
    print(D.groupby("tier").outlier.agg(["size", "sum", "mean"]).to_string(), file=sys.stderr)


if __name__ == "__main__":
    main()
