"""Tail asymmetry P_e/P_c of the human disordered proteome on the plane of charge decoration s and contact energy CE.

    python proteome/tail_asymmetry.py [--frames data/idp_rg_frames.npy]

P_e/P_c depends only on the shape of P(R_g): the scale S_m cancels, and the bounds (idprg.distribution.population_bounds) are
relative to each chain's own mean and follow its length. Panel a is the proteome from the network, the median in each cell.
Panel b is the validation proteins' MD frames on the same axes and colour scale, so the map can be checked. s = s_+ + s_-.

The script also tests what the map assumes: that at fixed s, CE and length, mean hydropathy changes size and not shape.

Writes to proteome/results/:
    tail_asymmetry.{png,pdf}   the figure
    tail_asymmetry.json        every number quoted
"""
import argparse
import json
import pathlib
import sys

import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE.parent / "paper"))

import idprg  # noqa: E402,F401  (float64 in JAX)
from idprg import closed_form as CF, data, distribution as D, network  # noqa: E402
from make_figures import COL2, INK, INK2, MODEL, ORANGE, panel_label, set_style  # noqa: E402
from populations import BASELINE, frame_tails, load_proteome, replicate_groups  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patheffects import withStroke  # noqa: E402
from scipy.ndimage import gaussian_filter  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402

RES = HERE / "results"
S_LIM, CE_LIM = (-0.4, 1.5), (-2.0, 3.0)
RATIO_RANGE = (0.7, 2.0)               # logarithmic; holds 98% of the proteome
RATIO_TICKS = (0.7, 0.8, 0.9, 1, 1.25, 1.5, 2)
MIN_CELL = 3                    # chains in a cell of panel a
SHAPE_COLORS = LinearSegmentedColormap.from_list("asymmetry", ["#123f7a", MODEL, "#cfccc4", ORANGE, "#7d3208"])
NAMED = (("P37840", "α-synuclein"), ("P35637_idr1", "FUS N-term"), ("P06454", "prothymosin-α"), ("P07199_idr2", "CENP-B IDR"), ("Q9UPN9_idr5", "TRIM33 IDR"))
LABEL_AT = {"α-synuclein": (0.07, -0.5), "FUS N-term": (0.09, 0.14), "TRIM33 IDR": (-0.38, 1.45), "prothymosin-α": (1.0, 0.42)}       # CENP-B lies beyond the right edge
REGIONS = dict(blocky_sticky=lambda s, c: (s < -0.1) & (c > 0.2), like_charged_weak=lambda s, c: (s > 0.2) & (c < -1.2), middle=lambda s, c: (np.abs(s) <= 0.1) & (c <= 0.2))


def rho(a, b):
    return float(spearmanr(a, b).statistic)


def within(table, x, y, by, bins):
    """Spearman correlation of x and y after removing the mean of each in cells of `by` (quantile bins)."""
    cell = sum(pd.qcut(table[c].rank(method="first"), bins, labels=False) * bins ** i for i, c in enumerate(by))
    return rho(table[x] - table.groupby(cell)[x].transform("mean"), table[y] - table.groupby(cell)[y].transform("mean"))


def tails(F, table, floor):
    bound = lambda i: (lambda n, c: D.population_bounds(n)[i])
    return frame_tails(F, table.frame_row.to_numpy(), table.N.to_numpy(float), bound(0), bound(1), "own", floor)


def density_levels(x, y, shares):
    """Grid, smoothed density and the density levels enclosing the given shares of the chains."""
    H, xe, ye = np.histogram2d(x, y, bins=(120, 120), range=(S_LIM, CE_LIM))
    H = gaussian_filter(H, 2.0)
    ordered = np.sort(H.ravel())[::-1]
    cum = np.cumsum(ordered) / ordered.sum()
    return 0.5 * (xe[1:] + xe[:-1]), 0.5 * (ye[1:] + ye[:-1]), H.T, sorted(float(ordered[np.searchsorted(cum, q)]) for q in shares)


def figure(Q, V):
    fig, axes = plt.subplots(1, 2, figsize=(COL2, 3.5), layout="constrained", sharey=True, sharex=True)
    norm = TwoSlopeNorm(vmin=np.log(RATIO_RANGE[0]), vcenter=0, vmax=np.log(RATIO_RANGE[1]))
    grid = density_levels(Q.s.to_numpy(), Q.contact_energy.to_numpy(), (0.5, 0.9))
    inside = Q[Q.s.between(*S_LIM) & Q.contact_energy.between(*CE_LIM)]
    cells = axes[0].hexbin(inside.s, inside.contact_energy, C=inside.ln_ratio, reduce_C_function=np.median, gridsize=(52, 26), extent=(*S_LIM, *CE_LIM), mincnt=MIN_CELL,
                           cmap=SHAPE_COLORS, norm=norm, linewidths=0.1, edgecolors="face", rasterized=True, zorder=1)
    order = V.md_ln_ratio.abs().to_numpy().argsort()
    axes[1].scatter(V.s.to_numpy()[order], V.contact_energy.to_numpy()[order], c=V.md_ln_ratio.to_numpy()[order], s=9, cmap=SHAPE_COLORS, norm=norm, linewidths=0.25, edgecolors="white", zorder=2)
    for ax in axes:
        ax.contour(*grid[:3], levels=grid[3], colors=INK2, linewidths=(0.5, 0.8), zorder=3)
        ax.axvline(0, color=INK2, lw=0.5, ls=":", zorder=1)
        ax.set_xlim(*S_LIM), ax.set_ylim(*CE_LIM)
        ax.set_xlabel("charge decoration $s$")
    named = Q.set_index("chain")
    for chain, name in NAMED:
        r = named.loc[chain]
        if name not in LABEL_AT:
            continue
        axes[0].plot(r.s, r.contact_energy, "o", ms=3.6, mfc=INK, mec="white", mew=0.5, zorder=5)
        axes[0].annotate(name, (r.s, r.contact_energy), LABEL_AT[name], fontsize=6, ha="left", va="center", color=INK, zorder=6, path_effects=[withStroke(linewidth=1.6, foreground="white")], arrowprops=dict(arrowstyle="-", color=INK, lw=0.5, shrinkA=0.5, shrinkB=1.5))
    axes[0].set_ylabel("contact energy CE")
    for ax in axes:
        ax.set_xticks([-0.4, 0, 0.4, 0.8, 1.2])
    axes[0].set_title("human proteome, network (median in each cell)", fontsize=7)
    axes[1].set_title(f"{len(V)} validation proteins, MD (one dot each)", fontsize=7)
    panel_label(axes[0], "a", x=-0.12, y=1.03), panel_label(axes[1], "b", x=-0.05, y=1.03)
    axes[1].legend(handles=[Line2D([], [], color=INK2, lw=0.8, label="50% of proteome chains"), Line2D([], [], color=INK2, lw=0.5, label="90% of proteome chains")], loc="upper right", frameon=False, fontsize=6, handlelength=1.4)
    bar = fig.colorbar(cells, ax=axes, shrink=0.9, pad=0.015, extend="both")
    bar.set_ticks([np.log(v) for v in RATIO_TICKS]), bar.set_ticklabels([f"{v:g}" for v in RATIO_TICKS])
    bar.set_label(r"tail asymmetry $P_e/P_c$", fontsize=7), bar.ax.tick_params(labelsize=6.5), bar.outline.set_linewidth(0.4)
    for ext in ("png", "pdf"):
        fig.savefig(RES / f"tail_asymmetry.{ext}", dpi=300)
    plt.close(fig)
    print(f"-> {RES / 'tail_asymmetry'}.png / .pdf", file=sys.stderr)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--frames", default=str(HERE.parent / "data" / "idp_rg_frames.npy"))
    args = ap.parse_args(argv)
    set_style()

    P, F = data.proteins(), data.frames(args.frames)
    floor = 1 / F.shape[1]
    tr, ev = data.split(P, "standard")
    T, V = P.iloc[tr].copy().reset_index(drop=True), P.iloc[ev].copy().reset_index(drop=True)
    for table in (T, V):
        table["s"] = table.scd_pos + table.scd_neg
        table["md_mean"] = np.asarray(F[table.frame_row.to_numpy()], float).mean(1)
        table["md_ln_size"] = np.log(table.md_mean / BASELINE(table.N.to_numpy(float)))
        pc, pe = tails(F, table, floor)
        table["md_ln_ratio"] = np.log(pe / pc)
    Nv = V.N.to_numpy(float)
    for tag, (lnA, lnS, u) in (("net", network.parameters(*network.load(), V)), ("cf", CF.parameters(CF.load(), V))):
        lnA, lnS, u = (np.asarray(a, float) for a in (lnA, lnS, u))
        nu = D.nu_of(u)
        pc, pe = D.populations(lnA, nu, Nv)
        V[f"{tag}_ln_ratio"] = np.log(pe / pc)
        V[f"{tag}_ln_size"] = np.log(D.moments(lnA, lnS, nu)[0] / BASELINE(Nv))

    Q = load_proteome()
    Q["s"] = Q.scd_pos + Q.scd_neg
    Q["ln_ratio"] = np.log(Q.net_pe / Q.net_pc)
    Q["cf_ln_ratio"] = np.log(Q.cf_pe / Q.cf_pc)
    Q["ln_size"] = np.log(Q.net_mean_rg / BASELINE(Q.N.to_numpy(float)))

    # replicate noise in ln(P_e/P_c)
    noise = []
    for g in replicate_groups(P):
        pc, pe = tails(F, g, floor)
        r = np.log(pe / pc)
        noise += [(r[i] - r[j]) ** 2 / 2 for i in range(len(r)) for j in range(i)]

    keys = ["s", "contact_energy", "log_seq_length"]
    hydropathy = dict(note="Spearman correlation with mean hydropathy after removing the mean in cells of s, CE and length", size=dict(), shape=dict())
    hydropathy["size"] = dict(md_training=within(T, "mean_hydro", "md_ln_size", keys, 5), md_validation=within(V, "mean_hydro", "md_ln_size", keys, 3), network_validation=within(V, "mean_hydro", "net_ln_size", keys, 3), closed_form_validation=within(V, "mean_hydro", "cf_ln_size", keys, 3), network_proteome=within(Q, "mean_hydro", "ln_size", keys, 7))
    hydropathy["shape"] = dict(md_training=within(T, "mean_hydro", "md_ln_ratio", keys, 5), md_validation=within(V, "mean_hydro", "md_ln_ratio", keys, 3), network_validation=within(V, "mean_hydro", "net_ln_ratio", keys, 3), closed_form_validation=within(V, "mean_hydro", "cf_ln_ratio", keys, 3), network_proteome=within(Q, "mean_hydro", "ln_ratio", keys, 7), closed_form_proteome=within(Q, "mean_hydro", "cf_ln_ratio", keys, 7))

    regions = {}
    for name, test in REGIONS.items():
        q, v = test(Q.s, Q.contact_energy), test(V.s, V.contact_energy)
        regions[name] = dict(proteome_share=float(q.mean()), proteome_network_median=float(np.exp(Q.ln_ratio[q].median())), validation_proteins=int(v.sum()), validation_md_median=float(np.exp(V.md_ln_ratio[v].median())), validation_network_median=float(np.exp(V.net_ln_ratio[v].median())))
    top = Q[Q.tier == "confident"].nlargest(8, "ln_ratio")
    m = dict(n_chains=len(Q), n_val=len(V), replicate_noise_sd_ln_ratio=float(np.sqrt(np.mean(noise))), replicate_pairs=len(noise),
             ratio_quantiles_proteome={str(q): float(np.exp(Q.ln_ratio.quantile(q))) for q in (.01, .05, .25, .5, .75, .95, .99)},
             ratio_quantiles_validation_md={str(q): float(np.exp(V.md_ln_ratio.quantile(q))) for q in (.01, .05, .25, .5, .75, .95, .99)},
             network_against_md=dict(rho=rho(V.net_ln_ratio, V.md_ln_ratio), sd_network=float(V.net_ln_ratio.std()), sd_md=float(V.md_ln_ratio.std())),
             closed_form_against_md=dict(rho=rho(V.cf_ln_ratio, V.md_ln_ratio), sd_closed_form=float(V.cf_ln_ratio.std())),
             rho_axes=dict(proteome_network=dict(s=rho(Q.s, Q.ln_ratio), CE=rho(Q.contact_energy, Q.ln_ratio), length=rho(Q.log_seq_length, Q.ln_ratio)), validation_md=dict(s=rho(V.s, V.md_ln_ratio), CE=rho(V.contact_energy, V.md_ln_ratio), length=rho(V.log_seq_length, V.md_ln_ratio))),
             rho_with_size=dict(proteome_network=rho(Q.ln_size, Q.ln_ratio), validation_md=rho(V.md_ln_size, V.md_ln_ratio)),
             hydropathy_at_fixed_s_CE_length=hydropathy, regions=regions,
             outside_axes=dict(proteome=int((~(Q.s.between(*S_LIM) & Q.contact_energy.between(*CE_LIM))).sum()), validation=int((~(V.s.between(*S_LIM) & V.contact_energy.between(*CE_LIM))).sum())),
             outside_colors=dict(proteome=int(((Q.ln_ratio < np.log(RATIO_RANGE[0])) | (Q.ln_ratio > np.log(RATIO_RANGE[1]))).sum()), validation=int(((V.md_ln_ratio < np.log(RATIO_RANGE[0])) | (V.md_ln_ratio > np.log(RATIO_RANGE[1]))).sum())),
             named={name: dict(chain=c, s=float(r.s), CE=float(r.contact_energy), ratio=float(np.exp(r.ln_ratio))) for c, name in NAMED for r in [Q.set_index("chain").loc[c]]},
             largest_ratio=[dict(chain=r.chain, gene=r.gene, N=int(r.N), s=float(r.s), CE=float(r.contact_energy), ratio=float(np.exp(r.ln_ratio))) for r in top.itertuples()])
    (RES / "tail_asymmetry.json").write_text(json.dumps(m, indent=1, default=float))
    figure(Q, V)
    print(json.dumps(m, indent=1, default=float))


if __name__ == "__main__":
    main()
