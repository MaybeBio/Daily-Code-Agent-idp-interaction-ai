"""Compact and extended populations of the human disordered proteome: size trend, equal-mean pairs, and sequence logic.

    python proteome/compact_states.py [--frames data/idp_rg_frames.npy] [--tol 0.01] [--keep-extrapolated]

Compact and extended probabilities are those of fig26, each chain against its own mean:
    P_c = P(Rg < 0.8 <Rg>),   P_e = P(Rg > 1.3 <Rg>).
Both depend only on the shape (A, nu), not on S_m. Chains are one per distinct sequence, without the extrapolated tier
(rg_groups.chains). The validation split supplies the calibration: network, closed form and MD frames, the same way.

Equal-mean pairs have <Rg> within --tol, each set paired on its own <Rg>. The replicate curve holds the proteins simulated
more than once (data/proteins.csv has duplicate sequences), so any fold difference there is sampling noise.

The residual is ln P_c minus the running median of ln P_c against ln N. The closed form attributes it by exact Shapley
values over its four sequence inputs (s+, s-, <lambda>, h), against a reference chain of the same N holding the proteome's
median inputs. Each term is then tested against MD: validation MD residual ~ sum of beta_t x term_t (slope 1 = calibrated,
0 = no signal, negative = wrong direction). Intervals are cluster bootstraps.

Writes to proteome/results/:
    compact_states_size_pairs.{png,pdf}   a, b P_c and P_e against N; c fold difference in equal-mean pairs; d free energy
    compact_states_sequence.{png,pdf}     a, b (s, h) maps; c closed-form terms of the outliers; d those terms against MD
    compact_states.json                   every number quoted
    compact_states_pairs.csv              fold-difference statistics of every pair set and threshold
    compact_states_chains.csv             one row per chain: P_c, P_e, residuals, Shapley terms
"""
import argparse
import itertools
import json
import math
import pathlib
import sys

import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "paper"))
import idprg  # noqa: E402,F401  (float64 in JAX)
from idprg import closed_form as CF, data, distribution as D, network  # noqa: E402
from make_figures import BURGUNDY, CF_STYLE, COL2, INK, INK2, MODEL, MOSS, ORANGE, panel_label, set_style  # noqa: E402
from rg_groups import chains as proteome_chains  # noqa: E402

import matplotlib  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402

RES = HERE / "results"
TC, TE = 0.8, 1.3                    # fig26's thresholds
KT_KCAL = 0.0019872 * 298            # kcal/mol at the simulations' 298 K
TERMS = ("s+", "s-", "lambda", "h")
LABEL = {"s+": r"$s_+$", "s-": r"$s_-$", "lambda": r"$\langle\lambda\rangle$", "h": r"$h$"}
INPUT = {"s+": "scd_pos", "s-": "scd_neg", "lambda": "mean_hydro", "h": "h"}
LIGHT = "#d9d8d3"
DIVERGING = LinearSegmentedColormap.from_list("resid", [MODEL, "#e8e7e3", BURGUNDY])
DENSITY = LinearSegmentedColormap.from_list("density", ["#f1f0ec", "#a3a29d", "#52514e"])


# ── probabilities ─────────────────────────────────────────────────────────────
def cdf_at(x, lnA, lnSm, nu, block=2000):
    """P(Rg <= x) per chain; D.cdf takes one chain at a time and this runs over the whole proteome."""
    x, lnA, lnSm, nu = [np.asarray(a, float) for a in (x, lnA, lnSm, nu)]
    out, T, dT = np.empty(len(x)), D.T, D.T[1] - D.T[0]
    for b in range(0, len(x), block):
        i = slice(b, b + block)
        _, w = D.grid(lnA[i], lnSm[i], nu[i])
        c = np.concatenate([np.zeros((w.shape[0], 1)), np.cumsum(0.5 * (w[:, 1:] + w[:, :-1]) * dT, axis=1)], axis=1)
        c /= c[:, -1:]
        t = np.log(x[i]) - lnSm[i]
        k = np.clip(((t - T[0]) / dT).astype(int), 0, len(T) - 2)
        f, n = (t - T[k]) / dT, np.arange(len(k))
        out[i] = c[n, k] * (1 - f) + c[n, k + 1] * f
    return out


def mean_of(lnA, lnSm, nu, block=2000):
    return np.concatenate([D.moments(lnA[b:b + block], lnSm[b:b + block], nu[b:b + block])[0]
                           for b in range(0, len(lnA), block)])


def tails(lnA, lnSm, nu, tc=TC, te=TE):
    m = mean_of(lnA, lnSm, nu)
    return m, cdf_at(tc * m, lnA, lnSm, nu), 1 - cdf_at(te * m, lnA, lnSm, nu)


# ── size trend, pairs ─────────────────────────────────────────────────────────
def trend(N, y, bins):
    """Running median of y against ln N, on quantile bins of N (ties at short N merge)."""
    e = np.unique(np.quantile(N, np.linspace(0, 1, bins + 1)))
    x, m = [], []
    for a, b, last in zip(e[:-1], e[1:], [False] * (len(e) - 2) + [True]):
        k = (N >= a) & ((N <= b) if last else (N < b))
        if k.any():
            x.append(np.median(np.log(N[k]))), m.append(np.median(y[k]))
    x, m = np.array(x), np.array(m)
    return (lambda n: np.interp(np.log(n), x, m)), x, m


def pairs(rg, tol):
    """All pairs (i, j) with <Rg> within tol."""
    o = np.argsort(rg)
    lr = np.log(rg[o])
    hi = np.searchsorted(lr, lr + np.log1p(tol), side="right")
    I = [np.full(h - i - 1, o[i]) for i, h in enumerate(hi)]
    J = [o[i + 1:h] for i, h in enumerate(hi)]
    return np.concatenate(I), np.concatenate(J)


def fold_stats(f):
    """Statistics of ln fold differences, reported as folds."""
    return dict(pairs=len(f), median=float(np.exp(np.median(f))), p90=float(np.exp(np.quantile(f, .9))),
                p99=float(np.exp(np.quantile(f, .99))), max=float(np.exp(f.max())),
                frac_ge_1p25=float((f >= np.log(1.25)).mean()), frac_ge_1p5=float((f >= np.log(1.5)).mean()),
                frac_ge_2=float((f >= np.log(2)).mean()))


# ── closed-form attribution ───────────────────────────────────────────────────
def cf_shape(coef, N, sp, sn, lam, h):
    """(ln A, nu) of the closed form: P_c and P_e need no ln S_m, so its coefficients play no part."""
    c = coef["u"]["coef"]
    u = (coef["u"]["intercept"] + c["scd_pos"] * sp + c["scd_pos/N"] * sp / N + c["scd_neg"] * sn
         + c["mean_hydro"] * lam + c["h"] * h)
    nu, a = D.nu_of(u), coef["ln_A"]
    return a["1"] + a["nu"] * nu + a["1/N"] / N + a["mean_hydro"] * lam, nu


def shapley(coef, T, ref):
    """Exact Shapley values of ln P_c over the four inputs (s+/N travels with s+). Returns ({term: phi}, ln P_c of the
    reference chain). `ref` holds each term's reference value per row. The values sum to ln P_c(x) - ln P_c(ref)."""
    N, n = T.N.to_numpy(float), len(TERMS)
    value = {}
    for r in range(n + 1):
        for on in itertools.combinations(range(n), r):
            x = {t: (T[INPUT[t]].to_numpy(float) if i in on else ref[t]) for i, t in enumerate(TERMS)}
            lnA, nu = cf_shape(coef, N, x["s+"], x["s-"], x["lambda"], x["h"])
            value[on] = np.log(tails(lnA, np.zeros(len(N)), nu)[1])
    phi = {t: np.zeros(len(N)) for t in TERMS}
    for i, t in enumerate(TERMS):
        rest = [j for j in range(n) if j != i]
        for r in range(n):
            w = math.factorial(r) * math.factorial(n - r - 1) / math.factorial(n)
            for on in itertools.combinations(rest, r):
                phi[t] += w * (value[tuple(sorted(on + (i,)))] - value[on])
    return phi, value[()]


# ── data ──────────────────────────────────────────────────────────────────────
def proteome(keep_extrapolated):
    Q = proteome_chains(keep_extrapolated, "none").copy().reset_index(drop=True)
    Q["s"] = Q.scd_pos + Q.scd_neg
    for m in ("net", "cf"):
        lnA, lnS, nu = Q[f"{m}_lnA"].to_numpy(), np.log(Q[f"{m}_S_m"].to_numpy()), Q[f"{m}_nu"].to_numpy()
        _, Q[f"{m}_pc"], Q[f"{m}_pe"] = tails(lnA, lnS, nu)
    return Q


def validation(frames_path):
    P = data.proteins()
    V = P.iloc[data.split(P, "standard")[1]].copy().reset_index(drop=True)
    V["h"] = CF.regressors(V, **CF.load()["ce_threshold"]).h
    net, spec = network.load()
    coef = CF.load()
    for tag, (lnA, lnS, u) in (("net", network.parameters(net, spec, V)), ("cf", CF.parameters(coef, V))):
        lnA, lnS, u = (np.asarray(a, float) for a in (lnA, lnS, u))
        V[f"{tag}_mean"], V[f"{tag}_pc"], V[f"{tag}_pe"] = tails(lnA, lnS, D.nu_of(u))
    F = data.frames(frames_path)
    floor = 1 / F.shape[1]                  # one frame: a tail with no stored frame still has a logarithm
    fr = np.asarray(F[V.frame_row.to_numpy()], float)
    m = fr.mean(1)
    V["md_mean"], V["md_pc"] = m, np.maximum((fr < TC * m[:, None]).mean(1), floor)
    V["md_pe"] = np.maximum((fr > TE * m[:, None]).mean(1), floor)
    return V, P, F, floor


def replicates(P, F, floor):
    """Every pair of simulations of one sequence: |ln fold| of P_c and P_e, which can only be sampling noise."""
    rows = []
    for _, g in P[P.sequence.duplicated(keep=False)].groupby("sequence"):
        fr = np.asarray(F[g.frame_row.to_numpy()], float)
        m = fr.mean(1)
        pc, pe = np.maximum((fr < TC * m[:, None]).mean(1), floor), np.maximum((fr > TE * m[:, None]).mean(1), floor)
        for i, j in itertools.combinations(range(len(g)), 2):
            rows.append(dict(N=int(g.N.iloc[0]), dpc=np.log(pc[i] / pc[j]), dpe=np.log(pe[i] / pe[j]),
                             dmean=np.log(m[i] / m[j])))
    return pd.DataFrame(rows)


# ── figures ───────────────────────────────────────────────────────────────────
def survival(ax, lnf, **kw):
    x = np.exp(np.sort(lnf))
    ax.step(x, 1 - np.arange(len(x)) / len(x), where="post", **kw)


def size_panel(ax, N, y, trend_net, trend_cf, md, ylabel, title, note, ylim):
    ax.hexbin(N, y, gridsize=(46, 34), xscale="log", extent=(np.log10(20), np.log10(490), *ylim), cmap=DENSITY, bins="log",
              mincnt=1, linewidths=0, rasterized=True)
    n = np.exp(np.linspace(np.log(20), np.log(490), 200))
    ax.plot(n, trend_net(n), label="proteome, network (median)", color=BURGUNDY, lw=1.5)
    ax.plot(n, trend_cf(n), label="proteome, closed form (median)", **CF_STYLE)
    ax.plot(md[0], md[1], "o", ms=3.4, mfc=MOSS, mec="white", mew=0.4, label="validation, MD (median)", zorder=5)
    ax.set_xscale("log"), ax.set_xlim(20, 490), ax.set_ylim(*ylim)
    ax.set_xticks([20, 50, 100, 200, 400]), ax.set_xticklabels(["20", "50", "100", "200", "400"])
    ax.minorticks_off()
    ax.set_xlabel("chain length N"), ax.set_ylabel(ylabel)
    ax.set_title(title, loc="left", pad=3, fontsize=7.5)
    ax.text(0.97, 0.05, note, transform=ax.transAxes, ha="right", va="bottom", fontsize=6.3, color=INK2, linespacing=1.3)


def figure_size_pairs(Q, V, trends, fold, rep, pick, m):
    fig = plt.figure(figsize=(COL2, 5.1), layout="constrained")
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.15])
    a, b, c, d = [fig.add_subplot(gs[i, j]) for i, j in ((0, 0), (0, 1), (1, 0), (1, 1))]
    N = Q.N.to_numpy(float)
    bins = np.unique(np.quantile(V.N, np.linspace(0, 1, 11)))
    grp = pd.cut(V.N, bins, include_lowest=True)
    md = lambda col: (V.groupby(grp, observed=True).N.median().to_numpy(), V.groupby(grp, observed=True)[col].median().to_numpy())
    size_panel(a, N, Q.net_pc, trends["net_pc"], trends["cf_pc"], md("md_pc"), r"$P_c=P(R_g<0.8\langle R_g\rangle)$", "compact population",
               f"size explains {100 * m['size_r2']['net_pc']:.0f}% of the variance of ln $P_c$", (0, 0.22))
    size_panel(b, N, Q.net_pe, trends["net_pe"], trends["cf_pe"], md("md_pe"), r"$P_e=P(R_g>1.3\langle R_g\rangle)$", "extended population",
               f"size explains {100 * m['size_r2']['net_pe']:.0f}% of the variance of ln $P_e$", (0, 0.13))
    a.legend(loc="upper left", borderaxespad=0.1, handlelength=1.8, labelspacing=0.2, fontsize=6)

    survival(c, fold["proteome_net"], color=BURGUNDY, lw=1.6, label="proteome, network", zorder=5)
    survival(c, fold["proteome_cf"], label="proteome, closed form", **CF_STYLE)
    survival(c, fold["val_net"], color=BURGUNDY, lw=1.0, ls=(0, (1.2, 1.2)), label="validation, network", zorder=4)
    survival(c, fold["val_md"], color=MOSS, lw=1.6, label="validation, MD", zorder=4)
    survival(c, np.abs(rep.dpc.to_numpy()), color=INK2, lw=0.9, label="same sequence simulated twice", zorder=3)
    c.set_xscale("log"), c.set_yscale("log"), c.set_xlim(1, 4), c.set_ylim(1e-5, 1.3)
    c.set_xticks([1, 1.25, 1.5, 2, 3, 4]), c.set_xticklabels(["1", "1.25", "1.5", "2", "3", "4"]), c.minorticks_off()
    c.set_yticks([1, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5])
    c.set_xlabel(r"fold difference in $P_c$ within a pair"), c.set_ylabel("fraction of pairs at or above")
    c.set_title(r"pairs with equal $\langle R_g\rangle$ (within 1%)", loc="left", pad=3, fontsize=7.5)
    c.legend(loc="upper right", borderaxespad=0.1, handlelength=1.8, labelspacing=0.25, fontsize=6)

    A, B = pick["compact"], pick["extended"]
    dF_pair = m["free_energy_pair"]["ddF_kT"]
    lo, hi = 0.45 * A.net_mean_rg, 1.9 * A.net_mean_rg
    x = np.linspace(lo, hi, 500)
    for row, color in ((A, ORANGE), (B, MODEL)):
        F = -D.logpdf(x, row.net_lnA, np.log(row.net_S_m), row.net_nu)
        F -= F.min()
        d.plot(x, F, color=color, lw=1.5)
    mean_rg = 0.5 * (A.net_mean_rg + B.net_mean_rg)
    d.axvline(mean_rg, color=INK2, lw=0.5, ls=":"), d.axvline(TC * mean_rg, color=INK2, lw=0.5, ls=":")
    d.axvspan(lo, TC * mean_rg, color=LIGHT, alpha=0.35, lw=0, zorder=0)
    d.set_xlim(lo, hi), d.set_ylim(0, 6)
    d.text(TC * mean_rg, 5.85, r"$0.8\langle R_g\rangle$ ", ha="right", va="top", fontsize=6.3, color=INK2)
    d.text(mean_rg, 5.85, r" $\langle R_g\rangle$", ha="left", va="top", fontsize=6.3, color=INK2)
    d.set_xlabel(r"$R_g$ (nm)"), d.set_ylabel(r"$-\ln P(R_g)$, minimum at 0 ($k_BT$)")
    d.set_title(rf"an equal-mean pair from c: $\Delta\Delta F_c$ = {dF_pair:.2f} $k_BT$", loc="left", pad=3, fontsize=7.5)
    dF = lambda p: -np.log(p / (1 - p))
    label = lambda r: f"{r.gene if isinstance(r.gene, str) and r.gene else r.accession}, N = {r.N}\n" + rf"$P_c$ = {r.net_pc:.3f},  $\Delta F_c$ = {dF(r.net_pc):.2f} $k_BT$"
    d.legend(handles=[Line2D([], [], color=ORANGE, lw=1.5, label=label(A)), Line2D([], [], color=MODEL, lw=1.5, label=label(B))],
             loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=2, borderaxespad=0, handlelength=1.3, columnspacing=1.2, fontsize=6.2)
    for ax, letter in zip((a, b, c, d), "abcd"):
        panel_label(ax, letter, x=-0.1, y=1.02)
    return fig


def figure_sequence(Q, out, slopes):
    fig = plt.figure(figsize=(COL2, 5.0), layout="constrained")
    gs = fig.add_gridspec(2, 2, height_ratios=[1.12, 1])
    a, b, c, d = [fig.add_subplot(gs[i, j]) for i, j in ((0, 0), (0, 1), (1, 0), (1, 1))]
    norm = TwoSlopeNorm(vmin=-0.25, vcenter=0, vmax=0.25)
    ext = (-0.6, 2.0, 0, 4)
    for ax, col, title in ((a, "net_rpc", "network"), (b, "cf_rpc", "closed form")):
        hb = ax.hexbin(Q.s, Q.h, C=Q[col], reduce_C_function=np.median, gridsize=(34, 22), extent=ext, mincnt=3, cmap=DIVERGING, norm=norm, linewidths=0, rasterized=True)
        ax.set_xlim(ext[0], ext[1]), ax.set_ylim(ext[2] - 0.08, ext[3])
        ax.axvline(0, color=INK2, lw=0.4, ls=":")
        ax.set_xlabel(r"$s=s_++s_-$")
        ax.text(0.98, 0.96, "s < 0: opposite charges in blocks\ns > 0: like charges dominate", transform=ax.transAxes, ha="right", va="top", fontsize=6.2, color=INK2, linespacing=1.3)
        ax.set_ylabel(r"$h$, excess contact energy")
        ax.set_title(title, loc="left", pad=3, fontsize=7.5)
    cb = fig.colorbar(hb, ax=[a, b], location="right", shrink=0.92, pad=0.01, aspect=26, extend="both")
    cb.set_label(r"$\ln P_c$ minus size trend (cell median)", fontsize=6.8)
    cb.ax.tick_params(labelsize=6.5), cb.outline.set_linewidth(0.4)

    y, h = np.arange(len(TERMS)), 0.36
    for sign, color, off, name in ((-1, MODEL, -h / 2 - 0.015, "depleted"), (1, BURGUNDY, h / 2 + 0.015, "enriched")):
        vals = [out[name]["mean_term"][t] for t in TERMS]
        dom = [out[name]["dominant_share"][t] for t in TERMS]
        c.barh(y + off, vals, h, color=color, linewidth=0, label=f"{out[name]['label']}: mean residual {out[name]['net_residual']:+.2f}\n(closed form {out[name]['cf_residual']:+.2f})")
        for yi, v, s in zip(y + off, vals, dom):
            if s > 0:
                c.text(v + (0.008 if v >= 0 else -0.008), yi, f"{100 * s:.0f}%", ha="left" if v >= 0 else "right", va="center", fontsize=5.8, color=INK2)
    c.axvline(0, color=INK2, lw=0.5)
    c.set_yticks(y), c.set_yticklabels([LABEL[t] for t in TERMS], fontsize=8), c.invert_yaxis()
    span = max(abs(v) for g in ("depleted", "enriched") for v in out[g]["mean_term"].values())
    c.set_xlim(-1.32 * span, 1.32 * max(span * 0.6, max(out["enriched"]["mean_term"].values()), max(out["depleted"]["mean_term"].values())))
    c.set_xlabel(r"mean closed-form term of $\ln P_c$")
    c.set_title("closed-form terms of the network's outliers", loc="left", pad=3, fontsize=7.5)
    c.text(0.99, 0.97, "%: share led by term", transform=c.transAxes, ha="right", va="top", fontsize=5.8, color=INK2)
    c.legend(loc="lower left", borderaxespad=0.1, handlelength=1.1, labelspacing=0.3, fontsize=5.9)

    d.axvspan(-10, 0, color=LIGHT, alpha=0.3, lw=0, zorder=0)
    d.axvline(0, color=INK2, lw=0.5), d.axvline(1, color=INK2, lw=0.5, ls=(0, (4, 3)))
    for i, t in enumerate(TERMS):
        s = slopes[t]
        d.plot(s["ci"], [i, i], color=INK, lw=1.0, solid_capstyle="butt")
        d.plot(s["beta"], i, "o", ms=4.2, mfc=INK, mec="white", mew=0.4, zorder=4)
    d.set_yticks(y), d.set_yticklabels([LABEL[t] for t in TERMS], fontsize=8), d.invert_yaxis()
    lo = min(s["ci"][0] for s in slopes.values())
    hi = max(s["ci"][1] for s in slopes.values())
    d.set_xlim(min(-0.4, lo - 0.25), max(2.6, hi + 0.25))
    d.set_xlabel(r"slope of the MD residual on the term (validation)")
    d.set_title("does MD support each term?", loc="left", pad=3, fontsize=7.5)
    d.text(1.08, -0.62, "1: calibrated", ha="left", va="center", fontsize=5.8, color=INK2)
    d.text(-0.08, -0.62, "0: no signal", ha="right", va="center", fontsize=5.8, color=INK2)
    d.text(0.02, 0.03, "wrong direction", transform=d.transAxes, ha="left", va="bottom", fontsize=5.8, color=INK2)
    d.set_ylim(len(TERMS) - 0.5, -0.72)
    for ax, letter in zip((a, b, c, d), "abcd"):
        panel_label(ax, letter, x=-0.1, y=1.02)
    return fig


def save(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(RES / f"{name}.{ext}", dpi=300)
    plt.close(fig)
    print(f"-> {RES / name}.png / .pdf", file=sys.stderr)


# ── main ──────────────────────────────────────────────────────────────────────
def boot_clusters(cluster, rng, n):
    idx = {k: np.flatnonzero(cluster == k) for k in np.unique(cluster)}
    keys = list(idx)
    for _ in range(n):
        yield np.concatenate([idx[keys[i]] for i in rng.integers(len(keys), size=len(keys))])


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--frames", default=str(HERE.parent / "data" / "idp_rg_frames.npy"))
    ap.add_argument("--tol", type=float, default=0.01)
    ap.add_argument("--keep-extrapolated", action="store_true")
    ap.add_argument("--boots", type=int, default=2000)
    args = ap.parse_args(argv)
    set_style()
    rng = np.random.default_rng(2026)
    coef = CF.load()

    Q = proteome(args.keep_extrapolated)
    V, P, F, floor = validation(args.frames)
    rep = replicates(P, F, floor)
    N, Nv = Q.N.to_numpy(float), V.N.to_numpy(float)
    m = dict(n_chains=len(Q), tiers=Q.tier.value_counts().to_dict(), n_val=len(V), tol=args.tol,
             thresholds=dict(compact=TC, extended=TE))

    # size trends and residuals
    trends, m["size_r2"] = {}, {}
    for tag in ("net", "cf"):
        for q in ("pc", "pe"):
            y = np.log(Q[f"{tag}_{q}"].to_numpy())
            trends[f"{tag}_{q}"] = lambda n, f=trend(N, y, 40)[0]: np.exp(f(n))
            r = y - np.log(trends[f"{tag}_{q}"](N))
            Q[f"{tag}_r{q}"] = r
            m["size_r2"][f"{tag}_{q}"] = float(1 - r.var() / y.var())
    for tag in ("net", "cf", "md"):
        for q in ("pc", "pe"):
            y = np.log(V[f"{tag}_{q}"].to_numpy())
            V[f"{tag}_r{q}"] = y - trend(Nv, y, 20)[0](Nv)
    m["population"] = {f"{t}_{q}": dict(p5=float(Q[f"{t}_{q}"].quantile(.05)), median=float(Q[f"{t}_{q}"].median()), p95=float(Q[f"{t}_{q}"].quantile(.95)),
                                        rho_N=float(spearmanr(Q[f"{t}_{q}"], Q.N).statistic)) for t in ("net", "cf") for q in ("pc", "pe")}
    m["population_vs_spread_rho"] = {q: float(spearmanr(Q[f"net_{q}"], Q.net_std_rg / Q.net_mean_rg).statistic) for q in ("pc", "pe")}

    # thresholds set by the spread instead: do they keep the signal?
    lnA, lnS, nu, mr = Q.net_lnA.to_numpy(), np.log(Q.net_S_m.to_numpy()), Q.net_nu.to_numpy(), Q.net_mean_rg.to_numpy()
    cv = (Q.net_std_rg / Q.net_mean_rg).to_numpy()
    span = lambda p: float(np.quantile(p, .95) / np.quantile(p, .05))
    m["threshold_choice"] = dict(ratio=dict(pc_p95_over_p5=span(Q.net_pc), pe_p95_over_p5=span(Q.net_pe)), spread={})
    for k in (1.33, 2.0):
        lo, up = cdf_at(mr * (1 - k * cv), lnA, lnS, nu), 1 - cdf_at(mr * (1 + k * cv), lnA, lnS, nu)
        m["threshold_choice"]["spread"][f"{k}sd"] = dict(low_p95_over_p5=span(lo), high_p95_over_p5=span(up), low_median=float(np.median(lo)), high_median=float(np.median(up)),
                                                         low_rho_nu=float(spearmanr(lo, nu).statistic), high_rho_nu=float(spearmanr(up, nu).statistic))

    # equal-mean pairs
    rows, fold = [], {}
    def add(name, model, quantity, f):
        rows.append(dict(set=name, model=model, quantity=quantity, **fold_stats(f)))
        return f
    Ip, Jp = pairs(mr, args.tol)
    Ic, Jc = pairs(Q.cf_mean_rg.to_numpy(), args.tol)
    Imd, Jmd = pairs(V.md_mean.to_numpy(), args.tol)
    Inet, Jnet = pairs(V.net_mean.to_numpy(), args.tol)
    Icf, Jcf = pairs(V.cf_mean.to_numpy(), args.tol)
    lnfold = lambda x, i, j: np.abs(np.log(x[i] / x[j]))
    for q, nm in (("pc", "P_c"), ("pe", "P_e")):
        fold[f"proteome_net_{q}"] = add("proteome", "network", nm, lnfold(Q[f"net_{q}"].to_numpy(), Ip, Jp))
        fold[f"proteome_cf_{q}"] = add("proteome", "closed form", nm, lnfold(Q[f"cf_{q}"].to_numpy(), Ic, Jc))
        fold[f"val_net_{q}"] = add("validation", "network", nm, lnfold(V[f"net_{q}"].to_numpy(), Inet, Jnet))
        fold[f"val_cf_{q}"] = add("validation", "closed form", nm, lnfold(V[f"cf_{q}"].to_numpy(), Icf, Jcf))
        fold[f"val_md_{q}"] = add("validation", "MD", nm, lnfold(V[f"md_{q}"].to_numpy(), Imd, Jmd))
        add("replicates", "MD", nm, np.abs(rep[f"d{q}"].to_numpy()))
    conf = Q.tier.to_numpy() == "confident"
    keep = conf[Ip] & conf[Jp]
    add("proteome, confident tier only", "network", "P_c", fold["proteome_net_pc"][keep])
    add("proteome, confident tier only", "network", "P_e", fold["proteome_net_pe"][keep])
    for t in (0.70, 0.75, 0.85):
        add(f"proteome, threshold {t}", "network", f"P(Rg<{t}<Rg>)", lnfold(cdf_at(t * mr, lnA, lnS, nu), Ip, Jp))
    for t in (1.2, 1.4):
        add(f"proteome, threshold {t}", "network", f"P(Rg>{t}<Rg>)", lnfold(1 - cdf_at(t * mr, lnA, lnS, nu), Ip, Jp))
    T = pd.DataFrame(rows)
    T.to_csv(RES / "compact_states_pairs.csv", index=False, float_format="%.5g")
    big = fold["proteome_net_pc"] >= np.log(1.5)
    m["chains_with_partner"] = dict(ge_1p5=float(len(np.unique(np.r_[Ip[big], Jp[big]])) / len(Q)))
    cau = ~conf
    m["pairs_ge_1p5_with_caution_chain"] = float((cau[Ip[big]] | cau[Jp[big]]).mean())
    m["pair_ln_N_ratio_median"] = float(np.median(np.abs(np.log(N[Ip] / N[Jp]))))
    rb = pd.cut(np.sqrt(mr[Ip] * mr[Jp]), [0, 2, 3, 4, 6, 20])
    m["proteome_by_rg"] = {str(k): dict(pairs=len(g), median=float(np.exp(np.median(g))), p99=float(np.exp(np.quantile(g, .99)))) for k, g in pd.Series(fold["proteome_net_pc"]).groupby(rb, observed=True)}

    # calibration: replicate noise, then compression on the same MD pairs and per protein
    sig = {q: float(rep[f"d{q}"].std() / np.sqrt(2)) for q in ("pc", "pe")}
    sig_rob = {q: float(1.4826 * np.median(np.abs(rep[f"d{q}"])) / np.sqrt(2)) for q in ("pc", "pe")}
    boot = pd.DataFrame([{q: rep.iloc[rng.integers(len(rep), size=len(rep))][f"d{q}"].std() / np.sqrt(2) for q in ("pc", "pe")} for _ in range(args.boots)])
    m["replicates"] = dict(pairs=len(rep), N_median=float(rep.N.median()), noise_sd_ln=sig, noise_sd_ln_robust=sig_rob, noise_sd_ln_ci={q: [float(v) for v in boot[q].quantile([.025, .975])] for q in ("pc", "pe")},
                           noise_sd_ln_N_ge_100={q: float(rep[rep.N >= 100][f"d{q}"].std() / np.sqrt(2)) for q in ("pc", "pe")}, mean_rg_noise_pct=float(100 * rep.dmean.std() / np.sqrt(2)))
    m["calibration"] = {}
    for q in ("pc", "pe"):
        md, net, cf_ = [np.log(V[f"{t}_{q}"].to_numpy()) for t in ("md", "net", "cf")]
        smd, snet, scf = md[Imd] - md[Jmd], net[Imd] - net[Jmd], cf_[Imd] - cf_[Jmd]
        scale = lambda x: 1.4826 * np.median(np.abs(x))        # outlier-resistant pair sd
        md_sd = np.sqrt(max(scale(smd) ** 2 - 2 * sig[q] ** 2, 0))
        rmd, rnet, rcf = V[f"md_r{q}"].to_numpy(), V[f"net_r{q}"].to_numpy(), V[f"cf_r{q}"].to_numpy()
        true_sd = lambda r: np.sqrt(max(r.var(ddof=1) - sig[q] ** 2, 1e-12))
        def per_protein(i):
            return dict(net_over_md_sd=rnet[i].std(ddof=1) / true_sd(rmd[i]), cf_over_md_sd=rcf[i].std(ddof=1) / true_sd(rmd[i]),
                        md_on_net=np.polyfit(rnet[i], rmd[i], 1)[0], md_on_cf=np.polyfit(rcf[i], rmd[i], 1)[0],
                        rho_net=spearmanr(rnet[i], rmd[i]).statistic, rho_cf=spearmanr(rcf[i], rmd[i]).statistic)
        est = per_protein(np.arange(len(V)))
        b = pd.DataFrame([per_protein(i) for i in boot_clusters(V.cluster.to_numpy(), rng, args.boots)])
        m["calibration"][q] = dict(pairs_on_md=len(smd), median_fold_md_raw=float(np.exp(np.median(np.abs(smd)))), median_fold_md_denoised=float(np.exp(0.6745 * md_sd)), median_fold_net_same_pairs=float(np.exp(np.median(np.abs(snet)))), median_fold_cf_same_pairs=float(np.exp(np.median(np.abs(scf)))),
                                   pair_rho_net=float(spearmanr(smd, snet).statistic), pair_rho_cf=float(spearmanr(smd, scf).statistic),
                                   per_protein={k: dict(value=float(v), ci=[float(x) for x in b[k].quantile([.025, .975])]) for k, v in est.items()})
    vb = pd.cut(V.md_mean, [0, 3, 4, 5, 20])
    m["calibration"]["pc"]["sd_by_rg"] = {str(k): dict(n=len(g), md=float(g.md_rpc.std()), net=float(g.net_rpc.std()), cf=float(g.cf_rpc.std())) for k, g in V.groupby(vb, observed=True)}

    # closed-form attribution, proteome and validation; reference inputs are the proteome medians at each N
    ref_fn = {t: trend(N, Q[INPUT[t]].to_numpy(float), 40)[0] for t in TERMS}
    phi, base = shapley(coef, Q, {t: ref_fn[t](N) for t in TERMS})
    for t in TERMS:
        Q[f"phi_{t}"] = phi[t]
    closure = np.log(Q.cf_pc.to_numpy()) - base - sum(phi.values())
    assert np.abs(closure).max() < 1e-6, f"Shapley terms do not sum to ln P_c ({np.abs(closure).max()})"
    phi_v, _ = shapley(coef, V, {t: ref_fn[t](Nv) for t in TERMS})
    m["cf_vs_net"] = dict(residual_rho_pc=float(spearmanr(Q.net_rpc, Q.cf_rpc).statistic), residual_rho_pe=float(spearmanr(Q.net_rpe, Q.cf_rpe).statistic))
    m["reference_vs_trend_ln_pc_p99_abs"] = float(np.quantile(np.abs(base - np.log(trends["cf_pc"](N))), .99))

    X = np.column_stack([phi_v[t] for t in TERMS] + [np.ones(len(V))])
    y = V.md_rpc.to_numpy()
    b_all = np.linalg.lstsq(X, y, rcond=None)[0]
    b_b = np.array([np.linalg.lstsq(X[i], y[i], rcond=None)[0] for i in boot_clusters(V.cluster.to_numpy(), rng, args.boots)])
    slopes = {t: dict(beta=float(b_all[i]), ci=[float(x) for x in np.percentile(b_b[:, i], [2.5, 97.5])], sd_term_val=float(np.std(phi_v[t])), sd_term_proteome=float(np.std(phi[t]))) for i, t in enumerate(TERMS)}
    m["term_vs_md"] = slopes

    out = {}
    lo_cut, hi_cut = Q.net_rpc.quantile([0.01, 0.99])
    for name, mask, sign, lab in (("depleted", Q.net_rpc <= lo_cut, -1, "lowest 1%"), ("enriched", Q.net_rpc >= hi_cut, 1, "highest 1%")):
        g = Q[mask]
        A_ = np.column_stack([g[f"phi_{t}"].to_numpy() for t in TERMS])
        dom = np.argmax(sign * A_, axis=1)
        out[name] = dict(label=lab, n=len(g), cut=float(lo_cut if sign < 0 else hi_cut), net_residual=float(g.net_rpc.mean()), cf_residual=float(g.cf_rpc.mean()), cf_same_sign=float((np.sign(g.cf_rpc) == sign).mean()),
                         N_median=float(g.N.median()), tiers=g.tier.value_counts().to_dict(), mean_term={t: float(A_[:, i].mean()) for i, t in enumerate(TERMS)},
                         dominant_share={t: float((dom == i).mean()) for i, t in enumerate(TERMS)}, medians={c: float(g[c].median()) for c in ("scd_pos", "scd_neg", "mean_hydro", "h", "fcr", "N")})
    m["outliers"] = out

    # the free-energy pair: the largest fold among confident, reviewed chains of similar length where the closed form agrees
    rev = Q.reviewed.fillna(False).to_numpy(bool)
    cf_fold = np.log(Q.cf_pc.to_numpy()[Ip] / Q.cf_pc.to_numpy()[Jp])
    net_fold = np.log(Q.net_pc.to_numpy()[Ip] / Q.net_pc.to_numpy()[Jp])
    ok = conf[Ip] & conf[Jp] & rev[Ip] & rev[Jp] & (np.abs(np.log(N[Ip] / N[Jp])) <= np.log(1.15)) & (np.sign(cf_fold) == np.sign(net_fold)) & (np.abs(cf_fold) >= 0.5 * np.abs(net_fold))
    k = np.flatnonzero(ok)[np.argmax(np.abs(net_fold[ok]))]
    hi_i, lo_i = (Ip[k], Jp[k]) if net_fold[k] > 0 else (Jp[k], Ip[k])
    A, B = Q.iloc[hi_i], Q.iloc[lo_i]
    dF = lambda p: float(-np.log(p / (1 - p)))
    rank = float((fold["proteome_net_pc"] < abs(net_fold[k])).mean())
    m["free_energy_pair"] = dict(rank_quantile=rank, fold_pc=float(np.exp(abs(net_fold[k]))), fold_pc_closed_form=float(np.exp(abs(cf_fold[k]))), ddF_kT=dF(B.net_pc) - dF(A.net_pc), ddF_kcal=(dF(B.net_pc) - dF(A.net_pc)) * KT_KCAL,
                                 proteins={r.chain: dict(gene=r.gene, description=r.description, N=int(r.N), mean_rg=float(r.net_mean_rg), P_c=float(r.net_pc), P_e=float(r.net_pe), dF_kT=dF(r.net_pc), P_c_closed_form=float(r.cf_pc), s=float(r.s), h=float(r.h), mean_hydro=float(r.mean_hydro), fcr=float(r.fcr), tier=r.tier) for r in (A, B)})

    chains_csv = ["chain", "gene", "N", "tier", "s", "scd_pos", "scd_neg", "h", "mean_hydro", "net_mean_rg", "net_pc", "net_pe", "net_rpc", "net_rpe", "cf_mean_rg", "cf_pc", "cf_pe", "cf_rpc", "cf_rpe"] + [f"phi_{t}" for t in TERMS]
    Q[chains_csv].to_csv(RES / "compact_states_chains.csv", index=False, float_format="%.5g")
    (RES / "compact_states.json").write_text(json.dumps(m, indent=1, default=float))

    fold_c = {k[:-3]: v for k, v in fold.items() if k.endswith("_pc")}
    save(figure_size_pairs(Q, V, trends, fold_c, rep, dict(compact=A, extended=B), m), "compact_states_size_pairs")
    save(figure_sequence(Q, out, slopes), "compact_states_sequence")
    print(T.to_string(index=False, float_format=lambda x: f"{x:.4g}"))
    print(json.dumps({k: m[k] for k in ("size_r2", "threshold_choice", "chains_with_partner", "pairs_ge_1p5_with_caution_chain", "pair_ln_N_ratio_median", "proteome_by_rg", "replicates", "calibration", "cf_vs_net", "reference_vs_trend_ln_pc_p99_abs", "term_vs_md", "outliers", "free_energy_pair", "population_vs_spread_rho")}, indent=1, default=float))


if __name__ == "__main__":
    main()
