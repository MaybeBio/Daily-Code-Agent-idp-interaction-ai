"""The human disordered proteome, one dot per chain: size against compact and extended populations, and the two populations together.

    python proteome/populations.py [--frames data/idp_rg_frames.npy] [--boots 1000]

Panel a is the distribution of the median chain with its two tails coloured. Panel b plots P_c against P_e, split at the 5% a
typical training chain has by construction, and coloured by size: <Rg>/R_Flory, where R_Flory is the training proteins' Flory fit
at the chain's length. Panels c and d plot P_c and P_e against the same ratio. Populations are in percent. P_c and P_e use
idprg.distribution.population_bounds, the 5th and 95th percentiles of a typical chain of the same length. Chains are one per
distinct sequence, without the extrapolated tier (rg_groups.chains). The lines are medians over the validation proteins, in MD
(green) and in the network (red), so the network can be checked.

This script also refits those bounds from the training frames and scores them against the alternatives on the validation
proteins, so the choice can be checked: replicate noise (proteins simulated more than once), reliability (1 - noise
variance / variance across proteins), agreement of the network with MD, and redundancy in MD with spread and compaction.

Writes to proteome/results/:
    populations.{png,pdf}              the figure
    populations.json                   every number quoted
    populations_thresholds.csv         the score of every threshold definition on the validation proteins
    populations_tails.csv              each tail at each depth, network against MD: validation and the benchmark
    populations_chains.csv             one row per chain
"""
import argparse
import itertools
import json
import pathlib
import pickle
import sys

import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "paper"))
import idprg  # noqa: E402,F401  (float64 in JAX)
from idprg import closed_form as CF, data, distribution as D, network  # noqa: E402
from make_figures import BURGUNDY, CACHE, COL2, INK, INK2, MODEL, MOSS, ORANGE, panel_label, set_style  # noqa: E402
from rg_groups import chains as proteome_chains  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm  # noqa: E402
from matplotlib.legend_handler import HandlerTuple  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from scipy.optimize import curve_fit  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402

RES = HERE / "results"
BASELINE = lambda N: np.exp(data.LN_A_F + data.NU_F * np.log(N))         # the training proteins' Flory fit
PERCENT = (0.02, 0.05, 0.10, 0.15, 0.20, 0.25)
METRICS = ("compaction", "compact", "extended")


def trend(N, y, bins):
    """Running median of y against ln N, on quantile bins of N (ties at short N merge). Returns a function of N."""
    e = np.unique(np.quantile(N, np.linspace(0, 1, bins + 1)))
    x, m = [], []
    for a, b, last in zip(e[:-1], e[1:], [False] * (len(e) - 2) + [True]):
        k = (N >= a) & ((N <= b) if last else (N < b))
        if k.any():
            x.append(np.median(np.log(N[k]))), m.append(np.median(y[k]))
    return lambda n: np.interp(np.log(n), x, m)


def residual(N, y, bins):
    return y - trend(N, y, bins)(N)


def fit_bounds(T, F):
    """a + b / sqrt(N) to the median per length bin of each training protein's percentiles of Rg / <Rg>."""
    fr = np.asarray(F[T.frame_row.to_numpy()], float)
    x, N = fr / fr.mean(1, keepdims=True), T.N.to_numpy(float)
    form = lambda n, a, b: a + b * n ** -0.5
    out = {}
    for q in PERCENT:
        for side, p in (("lo", q), ("hi", 1 - q)):
            d = pd.DataFrame(dict(n=N, y=np.quantile(x, p, axis=1)))
            b = d.groupby(pd.qcut(d.n, 20, duplicates="drop"), observed=True).median()
            out[side, q] = curve_fit(form, b.n, b.y)[0]
    return out


def definitions(bounds):
    """name -> (reference, t_c(N, cv), t_e(N, cv)). Reference `own` is the chain's own mean, `length` the Flory baseline."""
    form = lambda p: (lambda N, cv: p[0] + p[1] * N ** -0.5)
    d = {"fixed 0.8 / 1.3": ("own", lambda N, cv: 0.8 + 0 * N, lambda N, cv: 1.3 + 0 * N),
         "own mean -/+ 1.5 sd": ("own", lambda N, cv: 1 - 1.5 * cv, lambda N, cv: 1 + 1.5 * cv)}
    for q in PERCENT:
        d[f"length-matched {100 * q:g}%"] = ("own", form(bounds["lo", q]), form(bounds["hi", q]))
    d["absolute, 5% of own bounds"] = ("length", form(bounds["lo", 0.05]), form(bounds["hi", 0.05]))
    return d


def frame_tails(F, rows, N, tc, te, ref, floor):
    fr = np.asarray(F[rows], float)
    m = fr.mean(1)
    r = m if ref == "own" else BASELINE(N)
    cv = fr.std(1, ddof=1) / m
    return (np.maximum((fr < (tc(N, cv) * r)[:, None]).mean(1), floor),
            np.maximum((fr > (te(N, cv) * r)[:, None]).mean(1), floor))


def replicate_groups(P):
    return [g for _, g in P[P.sequence.duplicated(keep=False)].groupby("sequence")]


def score_thresholds(defs, V, net, Q, F, floor, groups):
    """One row per definition: noise, reliability and agreement on validation; redundancy in MD and in the network."""
    N, Nq = V.N.to_numpy(float), Q.N.to_numpy(float)
    res = lambda y: residual(N, y, 20)
    res_q = lambda y: residual(Nq, y, 40)
    r_cv, r_comp = res(np.log(V.md_cv.to_numpy())), res(np.log(V.md_mean.to_numpy()))
    n_cv, n_comp = res(np.log(net["cv"])), res(np.log(net["mean"]))
    rows = []
    for name, (ref, tc, te) in defs.items():
        pc, pe = frame_tails(F, V.frame_row.to_numpy(), N, tc, te, ref, floor)
        rv = 1.0 if ref == "own" else BASELINE(N) / net["mean"]
        rq = 1.0 if ref == "own" else BASELINE(Nq) / Q.net_mean_rg.to_numpy()
        cvq = (Q.net_std_rg / Q.net_mean_rg).to_numpy()
        nc, ne = D.tail_below(tc(N, net["cv"]) * rv, net["lnA"], net["nu"]), 1 - D.tail_below(te(N, net["cv"]) * rv, net["lnA"], net["nu"])
        qc, qe = D.tail_below(tc(Nq, cvq) * rq, Q.net_lnA.to_numpy(), Q.net_nu.to_numpy()), 1 - D.tail_below(te(Nq, cvq) * rq, Q.net_lnA.to_numpy(), Q.net_nu.to_numpy())
        dc, de = [], []
        for g in groups:
            a, b = frame_tails(F, g.frame_row.to_numpy(), g.N.to_numpy(float), tc, te, ref, floor)
            for i, j in itertools.combinations(range(len(g)), 2):
                dc.append(np.log(a[i] / a[j])), de.append(np.log(b[i] / b[j]))
        noise_c, noise_e = np.std(dc) / 2 ** 0.5, np.std(de) / 2 ** 0.5
        a, b, c, d = res(np.log(pc)), res(np.log(pe)), res(np.log(nc)), res(np.log(ne))
        rho = lambda x, y: float(spearmanr(x, y).statistic)
        size_r2 = lambda y: float(1 - res_q(y).var() / y.var())
        rows.append(dict(
            definition=name, median_Pc_MD=np.median(pc), median_Pe_MD=np.median(pe), noise_Pc=noise_c, noise_Pe=noise_e,
            reliability_Pc=1 - noise_c ** 2 / a.var(), reliability_Pe=1 - noise_e ** 2 / b.var(),
            rho_net_MD_Pc=rho(c, a), rho_net_MD_Pe=rho(d, b), slope_MD_on_net_Pc=np.polyfit(c, a, 1)[0], slope_MD_on_net_Pe=np.polyfit(d, b, 1)[0],
            MD_rho_Pc_Pe=rho(a, b), MD_rho_Pc_spread=rho(a, r_cv), MD_rho_Pe_spread=rho(b, r_cv), MD_rho_Pc_compaction=rho(a, r_comp), MD_rho_Pe_compaction=rho(b, r_comp),
            net_rho_Pc_Pe=rho(c, d), net_rho_Pc_spread=rho(c, n_cv), net_rho_Pe_spread=rho(d, n_cv), net_rho_Pc_compaction=rho(c, n_comp), net_rho_Pe_compaction=rho(d, n_comp),
            proteome_size_r2_Pc=size_r2(np.log(qc)), proteome_size_r2_Pe=size_r2(np.log(qe))))
    return pd.DataFrame(rows)


def tail_skill(bounds, V, params, F, floor, rng, boots):
    """Network against simulation for each tail at each depth, length removed: validation (all, low charge) and the
    benchmark proteins with N >= 20 (HPS1 trajectories and STARLING, 95% interval from resampling families, with and without
    the A1-LCD variants). The noise limit is what replicate noise alone would allow the benchmark's rank correlation."""
    blob = pickle.loads(CACHE.read_bytes())
    smp = np.load(HERE.parent / "runs" / "starling" / "benchmark" / "rg_samples.npz")
    real = [r for r in blob["real"] if r["length"] >= 20 and r["name"] in smp.files]
    Nb = np.array([r["length"] for r in real], float)
    fam = np.array([r["family"] for r in real])
    pb = np.array([r["p"] for r in real])
    hps = [r["frames"] / r["frames"].mean() for r in real]
    star = [smp[r["name"]] / smp[r["name"]].mean() for r in real]
    idx = {f: np.flatnonzero(fam == f) for f in np.unique(fam)}
    keys, rest = list(idx), fam != "A1_LCD"
    Nv = V.N.to_numpy(float)
    fr = np.asarray(F[V.frame_row.to_numpy()], float)
    xv = fr / fr.mean(1, keepdims=True)
    low = V.scd_pos.to_numpy() < 0.25
    rho = lambda a, b: float(spearmanr(a, b).statistic)
    rows = []
    for side, q in itertools.product(("compact", "extended"), PERCENT):
        p = bounds["lo" if side == "compact" else "hi", q]
        t = lambda n: p[0] + p[1] * n ** -0.5
        up = side == "extended"
        emp = lambda x, u: np.mean(x > u) if up else np.mean(x < u)
        model = lambda n, lnA, nu: 1 - D.tail_below(t(n), lnA, nu) if up else D.tail_below(t(n), lnA, nu)
        rv = residual(Nv, np.log(model(Nv, params["lnA"], params["nu"])), 20)
        tv = t(Nv)[:, None]
        rm = residual(Nv, np.log(np.maximum((xv > tv).mean(1) if up else (xv < tv).mean(1), floor)), 20)
        L = lambda v: residual(Nb, np.log(np.maximum(v, 1e-4)), 6)
        rn = L(model(Nb, pb[:, 0], pb[:, 2]))
        rh = L(np.array([emp(x, u) for x, u in zip(hps, t(Nb))]))
        rs = L(np.array([emp(x, u) for x, u in zip(star, t(Nb))]))
        ci = lambda a, b: [float(v) for v in np.percentile([rho(a[i], b[i]) for i in (np.concatenate([idx[keys[j]] for j in rng.integers(len(keys), size=len(keys))]) for _ in range(boots))], [2.5, 97.5])]
        rows.append(dict(tail=side, percent=100 * q, validation=rho(rv, rm), validation_low_charge=rho(rv[low], rm[low]), net_robust_sd_benchmark=1.4826 * np.median(np.abs(rn - np.median(rn))),
                         benchmark_hps1=rho(rn, rh), benchmark_hps1_ci=ci(rn, rh), benchmark_hps1_no_A1=rho(rn[rest], rh[rest]), benchmark_starling=rho(rn, rs), benchmark_starling_ci=ci(rn, rs), benchmark_starling_no_A1=rho(rn[rest], rs[rest])))
    return pd.DataFrame(rows)


def quantile_table(x, conf):
    """x is a ln fold. Quantiles as folds, and the share of chains at least 1.25-fold away on each side."""
    q = dict(zip(("p1", "p5", "p25", "p75", "p95", "p99"), np.exp(np.quantile(x, [.01, .05, .25, .75, .95, .99]))))
    lo, hi = np.exp(np.quantile(x[conf], [.05, .95]))
    return dict(**{k: float(v) for k, v in q.items()}, min=float(np.exp(x.min())), max=float(np.exp(x.max())),
                below_0p8=float((x <= -np.log(1.25)).mean()), above_1p25=float((x >= np.log(1.25)).mean()), confident_p5=float(lo), confident_p95=float(hi))


NAMED = (("P07199_idr2", "CENP-B IDR"), ("P06454", "prothymosin-α"), ("P37840", "α-synuclein"), ("P35637_idr1", "FUS N-term"), ("Q9UPN9_idr5", "TRIM33 IDR"))
PHASE_LABEL_AT = {"CENP-B IDR": (5.15, 1.6), "prothymosin-α": (5.15, 2.75), "α-synuclein": (6.2, 6.75), "FUS N-term": (6.45, 4.35), "TRIM33 IDR": (0.4, 4.55)}
QUADRANTS = dict(Lb="fluctuation-suppressed\nstiff coils", Lt="asymmetric chains", Rt="broadly fluctuating coils", Rb="")
MD_BINS = [0.4, 0.8, 0.9, 0.95, 1, 1.05, 1.1, 1.2, 1.6]
XLIM, YLIM, SIZE_RANGE = (0.6, 1.52), (0, 8), (0.7, 1.4)
FLORY_RATIO = r"$\langle R_g\rangle / R_{\mathrm{Flory}}$"
SIZE_COLORS = LinearSegmentedColormap.from_list("size", ["#123f7a", MODEL, "#cfccc4", "#d9604c", "#7a1420"])


def scatter_data(Q, V):
    """Every population in percent, and <Rg>/R_Flory: the proteome, and validation in MD and the network."""
    S = Q.assign(size=Q.net_mean_rg / BASELINE(Q.N.to_numpy(float)), compact=100 * Q.net_pc, extended=100 * Q.net_pe)
    base = BASELINE(V.N.to_numpy(float))
    M = pd.DataFrame(dict(size=V.md_mean / base, compact=100 * V.md_pc, extended=100 * V.md_pe,
                          net_size=V.net_mean / base, net_compact=100 * V.net_pc, net_extended=100 * V.net_pe))
    return S, M


def validation_line(M, x, y):
    """Median size and population in the bins of MD_BINS; each source is binned on its own size."""
    grp = pd.cut(M[x], MD_BINS)
    return M.groupby(grp, observed=True)[x].median().to_numpy(), M.groupby(grp, observed=True)[y].median().to_numpy()


def median_chain(S):
    """The confident chain nearest the proteome medians of length, size, P_c and P_e, in units of each one's spread."""
    C = S[S.tier == "confident"]
    Z = np.column_stack([np.log(C.N), C["size"], C.compact, C.extended])
    return C.iloc[(((Z - np.median(Z, 0)) / Z.std(0)) ** 2).sum(1).argmin()]


def quadrant(compact, extended):
    return np.where(np.asarray(compact) < 5, "L", "R").astype(object) + np.where(np.asarray(extended) < 5, "b", "t").astype(object)


def typical_panel(ax, r):
    tc, te = (float(v) for v in D.population_bounds(r.N))
    x = np.linspace(0.3, 2.0, 1200)
    y = D.pdf(x * r.net_mean_rg, r.net_lnA, np.log(r.net_S_m), r.net_nu) * r.net_mean_rg
    ax.fill_between(x[x <= tc], y[x <= tc], color=MODEL, alpha=0.75, lw=0)
    ax.fill_between(x[x >= te], y[x >= te], color=ORANGE, alpha=0.85, lw=0)
    ax.plot(x, y, color=INK, lw=1.3)
    ax.axvline(1, color=INK2, lw=0.5, ls=":")
    top = y.max()
    ax.text(tc - 0.03, y[np.searchsorted(x, tc)] + 0.05 * top, rf"compact, $P_c$" + f"\n{r.compact:.1f}%", fontsize=6.5, ha="right", va="bottom", linespacing=1.25)
    ax.text(te + 0.03, y[np.searchsorted(x, te)] + 0.05 * top, rf"extended, $P_e$" + f"\n{r.extended:.1f}%", fontsize=6.5, ha="left", va="bottom", linespacing=1.25)
    ax.set_xlim(x[0], x[-1]), ax.set_ylim(0, 1.08 * top), ax.set_yticks(np.arange(0, top, 0.5))
    ax.set_xlabel(r"$R_g/\langle R_g\rangle$"), ax.set_ylabel("probability density")
    ax.set_title(f"the median chain: {r.gene} IDR, {int(r.N)} residues", fontsize=7)


def phase_panel(fig, ax, S):
    order = (S["size"] - 1).abs().to_numpy().argsort()                       # the rare, extreme chains on top
    dots = ax.scatter(S.compact.to_numpy()[order], S.extended.to_numpy()[order], c=S["size"].to_numpy()[order], s=2, cmap=SIZE_COLORS, linewidths=0, rasterized=True, zorder=1,
                      norm=TwoSlopeNorm(vcenter=1, vmin=SIZE_RANGE[0], vmax=SIZE_RANGE[1]))
    ax.axvline(5, color=INK2, lw=0.7, ls="--", zorder=2), ax.axhline(5, color=INK2, lw=0.7, ls="--", zorder=2)
    named = S.set_index("chain")
    for chain, name in NAMED:
        r = named.loc[chain]
        ax.plot(r.compact, r.extended, "o", ms=3.8, mfc=INK, mec="white", mew=0.5, zorder=4)
        ax.annotate(name, (r.compact, r.extended), PHASE_LABEL_AT[name], fontsize=6, ha="left", va="center", color=INK, zorder=5, arrowprops=dict(arrowstyle="-", color=INK2, lw=0.4, shrinkA=0.5, shrinkB=1.5))
    share = pd.Series(quadrant(S.compact, S.extended)).value_counts(normalize=True)
    for key, xy, ha, va in (("Lb", (0.015, 0.53), "left", "top"), ("Lt", (0.015, 0.985), "left", "top"), ("Rt", (0.985, 0.985), "right", "top"), ("Rb", (0.985, 0.015), "right", "bottom")):
        head = f"{QUADRANTS[key]}\n" if QUADRANTS[key] else ""
        ax.text(*xy, f"{head}{100 * share[key]:.0f}% of chains", transform=ax.transAxes, ha=ha, va=va, fontsize=6.2, color=INK2, linespacing=1.25, style="italic", zorder=5)
    ax.set_xlim(*YLIM), ax.set_ylim(*YLIM), ax.set_aspect("equal")
    ticks = range(0, 9, 2)
    ax.set_xticks(list(ticks)), ax.set_xticklabels([f"{t}%" for t in ticks]), ax.set_yticks(list(ticks)), ax.set_yticklabels([f"{t}%" for t in ticks])
    ax.set_xlabel(r"compact population $P_c$"), ax.set_ylabel(r"extended population $P_e$")
    bar = fig.colorbar(dots, ax=ax, shrink=0.85, pad=0.02, extend="both")
    tick_values = [0.7, 0.85, 1, 1.2, 1.4]
    bar.set_ticks(tick_values), bar.set_ticklabels([f"{v:g}" for v in tick_values])
    bar.set_label(FLORY_RATIO, fontsize=7), bar.ax.tick_params(labelsize=6.5), bar.outline.set_linewidth(0.4)


def figure(S, M):
    fig, axs = plt.subplots(2, 2, figsize=(COL2, 6.1), layout="constrained")
    (typical, phase), (compact, extended) = axs
    typical_panel(typical, median_chain(S))
    phase_panel(fig, phase, S)
    for ax, col, color, symbol, name in ((compact, "compact", MODEL, "P_c", "compact"), (extended, "extended", ORANGE, "P_e", "extended")):
        ax.scatter(S["size"], S[col], s=1.8, color=color, alpha=0.16, linewidths=0, rasterized=True, zorder=1)
        ax.axvline(1, color=INK2, lw=0.5, ls=":", zorder=2), ax.axhline(5, color=INK2, lw=0.7, ls="--", zorder=2)
        ax.plot(*validation_line(M, "size", col), "-o", color=MOSS, lw=1.3, ms=3, mec="white", mew=0.4, zorder=3)
        ax.plot(*validation_line(M, "net_size", f"net_{col}"), "-o", color=BURGUNDY, lw=1.3, ms=3, mec="white", mew=0.4, zorder=4)
        ax.set_xlim(*XLIM), ax.set_ylim(*YLIM)
        ax.set_xticks([0.6, 0.8, 1, 1.2, 1.4])
        ax.set_yticks(range(0, 9, 2)), ax.set_yticklabels([f"{t}%" for t in range(0, 9, 2)])
        ax.set_xlabel(FLORY_RATIO)
        ax.set_ylabel(rf"{name} population ${symbol}$")
    for ax, letter, x in ((typical, "a", -0.13), (phase, "b", -0.17), (compact, "c", -0.13), (extended, "d", -0.13)):
        panel_label(ax, letter, x=x, y=1.02)
    fig.legend(handles=[(Line2D([], [], ls="", marker="o", ms=3, color=MODEL), Line2D([], [], ls="", marker="o", ms=3, color=ORANGE)),
                        Line2D([], [], color=MOSS, lw=1.3, marker="o", ms=3, mec="white", mew=0.4), Line2D([], [], color=BURGUNDY, lw=1.3, marker="o", ms=3, mec="white", mew=0.4)],
               labels=["one human chain (network)", "validation proteins, MD (median)", "validation proteins, network (median)"], handler_map={tuple: HandlerTuple(ndivide=None, pad=0.3)},
               loc="outside lower center", ncol=3, handletextpad=0.5, columnspacing=1.6, fontsize=6.5)
    for ext in ("png", "pdf"):
        fig.savefig(RES / f"populations.{ext}", dpi=300)
    plt.close(fig)
    print(f"-> {RES / 'populations'}.png / .pdf", file=sys.stderr)


def load_proteome():
    Q = proteome_chains(False, "none").copy().reset_index(drop=True)
    N = Q.N.to_numpy(float)
    for tag in ("net", "cf"):
        Q[f"{tag}_pc"], Q[f"{tag}_pe"] = D.populations(Q[f"{tag}_lnA"].to_numpy(), Q[f"{tag}_nu"].to_numpy(), N)
        Q[f"{tag}_compaction"] = residual(N, np.log(Q[f"{tag}_mean_rg"].to_numpy()), 40)
        Q[f"{tag}_compact"] = residual(N, np.log(Q[f"{tag}_pc"].to_numpy()), 40)
        Q[f"{tag}_extended"] = residual(N, np.log(Q[f"{tag}_pe"].to_numpy()), 40)
    return Q


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--frames", default=str(HERE.parent / "data" / "idp_rg_frames.npy"))
    ap.add_argument("--boots", type=int, default=1000)
    ap.add_argument("--figure-only", action="store_true", help="redraw the figure without the scoring")
    args = ap.parse_args(argv)
    set_style()
    rng = np.random.default_rng(2026)

    P, F = data.proteins(), data.frames(args.frames)
    floor = 1 / F.shape[1]                 # one frame: a tail with no stored frame still has a logarithm
    tr, ev = data.split(P, "standard")
    T, V = P.iloc[tr], P.iloc[ev].copy().reset_index(drop=True)

    # validation: network, closed form, MD
    fr = np.asarray(F[V.frame_row.to_numpy()], float)
    V["md_mean"] = fr.mean(1)
    V["md_cv"] = fr.std(1, ddof=1) / V.md_mean.to_numpy()
    Nv = V.N.to_numpy(float)
    V["md_pc"], V["md_pe"] = frame_tails(F, V.frame_row.to_numpy(), Nv, lambda n, c: D.population_bounds(n)[0], lambda n, c: D.population_bounds(n)[1], "own", floor)
    params = {}
    for tag, (lnA, lnS, u) in (("net", network.parameters(*network.load(), V)), ("cf", CF.parameters(CF.load(), V))):
        lnA, lnS, u = (np.asarray(a, float) for a in (lnA, lnS, u))
        nu = D.nu_of(u)
        mean, sd = D.moments(lnA, lnS, nu)[:2]
        V[f"{tag}_mean"] = mean
        V[f"{tag}_pc"], V[f"{tag}_pe"] = D.populations(lnA, nu, Nv)
        params[tag] = dict(lnA=lnA, lnS=lnS, nu=nu, mean=mean, cv=sd / mean)

    # the proteome
    Q = load_proteome()
    N = Q.N.to_numpy(float)
    S_, M_ = scatter_data(Q, V)
    if args.figure_only:
        figure(S_, M_)
        return

    # the bounds: refit, and confirm the constants in idprg.distribution
    bounds = fit_bounds(T, F)
    for got, want in ((bounds["lo", 0.05], D.COMPACT_BOUND), (bounds["hi", 0.05], D.EXTENDED_BOUND)):
        assert np.allclose(got, want, atol=1e-3), f"population bounds changed: {got} against {want}"

    conf = Q.tier.to_numpy() == "confident"
    m = dict(n_chains=len(Q), tiers=Q.tier.value_counts().to_dict(), n_val=len(V), bounds=dict(compact=list(D.COMPACT_BOUND), extended=list(D.EXTENDED_BOUND)),
             bounds_at={str(n): [float(v) for v in D.population_bounds(n)] for n in (20, 50, 100, 200, 400)})
    m["distributions"] = {k: {t: quantile_table(Q[f"{t}_{k}"].to_numpy(), conf) for t in ("net", "cf")} for k in METRICS}
    flory_ratio = Q.net_mean_rg / BASELINE(N)
    m["plain_units"] = {k: dict(zip(("p5", "median", "p95"), (float(x) for x in np.quantile(v, [.05, .5, .95])))) for k, v in (("flory_ratio", flory_ratio), ("P_c_pct", 100 * Q.net_pc), ("P_e_pct", 100 * Q.net_pe), ("P_c_pct_closed_form", 100 * Q.cf_pc), ("P_e_pct_closed_form", 100 * Q.cf_pe))}
    m["validation_medians_pct"] = {f"{t_}_{q}": float(100 * V[f"{t_}_{q}"].median()) for t_ in ("md", "net", "cf") for q in ("pc", "pe")}
    m["bias_validation_points"] = {q: float(100 * (V[f"net_{q}"] - V[f"md_{q}"]).median()) for q in ("pc", "pe")}
    rho_ = lambda a_, b_: float(spearmanr(a_, b_).statistic)
    edges = [0.4, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.05, 1.1, 1.2, 1.6]
    bin_medians = lambda D_, cols: {str(k): dict(chains=len(g), **{c: float(g[c].median()) for c in cols}) for k, g in D_.groupby(pd.cut(D_["size"], edges), observed=True)}
    m["scatter"] = dict(rho_size_compact=dict(proteome_network=rho_(S_["size"], S_.compact), validation_md=rho_(M_["size"], M_.compact), validation_network=rho_(M_.net_size, M_.net_compact)),
                        rho_size_extended=dict(proteome_network=rho_(S_["size"], S_.extended), validation_md=rho_(M_["size"], M_.extended), validation_network=rho_(M_.net_size, M_.net_extended)),
                        rho_size_compact_ratio_0p9_to_1p1=dict(proteome_network=rho_(*S_.loc[(S_["size"] - 1).abs() <= 0.1, ["size", "compact"]].to_numpy().T), validation_md=rho_(*M_.loc[(M_["size"] - 1).abs() <= 0.1, ["size", "compact"]].to_numpy().T)),
                        proteome_bins=bin_medians(S_, ("compact", "extended", "N")), validation_bins=bin_medians(M_, ("compact", "extended", "net_compact", "net_extended")),
                        outside_axes=dict(x=int(((S_["size"] < XLIM[0]) | (S_["size"] > XLIM[1])).sum()), compact=int(((S_.compact < YLIM[0]) | (S_.compact > YLIM[1])).sum()), extended=int(((S_.extended < YLIM[0]) | (S_.extended > YLIM[1])).sum())),
                        named={name: dict(chain=c, N=int(r.N), tier=r.tier, flory_ratio=float(r["size"]), compact=float(r.compact), extended=float(r.extended), fcr=float(r.fcr), mean_hydro=float(r.mean_hydro)) for c, name in NAMED for r in [S_.set_index("chain").loc[c]]})
    m["populations_at_length"] = {str(k): dict(chains=len(g), P_c=float(g.net_pc.median()), P_e=float(g.net_pe.median())) for k, g in Q.groupby(pd.cut(Q.N, [19, 30, 50, 100, 200, 300, 490]), observed=True)}
    m["validation_md_at_length"] = {str(k): dict(proteins=len(g), P_c=float(g.md_pc.median()), P_e=float(g.md_pe.median())) for k, g in V.groupby(pd.cut(V.N, [19, 50, 100, 200, 300, 500]), observed=True)}
    flory = Q.net_mean_rg / BASELINE(N)
    m["flory_ratio"] = dict(rho_N=float(spearmanr(flory, N).statistic), median_by_length={str(k): float(v) for k, v in flory.groupby(pd.cut(Q.N, [19, 30, 50, 100, 200, 300, 490]), observed=True).median().items()})
    rho = lambda a, b: float(spearmanr(a, b).statistic)
    m["redundancy_proteome_network"] = dict(compaction_compact=rho(Q.net_compaction, Q.net_compact), compaction_extended=rho(Q.net_compaction, Q.net_extended), compact_extended=rho(Q.net_compact, Q.net_extended))

    # replicate noise and calibration on validation
    groups = replicate_groups(P)
    d = {k: [] for k in METRICS}
    for g in groups:
        pc, pe = frame_tails(F, g.frame_row.to_numpy(), g.N.to_numpy(float), lambda n, c: D.population_bounds(n)[0], lambda n, c: D.population_bounds(n)[1], "own", floor)
        mean = np.asarray(F[g.frame_row.to_numpy()], float).mean(1)
        for i, j in itertools.combinations(range(len(g)), 2):
            d["compaction"].append(np.log(mean[i] / mean[j])), d["compact"].append(np.log(pc[i] / pc[j])), d["extended"].append(np.log(pe[i] / pe[j]))
    noise = {k: float(np.std(v) / 2 ** 0.5) for k, v in d.items()}
    m["replicates"] = dict(pairs=len(d["compact"]), noise_sd_ln=noise)
    col = dict(compaction="mean", compact="pc", extended="pe")
    resid_v = {(t, k): residual(Nv, np.log(V[f"{t}_{col[k]}"].to_numpy()), 20) for t in ("md", "net", "cf") for k in METRICS}
    idx = {c: np.flatnonzero(V.cluster.to_numpy() == c) for c in np.unique(V.cluster)}
    keys = list(idx)
    m["calibration"] = {}
    for k in METRICS:
        md = resid_v["md", k]
        m["calibration"][k] = dict(noise_sd=noise[k], md_sd=float(md.std(ddof=1)))
        for t in ("net", "cf"):
            pred = resid_v[t, k]
            def stats(i):
                true = np.sqrt(max(md[i].var(ddof=1) - noise[k] ** 2, 1e-12))
                return dict(rho=rho(pred[i], md[i]), width_ratio=pred[i].std(ddof=1) / true, slope=np.polyfit(pred[i], md[i], 1)[0])
            est = stats(np.arange(len(V)))
            b = pd.DataFrame([stats(np.concatenate([idx[keys[j]] for j in rng.integers(len(keys), size=len(keys))])) for _ in range(args.boots)])
            m["calibration"][k][t] = {s: dict(value=float(v), ci=[float(x) for x in b[s].quantile([.025, .975])]) for s, v in est.items()}
    m["redundancy_validation_md"] = dict(compaction_compact=rho(resid_v["md", "compaction"], resid_v["md", "compact"]), compaction_extended=rho(resid_v["md", "compaction"], resid_v["md", "extended"]), compact_extended=rho(resid_v["md", "compact"], resid_v["md", "extended"]))

    # the threshold definitions
    S = score_thresholds(definitions(bounds), V, params["net"], Q, F, floor, groups)
    S.to_csv(RES / "populations_thresholds.csv", index=False, float_format="%.4g")
    tails = tail_skill(bounds, V, params["net"], F, floor, rng, args.boots)
    tails.to_csv(RES / "populations_tails.csv", index=False, float_format="%.4g")
    m["fitted_bounds"] = {f"{s}_{100 * q:g}": [float(x) for x in v] for (s, q), v in bounds.items()}

    mid = median_chain(S_)
    m["median_chain"] = dict(chain=mid.chain, gene=mid.gene, N=int(mid.N), flory_ratio=float(mid["size"]), P_c_pct=float(mid.compact), P_e_pct=float(mid.extended))
    md_q, net_q = quadrant(M_.compact, M_.extended), quadrant(M_.net_compact, M_.net_extended)
    deep = (M_.net_compact < 3.5) & (M_.net_extended > 5.5)
    m["phase_diagram"] = dict(proteome_share=pd.Series(quadrant(S_.compact, S_.extended)).value_counts(normalize=True).to_dict(),
                              validation_md_share=pd.Series(md_q).value_counts(normalize=True).to_dict(), validation_network_share=pd.Series(net_q).value_counts(normalize=True).to_dict(),
                              md_quadrant_given_network_quadrant=pd.crosstab(net_q, md_q, normalize="index").to_dict("index"),
                              deep_top_left=dict(chains=int(deep.sum()), md_medians=M_.loc[deep, ["compact", "extended", "size"]].median().to_dict(), network_medians=M_.loc[deep, ["net_compact", "net_extended", "net_size"]].median().to_dict(), md_quadrants=pd.Series(md_q[deep.to_numpy()]).value_counts().to_dict()))
    m["scatter"]["lines"] = {f"{col}_{tag}": dict(zip(("flory_ratio", "population"), (v.tolist() for v in validation_line(M_, x, f"{pre}{col}")))) for col in ("compact", "extended") for tag, x, pre in (("md", "size", ""), ("network", "net_size", "net_"))}
    (RES / "populations.json").write_text(json.dumps(m, indent=1, default=float))
    cols = ["chain", "gene", "description", "N", "tier", "net_mean_rg", "net_pc", "net_pe", "net_compaction", "net_compact", "net_extended", "cf_mean_rg", "cf_pc", "cf_pe", "cf_compaction", "cf_compact", "cf_extended"]
    Q[cols].to_csv(RES / "populations_chains.csv", index=False, float_format="%.5g")
    figure(S_, M_)
    pd.set_option("display.width", 200), pd.set_option("display.max_columns", 40)
    print(S.set_index("definition").T.round(3).to_string())
    print(tails.round(3).to_string(index=False))
    print(json.dumps({k: m[k] for k in ("n_chains", "tiers", "bounds_at", "distributions", "populations_at_length", "validation_md_at_length", "flory_ratio", "redundancy_proteome_network", "redundancy_validation_md", "replicates", "calibration")}, indent=1, default=float))


if __name__ == "__main__":
    main()
