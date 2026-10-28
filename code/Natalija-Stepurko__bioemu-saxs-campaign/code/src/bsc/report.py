"""Stage `report`: figures from results/ and the campaign proposal (docs/CAMPAIGN.md) whose
numbers and candidate table are read from results/analysis.json and results/candidates.csv.
"""
from __future__ import annotations

import json

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from bsc import analyse, envelopes  # noqa: E402
from bsc import config as C  # noqa: E402

INK, MUTED, RULE = "#16191D", "#5B646E", "#DDE1E4"
CLS_COL = {"folded": "#2F5D8A", "partly disordered": "#8A6FA8", "disordered": "#C2681A"}
MODEL_NAME = {"bioemu": "BioEmu-1", "alphafold": "AlphaFold2", "esmfold": "ESMFold", "boltz2": "Boltz-2",
              "idpfold": "IDPFold", "peptron": "PepTron", "boltz1x": "Boltz-1x", "esmflow": "ESMFlow",
              "idpsam": "idpSAM", "idpgan": "idpGAN", "idp-o": "IDP-o"}
# Role palette, used in every figure and on the page: one colour per role, none of them ink. The three
# line roles are checked with the dataviz palette validator on the white figure surface (all pairs pass).
ROLE_COL = {"data": "#8C96A0",              # measured points: light neutral, drawn semi-transparent
            "bioemu": "#1F6FE5",            # BioEmu-1 raw ensemble (and BioEmu-1 wherever it is named)
            "reweighted": "#F05A0A",        # BioEmu-1 reweighted towards the profile, drawn dotted
            "single_structure": "#D6338C"}  # AlphaFold2 / ESMFold reference, drawn dashed
CURVE_LW = {"raw": 2.0, "operating": 2.2, "alphafold": 1.8}
PAGE_ACCENT = "#2F5D8A"                     # links and headings on the page; not a figure role
NEUTRAL_BAR = "#8A939B"
# the model under study in its role colour; comparators in desaturated tones
MODEL_COL = {"bioemu": ROLE_COL["bioemu"], "peptron": "#8A6FA8", "boltz2": "#7A9E7E", "idpfold": "#B08968",
             "idpsam": "#9A8FA3", "alphafold": "#8A939B", "esmfold": "#B9C0C6",
             "boltz1x": "#A9B8A3", "esmflow": "#C9CDD1", "idpgan": "#B5ADB9", "idp-o": "#C4BDC7"}
ACCENT = PAGE_ACCENT
SINGLE_STRUCTURE = ("alphafold", "esmfold")      # reference baselines, one structure each
CLASS_ORDER = [c[0] for c in C.DISORDER_CLASSES]
OUT = C.RESULTS / "figures"


PX_PER_IN = 80                   # figure inches -> CSS pixels at the display width
PIXEL_DENSITY = 3                # PNG pixels per CSS pixel; the page divides the PNG width by this
DPI = PIXEL_DENSITY * PX_PER_IN
LINE_ALPHA = 0.9                 # one transparency for every line plot on the page
RESIDUAL_LW = 0.8
KIND_GREYS = ["#C3C9CF", "#959DA5", "#66707A", "#3A4148"]


def style():
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10.5, "axes.titlesize": 10.5,
                         "axes.labelsize": 10.5, "xtick.labelsize": 9.5, "ytick.labelsize": 9.5,
                         "legend.fontsize": 9.5, "axes.edgecolor": RULE, "axes.labelcolor": INK,
                         "xtick.color": MUTED, "ytick.color": MUTED, "axes.spines.top": False,
                         "axes.spines.right": False, "figure.dpi": DPI, "savefig.dpi": DPI})


def figure(width_px: int, height_px: int, **kw):
    """A figure sized for its display width, so 10 pt text renders at about 11 px on the page."""
    return plt.subplots(figsize=(width_px / PX_PER_IN, height_px / PX_PER_IN), **kw)


def save(fig, name: str) -> None:
    """Full-colour PNG at PIXEL_DENSITY pixels per display pixel (no palette quantisation, which
    roughens anti-aliased lines)."""
    from PIL import Image
    p = OUT / f"{name}.png"
    fig.savefig(p, dpi=DPI, facecolor="white")
    plt.close(fig)
    Image.open(p).convert("RGB").save(p, optimize=True)


def _strip(ax, groups: list[tuple[str, np.ndarray]], seed: int = 0, s: float = 12):
    for i, (cls, v) in enumerate(groups):
        x = i + np.random.default_rng(seed + i).uniform(-0.22, 0.22, len(v))
        ax.scatter(x, v, s=s, color=CLS_COL[cls], alpha=0.6, lw=0)


def fig_error_map(d: pd.DataFrame, A: dict) -> None:
    """Result 1: raw chi2 by class (left) and the error-free NRMSD score by class (right)."""
    m = d[(d["model"] == C.MODEL) & d["clean"]]
    fig, axes = figure(720, 300, ncols=2)
    ax = axes[0]
    groups = [(c, np.log10(m.loc[m["disorder_class"] == c, "chi2_raw"].to_numpy())) for c in CLASS_ORDER]
    _strip(ax, groups)
    for i, (_, v) in enumerate(groups):
        med = np.median(v)
        ax.hlines(med, i - 0.3, i + 0.3, color=INK, lw=2)
        ax.text(i + 0.33, med, f"{10 ** med:.1f}", va="center", fontsize=9.5, color=INK)
    ax.axhline(np.log10(2), color=MUTED, lw=0.8, ls="--")
    ax.text(2.58, np.log10(2) - 0.06, "χ² = 2", fontsize=9.5, color=MUTED, ha="right", va="top")
    ax.set_xticks(range(3))
    ax.set_xticklabels([f"{c}\nn = {len(v)}" for c, v in groups])
    ax.set_xlim(-0.5, 2.6)
    ax.set_ylim(-0.6, 3.4)
    ax.set_yticks([0, 1, 2, 3])
    ax.set_yticklabels(["1", "10", "100", "1000"])
    ax.set_ylabel("raw reduced χ² (log scale)")
    ax.set_title("Reduced χ², weighted by the reported errors", loc="left")
    ax = axes[1]
    nr = A["headline_robustness"].get("nrmsd", {}).get("by_class")
    if nr and "nrmsd_raw" in m:
        groups = [(c, m.loc[m["disorder_class"] == c, "nrmsd_raw"].dropna().to_numpy()) for c in CLASS_ORDER]
        _strip(ax, groups, seed=10)
        for i, cls in enumerate(CLASS_ORDER):
            st = nr[cls]
            ax.hlines(st["nrmsd_raw_median"], i - 0.3, i + 0.3, color=INK, lw=2)
            ax.vlines(i + 0.36, *st["nrmsd_raw_ci95"], color=INK, lw=1.6)
        ax.set_yscale("log")
        ax.set_ylim(0.005, 0.35)
        ax.set_yticks([0.01, 0.03, 0.1, 0.3])
        ax.set_yticklabels(["0.01", "0.03", "0.1", "0.3"])
        ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        ax.set_xticks(range(3))
        ax.set_xticklabels([f"{c}\nn = {len(v)} · {nr[c]['nrmsd_raw_median']:.3f}" for c, v in groups])
        ax.set_xlim(-0.5, 2.55)
        ax.set_ylabel("NRMSD of ln I (log scale)")
        ax.set_title("NRMSD, no error weighting (bar: 95% interval)", loc="left")
    else:
        ax.text(0.5, 0.5, "NRMSD not computed", ha="center", va="center", transform=ax.transAxes, color=MUTED)
        ax.set_axis_off()
    fig.tight_layout(w_pad=2.5)
    save(fig, "fig_error_map")


def fig_length(d: pd.DataFrame, A: dict) -> None:
    """E2: raw chi2 against chain length, coloured by class."""
    m = d[(d["model"] == C.MODEL) & d["clean"]]
    fig, ax = figure(560, 300)
    for cls in CLASS_ORDER:
        g = m[m["disorder_class"] == cls]
        ax.scatter(g["length"], np.log10(g["chi2_raw"]), s=12, color=CLS_COL[cls], alpha=0.7, lw=0, label=cls)
    ax.axhline(np.log10(2), color=MUTED, lw=0.8, ls="--")
    ax.set_xscale("log")
    ax.set_xticks([20, 50, 100, 200, 500])
    ax.set_xticklabels(["20", "50", "100", "200", "500"])
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_yticks([0, 1, 2, 3])
    ax.set_yticklabels(["1", "10", "100", "1000"])
    ax.set_xlabel("chain length (residues)")
    ax.set_ylabel("raw reduced χ² (log scale)")
    ax.set_ylim(-0.6, 3.4)
    ax.legend(frameon=False, loc="upper left")
    e2 = A["expectations"]["E2"]
    ax.set_title(f"Folded proteins: Spearman ρ = {e2['spearman_rho']:+.2f}, p = {e2['p']:.2g}", loc="left")
    fig.tight_layout()
    save(fig, "fig_length")


def fig_rg(d: pd.DataFrame) -> None:
    """Ensemble Rg against experimental Guinier Rg."""
    m = d[(d["model"] == C.MODEL) & d["clean"]].dropna(subset=["rg_exp", "rg_model"])
    fig, ax = figure(480, 420)
    lim = [0.8 * min(m["rg_exp"].min(), m["rg_model"].min()), 1.1 * max(m["rg_exp"].max(), m["rg_model"].max())]
    ax.plot(lim, lim, color=MUTED, lw=0.8)
    for cls in CLASS_ORDER:
        g = m[m["disorder_class"] == cls]
        ax.scatter(g["rg_exp"], g["rg_model"], s=14, color=CLS_COL[cls], alpha=0.75, lw=0, label=cls)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ticks = [t for t in (10, 20, 30, 50, 100, 150) if lim[0] <= t <= lim[1]]
    for axis in (ax.xaxis, ax.yaxis):
        axis.set_major_locator(matplotlib.ticker.FixedLocator(ticks))
        axis.set_major_formatter(matplotlib.ticker.FixedFormatter([str(t) for t in ticks]))
        axis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_xlabel("experimental Rg (Å, Guinier)")
    ax.set_ylabel("BioEmu-1 ensemble Rg (Å)")
    ax.legend(frameon=False, loc="upper left")
    ax.set_title("Above the diagonal: ensemble larger than measured", loc="left")
    fig.tight_layout()
    save(fig, "fig_rg")


def fig_rg_robustness(d: pd.DataFrame, A: dict) -> None:
    """Result 2: Rg ratio per class from the curves (left) and curve against Cα coordinates (right)."""
    R = A["rg_robustness"]
    m = d[(d["model"] == C.MODEL) & d["clean"]].dropna(subset=["rg_ratio"])
    fig, axes = figure(720, 340, ncols=2, gridspec_kw={"width_ratios": [1.05, 1]})
    yt = [0.6, 0.8, 1.0, 1.2, 1.5, 2.0, 2.5]

    def frame(ax):
        ax.set_yscale("log")
        ax.set_yticks(yt)
        ax.set_yticklabels([str(t) for t in yt])
        ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        ax.set_ylim(0.58, 2.7)
        ax.axhline(1, color=MUTED, lw=0.8, ls="--")

    ax = axes[0]
    frame(ax)
    xl = []
    for i, cls in enumerate(CLASS_ORDER):
        g = m[m["disorder_class"] == cls]
        x = i + np.random.default_rng(i).uniform(-0.22, 0.22, len(g))
        ax.scatter(x, g["rg_ratio"], s=12 if cls != "disordered" else 18, color=CLS_COL[cls], alpha=0.55, lw=0)
        st = R["curve_rg_ratio_by_class"][cls]
        ax.plot([i - 0.3, i + 0.3], [st["median"]] * 2, color=INK, lw=2)
        ax.plot([i + 0.36, i + 0.36], st["ci95"], color=INK, lw=1.6)
        xl.append(f"{cls}\nn = {st['n']} · {st['median']:.2f}")
    ax.set_xticks(range(3))
    ax.set_xticklabels(xl, fontsize=9.5)
    ax.set_xlim(-0.5, 2.55)
    ax.plot([-0.42, -0.22], [2.42, 2.42], color=MUTED, lw=0.8, ls="--")
    ax.text(-0.18, 2.42, "ratio 1: predicted size matches experiment", fontsize=9.5, color=MUTED, va="center")
    ax.annotate("Expected: over-compaction\nObserved: extension", xy=(2, R["curve_rg_ratio_by_class"]["disordered"]["median"]),
                xytext=(0.62, 1.9), fontsize=9.5, color=INK, ha="left", va="center",
                arrowprops={"arrowstyle": "-|>", "color": MUTED, "lw": 0.8})
    ax.set_ylabel("ensemble Rg / experimental Rg")
    ax.set_title("From the scattering curves (bar: 95% interval)", loc="left")
    ax = axes[1]
    cr = R.get("coordinate_rg", {}).get("by_class", {})
    p = C.RESULTS / "coord_rg.csv"
    if cr and p.exists():
        frame(ax)
        tab = pd.read_csv(p)
        tab = tab[tab["label"].isin(m["label"]) & (tab.get("reason", "") != "worked example")]
        pos, ticks, labels = 0.0, [], []
        for cls in ("folded", "disordered"):
            if cls not in cr:
                continue
            g = tab[tab["disorder_class"] == cls]
            for col, nm, key in (("ratio_curve_exp", "curve", "curve_over_exp"),
                                 ("ratio_coord_exp", "Cα coords", "coord_over_exp")):
                x = pos + np.random.default_rng(int(10 * pos)).uniform(-0.22, 0.22, len(g))
                ax.scatter(x, g[col], s=14, color=CLS_COL[cls], alpha=0.55 if nm == "curve" else 0.35, lw=0,
                           marker="o" if nm == "curve" else "D")
                ax.plot([pos - 0.3, pos + 0.3], [cr[cls][f"{key}_median"]] * 2, color=INK, lw=2)
                ax.plot([pos + 0.36, pos + 0.36], cr[cls][f"{key}_ci95"], color=INK, lw=1.6)
                ticks.append(pos)
                labels.append(f"{cls}\n{nm}\n{cr[cls][f'{key}_median']:.2f}")
                pos += 1.0
            pos += 0.5
        ax.set_xticks(ticks)
        ax.set_xticklabels(labels, fontsize=9.5)
        ax.set_xlim(-0.5, pos - 0.6)
        ax.set_title(f"Curve and conformer coordinates ({cr['disordered']['n']} + {cr['folded']['n']})", loc="left")
    else:
        ax.text(0.5, 0.5, "coordinate Rg not computed\n(run `bsc coordrg`)", ha="center", va="center",
                transform=ax.transAxes, color=MUTED)
        ax.set_axis_off()
    fig.tight_layout(w_pad=2.0)
    save(fig, "fig_rg_robustness")


def fig_kinds(d: pd.DataFrame, paths: list[dict], A: dict) -> None:
    """Result 3: the four kinds as one stacked bar, and one representative reweighting path."""
    m = d[(d["model"] == C.MODEL) & d["clean"]]
    counts = m["resolvability"].value_counts().reindex(analyse.KINDS, fill_value=0)
    n = int(counts.sum())
    fig, axes = figure(720, 220, ncols=2, gridspec_kw={"width_ratios": [1.6, 1]})
    ax = axes[0]
    left = 0.0
    names = ["raw fit", "modestly\nreweightable", "strongly\nreweightable", "target fit\nnot reached"]
    for k, (kind, c, nm) in enumerate(zip(analyse.KINDS, KIND_GREYS, names, strict=True)):
        w = counts[kind] / n
        ax.barh(0, w, left=left, color=c, height=0.5, edgecolor="white", lw=2)
        ax.text(left + w / 2, 0, f"{counts[kind]}\n{100 * w:.0f}%", ha="center", va="center", fontsize=9.5,
                color=INK if k < 2 else "white")
        if k == 1:      # the narrow segment's name goes below the bar so neighbouring names do not collide
            ax.text(left + w / 2, -0.33, nm, ha="center", va="top", fontsize=9.5, color=INK)
        else:
            ax.text(left + w / 2, 0.33, nm, ha="center", va="bottom", fontsize=9.5, color=INK)
        left += w
    ax.set_xlim(0, 1)
    ax.set_ylim(-1.0, 0.75)
    ax.set_axis_off()
    ax.text(0, -0.82, f"{n} profiles · kinds read where the path first reaches χ² ≤ 2; φ threshold 0.5",
            fontsize=9.5, color=MUTED, va="top")
    ax = axes[1]
    lab = A.get("representative_path")
    P = pd.DataFrame(paths)
    g = P[P["label"] == lab].sort_values("theta", ascending=False) if lab else P.iloc[:0]
    if len(g):
        ax.plot(g["phi"], g["chi2"], color=ROLE_COL["reweighted"], lw=1.6, alpha=LINE_ALPHA,
                marker="o", ms=3)
        ax.axhline(2, color=MUTED, lw=0.8, ls="--")
        hit = g[g["chi2"] <= 2]
        if len(hit):
            ph = hit.iloc[0]
            ax.plot([ph["phi"]], [ph["chi2"]], "o", ms=7, mfc="white", mec=ROLE_COL["reweighted"], mew=1.4)
            ax.annotate(f"χ² ≤ 2 at φ = {ph['phi']:.2f}", xy=(ph["phi"], ph["chi2"]), xytext=(0.6, 0.62),
                        textcoords="axes fraction", fontsize=9.5, arrowprops={"arrowstyle": "-", "color": MUTED, "lw": 0.8})
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(1.05, max(0.005, 0.7 * g["phi"].min()))
        ax.set_xlabel("φ kept (← more reweighting)")
        ax.set_ylabel("reduced χ²")
        ax.set_title(f"{lab}, strongly reweightable", loc="left")
        ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
        yt = [t for t in (1, 1.5, 2, 3, 5, 10, 20, 50, 100) if g["chi2"].min() * 0.8 <= t <= g["chi2"].max() * 1.2]
        ax.yaxis.set_major_locator(matplotlib.ticker.FixedLocator(yt))
        ax.yaxis.set_major_formatter(matplotlib.ticker.FixedFormatter([f"{t:g}" for t in yt]))
        ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    else:
        ax.set_axis_off()
    fig.tight_layout(w_pad=1.5)
    save(fig, "fig_kinds")


def fig_resolvability(d: pd.DataFrame, paths: list[dict]) -> None:
    """Every reweighting path (chi2 against effective sample fraction) and the kinds by class."""
    m = d[(d["model"] == C.MODEL) & d["clean"]].set_index("label")
    fig, axes = figure(720, 300, ncols=2, gridspec_kw={"width_ratios": [1.5, 1]})
    ax = axes[0]
    P = pd.DataFrame(paths)
    for lab, g in P.groupby("label"):
        if lab not in m.index:
            continue
        g = g.sort_values("phi", ascending=False)
        ax.plot(g["phi"], g["chi2"], color=CLS_COL[m.loc[lab, "disorder_class"]], alpha=0.25, lw=0.8)
    ax.axhline(2, color=MUTED, lw=0.8, ls="--")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(1.05, 0.005)
    ax.set_ylim(0.5, 300)
    for axis in (ax.xaxis, ax.yaxis):
        axis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
        axis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_xlabel("φ kept (← more reweighting)")
    ax.set_ylabel("reduced χ²")
    ax.set_title("Every path, coloured by class", loc="left")
    ax = axes[1]
    tab = (m.groupby(["disorder_class", "resolvability"]).size().unstack(fill_value=0)
           .reindex(index=CLASS_ORDER, columns=analyse.KINDS, fill_value=0))
    share = tab.div(tab.sum(axis=1), axis=0)
    left = np.zeros(3)
    for col, c in zip(analyse.KINDS, KIND_GREYS, strict=True):
        ax.barh(range(3), share[col], left=left, color=c, label=col, height=0.6, edgecolor="white", lw=1.5)
        left += share[col].to_numpy()
    ax.set_yticks(range(3))
    ax.set_yticklabels(["folded", "partly\ndisordered", "disordered"])
    ax.set_xlim(0, 1)
    ax.set_xlabel("share of entries")
    ax.legend(frameon=False, fontsize=9.5, ncol=2, loc="lower center", bbox_to_anchor=(0.4, 1.0))
    ax.invert_yaxis()
    fig.tight_layout()
    save(fig, "fig_resolvability")


def fig_models(A: dict) -> None:
    """Model comparison: median raw chi2 by class per model (left) and the disordered Rg ratio (right);
    BioEmu-1 and PepTron emphasised, single-structure predictors grouped below."""
    comps = [r["model"] for r in A["model_comparison"]]
    gens = sorted([m for m in comps if m not in SINGLE_STRUCTURE and m != "peptron"])
    rows = [C.MODEL] + (["peptron"] if "peptron" in comps else []) + gens + [None] + \
        [m for m in SINGLE_STRUCTURE if m in comps]
    ypos, y = {}, 0.0
    for mo in rows:
        if mo is None:
            y += 0.6
            continue
        ypos[mo] = y
        y += 1
    bc = pd.DataFrame(A["by_class"])
    rgd = A["rg_ratio_disordered_by_model"]
    fig, axes = figure(720, 300, ncols=2, sharey=True, gridspec_kw={"width_ratios": [1.45, 1]})
    strong = {C.MODEL, "peptron"}
    marker = {"folded": "o", "partly disordered": "s", "disordered": "^"}
    ax = axes[0]
    for mo, yy in ypos.items():
        g = bc[bc["model"] == mo].set_index("disorder_class")
        vals = [g.loc[c, "chi2_raw_median"] for c in CLASS_ORDER if c in g.index]
        ax.plot([min(vals), max(vals)], [yy, yy], color=RULE, lw=1.2, zorder=1)
        for cls in CLASS_ORDER:
            if cls in g.index:
                ax.scatter(g.loc[cls, "chi2_raw_median"], yy, marker=marker[cls], s=42 if mo in strong else 28,
                           color=CLS_COL[cls], alpha=1.0 if mo in strong else 0.45, lw=0, zorder=2,
                           label=cls if mo == C.MODEL else None)
    ax.axvline(2, color=MUTED, lw=0.8, ls="--")
    ax.set_xscale("log")
    ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
    ax.set_xlabel("median raw reduced χ² (log scale)")
    ax.set_yticks(list(ypos.values()))
    ax.set_yticklabels([MODEL_NAME.get(mo, mo) for mo in ypos])
    for t, mo in zip(ax.get_yticklabels(), ypos, strict=True):
        t.set_color(INK if mo in strong else MUTED)
        t.set_fontweight("bold" if mo in strong else "normal")
    ax.invert_yaxis()
    ax.legend(frameon=False, fontsize=9.5, loc="lower right", handletextpad=0.2)
    if SINGLE_STRUCTURE[0] in ypos:
        ax.text(ax.get_xlim()[0], ypos[SINGLE_STRUCTURE[0]] - 0.8, " single-structure predictors", fontsize=9.5,
                color=MUTED, va="center")
    ax.set_title("Raw fit by class", loc="left")
    ax = axes[1]
    for mo, yy in ypos.items():
        if mo not in rgd:
            continue
        ax.plot([1, rgd[mo]], [yy, yy], color=RULE, lw=1.2, zorder=1)
        dot = (MODEL_COL.get(mo, MUTED) if mo in strong
               else ROLE_COL["single_structure"] if mo in SINGLE_STRUCTURE else MUTED)
        ax.scatter(rgd[mo], yy, s=42 if mo in strong else 28, color=dot, alpha=1.0 if mo in strong else 0.6,
                   lw=0, zorder=2)
        ax.text(rgd[mo] + (0.03 if rgd[mo] >= 1 else -0.03), yy, f"{rgd[mo]:.2f}", va="center", fontsize=9.5,
                ha="left" if rgd[mo] >= 1 else "right", color=INK if mo in strong else MUTED, zorder=3,
                bbox={"facecolor": "white", "edgecolor": "none", "pad": 0.3})
    ax.axvline(1, color=MUTED, lw=0.8)
    ax.axvline(1.1, color=MUTED, lw=0.8, ls="--")
    ax.text(1.12, -0.55, "1.1", fontsize=9.5, color=MUTED, va="center")
    ax.set_xlim(0.2, 1.5)
    ax.set_xlabel("disordered: median Rg ratio")
    ax.set_title("Size of disordered ensembles", loc="left")
    fig.tight_layout(w_pad=1.0)
    save(fig, "fig_models")




def fig_predict() -> None:
    """Result 4: held-out prediction of the raw error (A) and the high-discrepancy target yield of
    each acquisition policy among the top 10 (B)."""
    pp, sp = C.RESULTS / "predict.json", C.RESULTS / "select_eval.json"
    if not (pp.exists() and sp.exists()):
        print("  fig_predict: predict.json or select_eval.json missing; skipped", flush=True)
        return
    P, S = json.load(open(pp)), json.load(open(sp))
    oof = pd.read_csv(C.RESULTS / "predict_oof.csv")
    fig, axes = figure(720, 320, ncols=2, gridspec_kw={"width_ratios": [1, 1.15]})
    ax = axes[0]
    best = P["best_regressor"]
    for cls in CLASS_ORDER:
        g = oof[oof["disorder_class"] == cls]
        ax.scatter(g[f"pred_{best}"], g["y_log10_chi2"], s=11, color=CLS_COL[cls], alpha=0.65, lw=0, label=cls)
    lo, hi = oof["y_log10_chi2"].min() - 0.1, oof["y_log10_chi2"].max() + 0.1
    ax.plot([lo, hi], [lo, hi], color=MUTED, lw=0.8)
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_aspect("equal")
    ax.set_xlabel("predicted log₁₀ raw χ² (held out)")
    ax.set_ylabel("observed log₁₀ raw χ²")
    r2 = P["regression"][best]["r2"]["mean"]
    ax.set_title(f"A · the model explains {100 * r2:.0f}% of held-out variation", loc="left")
    ax.legend(frameon=False, fontsize=9.5, loc="lower right", handletextpad=0.1)
    ax = axes[1]
    k = str(S["top_k"][0])
    pol = S["policies"][k]
    rnd_sd = S["random_single_draw_sd"][k]["strong_or_notfit"]
    bars = [("random", pol["random"]["strong_or_notfit"]["mean"], rnd_sd, RULE),
            ("predicted\nerror", pol["predicted error"]["strong_or_notfit"]["mean"],
             pol["predicted error"]["strong_or_notfit"]["sd"], NEUTRAL_BAR),
            ("acquisition\npolicy fed\npredictions", pol["predicted priority"]["strong_or_notfit"]["mean"],
             pol["predicted priority"]["strong_or_notfit"]["sd"], ACCENT),
            ("oracle\nranking", pol["observed priority (ceiling)"]["strong_or_notfit"]["mean"],
             pol["observed priority (ceiling)"]["strong_or_notfit"]["sd"], "white")]
    xs = [0, 1, 2, 3.5]
    for x, (nm, v, e, c) in zip(xs, bars, strict=True):
        ax.bar(x, v, width=0.68, color=c, edgecolor=INK if nm.startswith("oracle") else "none", lw=1,
               hatch="///" if nm.startswith("oracle") else None)
        ax.errorbar(x, v, yerr=e, fmt="none", ecolor=INK, lw=0.9, capsize=3)
        ax.text(x, v + e + 0.03, f"{v:.2f}", ha="center", fontsize=9.5, color=INK)
    ax.axvline(2.75, color=MUTED, lw=0.8, ls=":")
    ax.text(3.5, 1.22, "uses measured\nSAXS; not\ndeployable", ha="center", fontsize=9.5, color=MUTED, va="center")
    ax.axhline(S["pool_mean"]["strong_or_notfit"], color=MUTED, lw=0.8, ls="--")
    ax.text(2.8, S["pool_mean"]["strong_or_notfit"] + 0.02, "pool", ha="left", va="bottom", fontsize=9.5, color=MUTED)
    ax.set_xticks(xs)
    ax.set_xticklabels([b[0] for b in bars], fontsize=9.5)
    ax.set_ylim(0, 1.42)
    ax.set_xlim(-0.6, 4.2)
    ax.set_ylabel(f"high-discrepancy target yield, top {k}")
    ax.set_title("B · what each policy picks, held-out folds", loc="left")
    fig.tight_layout(w_pad=1.5)
    save(fig, "fig_predict")


def fig_predict_detail() -> None:
    """Calibration of the tree classifier and permutation importance of the features (details)."""
    pp = C.RESULTS / "predict.json"
    if not pp.exists():
        return
    P = json.load(open(pp))
    fig, axes = figure(720, 300, ncols=2)
    ax = axes[0]
    cal = P["calibration_hgb"]
    ax.plot([0, 1], [0, 1], color=MUTED, lw=0.8)
    ax.plot([c["mean_predicted"] for c in cal], [c["observed_rate"] for c in cal], "o-", color=ACCENT, ms=5, alpha=LINE_ALPHA)
    for c in cal:
        ax.text(c["mean_predicted"], c["observed_rate"] + 0.04, f"n={c['n']}", ha="center", fontsize=9.5, color=MUTED)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel("predicted probability (quintile means)")
    ax.set_ylabel("observed high-discrepancy share")
    c = P["classification"]
    ax.set_title(f"Calibration · AUC {c['hgb']['auc']['mean']:.2f} ± {c['hgb']['auc']['sd']:.2f}", loc="left")
    ax = axes[1]
    imp = P["permutation_importance_hgb"][:8][::-1]
    ax.barh(range(len(imp)), [i["mean"] for i in imp], xerr=[i["sd"] for i in imp], color=NEUTRAL_BAR,
            height=0.6, error_kw={"lw": 0.8, "ecolor": INK})
    ax.set_yticks(range(len(imp)))
    ax.set_yticklabels([i["feature"].replace("_", " ") for i in imp], fontsize=9.5)
    ax.axvline(0, color=INK, lw=0.8)
    ax.set_xlabel("drop in held-out R² when permuted")
    ax.set_title("Permutation importance (trees)", loc="left")
    fig.tight_layout()
    save(fig, "fig_predict_detail")


SILHOUETTE = {"particle": "#DDE3E9", "protein": "#AAB7C5"}    # neutral grey-blue: lighter outside, darker inside
TRACE_COL = "#3A4148"


def _silhouette(ax, xy: np.ndarray, colour: str, lo: np.ndarray, hi: np.ndarray, step: float) -> None:
    """Filled projection of a set of map points (spaced `step` Å apart) onto the drawing plane, with
    a smoothed outline."""
    from scipy import ndimage
    bx = np.arange(lo[0], hi[0] + step, step)
    by = np.arange(lo[1], hi[1] + step, step)
    h, _, _ = np.histogram2d(xy[:, 0], xy[:, 1], bins=[bx, by])
    occ = ndimage.gaussian_filter(ndimage.binary_closing(h > 0, iterations=2).astype(float), 0.8)
    cx, cy = 0.5 * (bx[1:] + bx[:-1]), 0.5 * (by[1:] + by[:-1])
    ax.contourf(cx, cy, occ.T, levels=[0.5, 2.0], colors=[colour], zorder=1)


def fig_examples(A: dict, ent: pd.DataFrame) -> dict | None:
    """Figure 2: one column per class. The measured profile (every point, with its error as a band),
    the raw and reweighted BioEmu-1 curves and the AlphaFold2 curve under one shared legend, with the
    chi2 values in a second title line; residuals; and the DENSS envelope as projected silhouettes of
    two density levels with the C-alpha trace of the highest-weight conformer docked inside."""
    stats_path = C.RESULTS / "denss_stats.json"
    if not stats_path.exists():
        print("  fig_examples: no results/denss_stats.json (run `bsc envelopes`); tracked figure kept", flush=True)
        return None
    stats = json.load(open(stats_path))
    examples = A["examples"]
    need = [C.RESULTS / stats[lab]["avg_map"] for lab in examples.values() if lab in stats]
    if len(need) < len(examples) or not all(p.exists() for p in need) or not C.SAXS_DIR.exists():
        print("  fig_examples: DENSS maps or input curves missing; tracked figure kept", flush=True)
        return None
    e = ent.set_index("label")
    W, H = 860, 650
    fig = plt.figure(figsize=(W / PX_PER_IN, H / PX_PER_IN))
    gs = fig.add_gridspec(2, 3, height_ratios=[2.1, 0.85], hspace=0.07, wspace=0.2,
                          left=0.065, right=0.99, top=0.862, bottom=0.44)
    gs3 = fig.add_gridspec(1, 3, wspace=0.06, left=0.02, right=0.99, top=0.34, bottom=0.085)
    data_style = {"facecolors": "none", "edgecolors": ROLE_COL["data"], "alpha": 0.25, "linewidths": 0.5}
    col = {"raw": ROLE_COL["bioemu"], "operating": ROLE_COL["reweighted"], "alphafold": ROLE_COL["single_structure"]}
    style = {"raw": "-", "operating": ":", "alphafold": "--"}
    out = {}
    for j, (cls, label) in enumerate(examples.items()):
        f = envelopes.example_fit(label)
        st = stats[label]
        ax = fig.add_subplot(gs[0, j])
        axr = fig.add_subplot(gs[1, j], sharex=ax)
        q, I, s = f["q"], f["I"], f["sigma"]
        pos = I > 0
        y0 = np.log(I[pos]).min() - 0.08 * np.ptp(np.log(I[pos]))
        lo = np.where(I - s > 0, np.log(np.clip(I - s, 1e-30, None)), y0)
        hi = np.log(np.clip(I + s, 1e-30, None))
        ax.fill_between(q, np.maximum(lo, y0), np.maximum(hi, y0), color="#D3D9DF", lw=0, zorder=0)
        ax.scatter(q[pos], np.log(I[pos]), s=6, zorder=1, **data_style)
        for k, name in enumerate(("alphafold", "raw", "operating")):
            fit = f["curves"][name]["fit"]
            okf = fit > 0
            ax.plot(q[okf], np.log(fit[okf]), color=col[name], lw=CURVE_LW[name], ls=style[name], alpha=LINE_ALPHA,
                    zorder=3 + k)
            axr.plot(q, (fit - I) / s, color=col[name], lw=RESIDUAL_LW, ls="-", alpha=LINE_ALPHA, zorder=3 + k)
        c2 = f["curves"]
        rew = (f"reweighted {c2['operating']['chi2']:.2f} (φ {f['phi_operating']:.2f})" if f["operating_reached"]
               else f"χ² minimum {c2['operating']['chi2']:.2f} (φ {f['phi_operating']:.2f})")
        ax.set_title(f"{cls}: {label} · {int(e.loc[label, 'length'])} residues", loc="left", fontsize=10, pad=30)
        ax.text(0, 1.015, f"χ² raw {c2['raw']['chi2']:.2f} · AlphaFold2 {c2['alphafold']['chi2']:.1f}\n{rew}",
                transform=ax.transAxes, fontsize=8.5, color=MUTED, va="bottom", linespacing=1.25)
        ax.set_ylim(y0, np.log(I[pos]).max() + 0.2)
        plt.setp(ax.get_xticklabels(), visible=False)
        axr.axhline(0, color=MUTED, lw=0.7, zorder=1)
        for g in (-3, 3):
            axr.axhline(g, color="#C3C9CF", lw=0.7, ls="--", zorder=1)
        lim = np.nanmax(np.abs(np.concatenate([(c2[n]["fit"] - I) / s for n in c2])))
        axr.set_ylim(-min(lim * 1.05, 25), min(lim * 1.05, 25))
        axr.set_xlabel("q (Å⁻¹)")
        if j == 0:
            ax.set_ylabel("ln I(q)")
            axr.set_ylabel("residual / σ")
        # envelope: projected silhouettes of the two density levels, C-alpha trace on top
        volume = envelopes.POROD_A3_PER_DA * 1000.0 * float(e.loc[label, "mw_seq_kda"])
        support = float(st.get("support_volume_mean_A3", volume))
        dk = envelopes.docked_example(label, f["top_conformer"], C.RESULTS / st["avg_map"], volume, support)
        ax3 = fig.add_subplot(gs3[0, j])
        allxy = np.vstack([sf["points"][:, :2] for sf in dk["surfaces"]] + [dk["ca"][:, :2]])
        pad = 2 * dk["voxel_A"]
        lo2, hi2 = allxy.min(axis=0) - pad, allxy.max(axis=0) + pad
        for sf in dk["surfaces"]:
            _silhouette(ax3, sf["points"][:, :2], SILHOUETTE[sf["name"]], lo2, hi2, dk["voxel_A"])
        ax3.plot(dk["ca"][:, 0], dk["ca"][:, 1], color=TRACE_COL, lw=0.8, alpha=LINE_ALPHA, zorder=3,
                 solid_joinstyle="round")
        ax3.set_xlim(lo2[0], hi2[0])
        ax3.set_ylim(lo2[1], hi2[1])
        ax3.set_aspect("equal")
        ax3.set_anchor("S")
        ax3.set_axis_off()
        fig.text(gs3[0, j].get_position(fig).x0 + 0.01, 0.05, f"{st.get('resolution_A', 0):.0f} Å envelope · "
                 f"w = {f['top_weight']:.2f} · {100 * dk['inside']:.0f}% of Cα inside", fontsize=9.5, color=MUTED)
        out[label] = {"class": cls, "length": int(e.loc[label, "length"]),
                      "chi2": {k: float(v["chi2"]) for k, v in c2.items()},
                      "chi2_minimum": f["chi2_minimum"], "phi_operating": float(f["phi_operating"]),
                      "phi_minimum": float(f["phi_minimum"]), "operating_is_raw": f["operating_is_raw"],
                      "operating_reached": f["operating_reached"],
                      "rg_exp": float(f["rg_exp"]), "rg_raw": float(f["rg_raw"]),
                      "rg_operating": float(f["rg_operating"]), "rg_minimum": float(f["rg_minimum"]),
                      "top_conformer": f["top_conformer"], "top_weight": f["top_weight"],
                      "ca_inside_envelope": dk["inside"], "ca_inside_support_surface": dk["inside_outer"],
                      "density_levels": [{k: lv[k] for k in ("name", "volume_A3", "iso_over_max", "ca_inside")}
                                         for lv in dk["surfaces"]],
                      "chain_breaks": dk["chain_breaks"], "n_points": int(len(q)), "n_points_positive": int(pos.sum()),
                      "isovalue_volume_A3": volume, "support_volume_A3": support,
                      "denss": {k: st[k] for k in ("dmax_A", "n_maps", "chi2_median", "resolution_A",
                                                  "maps_accepted", "rg_per_map_mean") if k in st}}
    # one shared legend above the three columns
    handles = [matplotlib.lines.Line2D([], [], color=col["raw"], lw=CURVE_LW["raw"], ls=style["raw"], label="BioEmu-1 raw"),
               matplotlib.lines.Line2D([], [], color=col["operating"], lw=CURVE_LW["operating"], ls=style["operating"],
                                       label="reweighted"),
               matplotlib.lines.Line2D([], [], color=col["alphafold"], lw=CURVE_LW["alphafold"], ls=style["alphafold"],
                                       label="AlphaFold2"),
               matplotlib.lines.Line2D([], [], marker="o", ls="none", mfc="none", mec=ROLE_COL["data"], mew=0.7,
                                       ms=4, label="measured"),
               matplotlib.patches.Patch(color="#D3D9DF", label="± σ")]
    fig.legend(handles=handles, loc="upper center", ncol=5, frameon=False, fontsize=9.5,
               bbox_to_anchor=(0.53, 1.0), handlelength=2.6, columnspacing=1.6)
    # density key
    x = 0.03
    fig.text(x, 0.02, "envelope density", fontsize=9.5, color=MUTED, va="center")
    for name, label_txt in (("particle", "lower: particle volume"), ("protein", "higher: protein volume")):
        x += 0.13 if name == "particle" else 0.2
        fig.add_artist(matplotlib.patches.Rectangle((x, 0.008), 0.022, 0.025, transform=fig.transFigure,
                                                    facecolor=SILHOUETTE[name], edgecolor=RULE, lw=0.6))
        fig.text(x + 0.028, 0.02, label_txt, fontsize=9.5, color=INK, va="center")
    fig.add_artist(matplotlib.lines.Line2D([0.61, 0.645], [0.02, 0.02], transform=fig.transFigure, color=TRACE_COL,
                                           lw=0.8))
    fig.text(0.652, 0.02, "Cα trace of the top-weight conformer", fontsize=9.5, color=INK, va="center")
    save(fig, "fig_examples")
    json.dump(out, open(C.RESULTS / "examples.json", "w"), indent=1)
    return out


def pct(x: float) -> str:
    return f"{100 * x:.0f}%"


def f2(x) -> str:
    return "n/a" if x is None or not np.isfinite(x) else f"{x:.2f}"


ARM1_MAX_LENGTH = 350      # residues; constructs above this are harder to express and to purify monomeric
ARM1_LISTED = 12           # systems listed; the order is by observed raw chi2, a sort key, not a value of measuring


def arm1_eligible(cand: pd.DataFrame) -> pd.DataFrame:
    """Arm 1 eligibility over the clean candidates (monomers by measured mass): reweighting reaches chi2 <= 2,
    at most ARM1_MAX_LENGTH residues, folded or partly disordered."""
    return cand[cand["resolvability"].isin([analyse.KIND_STRONG, analyse.KIND_MODEST])
                & (cand["length"] <= ARM1_MAX_LENGTH)
                & cand["disorder_class"].isin(["folded", "partly disordered"])]


def campaign_arms(A: dict, cand: pd.DataFrame, d: pd.DataFrame, ent: pd.DataFrame) -> dict:
    """The two arms of the campaign, from results only.

    Arm 1 (model improvement): high-error monomers that reweighting can fit, re-measured under one
    standard condition. Arm 2 (hypothesis test): disordered proteins whose BioEmu-1 ensemble Rg exceeds
    the disordered scaling-law Rg by more than 10% (a pre-acquisition quantity), spanning length and
    composition, with matched controls whose ratio is within 10% of 1."""
    from bsc import features
    elig = arm1_eligible(cand)
    arm1 = elig.sort_values("chi2_raw", ascending=False).head(ARM1_LISTED)
    m = d[(d["model"] == C.MODEL) & d["clean"]].copy()
    m["kind"] = m["resolvability"]
    e = ent.set_index("label")
    m["disorder_mean"] = e.loc[m["label"], "disorder_mean"].to_numpy()
    sc = pd.read_csv(C.RESULTS / "scores.csv")
    sc = sc[sc["model"] == C.MODEL].set_index("label")
    for c in ("rg_conformer_iqr", "rg_conformer_median"):
        m[c] = sc.loc[m["label"], c].to_numpy()
    tab = features.build(m, ent).set_index("label")
    dis = tab[tab["disorder_class"] == "disordered"].copy()
    dis["ratio_pred"] = dis["rg_model_over_law_disordered"]
    dis["ratio_obs"] = m.set_index("label").loc[dis.index, "rg_ratio"]
    dis["length_bin"] = pd.cut(dis["length"], [0, 100, 200, 350, 10_000],
                               labels=["≤100", "101–200", "201–350", ">350"])
    test = dis[dis["ratio_pred"] > 1.1].sort_values(["length_bin", "ncpr"])
    ctrl = dis[(dis["ratio_pred"] - 1).abs() <= 0.1].sort_values(["length_bin", "ncpr"])
    # one test entry per length bin and charge tercile where available, up to twelve; controls matched by bin
    picks = []
    for _, g in test.groupby("length_bin", observed=True):
        g = g.sort_values("ncpr")
        idx = np.unique(np.linspace(0, len(g) - 1, min(3, len(g))).round().astype(int))
        picks += g.index[idx].tolist()
    test_sel = test.loc[picks].head(12)
    ctrl_picks = []
    for b in test_sel["length_bin"].unique():
        ctrl_picks += ctrl[ctrl["length_bin"] == b].index.tolist()[:2]
    ctrl_sel = ctrl.loc[ctrl_picks]

    def rows(t):
        return [{"label": lab, "length": int(r["length"]), "ncpr": round(float(r["ncpr"]), 3),
                 "fcr": round(float(r["fcr"]), 3), "frac_proline": round(float(r["frac_proline"]), 3),
                 "rg_model": round(float(r["rg_model"]), 1), "rg_law_disordered": round(float(r["rg_law_disordered"]), 1),
                 "ratio_model_over_law": round(float(r["ratio_pred"]), 2),
                 "ratio_model_over_exp_observed": round(float(r["ratio_obs"]), 2)} for lab, r in t.iterrows()]
    P = json.load(open(C.RESULTS / "predict.json")) if (C.RESULTS / "predict.json").exists() else None
    predictive = bool(P and max(P["regression"]["ridge"]["r2"]["mean"], P["regression"]["hgb"]["r2"]["mean"]) >= 0.2)
    out = {"arm1": {"replication": json.loads(arm1.to_json(orient="records")),
                    "n_eligible": int(len(elig)), "max_length": ARM1_MAX_LENGTH,
                    "novel_acquisition_possible": predictive,
                    "predictor_r2": None if P is None else {k: P["regression"][k]["r2"]["mean"]
                                                               for k in ("class_mean", "ridge", "hgb")}},
           "arm2": {"threshold_ratio": 1.1, "n_disordered_pool": int(len(dis)),
                    "n_test_candidates": int(len(test)), "n_control_candidates": int(len(ctrl)),
                    "test": rows(test_sel), "controls": rows(ctrl_sel),
                    "test_observed_median_ratio": float(test["ratio_obs"].median()),
                    "control_observed_median_ratio": float(ctrl["ratio_obs"].median()) if len(ctrl) else None,
                    "outcome": "median Rg_exp / Rg_model across the test arm",
                    "refutes_bias_if": "median Rg_exp / Rg_model >= 0.95 in the test arm (the ensembles are "
                                       "not systematically too extended)",
                    "confirms_bias_if": "median Rg_exp / Rg_model <= 0.9 in the test arm and >= 0.95 in the controls"}}
    json.dump(out, open(C.RESULTS / "campaign_arms.json", "w"), indent=1)
    return out


def write_campaign(A: dict, cand: pd.DataFrame, arms: dict) -> None:
    """docs/CAMPAIGN.md: the proposal, with every number read from the analysis."""
    b = A["bioemu"]
    E = A["expectations"]
    res = b["resolvability_counts"]
    n_bioemu = sum(res.values())
    cls = {r["disorder_class"]: r for r in A["by_class"] if r["model"] == C.MODEL}
    R = A["rg_robustness"]
    H = A["headline_robustness"]
    S = json.load(open(C.RESULTS / "select_eval.json")) if (C.RESULTS / "select_eval.json").exists() else None
    arm1 = arms["arm1"]["replication"]
    rows1 = "\n".join(
        f"| {r['label']} | {int(r['length'])} | {r['disorder_class']} | {r['chi2_raw']:.1f} | {r['chi2_best']:.2f} | "
        f"{r['resolvability']} | {r['rg_direction']} |" for r in arm1)
    rows2 = "\n".join(
        f"| {r['label']} | {r['length']} | {r['ncpr']:+.2f} | {r['fcr']:.2f} | {r['rg_model']} | {r['rg_law_disordered']} | "
        f"{r['ratio_model_over_law']:.2f} |" for r in arms["arm2"]["test"])
    rows2c = "\n".join(
        f"| {r['label']} | {r['length']} | {r['ncpr']:+.2f} | {r['fcr']:.2f} | {r['rg_model']} | {r['rg_law_disordered']} | "
        f"{r['ratio_model_over_law']:.2f} |" for r in arms["arm2"]["controls"])
    cr = R.get("coordinate_rg", {}).get("by_class", {}).get("disordered")
    coord_sentence = (f"The excess is in the conformers: the Cα radius of gyration of the sampled conformers is "
                      f"{cr['coord_over_exp_median']:.2f} times the measured value at the median "
                      f"(95% interval {cr['coord_over_exp_ci95'][0]:.2f}–{cr['coord_over_exp_ci95'][1]:.2f}, n = {cr['n']})."
                      if cr else "")
    pr2 = arms["arm1"]["predictor_r2"]
    if pr2 is None:
        novel = "The error predictor (`bsc predict`) has not been run, so no novel-acquisition list is offered."
    elif not arms["arm1"]["novel_acquisition_possible"]:
        novel = ("A novel-acquisition list (sequences without a SASBDB entry, chosen by the error predictor) is not "
                 f"offered: the predictor of the raw error from sequence-level features reaches R² = "
                 f"{max(pr2['ridge'], pr2['hgb']):.2f} in grouped cross-validation (class means alone "
                 f"{pr2['class_mean']:.2f}), so its ranking of unmeasured sequences would be close to random. "
                 "Only the replication list stands.")
    else:
        novel = ("A novel-acquisition list follows the error predictor's highest predicted error with the largest "
                 "disagreement between models.")
    if S:
        k = str(S["top_k"][0])
        pol = S["policies"][k]
        sel_sentence = (f"On held-out folds of the existing pool, the rule fed with predicted quantities picks "
                        f"{pol['predicted priority']['strong_or_notfit']['mean']:.2f} strongly reweightable or target-not-reached "
                        f"entries per entry chosen (top {k}) against {pol['random']['strong_or_notfit']['mean']:.2f} for "
                        f"random choice and {pol['observed priority (ceiling)']['strong_or_notfit']['mean']:.2f} when the "
                        f"observed quantities are used; the predicted-error policy alone reaches "
                        f"{pol['predicted error']['strong_or_notfit']['mean']:.2f}.")
    else:
        sel_sentence = ""
    text = f"""# Experimental campaign proposal: SAXS data for the next BioEmu

*Generated by `bsc report` from `results/`. The reasoning is fixed in `docs/DESIGN.md` §6; the two-arm
structure and the validation experiments are the post hoc additions of §7; the numbers are read from the
results.*

## 1. What the error map says

Scored without reweighting against {n_bioemu} quality-filtered SASBDB profiles with a BioEmu-1 ensemble, the
ensemble reaches reduced χ² ≤ 2 for {pct(b['share_raw_fit_chi2_2'])} of entries (95% interval
{pct(H['share_raw_fit']['ci95'][0])}–{pct(H['share_raw_fit']['ci95'][1])}; median χ² {b['chi2_raw_median']:.1f},
quartiles {b['chi2_raw_q25']:.1f}–{b['chi2_raw_q75']:.1f}). By class, the median raw χ² is
{cls['folded']['chi2_raw_median']:.1f} for folded proteins (n = {cls['folded']['n']}),
{cls['partly disordered']['chi2_raw_median']:.1f} for partly disordered (n = {cls['partly disordered']['n']}) and
{cls['disordered']['chi2_raw_median']:.1f} for disordered (n = {cls['disordered']['n']}).
Expectation E1 (disordered worse than folded) is {'met' if E['E1']['met'] else 'not met'}
(one-sided p = {E['E1']['p_one_sided']:.2g}); E2 (error rising with length among folded proteins) is
{'met' if E['E2']['met'] else 'not met'} (ρ = {E['E2']['spearman_rho']:+.2f}; folded proteins of up to 100 residues fit
{pct(E['E2']['share_fit_le_100'])} of the time, longer ones {pct(E['E2']['share_fit_gt_100'])}, with no trend above that); E4 (over-compaction of
disordered proteins) came out reversed: the median ensemble Rg is
{E['E4']['median_rg_ratio_disordered']:.2f} times the measured value (95% interval
{E['E4']['ci95_median_rg_ratio'][0]:.2f}–{E['E4']['ci95_median_rg_ratio'][1]:.2f}; two-sided Wilcoxon p = {E['E4']['p_two_sided']:.1e}), and only {pct(E['E4']['share_below_1'])} of disordered
entries are below 1, so the model's disordered ensembles are too extended. {coord_sentence}

## 2. What SAXS can and cannot resolve

Reweighting each ensemble towards its profile sorts the entries into four kinds:
{res.get(analyse.KIND_RAW, 0)} raw fits; {res.get(analyse.KIND_MODEST, 0)} modestly reweightable (χ² ≤ 2 while keeping at least half the
effective sample); {res.get(analyse.KIND_STRONG, 0)} strongly reweightable (χ² ≤ 2 only below that, which requires substantial population
changes); {res.get(analyse.KIND_NOTFIT, 0)} whose target fit is not reached by the tested reweighting procedure. Expectation E3 is
{'met' if E['E3']['met'] else 'not met'}: {pct(E['E3']['share_phi_above_0.3'])} of entries reach χ² ≈ 1 while
keeping more than 30% of the effective sample.

As interpretation: a **strongly reweightable** entry is one where the right conformers are present and
mis-weighted, so a profile of it carries a usable fine-tuning signal; an entry whose **target fit is not
reached** needs a different measurement or a different model, because SAXS says the ensemble is wrong without saying how
(NMR relaxation or chemical shifts, HDX-MS, single-molecule FRET on labelled constructs); a **modestly
reweightable** entry needs a small, consistent correction, often of the Rg scale, for which repeat
measurements at several concentrations matter more than new systems. φ measures how concentrated the
weights become, not information gain; SAXS is low-dimensional, so distinct ensembles can give the same
curve; and a target fit that is not reached does not establish that compatible conformations are absent from
BioEmu's distribution: it describes this finite ensemble, this forward model and this procedure.

## 3. Arm 1, replication: high-error, reweightable monomers under one standard condition

The design (§6) proposed a priority score, log₁₀ raw χ² × a reweighting weight × a tractability weight, to
rank systems for measurement. A retrospective test (end of this section) shows that the score does not pick high-discrepancy
profiles better than random choice when it is fed the quantities known before a measurement, and its weights were
set by judgement, so it is not used to recommend measurements. Arm 1 takes the systems that meet stated
eligibility rules: reweighting reaches χ² ≤ 2 (strongly or modestly reweightable), the measured mass matches a
monomer, the chain has at most {arms['arm1']['max_length']} residues, and the protein is folded or partly disordered.
{arms['arm1']['n_eligible']} entries qualify. The {len(arm1)} with the largest observed raw χ² are listed, which is the
size of the known discrepancy a replication can test and not an estimate of what the measurement is worth.
They are the same {len(arm1)} the score ranks first; what is dropped is its judgement-set weighting. The final choice among eligible systems needs a feasibility review (construct availability,
expression, deposited buffer).

| SASBDB | length | class | raw χ² | best χ² | kind | Rg |
|---|---|---|---|---|---|---|
{rows1}

Re-measuring deposited systems under standardised conditions (one buffer, one temperature, tags removed, a
concentration series) reduces condition heterogeneity and tests whether the discrepancy is reproducible; one
common buffer and temperature can itself shift some ensembles, so a changed profile is read against the
deposited one before it is read against the model. It adds no new region of sequence space. {novel}

{sel_sentence} This is
retrospective, and a profile that already exists carries no new information for the model; the evaluation
says only how far the rule can be trusted to rank unmeasured proteins from what is known about them.

## 4. Arm 2, hypothesis test: are BioEmu-1's disordered ensembles too extended?

The pre-registered expectation E4 predicted over-compaction; the archive shows the opposite. The hypothesis
came from this archive, so this arm is a prospective replication test with pre-specified model-based
selection and new SAXS measurements. The split into test and control proteins uses only the model and the
sequence; the pool is the archive's disordered entries whose deposited profiles pass the quality checks. The test arm takes disordered
proteins whose BioEmu-1 ensemble Rg exceeds the disordered scaling-law Rg (1.927 N^0.598 Å) by more than
{int(100 * (arms['arm2']['threshold_ratio'] - 1))}%, spanning chain length and net charge
({arms['arm2']['n_test_candidates']} of {arms['arm2']['n_disordered_pool']} disordered entries qualify); the
scaling-law-matched controls are disordered proteins of the same length bins whose ratio to the scaling law is
within 10% of 1 ({arms['arm2']['n_control_candidates']} qualify). In the archive their observed Rg_model / Rg_exp
is {f2(arms['arm2']['control_observed_median_ratio'])} at the median, against
{f2(arms['arm2']['test_observed_median_ratio'])} for the test candidates, so the controls are not known to agree
with experiment; they are matched to the test arm on the scaling-law ratio only.

Test arm:

| SASBDB | length | NCPR | FCR | ensemble Rg (Å) | scaling-law Rg (Å) | ratio |
|---|---|---|---|---|---|---|
{rows2}

Controls:

| SASBDB | length | NCPR | FCR | ensemble Rg (Å) | scaling-law Rg (Å) | ratio |
|---|---|---|---|---|---|---|
{rows2c}

Pre-specified outcome: the median of Rg_exp / Rg_model across the test arm, from Guinier fits of the new
profiles. The bias is refuted if that median is at least 0.95; it is confirmed if the median is at most 0.90
in the test arm and at least 0.95 in the controls. The campaign repeats the archive comparison under one
condition, with the selection made from the model alone.

## 5. Assay design

- **Primary assay.** Size-exclusion-coupled SAXS (SEC-SAXS) at a synchrotron beamline, so that
  aggregates and oligomers are separated from the monomer before scattering is recorded. Three
  concentrations across the elution peak are integrated and compared; Guinier Rg agreement across them is
  the first QC gate.
- **Constructs.** Expressed with a cleavable tag; the tag is removed before measurement, because a 2–3 kDa
  tag changes Rg of a small protein measurably. Disordered termini longer than ~15 residues are either
  kept (if present in the deposited entry) or trimmed in a second construct, so the model's handling of
  termini can be tested directly.
- **Buffer.** Matched buffer from the final SEC step, with the exact dialysate recorded; BioEmu does not
  model pH or ionic strength, so buffer is held at one standard condition (20 mM HEPES pH 7.4, 150 mM NaCl,
  2 mM TCEP) and recorded with each dataset.
- **Standards and controls.** A protein standard (bovine serum albumin or glucose isomerase) in every
  session for absolute scale; water for intensity calibration; a buffer-only frame bracketing each
  sample; a repeat of one previously measured SASBDB entry per batch as a cross-site control.
- **Companion measurement.** For entries whose target fit is not reached, HDX-MS on the same batch of protein: it reports
  per-segment exchange that distinguishes a locally unfolded region from a globally wrong fold.

### Feasibility judgement

- A single common buffer and temperature is the aim, but it may not hold for every protein: some may be
  insoluble, unstable or aggregate in it.
- Gate before any SAXS: a solubility test in the standard buffer; the A280 trace of the purification and of
  analytical SEC (a single symmetric peak at the expected elution volume) with the A260/A280 ratio for
  nucleic-acid contamination; SDS-PAGE (a single band at the expected mass, for purity and identity); and
  circular dichroism where the fold is in question, to confirm the expected folded or disordered content.
- Aggregation signs on the SAXS data: a low-angle upturn, a non-linear Guinier region, Rg or I(0)/c rising
  with concentration, a molecular weight from I(0) above the monomer, and frames across the SEC peak that do
  not give a constant Rg; all are among the scripted criteria of §6.
- Fallback: a protein that fails in the standard buffer is measured in the nearest buffer in which it is
  monodisperse (or in its deposited buffer), the deviation recorded, and analysed separately, because its
  replication no longer controls the condition; one that fails in every buffer is dropped and replaced by the
  next eligible system.

## 6. QC criteria, applied by script

Every dataset passes or fails on recorded values:

1. Guinier region: at least 12 points with q·Rg ≤ 1.3 and a linear fit whose residuals show no trend.
2. No aggregation: lowest-angle intensity within 5% of the Guinier line; Rg from the three concentrations
   agree within 3%.
3. Radiation damage: frame-to-frame χ² across the exposure series below 1.5; otherwise frames are
   discarded from the end until it is.
4. Molecular weight from I(0) against the standard within 20% of the sequence mass of the monomer.
5. Buffer subtraction sanity: no negative intensities within the usable range; the high-q plateau is
   flat.
6. Reproducibility: the cross-site control reproduces its deposited profile with χ² < 2 after scaling.
7. No concentration dependence: I(0)/c agrees within 5% across the three concentrations (Rg is covered by
   criterion 2).
8. SEC-SAXS peak homogeneity: Rg from the frames across the elution peak is constant within 3%.

Several criteria watch for aggregation on the SAXS data themselves: a low-angle upturn and a non-linear
Guinier region (1, 2), Rg or I(0)/c rising with concentration (2, 7), a molecular weight from I(0) above the
monomer (4), and frames across the SEC peak that do not give a constant Rg (8).

Datasets that fail are repeated once; a second failure retires the construct and triggers the
contingency below.

## 7. Success criteria and contingencies

- **Dataset-level success:** QC passed, and the profile constrains the BioEmu ensemble (reweighting to
  χ² ≤ 1 changes the effective sample fraction by at least 0.2). A profile that the raw ensemble already
  fits is still deposited, but counted as a confirmation, not as training signal.
- **Arm 1 success:** at least 60% of commissioned systems yield a usable profile; the set covers the three
  disorder classes in roughly the proportions of the error map; the strongly reweightable kind is
  over-represented, because that is where the model can learn.
- **Arm 2 success:** the pre-specified outcome is reached with at least eight test and six control profiles
  passing QC, whichever way it falls.
- **Contingencies.** Poor expression → switch to a homologue from the same family with a deposited
  SASBDB entry. A construct that fails the solubility and sample-quality gate, or aggregates, in the
  standard buffer → measured in the nearest buffer in which it is monodisperse (or in its deposited buffer),
  with the deviation recorded, and analysed separately, because its replication no longer controls the
  condition; one that fails in every buffer is dropped and replaced by the next eligible system.
  Beamtime loss → a laboratory SAXS instrument for the smallest, most concentrated samples, accepting
  the lower q range and recording it. Ambiguous SAXS (target fit not reached) → HDX-MS first, SAXS second.

## 8. Work packages for an external provider

| Package | Deliverable | Review point |
|---|---|---|
| WP1 Constructs | Expression plasmids for both arms, sequence-verified, with and without tags | Design review before synthesis; sequence files checked by script |
| WP2 Protein production | ≥ 2 mg per construct at ≥ 95% purity by SEC and SDS-PAGE, monodisperse by DLS | Purity and DLS reports per batch; a batch failing DLS does not proceed |
| WP3 SEC-SAXS | Three-concentration SEC-SAXS per sample, with standards and buffer frames; raw frames and reduced curves | Scripted QC on delivery; results joined to the candidate table |
| WP4 HDX-MS (target fit not reached only) | Deuterium uptake per peptide at four time points, with a fully deuterated control | Peptide coverage ≥ 85%; back-exchange reported |
| WP5 Deposition | SASBDB deposition of every QC-passed dataset with full metadata | Accession codes recorded against each candidate |

Milestones are set per batch of twelve constructs: constructs at week 2, protein at week 6, SAXS at
week 8, QC and model re-scoring at week 9. The re-scoring repeats this pipeline on the new profiles, so
every batch updates the error map before the next batch is chosen.

## 9. What this proposal rests on, and what it does not claim

The ensembles are from the PeptoneBench archive, generated with BioEmu at code commit ac7455d; the archive
does not record the checkpoint version, and re-sampling with the current v1.2 checkpoint is the first
follow-up. SASBDB profiles differ in buffer, temperature and construct, none of which the model sees, so
part of the raw error is condition mismatch, which a campaign under one standard condition reduces. The
scoring treats the Pepsi-SAXS forward model as exact; its hydration-shell parameters are a known source of
Rg bias of a few per cent, which is why modestly reweightable cases are read with care. The comparison with
other models is restricted to the six comparators on the shared entries.
"""
    (C.REPO / "docs" / "CAMPAIGN.md").write_text(text)


def main() -> None:
    style()
    OUT.mkdir(parents=True, exist_ok=True)
    A = json.load(open(C.RESULTS / "analysis.json"))
    d = pd.read_csv(C.RESULTS / "resolvability.csv")
    cand = pd.read_csv(C.RESULTS / "candidates.csv")
    paths = json.load(open(C.RESULTS / "paths" / f"{C.MODEL}.json"))
    ent = pd.read_csv(C.RESULTS / "entries.csv")
    fig_error_map(d, A)
    fig_length(d, A)
    fig_rg(d)
    fig_rg_robustness(d, A)
    fig_kinds(d, paths, A)
    fig_resolvability(d, paths)
    if A["model_comparison"]:
        fig_models(A)
    fig_predict()
    fig_predict_detail()
    fig_examples(A, ent)
    arms = campaign_arms(A, cand, d, ent)
    write_campaign(A, cand, arms)
    print(f"  figures -> {OUT}; proposal -> docs/CAMPAIGN.md", flush=True)


if __name__ == "__main__":
    main()
