#!/usr/bin/env python3
"""
Layer-wise probes on raw embeddings.

Usage:
  python scripts/make_probe_figure.py
"""
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({
    "font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
    "axes.titlesize": 10, "axes.titlepad": 11, "axes.titlepad": 11, "figure.dpi": 200, "savefig.bbox": "standard",
})

ROOT = Path(__file__).resolve().parent.parent
PROBE = ROOT / "outputs" / "probe"

MODELS = [("esm2_8m", "ESM-2-8M", "#7FA8D4", "o", "-"),
          ("esm2_35m", "ESM-2-35M", "#3B6FB6", "s", "-"),
          ("seqdance", "SeqDance", "#C0504D", "D", "-")]
AXES = [("disprot-molten_globule", "Molten globule", "n=83"),
        ("disprot-pre_molten_globule", "Pre-molten globule", "n=52")]


def load_profiles(axis):
    out, nulls = {}, []
    for key, label, *_ in MODELS:
        f = PROBE / f"{key}_embedding.json"
        if not f.exists():
            continue
        r = json.load(open(f))
        pts = []
        for tag, d in r["layers"].items():
            if axis in d:
                pts.append((int(tag[5:]), d[axis]["f1"],
                            d[axis]["ci95"][0], d[axis]["ci95"][1]))
                nn = d.get(axis + "__null")
                if nn:
                    nulls.append(nn["f1"])
        out[key] = sorted(pts)
    comp = None
    fc = PROBE / "esm2_8m_composition.json"
    if fc.exists():
        c = json.load(open(fc))["layers"].get("composition", {}).get(axis)
        comp = c["f1"] if c else None
    return out, (float(np.mean(nulls)) if nulls else None), comp


def panel_profile(ax, axis, title, ntxt, show_legend):
    prof, null, comp = load_profiles(axis)
    if null is not None:
        ax.axhspan(0, null, color="#BBBBBB", alpha=0.35, zorder=0)
        ax.text(1.04, null + 0.008, "shuffled-region\nnull", transform=ax.get_yaxis_transform(),
                ha="left", va="top", fontsize=8, color="#555", linespacing=1.2,
                clip_on=False, zorder=4)
    if comp is not None:
        ax.axhline(comp, color="#6E9C6E", ls=":", lw=1.2, zorder=1)
        ax.text(1.04, comp - 0.008, "sequence\ncomposition", transform=ax.get_yaxis_transform(),
                ha="left", va="bottom", fontsize=8, color="#4E7C4E", linespacing=1.2,
                clip_on=False, zorder=4)

    for key, label, col, mk, ls in MODELS:
        pts = prof.get(key)
        if not pts:
            continue
        x = [p[0] for p in pts]; y = [p[1] for p in pts]
        lo = [p[1] - p[2] for p in pts]; hi = [p[3] - p[1] for p in pts]
        ax.errorbar(x, y, yerr=[lo, hi], color=col, marker=mk, ls=ls, ms=4, lw=1.4,
                    capsize=2, elinewidth=0.7, label=label, zorder=3)

    ax.set(xlabel="layer", ylabel="region-aware F1",
           title=f"(a) Molten-globule detection by layer ({ntxt} regions)")
    ax.set_ylim(bottom=0)
    if show_legend:
        ax.legend(fontsize=8, loc="upper left", frameon=False)


def panel_selectivity(ax):
    """Layer-wise gap at MATCHED depth (ESM-2-35M vs SeqDance: same architecture, size and
    layer index). The sign of the gap tracks axis TYPE, not depth — conformational axes sit
    above zero at every layer, compositional axes below. This is what rules out 'SeqDance
    simply has better disorder representations', and it does so without selecting a layer."""
    E = PROBE / "esm2_35m_embedding.json"
    S = PROBE / "seqdance_embedding.json"
    PB = PROBE / "paired_by_layer.json"
    if not (E.exists() and S.exists()):
        ax.text(.5, .5, "no data", ha="center"); return
    e = json.load(open(E))["layers"]; sd = json.load(open(S))["layers"]
    paired = json.load(open(PB)) if PB.exists() else {}

    SERIES = [("disprot-molten_globule", "Molten globule", "#C0504D", "D", "-"),
              ("disprot-pre_molten_globule", "Pre-molten globule", "#E08E8C", "o", "-"),
              ("prediction-disorder-alphafold", "AlphaFold disorder", "#3B6FB6", "s", "--"),
              ("prediction-proline_rich-mobidb_lite_sub", "Proline-rich", "#7FA8D4", "^", "--")]
    for axis, lab, col, mk, ls in SERIES:
        pts = []
        for tag in e:
            if axis in e[tag] and tag in sd and axis in sd[tag]:
                pr = paired.get(axis, {}).get(tag)
                # prefer the paired gap and its CI; fall back to the difference of point
                # estimates when a layer has no dumps to pair over
                g = pr["gap"] if pr else sd[tag][axis]["f1"] - e[tag][axis]["f1"]
                ci = pr["ci95"] if pr else None
                pts.append((int(tag[5:]), g, ci, bool(pr and pr["sig"])))
        if not pts:
            continue
        pts.sort()
        x = [p[0] for p in pts]; y = [p[1] for p in pts]
        if any(p[2] for p in pts):
            lo = [p[1] - p[2][0] if p[2] else 0 for p in pts]
            hi = [p[2][1] - p[1] if p[2] else 0 for p in pts]
            ax.errorbar(x, y, yerr=[lo, hi], color=col, marker=mk, ls=ls, ms=4, lw=1.4,
                        capsize=2, elinewidth=0.7, label=lab, zorder=3)
        else:
            ax.plot(x, y, color=col, marker=mk, ls=ls, ms=4, lw=1.4, label=lab, zorder=3)
        # hollow markers where the paired test does not exclude zero
        ns = [(p[0], p[1]) for p in pts if p[2] and not p[3]]
        if ns:
            ax.plot([p[0] for p in ns], [p[1] for p in ns], marker=mk, ls="none", ms=4,
                    mfc="white", mec=col, mew=1.2, zorder=4)

    ax.axhline(0, color="#333", lw=1.0, zorder=2)
    lim = ax.get_ylim()
    ax.axhspan(0, lim[1], color="#C0504D", alpha=0.06, zorder=0)
    ax.axhspan(lim[0], 0, color="#3B6FB6", alpha=0.06, zorder=0)
    ax.set_ylim(lim)
    ax.text(0.02, 0.96, "favors SeqDance", transform=ax.transAxes, ha="left", va="top",
            fontsize=8, color="#8B3A38")
    ax.text(0.30, 0.04, "favors ESM-2", transform=ax.transAxes, ha="center", va="bottom",
            fontsize=8, color="#2E5A93")
    ax.set(xlabel="layer", ylabel="ΔF1  (SeqDance − ESM-2-35M)",
           title="(b) Conformational versus compositional axes")
    ax.legend(fontsize=8, loc="lower right", frameon=False)


def panel_paired(ax):
    files = [("paired_matched_L10.json", "ESM-2-35M L10 → SeqDance L10", "matched scale + depth"),
             ("paired_paper_layers.json", "ESM-2-8M L4 → SeqDance L10", "paper's SAE layers"),
             ("paired_best_layers.json", "ESM-2-8M L1 → SeqDance L12", "best layer each (post hoc)")]
    rows = []
    for fn, lab, sub in files:
        p = PROBE / fn
        if p.exists():
            r = json.load(open(p))
            rows.append((lab, sub, r["gap"], r["gap_ci95"], r["boot_p_cmp_le_ref"]))
    if not rows:
        ax.text(.5, .5, "no paired results", ha="center"); return

    y = np.arange(len(rows))[::-1]
    for yi, (lab, sub, g, ci, p) in zip(y, rows):
        ax.plot([ci[0], ci[1]], [yi, yi], color="#444", lw=1.6, solid_capstyle="round", zorder=2)
        ax.plot([g], [yi], marker="D", ms=6, color="#C0504D", zorder=3)
        ax.text(ci[1] + 0.012, yi, f"p={p:.4f}", va="center", fontsize=7, color="#333")
    ax.axvline(0, color="#999", lw=0.9, ls="--", zorder=1)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{l}\n{s}" for l, s, *_ in rows], fontsize=7)
    ax.set(xlabel="paired ΔF1 (SeqDance − ESM-2), 95% CI",
           title="Molten globule: paired region bootstrap")
    ax.set_xlim(left=min(0, min(r[3][0] for r in rows)) - 0.03,
                right=max(r[3][1] for r in rows) + 0.09)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)


def main():
    fig, ax = plt.subplots(1, 2, figsize=(13.0, 3.6))
    panel_profile(ax[0], *AXES[0], show_legend=True)
    panel_selectivity(ax[1])
    # panels keep their 11in geometry; the extra 2in of canvas is blank margin,
    # so the panels are not reshaped - only the overall aspect changes
    pad = (13.0 - 11.0) / 2 / 13.0
    fig.subplots_adjust(left=pad + 0.125 * 11 / 13, right=1 - (pad + 0.10 * 11 / 13),
                        top=0.86, bottom=0.17, wspace=0.55)
    out = ROOT / "outputs" / "figures" / "fig_probe.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=220, bbox_inches=None)
    fig.savefig(out.with_suffix(".pdf"), bbox_inches=None)
    print(f"saved → {out}")


if __name__ == "__main__":
    main()
