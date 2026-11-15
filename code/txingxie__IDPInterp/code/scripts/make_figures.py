#!/usr/bin/env python3
"""
Render figures and tables.

Produces (in outputs/figures/):
  fig2_atlas.png       ESM-2 compositional atlas: UMAP families, per-axis F1, held-out r
  fig3_dissociation.png  ESM vs SeqDance per-axis F1, molten-globule bootstrap, cross-model uniqueness
  tables.md            Table 1 (SAE quality), Table 2 (per-axis F1), Table 3 (exemplar features)

Usage:  python scripts/make_figures.py
"""
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "outputs"
FIG = OUT / "figures"
FIG.mkdir(exist_ok=True)
ESM = OUT / "atlas" / "esm2_8m_topk_32x_k64_layer4"
SD = OUT / "atlas" / "seqdance_topk_32x_k64_layer10"
CMP = OUT / "compare" / "esm2_8m_topk_32x_k64_layer4__vs__seqdance_topk_32x_k64_layer10"

plt.rcParams.update({
    "font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
    "axes.titlesize": 10, "figure.dpi": 200, "savefig.bbox": "tight",
})
ESM_C, SD_C = "#3B6FB6", "#C0504D"  # blue / red

AXIS_LABEL = {
    "prediction-disorder-alphafold": "AlphaFold disorder",
    "curated-disorder-merge": "Curated disorder",
    "curated-phase_separation-merge": "Phase separation",
    "prediction-proline_rich-mobidb_lite_sub": "Proline-rich",
    "prediction-extended-mobidb_lite_sub": "Extended coil",
    "derived-binding_mode_disorder_to_disorder-mobi": "Fuzzy (D-to-D)",
    "prediction-compact-mobidb_lite_sub": "Compact (ν proxy)",
    "derived-binding_mode_disorder_to_order-mobi": "MoRF (MobiDB)",
    "disprot-disorder_to_order": "MoRF (experimental)",
    "disprot-molten_globule": "Molten globule",
    "disprot-pre_molten_globule": "Pre-molten globule",
    "derived-fixed-th_90": "Ordered (control)",
}
# generic/compositional first, then conformational sub-states
AXIS_ORDER = [
    "prediction-disorder-alphafold", "curated-disorder-merge", "curated-phase_separation-merge",
    "prediction-proline_rich-mobidb_lite_sub", "prediction-extended-mobidb_lite_sub",
    "derived-binding_mode_disorder_to_disorder-mobi", "derived-fixed-th_90",
    "prediction-compact-mobidb_lite_sub", "derived-binding_mode_disorder_to_order-mobi",
    "disprot-disorder_to_order", "disprot-molten_globule", "disprot-pre_molten_globule",
]


def load():
    esm_sum = json.load(open(ESM / "atlas_summary.json"))
    sd_sum = json.load(open(SD / "atlas_summary.json"))
    esm_sp = json.load(open(OUT / "diagnostics/esm2_8m_topk_32x_k64_layer4/sparsity_report.json"))
    sd_sp = json.load(open(OUT / "diagnostics/seqdance_topk_32x_k64_layer10/sparsity_report.json"))
    # Panel (b) reports the primary (ESM-2-8M vs SeqDance) molten-globule test under
    # five-fold, protein-held-out cross-validation over all 83 regions (Section 3.5),
    # not the underpowered single-split bootstrap.
    boot = json.load(open(OUT / "bootstrap/disprot-molten_globule_cv5_esm2_8m.json"))
    cmp_sum = json.load(open(CMP / "summary.json"))
    esm_at = pd.read_csv(ESM / "atlas_table.csv")
    clusters = pd.read_csv(OUT / "diagnostics/esm2_8m_topk_32x_k64_layer4/feature_clusters.csv")
    cmp_uni = pd.read_csv(CMP / "cmp_disorder_unique.csv")
    return dict(esm_sum=esm_sum, sd_sum=sd_sum, esm_sp=esm_sp, sd_sp=sd_sp,
               boot=boot, cmp_sum=cmp_sum, esm_at=esm_at, clusters=clusters, cmp_uni=cmp_uni)


# ---------------------------------------------------------------- Figure 2
def fig2(d):
    fig, ax = plt.subplots(1, 3, figsize=(11, 3.3))

    # (a) UMAP of decoder dictionary, colored by family
    cl = d["clusters"]
    noise = cl["cluster"] == -1
    ax[0].scatter(cl.loc[noise, "umap_x"], cl.loc[noise, "umap_y"], s=3, c="#DDDDDD", linewidths=0)
    fams = cl.loc[~noise, "cluster"]
    ax[0].scatter(cl.loc[~noise, "umap_x"], cl.loc[~noise, "umap_y"], s=5,
                  c=fams, cmap="tab20", linewidths=0)
    ax[0].set(xticks=[], yticks=[], xlabel="UMAP 1", ylabel="UMAP 2",
              title=f"(a) ESM-2 decoder families (n={d['esm_sum']['n_families']})")

    # label a few validated families at their exemplar feature's UMAP position
    # (exemplars from Table 2/3); text pushed radially outward with a leader line
    LABELS = {8306: "acidic (D/E)", 4236: "arginine / basic", 7063: "proline",
              3350: "glycine linker", 3401: "Ser/Thr", 8708: "histidine"}
    # manual offset overrides (pts) for the crowded top-left pair so they don't collide
    OFFSETS = {7063: (46, 4), 3350: (-26, 22)}
    if "feature_idx" in cl.columns:
        pos = cl.drop_duplicates("feature_idx").set_index("feature_idx")
        cx, cy = cl["umap_x"].mean(), cl["umap_y"].mean()
        for fid, lab in LABELS.items():
            if fid not in pos.index:
                continue
            px, py = float(pos.loc[fid, "umap_x"]), float(pos.loc[fid, "umap_y"])
            dx, dy = px - cx, py - cy
            nrm = (dx * dx + dy * dy) ** 0.5 + 1e-9
            off = OFFSETS.get(fid, (34 * dx / nrm, 34 * dy / nrm))
            ax[0].scatter([px], [py], s=24, facecolors="none",
                          edgecolors="#222", linewidths=0.6, zorder=5)
            ax[0].annotate(lab, (px, py), xytext=off,
                           textcoords="offset points", fontsize=6, weight="bold",
                           ha="center", va="center", color="#222", zorder=6,
                           arrowprops=dict(arrowstyle="-", lw=0.5, color="#555"),
                           bbox=dict(boxstyle="round,pad=0.15", fc="white",
                                     ec="#bbbbbb", lw=0.3, alpha=0.9))

    # (b) per-axis best F1 (ESM)
    axes = [a for a in AXIS_ORDER if a in d["esm_sum"]["per_axis_best_f1"]]
    vals = [d["esm_sum"]["per_axis_best_f1"][a] for a in axes]
    y = np.arange(len(axes))[::-1]
    ax[1].barh(y, vals, color=ESM_C)
    ax[1].set_yticks(y)
    ax[1].set_yticklabels([AXIS_LABEL[a] for a in axes], fontsize=7)
    ax[1].set(xlabel="best-feature test F1", title="(b) ESM-2 F1 by axis", xlim=(0, 0.9))

    # (c) held-out description r distribution
    r = d["esm_at"]["pearson_r"].dropna().values
    med = float(np.median(r))
    ax[2].hist(r, bins=24, color=ESM_C, alpha=0.85)
    ax[2].axvline(med, color="black", ls="--", lw=1)
    ax[2].text(med + 0.02, ax[2].get_ylim()[1] * 0.9, f"median {med:.2f}", fontsize=8)
    ax[2].set(xlabel="held-out Pearson r", ylabel="features",
              title=f"(c) Description validation (n={len(r)})")
    fig.tight_layout()
    fig.savefig(FIG / "fig2_atlas.png"); fig.savefig(FIG / "fig2_atlas.pdf")
    plt.close(fig)


# ---------------------------------------------------------------- Figure 3
def fig3(d):
    fig, ax = plt.subplots(1, 3, figsize=(12, 3.5))

    # (a) per-axis F1, ESM vs SeqDance
    axes = [a for a in AXIS_ORDER if a in d["esm_sum"]["per_axis_best_f1"]]
    e = [d["esm_sum"]["per_axis_best_f1"][a] for a in axes]
    s = [d["sd_sum"]["per_axis_best_f1"].get(a, 0) for a in axes]
    y = np.arange(len(axes))[::-1]
    h = 0.38
    ax[0].barh(y + h / 2, e, h, color=ESM_C, label="ESM-2")
    ax[0].barh(y - h / 2, s, h, color=SD_C, label="SeqDance")
    ax[0].set_yticks(y)
    ax[0].set_yticklabels([AXIS_LABEL[a] for a in axes], fontsize=7)
    ax[0].axhline(4.5, color="#999", lw=0.6, ls=":")  # between generic and conformational
    ax[0].set(xlabel="best-feature test F1", title="(a) Disorder axes by model", xlim=(0, 0.9))
    ax[0].legend(fontsize=8, loc="lower right")

    # (b) molten-globule bootstrap
    b = d["boot"]
    labels = ["ESM-2\n(feat 2594)", "SeqDance\n(feat 188)"]
    pts = [b["ref"]["point_f1"], b["cmp"]["point_f1"]]
    ci = np.array([b["ref"]["ci95"], b["cmp"]["ci95"]])
    err = np.abs(ci.T - np.array(pts))
    ax[1].bar([0, 1], pts, color=[ESM_C, SD_C], width=0.55)
    ax[1].errorbar([0, 1], pts, yerr=err, fmt="none", ecolor="black", capsize=4, lw=1)
    ax[1].set_xticks([0, 1]); ax[1].set_xticklabels(labels, fontsize=8)
    ax[1].set(ylabel="molten-globule F1", title="(b) Molten-globule detection", ylim=(0, 0.5))
    # Annotate the paired CI: marginal per-model intervals cannot support a claim
    # about the difference between models.
    lo, hi = b["gap_ci95"]
    ax[1].text(0.5, 0.455,
               f"paired gap +{b['gap_cmp_minus_ref']:.2f}  "
               f"95% CI [{lo:+.2f}, {hi:+.2f}]\n"
               f"{b['k_folds']}-fold CV, n={b['n_regions_all_splits']} regions",
               ha="center", fontsize=7.5)

    # (c) cross-model uniqueness: max corr to any ESM-2 feature, grouped by how
    # compaction-specific the conformational state is (generic disorder -> pre-molten
    # -> molten globule). Generic disorder is a real distribution (n~999, violin);
    # the compaction states are small discrete sets shown as strips, not violins.
    # Molten-globule detectors are pinned by validated identity (Sec 4.4: 188, 14136,
    # 3966, 5229) rather than DisProt argmax, which files 188/14136 under pre-molten
    # by a 0.32-vs-0.30 hair. Pre-molten is the remaining pre-molten-globule set, and
    # is deliberately heterogeneous: its charge-based features overlap ESM-2 (up to
    # ~0.42) while its hydrophobic ones do not -- the "charge shared, hydrophobic
    # unique" thesis made visible.
    u = d["cmp_uni"]
    MG = {188, 14136, 3966, 5229}
    gen = u[u["triplet"] == "prediction-disorder-alphafold"]["best_abs_corr"].values
    mg = u[u["cmp_feature"].isin(MG)]["best_abs_corr"].values
    pm = u[(u["triplet"] == "disprot-pre_molten_globule")
           & (~u["cmp_feature"].isin(MG))]["best_abs_corr"].values
    np.random.seed(1)
    # Violins only where n supports a density estimate (generic n~999, pre-molten n~83);
    # molten globule (n=4) is too small for a KDE and is shown as points only.
    for vals, pos, vcol in [(gen, 1, ESM_C), (pm, 2, "#E4956B")]:
        pc = ax[2].violinplot([vals], positions=[pos], showmeans=False, showextrema=False, widths=0.8)
        pc["bodies"][0].set_facecolor(vcol); pc["bodies"][0].set_alpha(0.35)
    # Violins carry group identity, so dots are one neutral colour at one alpha and
    # density reads from overlap. Molten globule has no violin (n=4), so its dots take
    # the group colour instead.
    dot = dict(s=14, alpha=0.45, linewidths=0)
    ax[2].scatter(np.random.normal(1, 0.05, len(gen)), gen, c="#555", **dot)
    ax[2].scatter(np.random.normal(2, 0.06, len(pm)), pm, c="#555", **dot)
    ax[2].scatter(np.random.normal(3, 0.03, len(mg)), mg, c=SD_C, s=14, alpha=0.45, linewidths=0, zorder=3)
    r188 = u[u["cmp_feature"] == 188]["best_abs_corr"]
    if len(r188):
        ax[2].annotate("feat 188", (3, float(r188.iloc[0])), xytext=(3.2, 0.17), fontsize=7,
                       arrowprops=dict(arrowstyle="->", lw=0.6, color="#333"))
    ax[2].axhline(0.10, color="#999", ls="--", lw=0.8)
    ax[2].set_xticks([1, 2, 3]); ax[2].set_xlim(0.5, 3.75)
    ax[2].set_xticklabels([f"Generic\ndisorder (n={len(gen)})",
                           f"Pre-molten\nglobule (n={len(pm)})",
                           f"Molten\nglobule (n={len(mg)})"], fontsize=7.5)
    ax[2].set(ylabel="max corr. to any ESM-2 feature",
              title=f"(c) Cross-model uniqueness (CKA={d['cmp_sum']['linear_cka']:.2f})",
              ylim=(0, 0.65))
    fig.tight_layout()
    fig.savefig(FIG / "fig3_dissociation.png"); fig.savefig(FIG / "fig3_dissociation.pdf")
    plt.close(fig)


# ---------------------------------------------------------------- Tables
def tables(d):
    L = []
    # Table 1: SAE quality
    es, ss, esp, ssp = d["esm_sum"], d["sd_sum"], d["esm_sp"], d["sd_sp"]
    L += ["### Table 1. Sparse-autoencoder quality.\n",
          "| Model | Layer | Dict. | L0 | Dead % | #F1>0.5 | #F1>0.3 | Families | Median r (n) |",
          "|---|---|---|---|---|---|---|---|---|",
          f"| ESM-2-8M | 4 | 10,240 | {esp['L0_mean']:.0f} | {esp['dead_pct']:.2f} | "
          f"{es['n_F1_gt_0.5']} | {es['n_F1_gt_0.3']} | {es['n_families']} | "
          f"{es['median_heldout_pearson_r']:.2f} ({es['n_features_validated']}) |",
          f"| SeqDance | 10 | 15,360 | {ssp['L0_mean']:.0f} | {ssp['dead_pct']:.2f} | "
          f"{ss['n_F1_gt_0.5']} | {ss['n_F1_gt_0.3']} | {ss['n_families']} | "
          f"{ss['median_heldout_pearson_r']:.2f} ({ss['n_features_validated']}) |", ""]

    # Table 2: per-axis F1
    L += ["### Table 2. Best-feature test F1 by disorder axis.\n",
          "| Axis | ESM-2 | SeqDance |", "|---|---|---|"]
    for a in AXIS_ORDER:
        if a in es["per_axis_best_f1"]:
            L.append(f"| {AXIS_LABEL[a]} | {es['per_axis_best_f1'][a]:.2f} | "
                     f"{ss['per_axis_best_f1'].get(a, float('nan')):.2f} |")
    L.append("")

    # Table 3: curated exemplar features spanning compositional families
    at = d["esm_at"].set_index("feature_idx")
    CURATED = [
        (8306, "Acidic low-complexity (poly-D/E)"),
        (2474, "Polar / charged homopolymer repeats"),
        (4236, "Arginine-rich / NLS-like"),
        (7063, "Proline-rich (polyproline-II)"),
        (3350, "Glycine-rich linker"),
        (3401, "Ser/Thr-rich (phospho-acceptor)"),
        (8708, "Histidine-rich"),
    ]
    L += ["### Table 3. Representative validated ESM-2 features spanning compositional families.\n",
          "| Feature | Identity (from validated description) | Axis | F1 | Held-out r |",
          "|---|---|---|---|---|"]
    for fid, ident in CURATED:
        if fid not in at.index:
            continue
        r = at.loc[fid]
        rr = f"{r['pearson_r']:.2f}" if pd.notna(r["pearson_r"]) else "—"
        L.append(f"| {fid} | {ident} | {AXIS_LABEL.get(r['triplet'], r['triplet'])} | "
                 f"{r['f1']:.2f} | {rr} |")
    L.append("")
    (FIG / "tables.md").write_text("\n".join(L))


# ---------------------------------------------------------------- Figure 1
def fig1():
    from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle

    fig, ax = plt.subplots(figsize=(12, 3.6))
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 3.6)
    ax.axis("off")

    def box(x, y, w, h, text, fc="#EEEEEE", ec="#444444", tc="black", fs=9, lw=1.0, tx=None):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.01,rounding_size=0.06",
                                    fc=fc, ec=ec, lw=lw))
        cx = x + w / 2 if tx is None else tx
        ax.text(cx, y + h / 2, text, ha="center", va="center", fontsize=fs, color=tc)

    def arrow(x0, y0, x1, y1, color="#444444"):
        ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=11,
                                     lw=1.0, color=color, shrinkA=1, shrinkB=1))

    # --- models -> per-residue embeddings
    box(0.15, 2.16, 1.8, 0.72, "ESM-2\n(evolution)", fc="#DCE6F2", ec=ESM_C, tc=ESM_C)
    box(0.15, 0.54, 1.8, 0.72, "SeqDance\n(dynamics)", fc="#F3DEDD", ec=SD_C, tc=SD_C)
    box(2.45, 1.35, 1.55, 0.72, "per-residue\nembeddings")
    arrow(1.95, 2.42, 2.45, 1.9)
    arrow(1.95, 0.9, 2.45, 1.52)
    arrow(4.0, 1.71, 4.55, 1.71)

    # --- SAE container: dense embedding -> encoder fan-out -> sparse feature dictionary
    box(4.55, 0.7, 3.2, 1.95, "", fc="#F2F2F2", ec="#888888")
    ax.text(6.15, 2.42, "sparse autoencoder (TopK + AuxK)", ha="center", fontsize=8.5)
    dxx, dch = 5.05, 0.14
    ddy0 = 1.42 - 2.5 * dch
    for i in range(5):  # dense embedding vector
        ax.add_patch(Rectangle((dxx, ddy0 + i * dch), 0.2, dch, fc="#B9C4D2", ec="#8A93A0", lw=0.5))
    ax.text(dxx + 0.1, ddy0 - 0.12, "$d$", ha="center", va="top", fontsize=7, color="#555")
    sxx, sch = 6.55, 0.1
    ssy0 = 1.42 - 6.5 * sch
    active = {3, 7, 10}
    for i in range(13):  # overcomplete sparse dictionary
        fc = "#555555" if i in active else "#FFFFFF"
        ax.add_patch(Rectangle((sxx, ssy0 + i * sch), 0.22, sch, fc=fc, ec="#AAAAAA", lw=0.5))
    ax.text(sxx + 0.11, ssy0 + 13 * sch + 0.06, "features", ha="center", fontsize=6.5, color="#555")
    for i in range(5):  # encoder fan-out lines (few dense dims -> many sparse features)
        y_d = ddy0 + (i + 0.5) * dch
        for y_s in (ssy0 + (i * 2.6 + 0.5) * sch, ssy0 + (i * 2.6 + 2.0) * sch):
            ax.plot([dxx + 0.2, sxx], [y_d, y_s], color="#A9B2BC", lw=0.6, zorder=1.5)

    # --- three analysis outputs, each with a small glyph on the left
    bx, bw = 8.2, 3.6
    cys = [2.7, 1.71, 0.72]
    labels = ["region-aware F1\n(12 disorder axes)",
              "automated description\n+ held-out validation",
              "cross-model comparison\n(activation corr., CKA)"]
    for cy, lab in zip(cys, labels):
        box(bx, cy - 0.38, bw, 0.76, lab, fs=8.5, tx=10.55)
        arrow(7.8, 1.71, bx, cy)
    gx, gw = 8.36, 0.92

    # glyph 1: grouped per-axis F1 bars (ESM-2 vs SeqDance)
    cy = cys[0]
    base = cy - 0.24
    bh = [0.72, 0.34, 0.55, 0.84, 0.4, 0.62]
    bcol = [ESM_C, SD_C, ESM_C, ESM_C, SD_C, ESM_C]
    slot = gw / len(bh)
    for i, (hgt, col) in enumerate(zip(bh, bcol)):
        ax.add_patch(Rectangle((gx + i * slot + slot * 0.15, base), slot * 0.6, hgt * 0.46, fc=col, ec="none"))
    ax.plot([gx, gx + gw], [base, base], color="#999999", lw=0.6)

    # glyph 2: held-out validation scatter along the diagonal
    cy = cys[1]
    x0, x1, y0, y1 = gx + 0.06, gx + gw - 0.06, cy - 0.26, cy + 0.26
    ax.plot([x0, x1], [y0, y1], color="#CCCCCC", lw=0.8, zorder=1.5)
    sx = np.linspace(x0 + 0.03, x1 - 0.03, 9)
    frac = (sx - x0) / (x1 - x0)
    jit = np.array([0.03, -0.04, 0.05, -0.02, 0.04, -0.05, 0.02, -0.03, 0.03])
    ax.scatter(sx, y0 + frac * (y1 - y0) + jit, s=7, c=ESM_C, edgecolors="none", zorder=2)

    # glyph 3: two model feature-vectors linked by activation correlation
    cy = cys[2]
    vh = 0.09
    for j, (vx, tint, patt) in enumerate([(gx + 0.06, "#DCE6F2", {1, 3}), (gx + gw - 0.28, "#F3DEDD", {0, 3})]):
        vy0 = cy - 2.5 * vh
        base_c, accent = ("#DCE6F2", ESM_C) if j == 0 else ("#F3DEDD", SD_C)
        for i in range(5):
            ax.add_patch(Rectangle((vx, vy0 + i * vh), 0.22, vh, fc=accent if i in patt else base_c,
                                   ec="#AAAAAA", lw=0.4))
    ax.annotate("", xy=(gx + gw - 0.28, cy), xytext=(gx + 0.28, cy),
                arrowprops=dict(arrowstyle="<->", color="#777777", lw=0.8))
    ax.text(gx + gw / 2 - 0.02, cy + 0.14, "$r$", ha="center", fontsize=7.5, color="#555")

    fig.savefig(FIG / "fig1_framework.png"); fig.savefig(FIG / "fig1_framework.pdf")
    plt.close(fig)


def main():
    d = load()
    fig1()
    fig2(d)
    fig3(d)
    tables(d)
    print("Wrote:", *(p.name for p in sorted(FIG.glob("*"))))


if __name__ == "__main__":
    main()
