from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

RES = Path("results/depmap_essentiality/results.json")
OUT = Path("paper/fig_main_1page.pdf")
BLUE, ORANGE, GREEN, GREY = "#0173B2", "#DE8F05", "#029E73", "#9AA7A5"
plt.rcParams.update({"font.size": 8, "font.family": "DejaVu Sans", "axes.spines.top": False,
                     "axes.spines.right": False, "axes.linewidth": 0.9})

DISPLAY = [
    ("PTM density", "ptm_density"),
    ("PTM count", "ptm_count"),
    ("Phospho density", "ptm_phospho_density"),
    ("Acetyl density", "ptm_acetyl_density"),
    ("Ubiquitination", "ptm_ubiq_density"),
    ("Charged fraction", "frac_charged"),
]


def main() -> None:
    r = json.load(open(RES))
    m = r["models"]
    ca = r["confound_analysis"]
    so = r["sequence_only"]
    gb = m["gradient_boosting"]
    auc = {u["feature"]: u["auc"] for u in r["feature_separation_top"]}

    missing = [k for _, k in DISPLAY if k not in auc]
    if missing:
        print("WARNING missing univariate AUC for:", missing)

    fig, ax = plt.subplots(1, 2, figsize=(5.0, 1.28), gridspec_kw={"width_ratios": [1.1, 1]})

    labels = ["Base", "Seq", "LR", "RF", "GB"]
    vals = [ca["confound_only_auc"], so["cv_auc"], m["logistic"]["roc_auc"],
            m["random_forest"]["roc_auc"], gb["roc_auc"]]
    colors = [GREY, ORANGE, BLUE, BLUE, GREEN]
    err = [[0] * 4 + [gb["roc_auc"] - gb["roc_auc_ci"][0]], [0] * 4 + [gb["roc_auc_ci"][1] - gb["roc_auc"]]]
    x = np.arange(len(labels))
    ax[0].bar(x, vals, color=colors, width=0.7, yerr=err, capsize=2, error_kw={"lw": 0.9})
    for xi, v in zip(x, vals):
        ax[0].text(xi, v + 0.02, f"{v:.3f}", ha="center", va="bottom", fontsize=6, rotation=90)
    ax[0].set_xticks(x); ax[0].set_xticklabels(labels, fontsize=8)
    ax[0].set_ylim(0.5, 1.12); ax[0].set_yticks([0.5, 0.7, 0.9])
    ax[0].set_ylabel("ROC-AUC", fontsize=8.5)
    ax[0].tick_params(labelsize=8)
    ax[0].set_title("Models (5-fold CV)", fontsize=8.5)

    rows = [(name, auc.get(key, float("nan"))) for name, key in DISPLAY]
    rows.sort(key=lambda t: t[1])
    names = [t[0] for t in rows]
    fv = [t[1] for t in rows]
    y = np.arange(len(names))
    ax[1].barh(y, fv, color=BLUE, height=0.72)
    for yi, v in zip(y, fv):
        ax[1].text(v + 0.008, yi, f"{v:.2f}", va="center", fontsize=7)
    ax[1].set_yticks(y); ax[1].set_yticklabels(names, fontsize=8)
    ax[1].set_xlim(0.5, 1.0); ax[1].set_xticks([0.5, 0.7, 0.9])
    ax[1].set_xlabel("Univariate ROC-AUC", fontsize=8.5)
    ax[1].tick_params(labelsize=8)
    ax[1].set_title("Top single features", fontsize=8.5)

    fig.tight_layout(pad=0.4, w_pad=1.2)
    fig.savefig(OUT, bbox_inches="tight", pad_inches=0.02)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
