from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

OUT = Path("results/depmap_essentiality")
BLUE, ORANGE, GREEN, GREY = "#0173B2", "#DE8F05", "#029E73", "#9AA7A5"
plt.rcParams.update({"font.size": 8, "font.family": "DejaVu Sans", "axes.spines.top": False,
                     "axes.spines.right": False, "axes.linewidth": 0.8})


def main() -> None:
    r = json.load(open(OUT / "results.json"))
    m = r["models"]
    ca = r["confound_analysis"]
    so = r["sequence_only"]
    fig, ax = plt.subplots(1, 2, figsize=(7.0, 2.7))

    labels = ["Confounds\n(size, disorder)", "Sequence\nonly", "Logistic", "Random\nforest", "Gradient\nboosting"]
    vals = [ca["confound_only_auc"], so["cv_auc"], m["logistic"]["roc_auc"],
            m["random_forest"]["roc_auc"], m["gradient_boosting"]["roc_auc"]]
    colors = [GREY, ORANGE, BLUE, BLUE, GREEN]
    gb = m["gradient_boosting"]
    err = [[0] * 4 + [gb["roc_auc"] - gb["roc_auc_ci"][0]], [0] * 4 + [gb["roc_auc_ci"][1] - gb["roc_auc"]]]
    x = np.arange(len(labels))
    ax[0].bar(x, vals, color=colors, width=0.68, yerr=err, capsize=2, error_kw={"lw": 0.8})
    for xi, v in zip(x, vals):
        ax[0].text(xi, v + 0.012, f"{v:.3f}", ha="center", va="bottom", fontsize=7.5, fontweight="bold")
    ax[0].axhline(0.5, color="#555", lw=0.7, ls=(0, (4, 3)))
    ax[0].set_xticks(x); ax[0].set_xticklabels(labels, fontsize=7)
    ax[0].set_ylim(0.5, 1.0); ax[0].set_ylabel("ROC-AUC (5-fold CV)")
    ax[0].set_title("Predicting gene essentiality from IDR features", fontsize=8.5)

    feats = r["feature_separation_top"][:10][::-1]
    names = [f["feature"].replace("ptm_", "PTM ").replace("_", " ") for f in feats]
    fv = [f["auc"] for f in feats]
    fc = [ORANGE if f["is_confound"] else BLUE for f in feats]
    y = np.arange(len(names))
    ax[1].barh(y, fv, color=fc, height=0.7)
    for yi, v in zip(y, fv):
        ax[1].text(v + 0.004, yi, f"{v:.2f}", va="center", fontsize=7)
    ax[1].set_yticks(y); ax[1].set_yticklabels(names, fontsize=7)
    ax[1].set_xlim(0.5, 0.95); ax[1].set_xlabel("Univariate separation (ROC-AUC)")
    ax[1].set_title("Top discriminating IDR features", fontsize=8.5)
    fig.tight_layout()
    fig.savefig(OUT / "fig_main.pdf"); fig.savefig(OUT / "fig_main.png", dpi=220)
    print("wrote fig_main.pdf / fig_main.png")


if __name__ == "__main__":
    main()
