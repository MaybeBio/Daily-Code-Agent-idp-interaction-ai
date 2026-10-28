from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import squareform

DATA = Path("data/depmap")
SUPP = Path("results/depmap_essentiality/supp")
SUPP.mkdir(parents=True, exist_ok=True)

CONFOUNDS = [
    "protein_length", "total_idr_len", "n_idr_segments", "idr_fraction",
    "disorder_mean", "disorder_frac", "metapredict_mean", "metapredict_frac",
    "aiupred_disorder_mean", "plddt_mean",
]


def main() -> None:
    gl = pd.read_csv(DATA / "gene_level_idr_features.csv", index_col=0)
    feats = [c for c in gl.columns if gl[c].dtype != object]
    keep = [c for c in feats if gl[c].notna().mean() >= 0.5]

    X = gl[keep].copy()
    X = X.fillna(X.median(numeric_only=True))
    corr = X.corr(method="spearman")
    corr.to_csv(SUPP / "feature_correlation.csv")

    d = 1.0 - corr.abs().fillna(0.0).to_numpy()
    d = np.nan_to_num(d, nan=1.0, posinf=1.0, neginf=1.0)
    np.fill_diagonal(d, 0.0)
    d = (d + d.T) / 2.0
    Z = linkage(squareform(d, checks=False), method="average")
    order = leaves_list(Z)
    cols = [keep[i] for i in order]
    corr_o = corr.loc[cols, cols]

    pairs = []
    arr = corr.to_numpy()
    for i in range(len(keep)):
        for j in range(i + 1, len(keep)):
            r = arr[i, j]
            if abs(r) >= 0.8:
                pairs.append((keep[i], keep[j], round(float(r), 3)))
    pairs.sort(key=lambda t: -abs(t[2]))
    pd.DataFrame(pairs, columns=["feature_a", "feature_b", "spearman"]).to_csv(
        SUPP / "feature_correlation_strong_pairs.csv", index=False)

    from scipy.cluster.hierarchy import fcluster
    cl = fcluster(Z, t=0.2, criterion="distance")
    groups: dict[int, list[str]] = {}
    for f, c in zip(keep, cl):
        groups.setdefault(int(c), []).append(f)
    multi = {c: g for c, g in groups.items() if len(g) > 1}

    n = len(cols)
    fig, ax = plt.subplots(figsize=(11, 9.5))
    im = ax.imshow(corr_o.to_numpy(), vmin=-1, vmax=1, cmap="RdBu_r", aspect="equal")
    ax.set_xticks(range(n)); ax.set_xticklabels(cols, rotation=90, fontsize=5.3)
    ax.set_yticks(range(n)); ax.set_yticklabels(cols, fontsize=5.3)
    ax.set_title("Spearman correlation among the 61 features (hierarchically ordered)", fontsize=11)
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label("Spearman rho", fontsize=9)
    fig.tight_layout()
    fig.savefig(SUPP / "feature_correlation_heatmap.png", dpi=170)
    fig.savefig(SUPP / "feature_correlation_heatmap.pdf")

    off = arr[np.triu_indices(len(keep), k=1)]
    print(f"features {len(keep)} | mean |rho| off-diagonal {np.mean(np.abs(off)):.3f} | "
          f"pairs with |rho|>=0.8: {len(pairs)}")
    print(f"correlated clusters (|rho|>=0.8), {len(multi)} groups:")
    for c, g in sorted(multi.items(), key=lambda kv: -len(kv[1])):
        print("  -", ", ".join(g))
    print("\nTop strong pairs:")
    for a, b, r in pairs[:12]:
        print(f"  {a:26s} {b:26s} {r:+.3f}")
    print(f"\nwrote correlation CSV, strong-pairs CSV, and heatmap (png/pdf) to {SUPP}")


if __name__ == "__main__":
    main()
