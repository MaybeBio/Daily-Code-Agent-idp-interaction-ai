from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.ensemble import HistGradientBoostingClassifier, RandomForestClassifier
from sklearn.impute import SimpleImputer
from sklearn.inspection import permutation_importance
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold

DATA = Path("data/depmap")
SUPP = Path("results/depmap_essentiality/supp")
SUPP.mkdir(parents=True, exist_ok=True)
SEED = 1729

CONFOUNDS = [
    "protein_length", "total_idr_len", "n_idr_segments", "idr_fraction",
    "disorder_mean", "disorder_frac", "metapredict_mean", "metapredict_frac",
    "aiupred_disorder_mean", "plddt_mean",
]


def strip_entrez(name: str) -> str:
    return re.sub(r"\s*\(\d+\)$", "", str(name)).strip()


def load_control_set(path: Path) -> set[str]:
    g = set()
    with open(path) as f:
        next(f)
        for line in f:
            s = line.strip()
            if s:
                g.add(strip_entrez(s))
    return g


def main() -> None:
    gl = pd.read_csv(DATA / "gene_level_idr_features.csv", index_col=0)
    feats = [c for c in gl.columns if gl[c].dtype != object]
    keep = [c for c in feats if gl[c].notna().mean() >= 0.5]

    ess = load_control_set(DATA / "AchillesCommonEssentialControls.csv")
    non = load_control_set(DATA / "AchillesNonessentialControls.csv")
    label = pd.Series(index=gl.index, dtype="float64")
    label[gl.index.isin(ess)] = 1.0
    label[gl.index.isin(non)] = 0.0
    train = gl[label.notna()]
    y = label[label.notna()].to_numpy().astype(int)

    imp = SimpleImputer(strategy="median")
    X = imp.fit_transform(train[keep])

    perm = np.zeros((5, len(keep)))
    cv = StratifiedKFold(5, shuffle=True, random_state=SEED)
    for k, (tr, te) in enumerate(cv.split(X, y)):
        gb = HistGradientBoostingClassifier(random_state=SEED, max_iter=400).fit(X[tr], y[tr])
        r = permutation_importance(gb, X[te], y[te], scoring="roc_auc",
                                   n_repeats=10, random_state=SEED, n_jobs=2)
        perm[k] = r.importances_mean
    perm_mean = perm.mean(0)
    perm_std = perm.std(0)

    rf = RandomForestClassifier(n_estimators=400, random_state=SEED, n_jobs=3).fit(X, y)
    rf_imp = rf.feature_importances_

    uauc = []
    for c in keep:
        x = train[c].to_numpy(dtype=float)
        ok = np.isfinite(x)
        a = roc_auc_score(y[ok], x[ok]) if ok.sum() >= 50 and len(np.unique(x[ok])) >= 3 else np.nan
        uauc.append(max(a, 1 - a) if a == a else np.nan)

    df = pd.DataFrame({
        "feature": keep,
        "role": ["baseline" if c in CONFOUNDS else "character" for c in keep],
        "gb_permutation_importance": np.round(perm_mean, 5),
        "gb_permutation_std": np.round(perm_std, 5),
        "rf_impurity_importance": np.round(rf_imp, 5),
        "univariate_auc": np.round(uauc, 3),
    }).sort_values("gb_permutation_importance", ascending=False)
    df.to_csv(SUPP / "feature_importance.csv", index=False)

    top = df.head(15).iloc[::-1]
    fig, ax = plt.subplots(figsize=(7.2, 6))
    ax.barh(top["feature"], top["gb_permutation_importance"],
            xerr=top["gb_permutation_std"], color="#0f766e", capsize=2)
    ax.set_xlabel("Gradient-boosting permutation importance (drop in ROC-AUC)")
    ax.set_title("Top 15 features by permutation importance")
    fig.tight_layout()
    fig.savefig(SUPP / "feature_importance_top15.png", dpi=170)
    fig.savefig(SUPP / "feature_importance_top15.pdf")

    print("Top 12 by gradient-boosting permutation importance:")
    print(f"{'feature':26s} {'perm':>8s} {'rf_imp':>8s} {'uni_auc':>8s}")
    for _, r in df.head(12).iterrows():
        print(f"{r['feature']:26s} {r['gb_permutation_importance']:>8.4f} "
              f"{r['rf_impurity_importance']:>8.4f} {r['univariate_auc']:>8.3f}")
    print(f"\nwrote {SUPP}/feature_importance.csv and bar figure")


if __name__ == "__main__":
    main()
