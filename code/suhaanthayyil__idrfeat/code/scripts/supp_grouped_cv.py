from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier, RandomForestClassifier
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import StratifiedGroupKFold, cross_val_predict
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

DATA = Path("data/depmap")
OUT = Path("results/depmap_essentiality")
SUPP = OUT / "supp"
SEED = 1729
N_BOOT = 1000


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


def auc_ci(y, p, rng):
    boots = []
    n = len(y)
    for _ in range(N_BOOT):
        idx = rng.integers(0, n, n)
        if len(np.unique(y[idx])) < 2:
            continue
        boots.append(roc_auc_score(y[idx], p[idx]))
    lo, hi = np.percentile(boots, [2.5, 97.5])
    return float(roc_auc_score(y, p)), float(lo), float(hi)


def main() -> None:
    gl = pd.read_csv(DATA / "gene_level_idr_features.csv", index_col=0)
    feats = [c for c in gl.columns if gl[c].dtype != object]
    keep = [c for c in feats if gl[c].notna().mean() >= 0.5]

    ess = load_control_set(DATA / "AchillesCommonEssentialControls.csv")
    non = load_control_set(DATA / "AchillesNonessentialControls.csv")
    label = pd.Series(index=gl.index, dtype="float64")
    label[gl.index.isin(ess)] = 1.0
    label[gl.index.isin(non)] = 0.0
    train = gl[label.notna()].copy()
    y = label[label.notna()].to_numpy().astype(int)

    cl = pd.read_csv(SUPP / "gene_clusters.csv")
    gene2cluster = dict(zip(cl["gene"], cl["cluster"]))
    groups = []
    for i, g in enumerate(train.index):
        groups.append(gene2cluster.get(g, f"__singleton_{i}"))
    groups = np.array(groups)
    codes = pd.Categorical(groups).codes

    n_clusters = len(set(groups))
    vc = pd.Series(groups).value_counts()
    multi_genes = int(vc[vc > 1].sum())
    print(f"train genes {len(y)} | clusters {n_clusters} | "
          f"genes in multi-gene clusters {multi_genes} ({100*multi_genes/len(y):.0f}%)")

    imp = SimpleImputer(strategy="median")
    X = imp.fit_transform(train[keep])
    rng = np.random.default_rng(SEED)

    models = {
        "logistic": make_pipeline(StandardScaler(), LogisticRegression(max_iter=2000, C=1.0)),
        "random_forest": RandomForestClassifier(n_estimators=400, n_jobs=3, random_state=SEED),
        "gradient_boosting": HistGradientBoostingClassifier(random_state=SEED, max_iter=400),
    }
    random_fold = {"logistic": 0.936, "random_forest": 0.954, "gradient_boosting": 0.956}

    cv = StratifiedGroupKFold(5, shuffle=True, random_state=SEED)
    report = {"n_train": int(len(y)), "n_clusters": int(n_clusters),
              "genes_in_multigene_clusters": multi_genes, "min_seq_id": 0.30, "models": {}}
    print(f"\n{'model':20s} {'random-fold':>12s} {'grouped-fold':>13s} {'95% CI':>18s} {'delta':>8s}")
    for name, m in models.items():
        p = cross_val_predict(m, X, y, cv=cv, groups=codes, method="predict_proba")[:, 1]
        auc, lo, hi = auc_ci(y, p, rng)
        pr = average_precision_score(y, p)
        report["models"][name] = {"grouped_roc_auc": auc, "grouped_roc_auc_ci": [lo, hi],
                                  "grouped_pr_auc": float(pr), "random_fold_roc_auc": random_fold[name],
                                  "delta": round(auc - random_fold[name], 3)}
        print(f"{name:20s} {random_fold[name]:>12.3f} {auc:>13.3f} "
              f"{'['+format(lo,'.3f')+', '+format(hi,'.3f')+']':>18s} {auc-random_fold[name]:>+8.3f}")

    json.dump(report, open(SUPP / "grouped_cv.json", "w"), indent=2)
    print(f"\nwrote {SUPP}/grouped_cv.json")


if __name__ == "__main__":
    main()
