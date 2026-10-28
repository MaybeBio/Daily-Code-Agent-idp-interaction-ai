from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.ensemble import HistGradientBoostingClassifier, RandomForestClassifier
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

DATA = Path("data/depmap")
OUT = Path("results/depmap_essentiality")
SEED = 1729
N_BOOT = 1000
N_PERM = 300

CONFOUNDS = [
    "protein_length", "total_idr_len", "n_idr_segments", "idr_fraction",
    "disorder_mean", "disorder_frac", "metapredict_mean", "metapredict_frac",
    "aiupred_disorder_mean", "plddt_mean",
]

ANNOTATION_PREFIXES = ("elm_", "ptm_", "phylop", "phastcons")


def strip_entrez(name: str) -> str:
    return re.sub(r"\s*\(\d+\)$", "", str(name)).strip()


def load_control_set(path: Path) -> set[str]:
    genes = set()
    with open(path) as f:
        next(f)
        for line in f:
            s = line.strip()
            if s:
                genes.add(strip_entrez(s))
    return genes


def load_hart(path: Path) -> set[str]:
    d = pd.read_csv(path, sep="\t")
    return set(d.iloc[:, 0].astype(str))


def auc_ci(y: np.ndarray, p: np.ndarray, rng: np.random.Generator) -> tuple[float, float, float]:
    base = roc_auc_score(y, p)
    boots = []
    n = len(y)
    for _ in range(N_BOOT):
        idx = rng.integers(0, n, n)
        if len(np.unique(y[idx])) < 2:
            continue
        boots.append(roc_auc_score(y[idx], p[idx]))
    lo, hi = np.percentile(boots, [2.5, 97.5])
    return float(base), float(lo), float(hi)


def oof_proba(model, X: np.ndarray, y: np.ndarray) -> np.ndarray:
    cv = StratifiedKFold(5, shuffle=True, random_state=SEED)
    return cross_val_predict(model, X, y, cv=cv, method="predict_proba")[:, 1]


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    gl = pd.read_csv(DATA / "gene_level_idr_features.csv", index_col=0)
    eff = pd.read_csv(DATA / "mean_gene_effect.csv", index_col=0)["mean_gene_effect"]

    ess = load_control_set(DATA / "AchillesCommonEssentialControls.csv")
    non = load_control_set(DATA / "AchillesNonessentialControls.csv")
    feats = [c for c in gl.columns if gl[c].dtype != object]
    keep = [c for c in feats if gl[c].notna().mean() >= 0.5]
    confound_cols = [c for c in CONFOUNDS if c in keep]
    character_cols = [c for c in keep if c not in confound_cols]

    label = pd.Series(index=gl.index, dtype="float64")
    label[gl.index.isin(ess)] = 1.0
    label[gl.index.isin(non)] = 0.0
    train = gl[label.notna()].copy()
    ytr = label[label.notna()].to_numpy().astype(int)
    print(f"train genes {len(ytr)} (pos {ytr.sum()}, neg {(ytr==0).sum()})")

    imp = SimpleImputer(strategy="median")
    Xall = imp.fit_transform(train[keep])
    Xconf = imp.fit_transform(train[confound_cols])
    Xchar = imp.fit_transform(train[character_cols])
    rng = np.random.default_rng(SEED)

    models = {
        "logistic": make_pipeline(StandardScaler(), LogisticRegression(max_iter=2000, C=1.0)),
        "random_forest": RandomForestClassifier(n_estimators=400, n_jobs=-1, random_state=SEED),
        "gradient_boosting": HistGradientBoostingClassifier(random_state=SEED, max_iter=400),
    }
    report: dict = {"n_train": int(len(ytr)), "n_pos": int(ytr.sum()), "n_neg": int((ytr == 0).sum()),
                    "n_features": len(keep), "n_confounds": len(confound_cols),
                    "n_character": len(character_cols), "models": {}}

    oof = {}
    for name, m in models.items():
        p = oof_proba(m, Xall, ytr)
        oof[name] = p
        auc, lo, hi = auc_ci(ytr, p, rng)
        pr = average_precision_score(ytr, p)
        report["models"][name] = {"roc_auc": auc, "roc_auc_ci": [lo, hi], "pr_auc": float(pr)}
        print(f"  {name:18s} ROC-AUC {auc:.3f} [{lo:.3f}, {hi:.3f}]  PR-AUC {pr:.3f}")

    conf_model = make_pipeline(StandardScaler(), LogisticRegression(max_iter=2000))
    p_conf = oof_proba(HistGradientBoostingClassifier(random_state=SEED, max_iter=400), Xconf, ytr)
    p_full = oof["gradient_boosting"]
    p_char = oof_proba(HistGradientBoostingClassifier(random_state=SEED, max_iter=400), Xchar, ytr)
    auc_conf = roc_auc_score(ytr, p_conf)
    auc_full = roc_auc_score(ytr, p_full)
    auc_char = roc_auc_score(ytr, p_char)
    deltas = []
    n = len(ytr)
    for _ in range(N_BOOT):
        idx = rng.integers(0, n, n)
        if len(np.unique(ytr[idx])) < 2:
            continue
        deltas.append(roc_auc_score(ytr[idx], p_full[idx]) - roc_auc_score(ytr[idx], p_conf[idx]))
    dlo, dhi = np.percentile(deltas, [2.5, 97.5])
    report["confound_analysis"] = {
        "confound_only_auc": float(auc_conf), "character_only_auc": float(auc_char),
        "full_auc": float(auc_full), "delta_full_minus_confound": float(auc_full - auc_conf),
        "delta_ci": [float(dlo), float(dhi)],
    }
    print(f"  confound-only {auc_conf:.3f} | character-only {auc_char:.3f} | full {auc_full:.3f} "
          f"| delta {auc_full-auc_conf:+.3f} [{dlo:+.3f},{dhi:+.3f}]")

    perm = []
    gbt = HistGradientBoostingClassifier(random_state=SEED, max_iter=300)
    for _ in range(N_PERM):
        yp = rng.permutation(ytr)
        pp = oof_proba(gbt, Xall, yp)
        perm.append(roc_auc_score(yp, pp))
    perm = np.array(perm)
    pval = (np.sum(perm >= auc_full) + 1) / (N_PERM + 1)
    report["permutation_null"] = {"true_auc": float(auc_full), "null_max": float(perm.max()),
                                  "null_mean": float(perm.mean()), "p_value": float(pval), "n_perm": N_PERM}
    print(f"  permutation: true {auc_full:.3f} vs null max {perm.max():.3f} (mean {perm.mean():.3f}), p={pval:.4f}")

    others = gl.index.difference(train.index).intersection(eff.index)
    imp_gen = SimpleImputer(strategy="median")
    Xtr_ctrl = imp_gen.fit_transform(train[keep])
    gbt_gen = HistGradientBoostingClassifier(random_state=SEED, max_iter=400).fit(Xtr_ctrl, ytr)
    p_oth = gbt_gen.predict_proba(imp_gen.transform(gl.loc[others, keep]))[:, 1]
    eff_oth = eff.loc[others].to_numpy()
    rho_gen = spearmanr(p_oth, -eff_oth).statistic
    dep = eff_oth < -0.5
    nondep = eff_oth > -0.1
    mask = dep | nondep
    auc_gen = roc_auc_score(dep[mask].astype(int), p_oth[mask])
    report["heldout_generalization"] = {
        "n": int(len(others)), "spearman": float(rho_gen),
        "n_dependency": int(dep.sum()), "n_nondependency": int(nondep.sum()),
        "auc": float(auc_gen),
    }
    print(f"  held-out proteome: n={len(others)} Spearman {rho_gen:.3f} | "
          f"AUC {auc_gen:.3f} (dep {int(dep.sum())}, non {int(nondep.sum())})")

    seq_cols = [c for c in keep if not c.startswith(ANNOTATION_PREFIXES)]
    imp_seq = SimpleImputer(strategy="median")
    Xseq = imp_seq.fit_transform(train[seq_cols])
    p_seq = oof_proba(HistGradientBoostingClassifier(random_state=SEED, max_iter=400), Xseq, ytr)
    auc_seq_cv = roc_auc_score(ytr, p_seq)
    report["sequence_only"] = {"n_features": len(seq_cols), "cv_auc": float(auc_seq_cv)}
    print(f"  sequence-only (no annotations): n_feat={len(seq_cols)} CV {auc_seq_cv:.3f}")

    common = gl.index.intersection(eff.index)
    Xc = imp.fit_transform(gl.loc[common, keep])
    yc = eff.loc[common].to_numpy()
    from sklearn.ensemble import HistGradientBoostingRegressor
    from sklearn.model_selection import KFold, cross_val_predict as cvp
    pred = cvp(HistGradientBoostingRegressor(random_state=SEED, max_iter=400), Xc, yc,
               cv=KFold(5, shuffle=True, random_state=SEED))
    rho = spearmanr(pred, yc).statistic
    report["continuous_regression"] = {"n": int(len(yc)), "spearman": float(rho)}
    print(f"  continuous dependency regression: n={len(yc)} Spearman {rho:.3f}")

    univ = []
    for c in character_cols + confound_cols:
        x = train[c].to_numpy(dtype=float)
        ok = np.isfinite(x)
        if ok.sum() < 50 or len(np.unique(x[ok])) < 3:
            continue
        a = roc_auc_score(ytr[ok], x[ok])
        univ.append({"feature": c, "auc": float(max(a, 1 - a)), "direction": "up" if a >= 0.5 else "down",
                     "is_confound": c in confound_cols})
    univ.sort(key=lambda r: -r["auc"])
    report["feature_separation_top"] = univ[:20]
    print("  top separating features:", [(u["feature"], round(u["auc"], 3)) for u in univ[:8]])

    json.dump(report, open(OUT / "results.json", "w"), indent=2)
    print(f"\nwrote {OUT}/results.json")


if __name__ == "__main__":
    main()
