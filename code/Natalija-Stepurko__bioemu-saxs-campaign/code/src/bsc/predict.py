"""Stage `predict` (DESIGN §7, D2): can the model's error on a profile be predicted before the
profile is measured? Targets: log10 raw chi2 (regression) and whether the entry is strongly
reweightable or not fit along the path (classification). Features: features.build (sequence
descriptors, disorder, scaling-law Rg, the model's own ensemble Rg and spread). Grouped 5-fold
cross-validation by sequence cluster, repeated 5 times; baselines, linear models and small
gradient-boosted trees; permutation importance and a calibration check.

Output: results/predict.json, results/predict_oof.csv (out-of-fold predictions of the last repeat).
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier, HistGradientBoostingRegressor
from sklearn.inspection import permutation_importance
from sklearn.linear_model import LogisticRegression, RidgeCV
from sklearn.metrics import brier_score_loss, mean_absolute_error, r2_score, roc_auc_score
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from bsc import analyse, features
from bsc import config as C

N_SPLITS, N_REPEATS = 5, 5
HGB = dict(max_depth=3, max_iter=150, learning_rate=0.05, min_samples_leaf=10, l2_regularization=1.0,
           early_stopping=False)
POSITIVE_KINDS = (analyse.KIND_STRONG, analyse.KIND_NOTFIT)


def regressors(seed: int) -> dict:
    return {"ridge": make_pipeline(StandardScaler(), RidgeCV(alphas=np.logspace(-2, 3, 16))),
            "hgb": HistGradientBoostingRegressor(random_state=seed, **HGB)}


def classifiers(seed: int) -> dict:
    return {"logistic": make_pipeline(StandardScaler(), LogisticRegression(C=1.0, max_iter=2000)),
            "hgb": HistGradientBoostingClassifier(random_state=seed, **HGB)}


class ClassMean:
    """Baseline: the training mean (or rate) of the target within each disorder class."""

    def __init__(self, cls_cols):
        self.cls_cols = cls_cols

    def fit(self, X, y):
        self.means_ = {c: float(y[X[c] == 1].mean()) if (X[c] == 1).any() else float(y.mean())
                       for c in self.cls_cols}
        self.global_ = float(y.mean())
        return self

    def predict(self, X):
        out = np.full(len(X), self.global_)
        for c, v in self.means_.items():
            out[X[c].to_numpy() == 1] = v
        return out


def folds(groups: np.ndarray, repeat: int):
    gkf = GroupKFold(n_splits=N_SPLITS, shuffle=True, random_state=C.SEED + repeat)
    return list(gkf.split(np.zeros(len(groups)), groups=groups))


def _summ(values: list[float]) -> dict:
    v = np.asarray(values, float)
    return {"mean": float(v.mean()), "sd": float(v.std(ddof=1)) if len(v) > 1 else 0.0,
            "min": float(v.min()), "max": float(v.max()), "per_repeat": [float(x) for x in v]}


def cross_validate(table: pd.DataFrame, cols: list[str], groups: np.ndarray,
                   n_repeats: int = N_REPEATS) -> tuple[dict, pd.DataFrame]:
    X = table[cols]
    y_reg = np.log10(table["chi2_raw"].to_numpy())
    y_cls = table["kind"].isin(POSITIVE_KINDS).to_numpy().astype(int)
    cls_cols = [c for c in cols if c.startswith("class_")]
    reg_scores = {k: {"r2": [], "mae": []} for k in ("constant", "class_mean", "ridge", "hgb")}
    cls_scores = {k: {"auc": [], "brier": []} for k in ("rate", "class_rate", "logistic", "hgb")}
    importances = []
    oof_last = {}
    for rep in range(n_repeats):
        pred_reg = {k: np.zeros(len(table)) for k in reg_scores}
        pred_cls = {k: np.zeros(len(table)) for k in cls_scores}
        for tr, te in folds(groups, rep):
            Xtr, Xte = X.iloc[tr], X.iloc[te]
            pred_reg["constant"][te] = y_reg[tr].mean()
            pred_reg["class_mean"][te] = ClassMean(cls_cols).fit(Xtr, y_reg[tr]).predict(Xte)
            for k, mdl in regressors(C.SEED + rep).items():
                pred_reg[k][te] = mdl.fit(Xtr, y_reg[tr]).predict(Xte)
                if k == "hgb":
                    pi = permutation_importance(mdl, Xte, y_reg[te], n_repeats=10, random_state=C.SEED + rep,
                                                scoring="r2")
                    importances.append(pi.importances_mean)
            pred_cls["rate"][te] = y_cls[tr].mean()
            pred_cls["class_rate"][te] = ClassMean(cls_cols).fit(Xtr, y_cls[tr]).predict(Xte)
            for k, mdl in classifiers(C.SEED + rep).items():
                pred_cls[k][te] = mdl.fit(Xtr, y_cls[tr]).predict_proba(Xte)[:, 1]
        for k in reg_scores:
            reg_scores[k]["r2"].append(r2_score(y_reg, pred_reg[k]))
            reg_scores[k]["mae"].append(mean_absolute_error(y_reg, pred_reg[k]))
        for k in cls_scores:
            p = np.clip(pred_cls[k], 1e-6, 1 - 1e-6)
            cls_scores[k]["auc"].append(roc_auc_score(y_cls, p) if k not in ("rate",) else 0.5)
            cls_scores[k]["brier"].append(brier_score_loss(y_cls, p))
        oof_last = {"reg": pred_reg, "cls": pred_cls}
    imp = np.mean(importances, axis=0)
    imp_sd = np.std(importances, axis=0)
    order = np.argsort(-imp)
    # calibration of the tree classifier: quintiles of predicted probability against observed rate
    p = oof_last["cls"]["hgb"]
    edges = np.quantile(p, np.linspace(0, 1, 6))
    bins = np.clip(np.searchsorted(edges, p, side="right") - 1, 0, 4)
    calib = [{"bin": int(b), "n": int((bins == b).sum()), "mean_predicted": float(p[bins == b].mean()),
              "observed_rate": float(y_cls[bins == b].mean())} for b in range(5) if (bins == b).any()]
    out = {
        "n": int(len(table)), "n_groups": int(len(np.unique(groups))), "n_splits": N_SPLITS,
        "n_repeats": n_repeats, "features": cols,
        "target_regression": "log10 raw chi2", "target_classification": " or ".join(POSITIVE_KINDS),
        "positive_rate": float(y_cls.mean()),
        "regression": {k: {m: _summ(v) for m, v in s.items()} for k, s in reg_scores.items()},
        "classification": {k: {m: _summ(v) for m, v in s.items()} for k, s in cls_scores.items()},
        "permutation_importance_hgb": [{"feature": cols[i], "mean": float(imp[i]), "sd": float(imp_sd[i])}
                                       for i in order],
        "calibration_hgb": calib,
        "best_regressor": max(("ridge", "hgb"), key=lambda k: np.mean(reg_scores[k]["r2"])),
        "best_classifier": max(("logistic", "hgb"), key=lambda k: np.mean(cls_scores[k]["auc"])),
    }
    oof = pd.DataFrame({"label": table["label"], "disorder_class": table["disorder_class"],
                        "y_log10_chi2": y_reg, "y_positive": y_cls,
                        **{f"pred_{k}": v for k, v in oof_last["reg"].items()},
                        **{f"prob_{k}": v for k, v in oof_last["cls"].items()}})
    return out, oof


def sequence_groups(table: pd.DataFrame) -> np.ndarray:
    cl = pd.read_csv(C.RESULTS / "clusters.csv").set_index("label")
    return cl.loc[table["label"], "group"].to_numpy()


def main() -> None:
    ent, d = analyse.load()
    d["kind"] = d.apply(analyse.kind, axis=1)
    table = features.build(d, ent)
    groups = sequence_groups(table)
    cols = features.feature_columns(table)
    out, oof = cross_validate(table, cols, groups)
    json.dump(out, open(C.RESULTS / "predict.json", "w"), indent=1)
    oof.to_csv(C.RESULTS / "predict_oof.csv", index=False)
    r, c = out["regression"], out["classification"]
    r2 = {k: r[k]["r2"]["mean"] for k in r}
    auc = {k: c[k]["auc"]["mean"] for k in c}
    print(f"  n = {out['n']} in {out['n_groups']} sequence groups; R2 (log10 chi2): class-mean "
          f"{r2['class_mean']:.2f}, ridge {r2['ridge']:.2f}, trees {r2['hgb']:.2f}; AUC (strongly "
          f"reweightable or not fit): class-rate {auc['class_rate']:.2f}, logistic {auc['logistic']:.2f}, "
          f"trees {auc['hgb']:.2f}", flush=True)


if __name__ == "__main__":
    main()
