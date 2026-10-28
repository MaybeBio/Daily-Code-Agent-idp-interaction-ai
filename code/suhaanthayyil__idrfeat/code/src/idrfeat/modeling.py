from __future__ import annotations

import numpy as np
import pandas as pd

from .evaluate import average_precision, bootstrap_ci, group_kfold, roc_auc
from .features import FEATURE_COLUMNS


def default_features(df: pd.DataFrame) -> list[str]:
    return [c for c in FEATURE_COLUMNS if c in df.columns]


def prepare_xy(df: pd.DataFrame, label_col: str, feature_cols: list[str] | None = None):
    cols = feature_cols or default_features(df)
    X = df[cols].to_numpy(dtype=float)
    y = df[label_col].to_numpy(dtype=float)
    groups = df["accession"].to_numpy() if "accession" in df.columns else None
    return X, y, groups, cols


def logistic_factory(seed: int = 1729):
    from sklearn.impute import SimpleImputer
    from sklearn.linear_model import LogisticRegression
    from sklearn.pipeline import Pipeline
    from sklearn.preprocessing import StandardScaler

    return Pipeline(
        [
            ("impute", SimpleImputer(strategy="median")),
            ("scale", StandardScaler()),
            (
                "clf",
                LogisticRegression(
                    penalty="l2", class_weight="balanced", max_iter=2000, random_state=seed
                ),
            ),
        ]
    )


def gbt_factory(seed: int = 1729):
    from sklearn.ensemble import HistGradientBoostingClassifier

    return HistGradientBoostingClassifier(random_state=seed)


def cross_val_predict(
    df: pd.DataFrame,
    label_col: str,
    model_factory,
    feature_cols: list[str] | None = None,
    n_splits: int = 5,
    seed: int = 1729,
):
    X, y, groups, _ = prepare_xy(df, label_col, feature_cols)
    if groups is None:
        raise ValueError("accession column required for group cross-validation")
    oof = np.full(len(y), np.nan)
    fold_rows = []
    for k, (train, test) in enumerate(group_kfold(groups, n_splits, seed)):
        model = model_factory()
        model.fit(X[train], y[train])
        proba = model.predict_proba(X[test])[:, 1]
        oof[test] = proba
        fold_rows.append(
            {
                "fold": k,
                "roc_auc": roc_auc(y[test], proba),
                "pr_auc": average_precision(y[test], proba),
                "n_test": int(len(test)),
            }
        )
    return oof, pd.DataFrame(fold_rows)


def evaluate_model(
    df: pd.DataFrame,
    label_col: str,
    model_factory,
    feature_cols: list[str] | None = None,
    n_splits: int = 5,
    seed: int = 1729,
    n_boot: int = 1000,
) -> dict:
    _, y, groups, _ = prepare_xy(df, label_col, feature_cols)
    oof, folds = cross_val_predict(df, label_col, model_factory, feature_cols, n_splits, seed)
    mask = ~np.isnan(oof)
    g = groups[mask] if groups is not None else None
    return {
        "roc_auc": bootstrap_ci(y[mask], oof[mask], roc_auc, n_boot, seed, g),
        "pr_auc": bootstrap_ci(y[mask], oof[mask], average_precision, n_boot, seed, g),
        "folds": folds,
        "oof": oof,
    }


def permutation_importance_scores(
    df: pd.DataFrame,
    label_col: str,
    model_factory,
    feature_cols: list[str] | None = None,
    seed: int = 1729,
    n_repeats: int = 10,
) -> pd.DataFrame:
    from sklearn.inspection import permutation_importance

    X, y, _, cols = prepare_xy(df, label_col, feature_cols)
    model = model_factory()
    model.fit(X, y)
    result = permutation_importance(
        model, X, y, scoring="roc_auc", n_repeats=n_repeats, random_state=seed
    )
    return (
        pd.DataFrame(
            {
                "feature": cols,
                "importance_mean": result.importances_mean,
                "importance_std": result.importances_std,
            }
        )
        .sort_values("importance_mean", ascending=False)
        .reset_index(drop=True)
    )


def shap_values(df: pd.DataFrame, label_col: str, model_factory, feature_cols=None):
    import shap

    X, y, _, cols = prepare_xy(df, label_col, feature_cols)
    model = model_factory()
    model.fit(X, y)
    explainer = shap.Explainer(model, X)
    return explainer(X), cols
