from __future__ import annotations

import numpy as np
import pandas as pd

from idrfeat.modeling import (
    evaluate_model,
    gbt_factory,
    logistic_factory,
    permutation_importance_scores,
)


def _synthetic(n_groups: int = 12, per_group: int = 6, seed: int = 0) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    rows = []
    for g in range(n_groups):
        y = g % 2
        for _ in range(per_group):
            rows.append(
                {
                    "accession": f"P{g}",
                    "aa_frac_A": float(y) + rng.normal(0, 0.2),
                    "fcr": rng.normal(0, 1),
                    "label": float(y),
                }
            )
    return pd.DataFrame(rows)


def test_logistic_learns_separable_signal() -> None:
    df = _synthetic()
    out = evaluate_model(
        df, "label", logistic_factory, feature_cols=["aa_frac_A", "fcr"], n_splits=4, n_boot=100
    )
    assert out["roc_auc"][0] > 0.7
    assert len(out["oof"]) == len(df)
    assert set(out["folds"].columns) >= {"fold", "roc_auc", "pr_auc", "n_test"}


def test_gbt_runs_and_scores() -> None:
    df = _synthetic()
    out = evaluate_model(
        df, "label", gbt_factory, feature_cols=["aa_frac_A", "fcr"], n_splits=4, n_boot=50
    )
    assert 0.0 <= out["roc_auc"][0] <= 1.0


def test_permutation_importance_ranks_informative_feature_first() -> None:
    df = _synthetic()
    imp = permutation_importance_scores(
        df, "label", logistic_factory, feature_cols=["aa_frac_A", "fcr"], n_repeats=5
    )
    assert imp.iloc[0]["feature"] == "aa_frac_A"
