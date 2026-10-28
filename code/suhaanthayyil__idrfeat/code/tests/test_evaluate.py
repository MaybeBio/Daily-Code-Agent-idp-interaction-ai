from __future__ import annotations

import math

import numpy as np
import pandas as pd
import pytest

from idrfeat.evaluate import (
    apply_holdout,
    average_precision,
    bootstrap_ci,
    brier_score,
    calibration_table,
    evaluate_baselines,
    freeze_holdout,
    group_kfold,
    load_holdout,
    roc_auc,
)


def test_roc_auc_perfect_and_reversed() -> None:
    y = [0, 0, 1, 1]
    s = [0.1, 0.2, 0.8, 0.9]
    assert roc_auc(y, s) == 1.0
    assert roc_auc(y, [-v for v in s]) == 0.0


def test_roc_auc_handles_ties() -> None:
    y = [0, 1, 0, 1]
    s = [0.5, 0.5, 0.5, 0.5]
    assert math.isclose(roc_auc(y, s), 0.5)


def test_roc_auc_nan_dropped_pairwise() -> None:
    y = [0, 1, 1, 0]
    s = [0.1, np.nan, 0.9, 0.2]
    assert roc_auc(y, s) == 1.0


def test_average_precision_perfect() -> None:
    assert math.isclose(average_precision([0, 0, 1, 1], [0.1, 0.2, 0.8, 0.9]), 1.0)


def test_average_precision_all_negative_is_nan() -> None:
    assert math.isnan(average_precision([0, 0, 0], [0.1, 0.2, 0.3]))


def test_brier_and_calibration() -> None:
    y = [0, 0, 1, 1]
    p = [0.0, 0.0, 1.0, 1.0]
    assert brier_score(y, p) == 0.0
    table = calibration_table(y, p, bins=5)
    assert int(table["count"].sum()) == 4


def test_group_kfold_no_group_leak_and_full_cover() -> None:
    groups = np.array([f"P{i // 3}" for i in range(30)])
    folds = group_kfold(groups, n_splits=5, seed=1729)
    covered = set()
    for train, test in folds:
        assert set(groups[train]).isdisjoint(set(groups[test]))
        covered.update(test.tolist())
    assert covered == set(range(30))


def test_group_kfold_rejects_too_many_splits() -> None:
    with pytest.raises(ValueError):
        group_kfold(["a", "a", "b"], n_splits=5)


def test_bootstrap_ci_brackets_point() -> None:
    rng = np.random.default_rng(0)
    y = rng.integers(0, 2, size=200)
    s = y + rng.normal(0, 1, size=200)
    point, lo, hi = bootstrap_ci(y, s, roc_auc, n=200, seed=1729)
    assert lo <= point <= hi


def test_bootstrap_ci_grouped_runs() -> None:
    rng = np.random.default_rng(0)
    groups = np.array([f"P{i // 4}" for i in range(80)])
    y = rng.integers(0, 2, size=80)
    s = rng.normal(0, 1, size=80)
    point, lo, hi = bootstrap_ci(y, s, roc_auc, n=100, seed=1, groups=groups)
    assert lo <= point <= hi


def test_evaluate_baselines_shape() -> None:
    df = pd.DataFrame(
        {
            "accession": ["A", "A", "B", "B", "C", "C"],
            "seg_length": [10, 20, 30, 40, 50, 60],
            "disorder_mean": [0.5, 0.6, 0.7, 0.8, 0.9, 0.95],
        }
    )
    y = [0, 0, 0, 1, 1, 1]
    out = evaluate_baselines(df, y, features=["seg_length", "disorder_mean"], n_boot=50)
    assert list(out["feature"]) and set(out["feature"]) == {"seg_length", "disorder_mean"}
    assert {"roc_auc", "pr_auc", "roc_auc_lo", "roc_auc_hi"}.issubset(out.columns)


def test_freeze_holdout_deterministic_disjoint_and_frozen(tmp_path) -> None:
    accs = [f"P{i}" for i in range(100)]
    path = tmp_path / "holdout.txt"
    first = freeze_holdout(accs, path, frac=0.2, seed=1729)
    assert len(first) == 20
    with pytest.raises(FileExistsError):
        freeze_holdout(accs, path, frac=0.2, seed=1729)
    again = freeze_holdout(accs, path, frac=0.2, seed=1729, force=True)
    assert first == again
    assert load_holdout(path) == set(first)


def test_apply_holdout_splits_by_accession(tmp_path) -> None:
    df = pd.DataFrame({"accession": ["A", "A", "B", "C"], "x": [1, 2, 3, 4]})
    dev, test = apply_holdout(df, {"B"})
    assert set(dev["accession"]) == {"A", "C"}
    assert set(test["accession"]) == {"B"}
    assert len(dev) + len(test) == len(df)
