from __future__ import annotations

import numpy as np
import pandas as pd

from idrfeat.qc import feature_summary, high_correlation_pairs, near_constant

COLS = ["a", "b", "c", "d"]


def _df() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "a": [1.0, 2.0, 3.0, 4.0],
            "b": [2.0, 4.0, 6.0, 8.0],
            "c": [5.0, 5.0, 5.0, 5.0],
            "d": [1.0, np.nan, 3.0, np.nan],
        }
    )


def test_feature_summary_missing_and_stats() -> None:
    s = feature_summary(_df(), feature_cols=COLS).set_index("feature")
    assert s.loc["d", "missing_frac"] == 0.5
    assert s.loc["a", "min"] == 1.0 and s.loc["a", "max"] == 4.0
    assert s.loc["c", "distinct"] == 1


def test_near_constant_flags_constant_column() -> None:
    assert near_constant(_df(), feature_cols=COLS) == ["c"]


def test_high_correlation_pairs_flags_collinear() -> None:
    pairs = high_correlation_pairs(_df(), feature_cols=COLS, threshold=0.95)
    found = {(r.a, r.b) for r in pairs.itertuples()}
    assert ("a", "b") in found
