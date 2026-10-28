from __future__ import annotations

import numpy as np
import pandas as pd

from .features import FEATURE_COLUMNS


def _cols(df: pd.DataFrame, feature_cols: list[str] | None) -> list[str]:
    return feature_cols or [c for c in FEATURE_COLUMNS if c in df.columns]


def feature_summary(df: pd.DataFrame, feature_cols: list[str] | None = None) -> pd.DataFrame:
    cols = _cols(df, feature_cols)
    n = len(df)
    rows = []
    for c in cols:
        v = df[c].to_numpy(dtype=float)
        finite = np.isfinite(v)
        vals = v[finite]
        rows.append(
            {
                "feature": c,
                "n": int(finite.sum()),
                "missing_frac": float(1.0 - finite.mean()) if n else float("nan"),
                "mean": float(vals.mean()) if vals.size else float("nan"),
                "std": float(vals.std()) if vals.size else float("nan"),
                "min": float(vals.min()) if vals.size else float("nan"),
                "max": float(vals.max()) if vals.size else float("nan"),
                "distinct": int(np.unique(vals).size),
            }
        )
    return pd.DataFrame(rows)


def near_constant(
    df: pd.DataFrame, feature_cols: list[str] | None = None, std_tol: float = 1e-9
) -> list[str]:
    s = feature_summary(df, feature_cols)
    hit = s[(s["std"].fillna(0.0) <= std_tol) | (s["distinct"] <= 1)]
    return hit["feature"].tolist()


def high_correlation_pairs(
    df: pd.DataFrame, feature_cols: list[str] | None = None, threshold: float = 0.95
) -> pd.DataFrame:
    cols = _cols(df, feature_cols)
    corr = df[cols].corr(numeric_only=True).abs()
    names = list(corr.columns)
    pairs = []
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            r = corr.iloc[i, j]
            if pd.notna(r) and r >= threshold:
                pairs.append({"a": names[i], "b": names[j], "abs_corr": float(r)})
    if not pairs:
        return pd.DataFrame(columns=["a", "b", "abs_corr"])
    return pd.DataFrame(pairs).sort_values("abs_corr", ascending=False).reset_index(drop=True)
