from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

BASELINE_FEATURES = ["seg_length", "disorder_mean", "disorder_frac", "aiupred_binding_mean"]


def _clean(y, s):
    y = np.asarray(y, dtype=float)
    s = np.asarray(s, dtype=float)
    m = ~(np.isnan(y) | np.isnan(s))
    return y[m], s[m]


def _avg_ranks(sorted_vals: np.ndarray) -> np.ndarray:
    n = len(sorted_vals)
    ranks = np.empty(n, dtype=float)
    i = 0
    while i < n:
        j = i
        while j + 1 < n and sorted_vals[j + 1] == sorted_vals[i]:
            j += 1
        ranks[i : j + 1] = (i + j) / 2.0 + 1.0
        i = j + 1
    return ranks


def roc_auc(y, s) -> float:
    y, s = _clean(y, s)
    n_pos = int((y == 1).sum())
    n_neg = int((y == 0).sum())
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    order = np.argsort(s, kind="mergesort")
    ranks = np.empty(len(s), dtype=float)
    ranks[order] = _avg_ranks(s[order])
    return float((ranks[y == 1].sum() - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg))


def average_precision(y, s) -> float:
    y, s = _clean(y, s)
    n_pos = int((y == 1).sum())
    if n_pos == 0:
        return float("nan")
    order = np.argsort(-s, kind="mergesort")
    y = y[order]
    tp = np.cumsum(y)
    fp = np.cumsum(1.0 - y)
    precision = tp / (tp + fp)
    recall = tp / n_pos
    prev = np.concatenate([[0.0], recall[:-1]])
    return float(np.sum((recall - prev) * precision))


def brier_score(y, p) -> float:
    y, p = _clean(y, p)
    if len(y) == 0:
        return float("nan")
    return float(np.mean((p - y) ** 2))


def calibration_table(y, p, bins: int = 10) -> pd.DataFrame:
    y, p = _clean(y, p)
    edges = np.linspace(0.0, 1.0, bins + 1)
    idx = np.clip(np.digitize(p, edges[1:-1]), 0, bins - 1)
    rows = []
    for b in range(bins):
        m = idx == b
        count = int(m.sum())
        rows.append(
            {
                "bin_lo": float(edges[b]),
                "bin_hi": float(edges[b + 1]),
                "mean_pred": float(p[m].mean()) if count else float("nan"),
                "obs_rate": float(y[m].mean()) if count else float("nan"),
                "count": count,
            }
        )
    return pd.DataFrame(rows)


def group_kfold(groups, n_splits: int = 5, seed: int = 1729) -> list[tuple[np.ndarray, np.ndarray]]:
    groups = np.asarray(groups)
    uniq = np.array(sorted(set(groups.tolist())))
    if n_splits > len(uniq):
        raise ValueError("n_splits exceeds number of groups")
    perm = np.random.default_rng(seed).permutation(len(uniq))
    fold_of = {}
    for order_pos, gi in enumerate(perm):
        fold_of[uniq[gi]] = order_pos % n_splits
    fold = np.array([fold_of[g] for g in groups])
    folds = []
    for k in range(n_splits):
        test = np.where(fold == k)[0]
        train = np.where(fold != k)[0]
        folds.append((train, test))
    return folds


def bootstrap_ci(
    y,
    s,
    metric=roc_auc,
    n: int = 1000,
    seed: int = 1729,
    groups=None,
    alpha: float = 0.05,
) -> tuple[float, float, float]:
    y = np.asarray(y, dtype=float)
    s = np.asarray(s, dtype=float)
    point = metric(y, s)
    rng = np.random.default_rng(seed)
    stats = []
    if groups is None:
        idx_all = np.arange(len(y))
        for _ in range(n):
            idx = rng.choice(idx_all, size=len(idx_all), replace=True)
            stats.append(metric(y[idx], s[idx]))
    else:
        groups = np.asarray(groups)
        uniq = np.array(sorted(set(groups.tolist())))
        by = {g: np.where(groups == g)[0] for g in uniq}
        for _ in range(n):
            pick = rng.choice(uniq, size=len(uniq), replace=True)
            idx = np.concatenate([by[g] for g in pick])
            stats.append(metric(y[idx], s[idx]))
    stats = np.asarray([v for v in stats if not np.isnan(v)], dtype=float)
    if stats.size == 0:
        return float(point), float("nan"), float("nan")
    lo, hi = np.quantile(stats, [alpha / 2.0, 1.0 - alpha / 2.0])
    return float(point), float(lo), float(hi)


def evaluate_baselines(
    df: pd.DataFrame,
    y,
    features: list[str] | None = None,
    groups=None,
    n_boot: int = 1000,
    seed: int = 1729,
) -> pd.DataFrame:
    y = np.asarray(y, dtype=float)
    if groups is None and "accession" in df.columns:
        groups = df["accession"].to_numpy()
    cols = features if features is not None else [c for c in BASELINE_FEATURES if c in df.columns]
    rows = []
    for col in cols:
        s = df[col].to_numpy(dtype=float)
        auc, auc_lo, auc_hi = bootstrap_ci(y, s, roc_auc, n_boot, seed, groups)
        ap, ap_lo, ap_hi = bootstrap_ci(y, s, average_precision, n_boot, seed, groups)
        rows.append(
            {
                "feature": col,
                "roc_auc": auc,
                "roc_auc_lo": auc_lo,
                "roc_auc_hi": auc_hi,
                "pr_auc": ap,
                "pr_auc_lo": ap_lo,
                "pr_auc_hi": ap_hi,
                "n": int(np.isfinite(s).sum()),
            }
        )
    return pd.DataFrame(rows).sort_values("roc_auc", ascending=False).reset_index(drop=True)


def freeze_holdout(
    accessions,
    path: str | Path,
    frac: float = 0.15,
    seed: int = 1729,
    force: bool = False,
) -> list[str]:
    path = Path(path)
    if path.exists() and not force:
        raise FileExistsError(f"{path} exists; holdout is frozen, pass force=True to overwrite")
    accs = sorted(set(accessions))
    k = int(round(len(accs) * frac))
    idx = np.random.default_rng(seed).permutation(len(accs))[:k]
    hold = sorted(accs[i] for i in idx)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(hold) + "\n")
    return hold


def load_holdout(path: str | Path) -> set[str]:
    return {line.strip() for line in Path(path).read_text().splitlines() if line.strip()}


def apply_holdout(df: pd.DataFrame, holdout: set[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    mask = df["accession"].isin(holdout)
    dev = df[~mask].reset_index(drop=True)
    test = df[mask].reset_index(drop=True)
    return dev, test
