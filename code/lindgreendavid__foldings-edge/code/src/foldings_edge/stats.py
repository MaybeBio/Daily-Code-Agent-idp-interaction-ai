"""Statistical tests for the pLDDT vs. DisProt-curated-disorder comparison.

See `docs/research-protocol.md` for the preregistered method: two-sample
Mann-Whitney U as the primary distributional test (H1), two-sample
Kolmogorov-Smirnov as a robustness check, a percentile bootstrap CI on the
median difference, and precision/recall/F1/MCC (with Wilson and bootstrap
confidence intervals) for the fixed pLDDT<70 classifier (H2).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy import stats as scipy_stats

DEFAULT_SEED = 20260813
DEFAULT_BOOTSTRAP_RESAMPLES = 10_000
ALPHA = 0.05
PLDDT_CONFIDENT_THRESHOLD = 70.0


@dataclass(frozen=True)
class TestResult:
    statistic: float
    p_value: float
    significant: bool


@dataclass(frozen=True)
class BootstrapResult:
    median_a: float
    median_b: float
    median_diff: float
    ci_low: float
    ci_high: float
    resamples: int
    seed: int


@dataclass(frozen=True)
class ProportionCi:
    point: float
    ci_low: float
    ci_high: float


@dataclass(frozen=True)
class ClassifierMetrics:
    n: int
    true_positive: int
    false_positive: int
    true_negative: int
    false_negative: int
    precision: ProportionCi
    recall: ProportionCi
    f1: ProportionCi
    mcc: ProportionCi


def mann_whitney_two_sample(a: list[float], b: list[float]) -> TestResult:
    """Two-sided two-sample Mann-Whitney U test."""
    result = scipy_stats.mannwhitneyu(a, b, alternative="two-sided")
    return TestResult(
        statistic=float(result.statistic),
        p_value=float(result.pvalue),
        significant=bool(result.pvalue < ALPHA),
    )


def ks_two_sample(a: list[float], b: list[float]) -> TestResult:
    """Two-sided two-sample Kolmogorov-Smirnov test (robustness check for H1)."""
    result = scipy_stats.ks_2samp(a, b)
    return TestResult(
        statistic=float(result.statistic),
        p_value=float(result.pvalue),
        significant=bool(result.pvalue < ALPHA),
    )


def bootstrap_median_diff(
    a: list[float],
    b: list[float],
    *,
    resamples: int = DEFAULT_BOOTSTRAP_RESAMPLES,
    seed: int = DEFAULT_SEED,
) -> BootstrapResult:
    """Percentile bootstrap 95% CI for (median(b) - median(a)).

    Each group is resampled independently, with replacement, at its own
    observed size, using a seeded generator for exact reproducibility.
    """
    rng = np.random.default_rng(seed)
    arr_a = np.asarray(a, dtype=float)
    arr_b = np.asarray(b, dtype=float)
    # Vectorized, but chunked, resampling: draw (chunk x n) index matrices in
    # bounded-size batches rather than either looping resample-by-resample in
    # Python (too slow at n in the hundreds of thousands) or materializing the
    # full (resamples x n) matrix at once (too much memory at the same scale).
    max_cells = 20_000_000  # bounds a chunk's index matrix to ~160 MB (int64)
    chunk = max(1, min(resamples, max_cells // max(arr_a.size, arr_b.size, 1)))
    diffs = np.empty(resamples, dtype=float)
    start = 0
    while start < resamples:
        size = min(chunk, resamples - start)
        idx_a = rng.integers(0, arr_a.size, size=(size, arr_a.size))
        idx_b = rng.integers(0, arr_b.size, size=(size, arr_b.size))
        medians_a = np.median(arr_a[idx_a], axis=1)
        medians_b = np.median(arr_b[idx_b], axis=1)
        diffs[start : start + size] = medians_b - medians_a
        start += size
    ci_low, ci_high = np.percentile(diffs, [2.5, 97.5])
    return BootstrapResult(
        median_a=float(np.median(arr_a)),
        median_b=float(np.median(arr_b)),
        median_diff=float(np.median(arr_b) - np.median(arr_a)),
        ci_low=float(ci_low),
        ci_high=float(ci_high),
        resamples=resamples,
        seed=seed,
    )


def wilson_ci(successes: int, n: int, *, z: float = 1.959963984540054) -> ProportionCi:
    """Wilson score interval for a binomial proportion."""
    if n == 0:
        return ProportionCi(point=float("nan"), ci_low=float("nan"), ci_high=float("nan"))
    p_hat = successes / n
    denom = 1 + z**2 / n
    center = (p_hat + z**2 / (2 * n)) / denom
    half_width = (z * math.sqrt((p_hat * (1 - p_hat) + z**2 / (4 * n)) / n)) / denom
    return ProportionCi(
        point=p_hat, ci_low=max(0.0, center - half_width), ci_high=min(1.0, center + half_width)
    )


def _f1(precision: float, recall: float) -> float:
    if precision + recall == 0:
        return 0.0
    return 2 * precision * recall / (precision + recall)


def _mcc(tp: int, fp: int, tn: int, fn: int) -> float:
    numerator = tp * tn - fp * fn
    denom_sq = (tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)
    if denom_sq == 0:
        return 0.0
    return numerator / math.sqrt(denom_sq)


def classifier_metrics(
    is_disorder: list[bool],
    plddt: list[float],
    *,
    threshold: float = PLDDT_CONFIDENT_THRESHOLD,
    resamples: int = DEFAULT_BOOTSTRAP_RESAMPLES,
    seed: int = DEFAULT_SEED,
) -> ClassifierMetrics:
    """Precision/recall/F1/MCC for the fixed pLDDT<threshold "predicts disorder" classifier.

    Positive class: residue is inside a DisProt-curated disorder region.
    Predicted positive: pLDDT < threshold (i.e. "not confident").
    """
    actual = np.asarray(is_disorder, dtype=bool)
    predicted = np.asarray(plddt, dtype=float) < threshold

    tp = int(np.sum(predicted & actual))
    fp = int(np.sum(predicted & ~actual))
    tn = int(np.sum(~predicted & ~actual))
    fn = int(np.sum(~predicted & actual))
    n = len(actual)

    precision_point = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    recall_point = tp / (tp + fn) if (tp + fn) > 0 else 0.0

    precision_ci = wilson_ci(tp, tp + fp) if (tp + fp) > 0 else ProportionCi(0.0, 0.0, 0.0)
    recall_ci = wilson_ci(tp, tp + fn) if (tp + fn) > 0 else ProportionCi(0.0, 0.0, 0.0)

    # Bootstrap resampling of n residues with replacement is exactly equivalent
    # (each residue has a fixed, deterministic tp/fp/tn/fn category) to drawing
    # category counts from Multinomial(n, [p_tp, p_fp, p_tn, p_fn]), which is
    # what is done here. This avoids materializing an (resamples x n) index or
    # boolean matrix, which is intractable at n in the hundreds of thousands
    # and resamples=10,000.
    rng = np.random.default_rng(seed)
    probs = np.array([tp, fp, tn, fn], dtype=float) / n if n > 0 else np.zeros(4)
    counts = rng.multinomial(n, probs, size=resamples) if n > 0 else np.zeros((resamples, 4))
    # float64, not int64: the MCC denominator multiplies four counts together,
    # which overflows int64 once n is in the hundreds of thousands.
    counts = counts.astype(np.float64)
    s_tp, s_fp, s_tn, s_fn = counts[:, 0], counts[:, 1], counts[:, 2], counts[:, 3]

    with np.errstate(divide="ignore", invalid="ignore"):
        s_precision = np.where((s_tp + s_fp) > 0, s_tp / np.maximum(s_tp + s_fp, 1), 0.0)
        s_recall = np.where((s_tp + s_fn) > 0, s_tp / np.maximum(s_tp + s_fn, 1), 0.0)
        f1_denom = s_precision + s_recall
        f1_numerator = 2 * s_precision * s_recall
        f1_samples = np.where(f1_denom > 0, f1_numerator / np.maximum(f1_denom, 1e-12), 0.0)
        mcc_numerator = s_tp * s_tn - s_fp * s_fn
        mcc_denom_sq = (s_tp + s_fp) * (s_tp + s_fn) * (s_tn + s_fp) * (s_tn + s_fn)
        mcc_samples = np.where(
            mcc_denom_sq > 0, mcc_numerator / np.sqrt(np.maximum(mcc_denom_sq, 1)), 0.0
        )

    f1_point = _f1(precision_point, recall_point)
    mcc_point = _mcc(tp, fp, tn, fn)
    f1_low, f1_high = np.percentile(f1_samples, [2.5, 97.5])
    mcc_low, mcc_high = np.percentile(mcc_samples, [2.5, 97.5])

    return ClassifierMetrics(
        n=n,
        true_positive=tp,
        false_positive=fp,
        true_negative=tn,
        false_negative=fn,
        precision=ProportionCi(precision_point, precision_ci.ci_low, precision_ci.ci_high),
        recall=ProportionCi(recall_point, recall_ci.ci_low, recall_ci.ci_high),
        f1=ProportionCi(f1_point, float(f1_low), float(f1_high)),
        mcc=ProportionCi(mcc_point, float(mcc_low), float(mcc_high)),
    )
