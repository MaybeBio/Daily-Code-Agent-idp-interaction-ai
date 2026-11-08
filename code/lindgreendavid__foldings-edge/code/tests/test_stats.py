from __future__ import annotations

from foldings_edge.stats import (
    bootstrap_median_diff,
    classifier_metrics,
    ks_two_sample,
    mann_whitney_two_sample,
    wilson_ci,
)


def test_mann_whitney_identical_distributions_not_significant() -> None:
    a = [1.0, 2.0, 3.0, 4.0, 5.0] * 10
    b = [1.0, 2.0, 3.0, 4.0, 5.0] * 10
    result = mann_whitney_two_sample(a, b)
    assert result.p_value == 1.0
    assert result.significant is False


def test_mann_whitney_clearly_different_distributions_significant() -> None:
    a = list(range(100))
    b = [x + 1000 for x in range(100)]
    result = mann_whitney_two_sample(a, b)
    assert result.significant is True


def test_ks_two_sample_clearly_different_significant() -> None:
    a = list(range(100))
    b = [x + 1000 for x in range(100)]
    result = ks_two_sample(a, b)
    assert result.statistic == 1.0
    assert result.significant is True


def test_bootstrap_median_diff_deterministic_given_seed() -> None:
    a = [1.0, 2.0, 3.0, 4.0, 5.0]
    b = [10.0, 20.0, 30.0, 40.0, 50.0]
    first = bootstrap_median_diff(a, b, resamples=200, seed=42)
    second = bootstrap_median_diff(a, b, resamples=200, seed=42)
    assert first == second
    assert first.median_a == 3.0
    assert first.median_b == 30.0
    assert first.median_diff == 27.0
    assert first.ci_low <= first.median_diff <= first.ci_high


def test_bootstrap_median_diff_different_seed_can_differ() -> None:
    a = [1.0, 2.0, 3.0, 4.0, 5.0, 100.0]
    b = [10.0, 20.0, 30.0, 40.0, 50.0, 5.0]
    first = bootstrap_median_diff(a, b, resamples=200, seed=1)
    second = bootstrap_median_diff(a, b, resamples=200, seed=2)
    assert first.seed != second.seed


def test_wilson_ci_basic_properties() -> None:
    result = wilson_ci(50, 100)
    assert result.point == 0.5
    assert result.ci_low < 0.5 < result.ci_high


def test_wilson_ci_zero_n() -> None:
    result = wilson_ci(0, 0)
    assert result.point != result.point  # NaN


def test_classifier_metrics_perfect_classifier() -> None:
    # pLDDT < 70 exactly predicts is_disorder=True everywhere.
    is_disorder = [True, True, False, False]
    plddt = [10.0, 20.0, 90.0, 95.0]
    metrics = classifier_metrics(is_disorder, plddt, resamples=100, seed=1)
    assert metrics.precision.point == 1.0
    assert metrics.recall.point == 1.0
    assert metrics.f1.point == 1.0
    assert metrics.mcc.point == 1.0
    assert metrics.true_positive == 2
    assert metrics.false_positive == 0
    assert metrics.false_negative == 0
    assert metrics.true_negative == 2


def test_classifier_metrics_worst_case_mcc_zero_when_denominator_zero() -> None:
    # All predicted negative (never below threshold): tp=fp=0.
    is_disorder = [True, False]
    plddt = [90.0, 95.0]
    metrics = classifier_metrics(is_disorder, plddt, resamples=50, seed=1)
    assert metrics.true_positive == 0
    assert metrics.false_positive == 0
    assert metrics.precision.point == 0.0
    assert metrics.mcc.point == 0.0


def test_classifier_metrics_mixed_case_reasonable_bounds() -> None:
    is_disorder = [True, True, True, False, False, False, False, False]
    plddt = [40.0, 85.0, 30.0, 90.0, 20.0, 88.0, 91.0, 95.0]
    metrics = classifier_metrics(is_disorder, plddt, resamples=300, seed=7)
    assert 0.0 <= metrics.precision.point <= 1.0
    assert 0.0 <= metrics.recall.point <= 1.0
    assert -1.0 <= metrics.mcc.point <= 1.0
    assert metrics.f1.ci_low <= metrics.f1.point <= metrics.f1.ci_high
