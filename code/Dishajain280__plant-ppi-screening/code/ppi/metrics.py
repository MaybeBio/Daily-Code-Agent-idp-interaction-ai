"""Ranking metrics for screening: AUC, recall in the top slice, precision at a cutoff.

With a few dozen pairs, a single AUC hides how unstable it is, so every metric
comes with a bootstrap interval. Ties are handled explicitly: a score tied with
a positive counts as half a win.
"""

from __future__ import annotations

import math
import random
from dataclasses import dataclass


@dataclass
class Metric:
    value: float
    low: float | None = None
    high: float | None = None
    kind: str = "rate"  # "rate" keeps three decimals; "count" prints an integer

    def __str__(self) -> str:
        if math.isnan(self.value):
            return "n/a"
        if self.kind == "count":
            return f"{int(round(self.value))}"
        has_interval = (self.low is not None and self.high is not None
                        and not math.isnan(self.low) and not math.isnan(self.high))
        if not has_interval:
            return f"{self.value:.3f}"
        return f"{self.value:.3f} (95% CI {self.low:.3f}-{self.high:.3f})"


def auc(scores: list[float], labels: list[int]) -> float:
    """Area under the ROC curve, computed by pairwise comparison."""
    positives = [s for s, y in zip(scores, labels) if y == 1]
    negatives = [s for s, y in zip(scores, labels) if y == 0]
    if not positives or not negatives:
        return float("nan")
    wins = 0.0
    for p in positives:
        for n in negatives:
            if p > n:
                wins += 1.0
            elif p == n:
                wins += 0.5
    return wins / (len(positives) * len(negatives))


def recall_in_top_fraction(scores: list[float], labels: list[int],
                           fraction: float = 0.1) -> tuple[float, int]:
    """Fraction of positives that land in the top ``fraction`` of ranked pairs."""
    if not scores:
        return float("nan"), 0
    total_positive = sum(labels)
    if total_positive == 0:
        return float("nan"), 0
    k = max(1, math.ceil(len(scores) * fraction))
    order = sorted(range(len(scores)), key=lambda i: (-scores[i], i))
    hits = sum(labels[i] for i in order[:k])
    return hits / total_positive, k


def precision_at(scores: list[float], labels: list[int], threshold: float) -> tuple[float, int]:
    """Precision among pairs scoring at or above ``threshold``."""
    selected = [y for s, y in zip(scores, labels) if s >= threshold]
    if not selected:
        return float("nan"), 0
    return sum(selected) / len(selected), len(selected)


def bootstrap_auc(scores: list[float], labels: list[int], n: int = 2000,
                  seed: int = 0) -> tuple[float, float]:
    """Percentile bootstrap interval for AUC."""
    if len(scores) < 4:
        return float("nan"), float("nan")
    rng = random.Random(seed)
    estimates = []
    for _ in range(n):
        idx = [rng.randrange(len(scores)) for _ in scores]
        value = auc([scores[i] for i in idx], [labels[i] for i in idx])
        if not math.isnan(value):
            estimates.append(value)
    if not estimates:
        return float("nan"), float("nan")
    estimates.sort()
    return estimates[int(0.025 * len(estimates))], estimates[int(0.975 * len(estimates)) - 1]


def summarise(scores: list[float], labels: list[int], fraction: float = 0.1,
              threshold: float | None = None) -> dict[str, Metric]:
    """All headline metrics in one call."""
    point = auc(scores, labels)
    low, high = bootstrap_auc(scores, labels)
    recall, k = recall_in_top_fraction(scores, labels, fraction)
    out = {
        "auc": Metric(point, low, high),
        "recall_top": Metric(recall),
        "top_k": Metric(float(k), kind="count"),
        "n_pairs": Metric(float(len(scores)), kind="count"),
        "n_positive": Metric(float(sum(labels)), kind="count"),
    }
    if threshold is not None:
        precision, selected = precision_at(scores, labels, threshold)
        out["precision_at_threshold"] = Metric(precision)
        out["selected_at_threshold"] = Metric(float(selected), kind="count")
    return out
