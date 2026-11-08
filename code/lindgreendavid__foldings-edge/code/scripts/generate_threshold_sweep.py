#!/usr/bin/env python3
"""Precompute classifier performance across a grid of pLDDT thresholds.

Usage: python scripts/generate_threshold_sweep.py \
    --output site/app/data/threshold-sweep.json \
    [--residues data/external/joined_residues.csv]

This is a presentation-layer helper for the site's interactive threshold
explorer, not part of the frozen v0.1.0 analysis: it does not modify
`src/foldings_edge/`, `data/external/joined_residues.csv`, or
`reports/v0.1-foldings-edge-registry.json`. It reads the same, real joined
per-residue CSV the frozen registry is built from and recomputes
precision/recall/F1/MCC/confusion-counts at an integer grid of thresholds
(0-100 inclusive, i.e. every whole pLDDT point the site's slider can land
on), overall and split by the same protein-level conditional-folding flag
`docs/research-report.md` already reports on. It deliberately reports point
estimates only (no bootstrap CIs) at every grid threshold — the existing,
CI-bearing frozen numbers at threshold=70 remain the ones displayed as
"the reported result" elsewhere on the site.

The classifier convention matches `foldings_edge.stats.classifier_metrics`
exactly: positive class = residue is inside a DisProt-curated disorder
region; predicted positive = pLDDT < threshold (strict "<", not "<=").

Sorted-array + binary-search implementation (not a per-threshold full scan)
because a naive O(thresholds x residues) scan over 228,662 residues x 101
thresholds is ~23M redundant comparisons; sorting once and using bisect
against the fixed threshold grid is O(n log n) total.
"""

from __future__ import annotations

import argparse
import bisect
import csv
import json
import math
from pathlib import Path
from typing import Any

THRESHOLD_MIN = 0
THRESHOLD_MAX = 100
PREREGISTERED_THRESHOLD = 70


def _mcc(tp: int, fp: int, tn: int, fn: int) -> float:
    numerator = tp * tn - fp * fn
    denom_sq = (tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)
    if denom_sq == 0:
        return 0.0
    return numerator / math.sqrt(denom_sq)


def _f1(precision: float, recall: float) -> float:
    if precision + recall == 0:
        return 0.0
    return 2 * precision * recall / (precision + recall)


def _metrics_at(
    threshold: int,
    disorder_sorted: list[float],
    non_disorder_sorted: list[float],
) -> dict[str, Any]:
    """Confusion counts + precision/recall/F1/MCC at one threshold.

    Predicted positive iff pLDDT < threshold, matching
    `foldings_edge.stats.classifier_metrics`. Because both input lists are
    sorted ascending, `bisect_left(list, threshold)` gives exactly the count
    of values strictly less than `threshold`.
    """
    tp = bisect.bisect_left(disorder_sorted, threshold)
    fn = len(disorder_sorted) - tp
    fp = bisect.bisect_left(non_disorder_sorted, threshold)
    tn = len(non_disorder_sorted) - fp

    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0

    return {
        "threshold": threshold,
        "confusion": {
            "true_positive": tp,
            "false_positive": fp,
            "true_negative": tn,
            "false_negative": fn,
        },
        "precision": precision,
        "recall": recall,
        "f1": _f1(precision, recall),
        "mcc": _mcc(tp, fp, tn, fn),
    }


def _sweep(disorder_sorted: list[float], non_disorder_sorted: list[float]) -> list[dict[str, Any]]:
    return [
        _metrics_at(t, disorder_sorted, non_disorder_sorted)
        for t in range(THRESHOLD_MIN, THRESHOLD_MAX + 1)
    ]


def build_sweep(residues_path: Path) -> dict[str, Any]:
    disorder_plddt: list[float] = []
    non_disorder_plddt: list[float] = []
    conditional_disorder_plddt: list[float] = []
    non_conditional_disorder_plddt: list[float] = []

    with residues_path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            plddt = float(row["plddt"])
            is_disorder = row["is_disorder"].strip().lower() == "true"
            is_conditional = row["is_conditional"].strip().lower() == "true"
            if is_disorder:
                disorder_plddt.append(plddt)
                if is_conditional:
                    conditional_disorder_plddt.append(plddt)
                else:
                    non_conditional_disorder_plddt.append(plddt)
            else:
                non_disorder_plddt.append(plddt)

    disorder_plddt.sort()
    non_disorder_plddt.sort()
    conditional_disorder_plddt.sort()
    non_conditional_disorder_plddt.sort()

    return {
        "schema_version": 1,
        "threshold_min": THRESHOLD_MIN,
        "threshold_max": THRESHOLD_MAX,
        "threshold_step": 1,
        "preregistered_threshold": PREREGISTERED_THRESHOLD,
        "overall": _sweep(disorder_plddt, non_disorder_plddt),
        "by_conditional_flag": {
            "conditional_disorder_regions": _sweep(conditional_disorder_plddt, non_disorder_plddt),
            "non_conditional_disorder_regions": _sweep(
                non_conditional_disorder_plddt, non_disorder_plddt
            ),
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    parser.add_argument(
        "--residues",
        default="data/external/joined_residues.csv",
        help="Path to the joined per-residue CSV produced by scripts/fetch_data.py",
    )
    args = parser.parse_args()

    sweep = build_sweep(Path(args.residues))
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(sweep, indent=2, sort_keys=True) + "\n")
    print(f"wrote {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
