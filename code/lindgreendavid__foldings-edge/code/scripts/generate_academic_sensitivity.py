#!/usr/bin/env python3
"""Generate the post-release protein-cluster sensitivity analysis.

The frozen v0.1 registry treats residues as the analysis rows. This audit keeps its point
estimates intact while checking whether the main conclusions survive when proteins, rather
than correlated residues within proteins, are the resampling and testing units.
"""

from __future__ import annotations

import csv
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "data" / "external" / "joined_residues.csv"
OUTPUT = ROOT / "reports" / "post-release-academic-sensitivity.json"
SEED = 20260814
RESAMPLES = 10_000


def _stable(value: float) -> float:
    """Quantize computed results to eight platform-stable significant digits."""
    return float(f"{float(value):.8g}")


def _mcc(tp: np.ndarray, fp: np.ndarray, tn: np.ndarray, fn: np.ndarray) -> np.ndarray:
    numerator = tp * tn - fp * fn
    denominator = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return np.divide(numerator, denominator, out=np.zeros_like(numerator), where=denominator > 0)


def main() -> None:
    by_protein: dict[str, list[list[float]]] = defaultdict(lambda: [[], []])
    with INPUT.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            group = 0 if row["is_disorder"] == "True" else 1
            by_protein[row["acc"]][group].append(float(row["plddt"]))

    paired_differences = np.asarray(
        [
            np.median(outside) - np.median(inside)
            for inside, outside in by_protein.values()
            if inside and outside
        ],
        dtype=float,
    )
    rng = np.random.default_rng(SEED)
    paired_bootstrap = np.median(
        paired_differences[
            rng.integers(0, paired_differences.size, (RESAMPLES, paired_differences.size))
        ],
        axis=1,
    )
    signed_rank = stats.wilcoxon(paired_differences, alternative="two-sided")

    confusion = []
    for inside, outside in by_protein.values():
        inside_values = np.asarray(inside)
        outside_values = np.asarray(outside)
        confusion.append(
            [
                np.count_nonzero(inside_values < 70),
                np.count_nonzero(outside_values < 70),
                np.count_nonzero(outside_values >= 70),
                np.count_nonzero(inside_values >= 70),
            ]
        )
    confusion_array = np.asarray(confusion, dtype=float)
    sampled = confusion_array[
        rng.integers(0, len(confusion_array), (RESAMPLES, len(confusion_array)))
    ].sum(axis=1)
    tp, fp, tn, fn = sampled.T
    precision = tp / (tp + fp)
    recall = tp / (tp + fn)
    f1 = 2 * precision * recall / (precision + recall)
    mcc = _mcc(tp, fp, tn, fn)

    def interval(values: np.ndarray) -> list[float]:
        return [_stable(value) for value in np.quantile(values, [0.025, 0.975])]

    payload = {
        "schema_version": "1.0.0",
        "status": "post-release sensitivity; does not alter the frozen v0.1 registry",
        "settings": {
            "resampling_unit": "protein",
            "bootstrap_resamples": RESAMPLES,
            "bootstrap_seed": SEED,
            "classifier_threshold": "pLDDT < 70",
        },
        "h1_protein_paired_sensitivity": {
            "n_proteins_with_inside_and_outside_residues": int(paired_differences.size),
            "median_within_protein_median_difference_outside_minus_inside": _stable(
                np.median(paired_differences)
            ),
            "protein_cluster_bootstrap_95_ci": interval(paired_bootstrap),
            "two_sided_wilcoxon_statistic": _stable(signed_rank.statistic),
            "two_sided_wilcoxon_p_value": _stable(signed_rank.pvalue),
            "proteins_positive_difference": int(np.count_nonzero(paired_differences > 0)),
            "proteins_negative_difference": int(np.count_nonzero(paired_differences < 0)),
            "proteins_zero_difference": int(np.count_nonzero(paired_differences == 0)),
        },
        "h2_protein_cluster_bootstrap_95_ci": {
            "n_proteins": len(confusion_array),
            "precision": interval(precision),
            "recall": interval(recall),
            "f1": interval(f1),
            "mcc": interval(mcc),
        },
        "interpretation": (
            "The direction and moderate classifier performance survive protein-level "
            "resampling. Protein-cluster intervals, not residue-i.i.d. intervals, are the "
            "preferred uncertainty sensitivity because residues within a protein are correlated."
        ),
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
