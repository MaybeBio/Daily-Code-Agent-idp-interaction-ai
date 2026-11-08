"""Build the frozen v0.1.0 result registry from the joined per-residue CSV.

Implements exactly the preregistered plan in `docs/research-protocol.md`:
load the joined residues, run the primary H1 distributional comparison and
the H2 fixed-threshold classifier evaluation (overall and broken down by
conditional-folding flag and by evidence code), and package a downsampled,
reproducible visualization subset for the site.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np

from foldings_edge.dataset import Residue, evidence_code_counts, load_residues, split_by_disorder
from foldings_edge.stats import (
    ALPHA,
    DEFAULT_BOOTSTRAP_RESAMPLES,
    DEFAULT_SEED,
    PLDDT_CONFIDENT_THRESHOLD,
    bootstrap_median_diff,
    classifier_metrics,
    ks_two_sample,
    mann_whitney_two_sample,
)

SCHEMA_VERSION = 1
MIN_EVIDENCE_CODE_RESIDUES = 200
DISTRIBUTION_SAMPLE_SIZE = 3000
MIN_PROTEIN_RESIDUES_FOR_BREAKDOWN = 20


def _proportion_ci_dict(pc: Any) -> dict[str, float]:
    return {"point": pc.point, "ci_95_low": pc.ci_low, "ci_95_high": pc.ci_high}


def _classifier_dict(residues: list[Residue]) -> dict[str, Any]:
    is_disorder = [r.is_disorder for r in residues]
    plddt = [r.plddt for r in residues]
    metrics = classifier_metrics(is_disorder, plddt)
    return {
        "n": metrics.n,
        "threshold": PLDDT_CONFIDENT_THRESHOLD,
        "confusion": {
            "true_positive": metrics.true_positive,
            "false_positive": metrics.false_positive,
            "true_negative": metrics.true_negative,
            "false_negative": metrics.false_negative,
        },
        "precision": _proportion_ci_dict(metrics.precision),
        "recall": _proportion_ci_dict(metrics.recall),
        "f1": _proportion_ci_dict(metrics.f1),
        "mcc": _proportion_ci_dict(metrics.mcc),
    }


def _protein_breakdown(residues: list[Residue]) -> list[dict[str, Any]]:
    """Per-protein false-positive/false-negative counts, for the report's failure analysis.

    Field names follow the same convention as `foldings_edge.stats.classifier_metrics`:
    positive class = "residue is inside a DisProt disorder region", predicted positive =
    pLDDT < threshold. So a false negative here is a curated-disorder residue that
    nonetheless received a confident pLDDT (the "conditionally folded"-like failure
    mode), and a false positive is a non-disorder residue that received a low pLDDT.
    """
    by_protein: dict[str, list[Residue]] = {}
    for r in residues:
        by_protein.setdefault(r.acc, []).append(r)

    rows: list[dict[str, Any]] = []
    for acc, prows in by_protein.items():
        n = len(prows)
        if n < MIN_PROTEIN_RESIDUES_FOR_BREAKDOWN:
            continue
        disorder_n = sum(1 for r in prows if r.is_disorder)
        if disorder_n == 0:
            continue
        false_negative = sum(
            1 for r in prows if r.is_disorder and r.plddt >= PLDDT_CONFIDENT_THRESHOLD
        )
        false_positive = sum(
            1 for r in prows if not r.is_disorder and r.plddt < PLDDT_CONFIDENT_THRESHOLD
        )
        rows.append(
            {
                "acc": acc,
                "disprot_id": prows[0].disprot_id,
                "n_residues": n,
                "n_disorder_residues": disorder_n,
                "false_negative_residues": false_negative,
                "false_negative_rate_of_disorder": (
                    false_negative / disorder_n if disorder_n else 0.0
                ),
                "false_positive_residues": false_positive,
                "is_conditional": prows[0].is_conditional,
            }
        )
    rows.sort(key=lambda r: float(r["false_negative_rate_of_disorder"]), reverse=True)
    return rows


def _downsample(values: list[float], size: int, seed: int) -> list[float]:
    if len(values) <= size:
        return values
    rng = np.random.default_rng(seed)
    idx = rng.choice(len(values), size=size, replace=False)
    return [values[i] for i in sorted(idx.tolist())]


def build_registry(residues_path: Path) -> dict[str, Any]:
    residues = load_residues(residues_path)
    inside, outside = split_by_disorder(residues)

    inside_plddt = [r.plddt for r in inside]
    outside_plddt = [r.plddt for r in outside]

    mann_whitney = mann_whitney_two_sample(outside_plddt, inside_plddt)
    ks = ks_two_sample(outside_plddt, inside_plddt)
    boot = bootstrap_median_diff(inside_plddt, outside_plddt)

    conditional = [r for r in inside if r.is_conditional]
    non_conditional = [r for r in inside if not r.is_conditional]

    ec_counts = evidence_code_counts(residues)
    ec_breakdown = {}
    for ec_id, count in sorted(ec_counts.items()):
        if count < MIN_EVIDENCE_CODE_RESIDUES:
            continue
        ec_rows = [r for r in residues if (r.is_disorder and r.ec_id == ec_id) or not r.is_disorder]
        ec_name = next((r.ec_name for r in residues if r.ec_id == ec_id), "")
        ec_breakdown[ec_id] = {
            "ec_name": ec_name,
            "n_disorder_residues": count,
            "classifier": _classifier_dict(ec_rows),
        }

    proteins = {r.acc for r in residues}

    registry: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "dataset": {
            "total_proteins": len(proteins),
            "total_residues": len(residues),
            "disorder_annotated_residues": len(inside),
            "non_disorder_residues": len(outside),
            "conditional_disorder_residues": len(conditional),
            "non_conditional_disorder_residues": len(non_conditional),
        },
        "settings": {
            "alpha": ALPHA,
            "bootstrap_resamples": DEFAULT_BOOTSTRAP_RESAMPLES,
            "bootstrap_seed": DEFAULT_SEED,
            "plddt_confident_threshold": PLDDT_CONFIDENT_THRESHOLD,
        },
        "h1_distribution": {
            "mann_whitney_u": {
                "statistic": mann_whitney.statistic,
                "p_value": mann_whitney.p_value,
                "significant": mann_whitney.significant,
                "alpha": ALPHA,
            },
            "ks_robustness": {
                "statistic": ks.statistic,
                "p_value": ks.p_value,
                "significant": ks.significant,
                "alpha": ALPHA,
            },
            "bootstrap_median_diff_outside_minus_inside": {
                "median_inside_disorder": boot.median_a,
                "median_outside_disorder": boot.median_b,
                "median_diff": boot.median_diff,
                "ci_95_low": boot.ci_low,
                "ci_95_high": boot.ci_high,
                "resamples": boot.resamples,
                "seed": boot.seed,
            },
        },
        "h2_classifier": {
            "overall": _classifier_dict(residues),
            "by_conditional_flag": {
                "conditional_disorder_regions": _classifier_dict(conditional + outside),
                "non_conditional_disorder_regions": _classifier_dict(non_conditional + outside),
            },
            "by_evidence_code": ec_breakdown,
        },
        "protein_breakdown": _protein_breakdown(residues),
        "distributions": {
            "note": (
                f"Downsampled to at most {DISTRIBUTION_SAMPLE_SIZE} values per group "
                f"(seed={DEFAULT_SEED}) for visualization only; all statistics above are "
                "computed on the full joined dataset, not this subsample."
            ),
            "inside_disorder": _downsample(inside_plddt, DISTRIBUTION_SAMPLE_SIZE, DEFAULT_SEED),
            "outside_disorder": _downsample(
                outside_plddt, DISTRIBUTION_SAMPLE_SIZE, DEFAULT_SEED + 1
            ),
        },
    }
    return registry
