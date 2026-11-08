"""Pool-independent subgroup audit of the pLDDT<70 classifier (POST-HOC, not preregistered).

The frozen registry reports precision/F1/MCC for subgroups of disorder residues (by
conditional-folding flag, by evidence code), each time scored against the *same* shared pool of
non-disorder residues. Precision then depends mechanically on subgroup size: with a fixed
false-positive pool, a small subgroup always has lower precision even if the classifier treats
it identically. This module therefore compares subgroups with quantities that do not depend on
subgroup size:

* recall (sensitivity) within the subgroup,
* the area under the ROC curve of (100 - pLDDT) separating the subgroup's disorder residues
  from the shared non-disorder pool,
* both with 95% intervals from a protein-cluster bootstrap (residues within a protein are not
  independent, and each protein belongs to exactly one conditional-folding group),

and reports the counterfactual precision each subgroup would have *if its recall equalled the
overall recall*, which isolates the size artifact.
"""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path
from typing import Any

import numpy as np
from numpy.typing import NDArray
from scipy.stats import rankdata

from foldings_edge.dataset import Residue, load_residues
from foldings_edge.stats import PLDDT_CONFIDENT_THRESHOLD

AUDIT_SEED = 20261009
RECALL_RESAMPLES = 10_000
AUC_RESAMPLES = 1_000
MIN_EVIDENCE_CODE_RESIDUES = 200


def _by_protein(residues: list[Residue]) -> dict[str, list[Residue]]:
    groups: dict[str, list[Residue]] = defaultdict(list)
    for residue in residues:
        groups[residue.acc].append(residue)
    return dict(sorted(groups.items()))


def _auc(disorder_plddt: NDArray[np.float64], background_plddt: NDArray[np.float64]) -> float:
    """P(a disorder residue has lower pLDDT than a background residue), ties counted 1/2."""
    n1, n2 = len(disorder_plddt), len(background_plddt)
    if n1 == 0 or n2 == 0:
        return float("nan")
    ranks = rankdata(np.concatenate([-disorder_plddt, -background_plddt]))
    return float((ranks[:n1].sum() - n1 * (n1 + 1) / 2) / (n1 * n2))


def _quantiles(values: NDArray[np.float64]) -> tuple[float, float]:
    finite = values[np.isfinite(values)]
    lo, hi = np.quantile(finite, [0.025, 0.975])
    return float(lo), float(hi)


def precision_size_artifact(residues: list[Residue]) -> dict[str, Any]:
    """Observed vs counterfactual (equal-recall) precision for the conditional-folding groups."""
    disorder = [r for r in residues if r.is_disorder]
    background = [r for r in residues if not r.is_disorder]
    fp_pool = sum(1 for r in background if r.plddt < PLDDT_CONFIDENT_THRESHOLD)
    overall_recall = sum(1 for r in disorder if r.plddt < PLDDT_CONFIDENT_THRESHOLD) / len(disorder)
    out: dict[str, Any] = {"shared_false_positives": fp_pool, "overall_recall": overall_recall}
    for name, flag in (("conditional", True), ("non_conditional", False)):
        group = [r for r in disorder if r.is_conditional is flag]
        tp = sum(1 for r in group if r.plddt < PLDDT_CONFIDENT_THRESHOLD)
        expected_tp = overall_recall * len(group)
        out[name] = {
            "n_disorder_residues": len(group),
            "observed_precision": tp / (tp + fp_pool),
            "equal_recall_precision": expected_tp / (expected_tp + fp_pool),
        }
    return out


def conditional_flag_audit(
    residues: list[Residue],
    resamples_recall: int = RECALL_RESAMPLES,
    resamples_auc: int = AUC_RESAMPLES,
    seed: int = AUDIT_SEED,
) -> dict[str, Any]:
    proteins = _by_protein(residues)
    names = list(proteins)
    cond_flag = np.array([any(r.is_conditional for r in proteins[n]) for n in names])
    n_dis = {True: np.zeros(len(names)), False: np.zeros(len(names))}
    n_tp = {True: np.zeros(len(names)), False: np.zeros(len(names))}
    dis_plddt: dict[bool, list[NDArray[np.float64]]] = {True: [], False: []}
    bg_plddt: list[NDArray[np.float64]] = []
    for i, name in enumerate(names):
        rows = proteins[name]
        flag = bool(cond_flag[i])
        d = np.array([r.plddt for r in rows if r.is_disorder])
        b = np.array([r.plddt for r in rows if not r.is_disorder])
        n_dis[flag][i] = len(d)
        n_tp[flag][i] = float(np.sum(d < PLDDT_CONFIDENT_THRESHOLD))
        dis_plddt[flag].append(d)
        dis_plddt[not flag].append(np.array([]))
        bg_plddt.append(b)

    def recall(flag: bool, idx: NDArray[np.intp]) -> float:
        total = n_dis[flag][idx].sum()
        return float(n_tp[flag][idx].sum() / total) if total else float("nan")

    def auc(flag: bool, idx: NDArray[np.intp]) -> float:
        d = np.concatenate([dis_plddt[flag][i] for i in idx])
        b = np.concatenate([bg_plddt[i] for i in idx])
        return _auc(d, b)

    everyone = np.arange(len(names))
    point = {
        "recall_conditional": recall(True, everyone),
        "recall_non_conditional": recall(False, everyone),
        "auc_conditional": auc(True, everyone),
        "auc_non_conditional": auc(False, everyone),
    }

    rng = np.random.default_rng(seed)
    recall_diff = np.empty(resamples_recall)
    for k in range(resamples_recall):
        idx = rng.integers(0, len(names), len(names))
        recall_diff[k] = recall(False, idx) - recall(True, idx)
    auc_diff = np.empty(resamples_auc)
    for k in range(resamples_auc):
        idx = rng.integers(0, len(names), len(names))
        auc_diff[k] = auc(False, idx) - auc(True, idx)

    rd_lo, rd_hi = _quantiles(recall_diff)
    ad_lo, ad_hi = _quantiles(auc_diff)
    return {
        "n_proteins": len(names),
        "n_conditional_proteins": int(cond_flag.sum()),
        "point": point,
        "recall_difference_non_minus_conditional": {
            "point": point["recall_non_conditional"] - point["recall_conditional"],
            "ci_95": [rd_lo, rd_hi],
            "resamples": resamples_recall,
        },
        "auc_difference_non_minus_conditional": {
            "point": point["auc_non_conditional"] - point["auc_conditional"],
            "ci_95": [ad_lo, ad_hi],
            "resamples": resamples_auc,
        },
    }


def evidence_code_audit(
    residues: list[Residue], resamples: int = RECALL_RESAMPLES, seed: int = AUDIT_SEED
) -> dict[str, Any]:
    """Recall by evidence code with protein-cluster 95% intervals and protein counts."""
    disorder = [r for r in residues if r.is_disorder]
    codes: dict[str, list[Residue]] = defaultdict(list)
    for r in disorder:
        codes[r.ec_id].append(r)
    rng = np.random.default_rng(seed)
    out: dict[str, Any] = {}
    for ec_id, rows in sorted(codes.items()):
        if len(rows) < MIN_EVIDENCE_CODE_RESIDUES:
            continue
        per_protein = _by_protein(rows)
        n = np.array([len(v) for v in per_protein.values()], dtype=float)
        tp = np.array(
            [
                sum(1 for r in v if r.plddt < PLDDT_CONFIDENT_THRESHOLD)
                for v in per_protein.values()
            ],
            dtype=float,
        )
        boot = np.empty(resamples)
        for k in range(resamples):
            idx = rng.integers(0, len(n), len(n))
            boot[k] = tp[idx].sum() / n[idx].sum()
        lo, hi = _quantiles(boot)
        out[ec_id] = {
            "ec_name": rows[0].ec_name,
            "n_disorder_residues": len(rows),
            "n_proteins": len(per_protein),
            "recall": float(tp.sum() / n.sum()),
            "recall_ci_95": [lo, hi],
        }
    return out


def build_subgroup_audit(residues_path: Path) -> dict[str, Any]:
    residues = load_residues(residues_path)
    return {
        "schema_version": 1,
        "label": "POST-HOC (not preregistered): pool-independent subgroup audit",
        "seed": AUDIT_SEED,
        "precision_size_artifact": precision_size_artifact(residues),
        "conditional_flag": conditional_flag_audit(residues),
        "evidence_codes": evidence_code_audit(residues),
    }
