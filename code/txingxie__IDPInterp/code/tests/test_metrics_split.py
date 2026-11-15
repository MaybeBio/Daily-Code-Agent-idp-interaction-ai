"""
Unit tests for the leakage-free, split-aware evaluation path in
``RegionAwareEvaluator`` (idpinterp/eval/metrics.py).

These tests verify that the per-feature activation threshold is selected on the
validation residues/regions only, and that the reported precision / region
recall / F1 / average_precision are computed on the held-out test split.

Requires torch. Run with:
    pytest tests/test_metrics_split.py -v
"""
import numpy as np
import pytest

torch = pytest.importorskip("torch")

from idpinterp.eval.metrics import RegionAwareEvaluator


def _hand_scenario():
    """
    8 residues, feature 0 is informative, feature 1 is dead (all zeros).
    Each positive residue is its own singleton region.

        idx :  0    1    2    3    4    5    6    7
        act : 0.9  0.8  0.7  0.6  0.5  0.4  0.3  0.2
        lab :  1    1    0    1    1    0    0    0
        spl :  v    t    v    t    v    t    v    t

    Hand-computed (see scratchpad mirror): the best validation F1 (=0.8)
    occurs at threshold 0.5; applied to the test split this yields
    precision=1.0, region_recall=1.0, F1=1.0, average_precision=1.0.
    """
    act = np.array([0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2], dtype=np.float32)
    lab = np.array([1, 1, 0, 1, 1, 0, 0, 0], dtype=np.float32)
    val_res = np.array([True, False, True, False, True, False, True, False])

    all_acts = np.zeros((8, 2), dtype=np.float16)
    all_acts[:, 0] = act  # feature 1 stays all-zero (dead)

    pos = np.where(lab > 0.5)[0]                       # [0, 1, 3, 4]
    regions = [(int(p), int(p) + 1) for p in pos]      # singleton regions
    region_val_mask = val_res[pos]                     # [T, F, F, T]
    return all_acts, lab, regions, val_res, region_val_mask


def test_split_thresholds_on_val_reports_on_test():
    all_acts, labels, regions, val_res, region_val_mask = _hand_scenario()
    ev = RegionAwareEvaluator(chunk_size=256)
    results = ev.evaluate_all_features(
        all_acts,
        labels,
        regions,
        triplet="unit-test",
        region_maxes_matrix=None,           # exercise internal region-max scan
        device="cpu",
        val_residue_mask=val_res,
        val_region_mask=region_val_mask,
    )

    by_feat = {r["feature_idx"]: r for r in results}

    # Feature 1 is dead -> must be absent.
    assert 1 not in by_feat, "dead feature should be skipped"

    # Feature 0 matches the hand-computed, leakage-free expectation.
    assert 0 in by_feat
    f0 = by_feat[0]
    assert f0["best_threshold"] == pytest.approx(0.5, abs=1e-3)
    assert f0["precision"] == pytest.approx(1.0, abs=1e-3)
    assert f0["region_recall"] == pytest.approx(1.0, abs=1e-3)
    assert f0["f1"] == pytest.approx(1.0, abs=1e-3)
    assert f0["average_precision"] == pytest.approx(1.0, abs=1e-3)


def test_split_residue_sets_are_disjoint():
    """The val mask and its complement must partition all residues exactly once."""
    _, _, _, val_res, _ = _hand_scenario()
    test_res = ~val_res
    assert (val_res & test_res).sum() == 0
    assert (val_res | test_res).all()


def test_legacy_path_unaffected_when_no_split():
    """Without split masks, the legacy full-corpus path still returns results."""
    all_acts, labels, regions, _, _ = _hand_scenario()
    ev = RegionAwareEvaluator(chunk_size=256)
    results = ev.evaluate_all_features(
        all_acts, labels, regions, triplet="unit-test", device="cpu"
    )
    feats = {r["feature_idx"] for r in results}
    assert 0 in feats and 1 not in feats
