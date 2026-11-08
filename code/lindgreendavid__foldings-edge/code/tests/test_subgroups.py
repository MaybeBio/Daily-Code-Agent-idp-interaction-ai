import numpy as np
import pytest

from foldings_edge.dataset import Residue
from foldings_edge.subgroups import (
    _auc,
    build_subgroup_audit,
    conditional_flag_audit,
    evidence_code_audit,
    precision_size_artifact,
)


def res(acc, plddt, disorder, cond=False, ec="ECO:1", n=0):
    return Residue(
        acc, "DP", n, plddt, disorder, cond, ec if disorder else "", "name" if disorder else ""
    )


def make(n_cond=4, n_non=8, per=60):
    """Each protein: `per` disorder residues (pLDDT ~50) and `per` background residues (~90)."""
    rows = []
    for i in range(n_cond + n_non):
        cond = i < n_cond
        acc = f"P{i}"
        for k in range(per):
            rows.append(res(acc, 50.0 + (k % 3), True, cond, n=k))
        for k in range(per):
            rows.append(res(acc, 90.0 - (k % 3) if k >= 5 else 40.0, False, cond, n=1000 + k))
    return rows


def test_auc_extremes_and_ties():
    assert _auc(np.array([10.0, 20.0]), np.array([80.0, 90.0])) == 1.0
    assert _auc(np.array([80.0, 90.0]), np.array([10.0, 20.0])) == 0.0
    assert _auc(np.array([50.0]), np.array([50.0])) == 0.5
    assert np.isnan(_auc(np.array([]), np.array([1.0])))


def test_precision_artifact_shows_size_effect_even_with_identical_classifier():
    rows = make(n_cond=2, n_non=10)
    out = precision_size_artifact(rows)
    assert out["conditional"]["n_disorder_residues"] < out["non_conditional"]["n_disorder_residues"]
    # identical behaviour in both groups, yet precision differs purely because of group size
    assert out["conditional"]["observed_precision"] < out["non_conditional"]["observed_precision"]
    assert out["conditional"]["observed_precision"] == pytest.approx(
        out["conditional"]["equal_recall_precision"]
    )


def test_conditional_audit_has_no_recall_gap_for_identical_groups():
    out = conditional_flag_audit(make(), resamples_recall=200, resamples_auc=20, seed=3)
    assert out["n_conditional_proteins"] == 4
    assert out["recall_difference_non_minus_conditional"]["point"] == pytest.approx(0.0)
    lo, hi = out["recall_difference_non_minus_conditional"]["ci_95"]
    assert lo <= 0.0 <= hi
    assert out["auc_difference_non_minus_conditional"]["point"] == pytest.approx(0.0, abs=1e-9)


def test_conditional_audit_detects_a_real_recall_gap():
    rows = []
    for i in range(30):
        cond = i < 10
        for k in range(80):
            # conditional proteins: half the disorder residues look confident (pLDDT 80)
            plddt = 80.0 if cond and k % 2 == 0 else 50.0
            rows.append(res(f"Q{i}", plddt, True, cond, n=k))
            rows.append(res(f"Q{i}", 92.0, False, cond, n=1000 + k))
    out = conditional_flag_audit(rows, resamples_recall=300, resamples_auc=30, seed=4)
    assert out["recall_difference_non_minus_conditional"]["point"] == pytest.approx(0.5)
    assert out["recall_difference_non_minus_conditional"]["ci_95"][0] > 0.2


def test_evidence_audit_filters_small_codes_and_counts_proteins():
    rows = make()
    rows += [res("Z1", 30.0, True, False, ec="ECO:RARE", n=i) for i in range(10)]
    out = evidence_code_audit(rows, resamples=100, seed=5)
    assert "ECO:RARE" not in out
    assert out["ECO:1"]["n_proteins"] == 12
    assert out["ECO:1"]["recall"] == pytest.approx(1.0)


def test_build_subgroup_audit_roundtrip(tmp_path):
    path = tmp_path / "joined.csv"
    header = "acc,disprot_id,residue_number,plddt,is_disorder,is_conditional,ec_id,ec_name\n"
    lines = [
        f"{r.acc},{r.disprot_id},{r.residue_number},{r.plddt},{r.is_disorder},{r.is_conditional},{r.ec_id},{r.ec_name}"
        for r in make(per=210)
    ]
    path.write_text(header + "\n".join(lines) + "\n")
    out = build_subgroup_audit(path)
    assert out["label"].startswith("POST-HOC")
    assert set(out) >= {"precision_size_artifact", "conditional_flag", "evidence_codes"}
