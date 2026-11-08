import json
from pathlib import Path

ROOT = Path(__file__).parent.parent
TEX = (ROOT / "paper" / "paper.tex").read_text(encoding="utf-8")
AUDIT = json.loads((ROOT / "reports" / "post-release-subgroup-audit.json").read_text())
REGISTRY = json.loads((ROOT / "reports" / "v0.1-foldings-edge-registry.json").read_text())


def test_audit_numbers_in_paper_match_the_audit_file():
    cond = AUDIT["conditional_flag"]
    p = cond["point"]
    assert f"{p['recall_conditional']:.3f}" in TEX and f"{p['recall_non_conditional']:.3f}" in TEX
    assert f"{p['auc_conditional']:.3f}" in TEX and f"{p['auc_non_conditional']:.3f}" in TEX
    assert f"{cond['recall_difference_non_minus_conditional']['point']:.3f}" in TEX
    assert f"{cond['auc_difference_non_minus_conditional']['point']:.3f}" in TEX
    art = AUDIT["precision_size_artifact"]
    for group in ("conditional", "non_conditional"):
        assert f"{art[group]['equal_recall_precision']:.3f}" in TEX
        assert f"{art[group]['observed_precision']:.3f}" in TEX


def test_evidence_code_claims_in_paper_match_the_audit_file():
    hdx = AUDIT["evidence_codes"]["ECO:0006236"]
    assert hdx["n_proteins"] == 4 and "four proteins" in TEX
    assert f"{hdx['recall']:.3f}" in TEX


def test_registry_headline_numbers_in_paper():
    o = REGISTRY["h2_classifier"]["overall"]
    for key in ("precision", "recall", "f1", "mcc"):
        assert f"{o[key]['point']:.3f}" in TEX
    c = o["confusion"]
    assert f"{c['true_positive']:,}".replace(",", "{,}") in TEX
