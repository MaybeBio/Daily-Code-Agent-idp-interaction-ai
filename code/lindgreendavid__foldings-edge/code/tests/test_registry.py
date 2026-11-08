from __future__ import annotations

import subprocess
import sys
from pathlib import Path

from foldings_edge.registry import build_registry

REPO_ROOT = Path(__file__).resolve().parent.parent


def test_build_registry_schema(sample_residues_csv: Path) -> None:
    registry = build_registry(sample_residues_csv)
    assert registry["schema_version"] == 1
    dataset = registry["dataset"]
    assert dataset["total_residues"] == 35
    assert dataset["disorder_annotated_residues"] == 15
    assert dataset["non_disorder_residues"] == 20
    assert dataset["conditional_disorder_residues"] == 5
    assert dataset["total_proteins"] == 2

    h1 = registry["h1_distribution"]
    assert 0.0 <= h1["mann_whitney_u"]["p_value"] <= 1.0
    assert 0.0 <= h1["ks_robustness"]["p_value"] <= 1.0

    h2 = registry["h2_classifier"]
    assert 0.0 <= h2["overall"]["precision"]["point"] <= 1.0
    assert 0.0 <= h2["overall"]["recall"]["point"] <= 1.0
    assert "conditional_disorder_regions" in h2["by_conditional_flag"]
    assert "non_conditional_disorder_regions" in h2["by_conditional_flag"]

    # Protein B (15 residues) is under the min-size cutoff (20); protein A (20
    # residues, exactly at the cutoff) is included with a perfect classifier.
    breakdown = registry["protein_breakdown"]
    assert len(breakdown) == 1
    assert breakdown[0]["acc"] == "A00001"
    assert breakdown[0]["false_positive_residues"] == 0
    assert breakdown[0]["false_negative_residues"] == 0

    dist = registry["distributions"]
    assert len(dist["inside_disorder"]) == 15
    assert len(dist["outside_disorder"]) == 20


def test_build_registry_evidence_code_breakdown_respects_minimum(
    sample_residues_csv: Path,
) -> None:
    registry = build_registry(sample_residues_csv)
    # Both evidence codes in the fixture have far fewer than
    # MIN_EVIDENCE_CODE_RESIDUES (200) disorder residues, so the breakdown
    # should be empty on this tiny synthetic sample.
    assert registry["h2_classifier"]["by_evidence_code"] == {}


def test_registry_generation_script_matches_frozen_registry() -> None:
    """The frozen registry must match within cross-platform float tolerance."""
    residues = REPO_ROOT / "data" / "external" / "joined_residues.csv"
    frozen = REPO_ROOT / "reports" / "v0.1-foldings-edge-registry.json"
    if not residues.exists():
        return  # joined residues CSV may be absent in some environments; only runs if fetched
    result = subprocess.run(
        [
            sys.executable,
            str(REPO_ROOT / "scripts" / "generate_registry.py"),
            "--output",
            "/tmp/foldings-edge-registry-check.json",
            "--residues",
            str(residues),
        ],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        check=True,
    )
    assert result.returncode == 0
    comparison = subprocess.run(
        [
            sys.executable,
            str(REPO_ROOT / "scripts" / "compare_registry.py"),
            str(frozen),
            "/tmp/foldings-edge-registry-check.json",
        ],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        check=True,
    )
    assert comparison.returncode == 0
