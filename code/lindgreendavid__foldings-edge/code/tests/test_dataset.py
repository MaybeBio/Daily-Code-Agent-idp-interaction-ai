from __future__ import annotations

from pathlib import Path

from foldings_edge.dataset import (
    evidence_code_counts,
    load_residues,
    split_by_disorder,
    unique_proteins,
)


def test_load_residues_parses_all_rows(sample_residues_csv: Path) -> None:
    residues = load_residues(sample_residues_csv)
    assert len(residues) == 35  # 20 (protein A) + 15 (protein B)


def test_load_residues_parses_booleans(sample_residues_csv: Path) -> None:
    residues = load_residues(sample_residues_csv)
    conditional = [r for r in residues if r.is_conditional]
    assert len(conditional) == 5
    assert all(r.is_disorder for r in conditional)


def test_split_by_disorder(sample_residues_csv: Path) -> None:
    residues = load_residues(sample_residues_csv)
    inside, outside = split_by_disorder(residues)
    assert len(inside) == 15  # 10 (A) + 5 (B)
    assert len(outside) == 20  # 10 (A) + 10 (B)
    assert all(r.is_disorder for r in inside)
    assert all(not r.is_disorder for r in outside)


def test_unique_proteins(sample_residues_csv: Path) -> None:
    residues = load_residues(sample_residues_csv)
    assert unique_proteins(residues) == {"A00001", "B00002"}


def test_evidence_code_counts(sample_residues_csv: Path) -> None:
    residues = load_residues(sample_residues_csv)
    counts = evidence_code_counts(residues)
    assert counts == {"ECO:0000006": 10, "ECO:0000024": 5}
