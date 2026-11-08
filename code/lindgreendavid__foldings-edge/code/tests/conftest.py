"""Shared pytest fixtures: a small synthetic joined-residue CSV for unit tests.

This fixture data is NOT real AlphaFold DB / DisProt data — it exists only to
exercise the parsing, splitting, and statistics code paths deterministically
and quickly. Real-data verification lives in scripts/fetch_data.py and is run
against the live DisProt/AlphaFold DB APIs, not in this offline test suite.
"""

from __future__ import annotations

import csv
from pathlib import Path

import pytest

FIELDNAMES = [
    "acc",
    "disprot_id",
    "residue_number",
    "plddt",
    "is_disorder",
    "is_conditional",
    "ec_id",
    "ec_name",
]


def _row(**overrides: object) -> dict[str, object]:
    base: dict[str, object] = {
        "acc": "P00000",
        "disprot_id": "DP00000",
        "residue_number": 1,
        "plddt": 50.0,
        "is_disorder": False,
        "is_conditional": False,
        "ec_id": "",
        "ec_name": "",
    }
    base.update(overrides)
    return base


@pytest.fixture
def sample_rows() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    # Protein A: 20 residues, 1-10 disordered with low pLDDT (true positives),
    # 11-20 ordered with high pLDDT (true negatives). Not conditional.
    for i in range(1, 11):
        rows.append(
            _row(
                acc="A00001",
                disprot_id="DP00001",
                residue_number=i,
                plddt=30.0 + i,
                is_disorder=True,
                ec_id="ECO:0000006",
                ec_name="experimental evidence",
            )
        )
    for i in range(11, 21):
        rows.append(_row(acc="A00001", disprot_id="DP00001", residue_number=i, plddt=90.0 - i / 2))
    # Protein B: 15 residues, 1-5 disordered but CONFIDENT pLDDT (false negatives
    # for the classifier / "conditionally folded"-like), 6-15 ordered, low pLDDT
    # for a couple (false positives).
    for i in range(1, 6):
        rows.append(
            _row(
                acc="B00002",
                disprot_id="DP00002",
                residue_number=i,
                plddt=85.0,
                is_disorder=True,
                is_conditional=True,
                ec_id="ECO:0000024",
                ec_name="X-ray crystallography evidence",
            )
        )
    for i in range(6, 16):
        plddt = 40.0 if i < 8 else 92.0
        rows.append(_row(acc="B00002", disprot_id="DP00002", residue_number=i, plddt=plddt))
    return rows


@pytest.fixture
def sample_residues_csv(tmp_path: Path, sample_rows: list[dict[str, object]]) -> Path:
    path = tmp_path / "joined_residues.csv"
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES)
        writer.writeheader()
        for row in sample_rows:
            out = dict(row)
            out["is_disorder"] = str(row["is_disorder"])
            out["is_conditional"] = str(row["is_conditional"])
            writer.writerow(out)
    return path
