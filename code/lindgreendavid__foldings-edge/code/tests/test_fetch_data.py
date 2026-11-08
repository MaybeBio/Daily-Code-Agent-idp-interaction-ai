from __future__ import annotations

import csv
import importlib.util
import json
import sys
from pathlib import Path
from typing import Any

import pytest

SCRIPT_PATH = Path(__file__).resolve().parent.parent / "scripts" / "fetch_data.py"
_spec = importlib.util.spec_from_file_location("fetch_data", SCRIPT_PATH)
assert _spec is not None and _spec.loader is not None
fetch_data = importlib.util.module_from_spec(_spec)
sys.modules["fetch_data"] = fetch_data
_spec.loader.exec_module(fetch_data)


def _entry(
    disprot_id: str, acc: str, sequence: str, regions: list[dict[str, Any]]
) -> dict[str, Any]:
    return {"disprot_id": disprot_id, "acc": acc, "sequence": sequence, "regions": regions}


def test_sample_entries_deterministic_given_seed() -> None:
    entries = [_entry(f"DP{i:05d}", f"ACC{i}", "M" * 10, []) for i in range(50)]
    first = fetch_data.sample_entries(entries, 10, seed=42)
    second = fetch_data.sample_entries(entries, 10, seed=42)
    assert [e["disprot_id"] for e in first] == [e["disprot_id"] for e in second]
    assert len(first) == 10


def test_sample_entries_different_seed_can_differ() -> None:
    entries = [_entry(f"DP{i:05d}", f"ACC{i}", "M" * 10, []) for i in range(50)]
    a = fetch_data.sample_entries(entries, 10, seed=1)
    b = fetch_data.sample_entries(entries, 10, seed=2)
    assert [e["disprot_id"] for e in a] != [e["disprot_id"] for e in b]


def test_sample_entries_caps_at_available_size() -> None:
    entries = [_entry(f"DP{i:05d}", f"ACC{i}", "M" * 10, []) for i in range(5)]
    sampled = fetch_data.sample_entries(entries, 10, seed=1)
    assert len(sampled) == 5


def test_disorder_regions_filters_to_structural_state_disorder() -> None:
    entry = _entry(
        "DP00001",
        "P00001",
        "M" * 20,
        [
            {
                "term_namespace": "Structural state",
                "term_name": "disorder",
                "start": 1,
                "end": 5,
                "ec_id": "ECO:1",
                "ec_name": "evidence one",
            },
            {
                "term_namespace": "Structural transition",
                "term_name": "disorder to order",
                "start": 6,
                "end": 10,
                "ec_id": "ECO:2",
                "ec_name": "evidence two",
            },
            {"term_namespace": "Molecular function", "term_name": "binding", "start": 1, "end": 3},
        ],
    )
    regions = fetch_data._disorder_regions(entry)
    assert regions == [(1, 5, "ECO:1", "evidence one")]


def test_has_transition_region() -> None:
    with_transition = _entry(
        "DP00001",
        "P00001",
        "MMM",
        [
            {
                "term_namespace": "Structural transition",
                "term_name": "disorder to order",
                "start": 1,
                "end": 2,
            }
        ],
    )
    without_transition = _entry(
        "DP00002",
        "P00002",
        "MMM",
        [{"term_namespace": "Structural state", "term_name": "disorder", "start": 1, "end": 2}],
    )
    assert fetch_data._has_transition_region(with_transition) is True
    assert fetch_data._has_transition_region(without_transition) is False


def test_process_entry_no_alphafold_prediction(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(fetch_data, "fetch_alphafold_metadata", lambda acc: None)
    entry = _entry("DP00001", "P00001", "MMM", [])
    result = fetch_data.process_entry(entry)
    assert result == {
        "status": "excluded",
        "reason": "no_alphafold_prediction",
        "acc": "P00001",
        "disprot_id": "DP00001",
    }


def test_process_entry_fragmented(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(
        fetch_data,
        "fetch_alphafold_metadata",
        lambda acc: [{"uniprotAccession": acc}, {"uniprotAccession": acc}],
    )
    entry = _entry("DP00001", "P00001", "MMM", [])
    result = fetch_data.process_entry(entry)
    assert result["status"] == "excluded"
    assert result["reason"] == "fragmented_alphafold_entry"
    assert result["n_fragments"] == 2


def test_process_entry_sequence_length_mismatch(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(
        fetch_data,
        "fetch_alphafold_metadata",
        lambda acc: [{"uniprotAccession": acc, "uniprotSequence": "MM"}],
    )
    entry = _entry("DP00001", "P00001", "MMM", [])  # length 3 vs alphafold length 2
    result = fetch_data.process_entry(entry)
    assert result["status"] == "excluded"
    assert result["reason"] == "sequence_length_mismatch"
    assert result["disprot_length"] == 3
    assert result["alphafold_length"] == 2


def test_process_entry_included_joins_plddt_and_disorder(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(
        fetch_data,
        "fetch_alphafold_metadata",
        lambda acc: [
            {"uniprotAccession": acc, "uniprotSequence": "MMMMM", "plddtDocUrl": "http://x"}
        ],
    )
    monkeypatch.setattr(
        fetch_data,
        "fetch_confidence",
        lambda url: {
            "residueNumber": [1, 2, 3, 4, 5],
            "confidenceScore": [10.0, 20.0, 90.0, 95.0, 30.0],
        },
    )
    entry = _entry(
        "DP00001",
        "P00001",
        "MMMMM",
        [
            {
                "term_namespace": "Structural state",
                "term_name": "disorder",
                "start": 1,
                "end": 2,
                "ec_id": "ECO:1",
                "ec_name": "evidence one",
            }
        ],
    )
    result = fetch_data.process_entry(entry)
    assert result["status"] == "included"
    rows = result["rows"]
    assert len(rows) == 5
    assert rows[0]["is_disorder"] is True
    assert rows[0]["ec_id"] == "ECO:1"
    assert rows[2]["is_disorder"] is False
    assert rows[2]["ec_id"] == ""


def test_run_writes_csv_and_exclusions(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    entries = [
        _entry(
            "DP00001",
            "P00001",
            "MMM",
            [
                {
                    "term_namespace": "Structural state",
                    "term_name": "disorder",
                    "start": 1,
                    "end": 1,
                    "ec_id": "",
                    "ec_name": "",
                }
            ],
        ),
        _entry("DP00002", "P00002", "MM", []),
    ]
    monkeypatch.setattr(fetch_data, "fetch_all_human_disprot", lambda cache_path: entries)
    monkeypatch.setattr(
        fetch_data,
        "fetch_alphafold_metadata",
        lambda acc: (
            [{"uniprotAccession": acc, "uniprotSequence": "MMM", "plddtDocUrl": "http://x"}]
            if acc == "P00001"
            else None
        ),
    )
    monkeypatch.setattr(
        fetch_data,
        "fetch_confidence",
        lambda url: {"residueNumber": [1, 2, 3], "confidenceScore": [10.0, 80.0, 90.0]},
    )

    output_csv = tmp_path / "joined.csv"
    exclusions_path = tmp_path / "exclusions.json"
    raw_cache = tmp_path / "raw.json"

    fetch_data.run(2, output_csv, exclusions_path, raw_cache, max_workers=2)

    with output_csv.open() as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 3
    assert all(row["acc"] == "P00001" for row in rows)

    exclusions = json.loads(exclusions_path.read_text())
    assert exclusions["included_proteins"] == 1
    assert exclusions["excluded_proteins"] == 1
    assert exclusions["exclusions_by_reason"] == {"no_alphafold_prediction": 1}
