from __future__ import annotations

import json
from pathlib import Path

import pytest

from foldings_edge.cli import main


def test_summarize_command(sample_residues_csv: Path, capsys: pytest.CaptureFixture[str]) -> None:
    exit_code = main(["summarize", str(sample_residues_csv)])
    assert exit_code == 0
    out = capsys.readouterr().out
    assert "residues from" in out
    assert "H1 Mann-Whitney U p-value" in out
    assert "H2 pLDDT<70 classifier" in out


def test_registry_command_stdout(
    sample_residues_csv: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    exit_code = main(["registry", str(sample_residues_csv)])
    assert exit_code == 0
    out = capsys.readouterr().out
    data = json.loads(out)
    assert data["schema_version"] == 1


def test_registry_command_writes_file(sample_residues_csv: Path, tmp_path: Path) -> None:
    output = tmp_path / "out.json"
    exit_code = main(["registry", str(sample_residues_csv), "--output", str(output)])
    assert exit_code == 0
    data = json.loads(output.read_text())
    assert data["schema_version"] == 1


def test_no_command_raises_system_exit() -> None:
    with pytest.raises(SystemExit):
        main([])
