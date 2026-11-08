import textwrap

import pytest

from zpoly.cli import main


@pytest.fixture
def run_inp(tmp_path, monkeypatch):
    """run_inp("KEY = value ...") writes job.inp in tmp_path (the working directory) and runs it."""
    monkeypatch.chdir(tmp_path)

    def run(text: str, name: str = "job.inp") -> int:
        path = tmp_path / name
        path.write_text(textwrap.dedent(text).strip() + "\n")
        return main([str(path)])

    return run
