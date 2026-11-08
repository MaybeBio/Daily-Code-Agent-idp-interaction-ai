from __future__ import annotations

import subprocess
import sys
from pathlib import Path


def test_academic_sensitivity_regenerates_deterministically() -> None:
    root = Path(__file__).resolve().parents[1]
    output = root / "reports" / "post-release-academic-sensitivity.json"
    expected = output.read_bytes()

    subprocess.run(
        [sys.executable, str(root / "scripts" / "generate_academic_sensitivity.py")],
        cwd=root,
        check=True,
    )

    assert output.read_bytes() == expected
