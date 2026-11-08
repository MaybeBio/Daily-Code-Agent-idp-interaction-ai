#!/usr/bin/env python3
"""Generate the frozen v0.1.0 result registry deterministically.

Usage: python scripts/generate_registry.py --output reports/v0.1-foldings-edge-registry.json
       [--residues data/external/joined_residues.csv]

Requires the joined per-residue CSV to already exist locally (see
``scripts/fetch_data.py``). Regenerates byte-identically given a fixed
residues CSV and the fixed statistics seed in ``foldings_edge.stats``, and is
checked against the committed registry by ``tests/test_registry.py`` and the
CI ``research-registry`` job.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from foldings_edge.registry import build_registry


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    parser.add_argument(
        "--residues",
        default="data/external/joined_residues.csv",
        help="Path to the joined per-residue CSV produced by scripts/fetch_data.py",
    )
    args = parser.parse_args()

    registry = build_registry(Path(args.residues))
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(registry, indent=2, sort_keys=True) + "\n")
    print(f"wrote {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
