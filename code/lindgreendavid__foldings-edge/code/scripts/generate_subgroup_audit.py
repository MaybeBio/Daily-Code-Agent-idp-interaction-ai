#!/usr/bin/env python3
"""Generate the POST-HOC subgroup audit (not part of the frozen v0.1 registry).

Usage: python scripts/generate_subgroup_audit.py --output reports/post-release-subgroup-audit.json
       [--residues data/external/joined_residues.csv]
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from foldings_edge.subgroups import build_subgroup_audit


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    parser.add_argument("--residues", default="data/external/joined_residues.csv")
    args = parser.parse_args()
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w") as handle:
        json.dump(build_subgroup_audit(Path(args.residues)), handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
