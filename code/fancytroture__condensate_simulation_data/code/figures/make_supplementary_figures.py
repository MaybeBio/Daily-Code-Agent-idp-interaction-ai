#!/usr/bin/env python3
"""Build all or selected manuscript supporting figures."""

from __future__ import annotations

import argparse
from pathlib import Path

from supplementary import DEFAULT_OUT, FIGURE_BUILDERS
from supplementary.shared import setup_style, write_readme


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--only", default="", help="Comma-separated figure numbers, e.g. 1,5,8")
    args = parser.parse_args()
    requested = (
        {int(value) for value in args.only.split(",") if value.strip()}
        if args.only
        else set(FIGURE_BUILDERS)
    )
    unknown = requested.difference(FIGURE_BUILDERS)
    if unknown:
        parser.error(f"unknown supplementary figure number(s): {sorted(unknown)}")

    setup_style()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for number in sorted(requested):
        FIGURE_BUILDERS[number](args.output_dir)
    write_readme(args.output_dir)


if __name__ == "__main__":
    main()
