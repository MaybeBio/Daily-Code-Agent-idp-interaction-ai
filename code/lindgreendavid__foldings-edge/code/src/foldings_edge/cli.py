"""Command-line entry point for Foldings Edge."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from foldings_edge.registry import build_registry


def _cmd_summarize(args: argparse.Namespace) -> int:
    registry = build_registry(Path(args.residues))
    dataset = registry["dataset"]
    print(
        f"{dataset['total_residues']} residues from {dataset['total_proteins']} proteins "
        f"({dataset['disorder_annotated_residues']} DisProt-disorder-annotated, "
        f"{dataset['non_disorder_residues']} not)"
    )
    h1 = registry["h1_distribution"]["mann_whitney_u"]
    h2 = registry["h2_classifier"]["overall"]
    print(f"H1 Mann-Whitney U p-value: {h1['p_value']:.4g} (significant={h1['significant']})")
    print(
        f"H2 pLDDT<70 classifier: precision={h2['precision']['point']:.3f}, "
        f"recall={h2['recall']['point']:.3f}, f1={h2['f1']['point']:.3f}, "
        f"mcc={h2['mcc']['point']:.3f}"
    )
    return 0


def _cmd_registry(args: argparse.Namespace) -> int:
    registry = build_registry(Path(args.residues))
    text = json.dumps(registry, indent=2, sort_keys=True) + "\n"
    if args.output:
        Path(args.output).write_text(text)
    else:
        sys.stdout.write(text)
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="foldings-edge", description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    summarize = subparsers.add_parser("summarize", help="Print a short summary of the dataset")
    summarize.add_argument("residues", help="Path to the joined per-residue CSV")
    summarize.set_defaults(func=_cmd_summarize)

    registry = subparsers.add_parser("registry", help="Build and print the full result registry")
    registry.add_argument("residues", help="Path to the joined per-residue CSV")
    registry.add_argument("--output", help="Write JSON to this path instead of stdout")
    registry.set_defaults(func=_cmd_registry)

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    result: int = args.func(args)
    return result


if __name__ == "__main__":
    sys.exit(main())
