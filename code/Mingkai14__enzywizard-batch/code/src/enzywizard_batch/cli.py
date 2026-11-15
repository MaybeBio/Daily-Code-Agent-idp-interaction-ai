from __future__ import annotations
import argparse

from .commands.batch import add_batch_parser


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="enzywizard-batch",
        description="EnzyWizard-Batch: Run a complete EnzyWizard analysis workflow."
    )

    add_batch_parser(parser)
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    args.func(args)
