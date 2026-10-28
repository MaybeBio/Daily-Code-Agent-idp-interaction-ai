from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

from idrfeat.evaluate import freeze_holdout
from idrfeat.io import read_fasta


def _accessions_from_table(path: Path) -> list[str]:
    df = pd.read_parquet(path) if path.suffix == ".parquet" else pd.read_csv(path)
    return df["accession"].astype(str).tolist()


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="freeze_holdout")
    src = parser.add_mutually_exclusive_group(required=True)
    src.add_argument("--table", help="feature table parquet or csv with an accession column")
    src.add_argument("--fasta", help="protein FASTA (plain or .gz)")
    parser.add_argument("--out", default="config/holdout_proteins.txt")
    parser.add_argument("--frac", type=float, default=0.15)
    parser.add_argument("--seed", type=int, default=1729)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args(argv if argv is not None else sys.argv[1:])

    if args.table:
        accessions = _accessions_from_table(Path(args.table))
    else:
        accessions = list(read_fasta(args.fasta))
    if not accessions:
        print("no accessions found", file=sys.stderr)
        return 1

    hold = freeze_holdout(accessions, args.out, frac=args.frac, seed=args.seed, force=args.force)
    print(f"froze {len(hold)} of {len(set(accessions))} proteins to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
