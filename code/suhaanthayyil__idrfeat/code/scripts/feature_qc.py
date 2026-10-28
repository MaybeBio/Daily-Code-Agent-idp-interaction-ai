from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

from idrfeat.qc import feature_summary, high_correlation_pairs, near_constant


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="feature_qc")
    parser.add_argument("--table", required=True, help="feature table parquet or csv")
    parser.add_argument("--outdir", default="data/processed/qc")
    parser.add_argument("--corr-threshold", dest="corr_threshold", type=float, default=0.95)
    args = parser.parse_args(argv if argv is not None else sys.argv[1:])

    path = Path(args.table)
    df = pd.read_parquet(path) if path.suffix == ".parquet" else pd.read_csv(path)
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    summary = feature_summary(df)
    summary.to_csv(out / "feature_summary.csv", index=False)
    pairs = high_correlation_pairs(df, threshold=args.corr_threshold)
    pairs.to_csv(out / "high_correlation_pairs.csv", index=False)
    nc = near_constant(df)

    print(f"features summarized: {len(summary)}")
    print(f"near-constant features: {nc}")
    print(f"high-correlation pairs (>= {args.corr_threshold}): {len(pairs)}")
    print(f"wrote {out}/feature_summary.csv and {out}/high_correlation_pairs.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
