"""Generate NARDINI sequence-patterning statistics."""
import os
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
                      

import argparse
import contextlib
import io
import warnings
from multiprocessing import Pool

import numpy as np
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from nardini.constants import TYPEALL
from nardini.score_and_plot import calculate_zscore

AA_GROUP_NAMES = ["".join(group) for group in TYPEALL]
NARDINI_FEATURE_NAMES = [
    f"nardini_{left}__{right}"
    for left_index, left in enumerate(AA_GROUP_NAMES)
    for right_index, right in enumerate(AA_GROUP_NAMES)
    if left_index <= right_index
]

def calculate_one(task):
    row_index, uid, start, end, sequence, num_scrambles, base_seed = task
    record = SeqRecord(Seq(sequence), id=str(row_index), description=f"{uid}:{start}-{end}")
    with contextlib.redirect_stdout(io.StringIO()), warnings.catch_warnings():
        warnings.simplefilter("ignore")
        result = calculate_zscore(
            [record], TYPEALL, num_scrambles, base_seed + row_index
        )
    matrix = np.asarray(result[str(row_index)][3], dtype=float)
    vector = matrix[np.triu_indices(len(TYPEALL))]
    if vector.shape != (36,) or not np.isfinite(vector).all():
        raise ValueError(f"Invalid NARDINI result for {uid}:{start}-{end}: {vector.shape}")
    return [uid, start, end, *vector.tolist()]

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--idr-manifest", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--num-scrambles", type=int, default=1000)
    parser.add_argument("--random-seed", type=int, default=1)
    parser.add_argument("--workers", type=int, default=1)
    return parser.parse_args()

def main():
    args = parse_args()
    frame = pd.read_csv(args.idr_manifest, sep="\t")
    required = {"uniprot", "idr_start", "idr_end", "sequence"}
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f"IDR manifest is missing columns: {sorted(missing)}")
    tasks = [
        (
            index,
            str(row.uniprot),
            int(row.idr_start),
            int(row.idr_end),
            str(row.sequence),
            args.num_scrambles,
            args.random_seed,
        )
        for index, row in enumerate(frame.itertuples(index=False))
    ]
    if args.workers == 1:
        rows = [calculate_one(task) for task in tasks]
    else:
        with Pool(args.workers) as pool:
            rows = list(pool.imap(calculate_one, tasks, chunksize=1))
    output = pd.DataFrame(
        rows,
        columns=["uniprot", "idr_start", "idr_end", *NARDINI_FEATURE_NAMES],
    )
    output.to_csv(args.output, sep="\t", index=False)
    print(f"[OK] NARDINI cache rows={len(output)} output={args.output}")

if __name__ == "__main__":
    main()
