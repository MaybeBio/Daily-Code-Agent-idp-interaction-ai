#!/usr/bin/env python3
"""
prepare_af3_inputs_gencode.py
=============================

Bridge between Phase 5's predict-queue (data/gencode/v47/predict_queue_v47/)
and the Leonardo AF3 SLURM pipeline (script.02a / script.02b).

Reads `manifest.tsv` produced by build_predict_queue.py (already sorted
priority_tier ASC → length ASC) and emits per-sequence AF3 JSONs as
`af3_inputs/seq_NNNNN.json`, in the same on-disk format the existing 02a /
02b scripts consume.

The tier-then-length ordering is the whole point of running this rather
than just concatenating the FASTAs and reusing script.01.prep_af3_inputs:
it lets you submit a MANE-only first wave with a clean SLURM array range,
then alt-GENCODE, then novels, instead of every tier interleaved.

JSON schema (AF3 v3.0+, matches the 02a/02b `--json_path=` expectation):

    {
      "name": "<md5>",                       ← lowercased, becomes
                                                af3_outputs/<md5>/ on Leonardo
      "modelSeeds": [1],
      "sequences": [{"protein": {"id": "A", "sequence": "MVLS..."}}],
      "dialect": "alphafold3",
      "version": 2
    }

Inputs
------
  --queue-dir   path to predict_queue_v47/  (contains manifest.tsv + fastas/)
  --out-dir     where to write seq_NNNNN.json files (typically ./af3_inputs)
  --start-idx   first seq index (default 1).  Use this if you're appending
                onto an existing af3_inputs/ directory.

Outputs
-------
  <out-dir>/seq_NNNNN.json                  one per md5
  <out-dir>/prep_manifest.tsv               seq_idx ↔ md5 ↔ target_id ↔
                                            tier ↔ queued_for ↔ seq_len

Usage
-----
    python prepare_af3_inputs_gencode.py \\
        --queue-dir gencode/v47/predict_queue_v47 \\
        --out-dir   af3_inputs
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from pathlib import Path
from typing import List, Tuple


# ── FASTA reader ─────────────────────────────────────────────────────────────

def read_fasta_sequence(path: Path) -> str:
    """Single-record FASTA → bare protein string, uppercased, '*' stripped."""
    seq_parts: List[str] = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                continue
            seq_parts.append(line.strip())
    return "".join(seq_parts).upper().rstrip("*")


# ── Main ─────────────────────────────────────────────────────────────────────

def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--queue-dir", required=True, type=Path,
                    help="predict_queue_v47/ produced by build_predict_queue.py")
    ap.add_argument("--out-dir", required=True, type=Path,
                    help="Where to write seq_NNNNN.json files (e.g. af3_inputs)")
    ap.add_argument("--start-idx", type=int, default=1,
                    help="First seq index (default: 1 → seq_00001.json)")
    ap.add_argument("--model-seeds", type=int, nargs="+", default=[1],
                    help="AF3 diffusion seeds (default: [1])")
    args = ap.parse_args()

    manifest_in = args.queue_dir / "manifest.tsv"
    fasta_dir   = args.queue_dir / "fastas"
    if not manifest_in.exists():
        sys.exit(f"error: {manifest_in} not found")
    if not fasta_dir.exists():
        sys.exit(f"error: {fasta_dir} not found")

    args.out_dir.mkdir(parents=True, exist_ok=True)

    # ── Read queue manifest (already tier-then-length sorted) ───────────────
    with open(manifest_in, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    print(f"queue manifest: {len(rows):,} rows", file=sys.stderr)

    # ── Write one JSON per row, in input order ──────────────────────────────
    prep_rows: List[Tuple[int, str, str, str, int, int, str]] = []
    tier_bounds = {1: [None, None], 2: [None, None], 3: [None, None]}

    written = 0
    for offset, row in enumerate(rows):
        seq_idx     = args.start_idx + offset
        md5         = row["protein_md5"].lower()
        target_id   = row["target_id"]
        tier        = int(row["priority_tier"])
        queued_for  = row["queued_for"]
        length      = int(row["protein_length"])
        fasta_name  = row["fasta_filename"]

        fa_path = fasta_dir / fasta_name
        if not fa_path.exists():
            sys.exit(f"error: FASTA missing: {fa_path}")
        sequence = read_fasta_sequence(fa_path)
        if len(sequence) != length:
            print(f"warn: length mismatch for {target_id} "
                  f"({len(sequence)} vs manifest {length})", file=sys.stderr)

        # Sanitise like the existing script.01 does (md5 is already
        # [0-9a-f]+, but be defensive in case future queues use ids).
        name = re.sub(r"[^a-z0-9_]", "_", md5.lower())

        job = {
            "name":       name,
            "modelSeeds": list(args.model_seeds),
            "sequences":  [{"protein": {"id": "A", "sequence": sequence}}],
            "dialect":    "alphafold3",
            "version":    2,
        }
        json_name = f"seq_{seq_idx:05d}.json"
        with open(args.out_dir / json_name, "w") as fh:
            json.dump(job, fh, indent=2)

        prep_rows.append((seq_idx, json_name, name, target_id, tier,
                          length, queued_for))

        # Track [start, end] index per tier
        if tier_bounds[tier][0] is None:
            tier_bounds[tier][0] = seq_idx
        tier_bounds[tier][1] = seq_idx
        written += 1

    # ── Write prep manifest (bridges seq_idx → DB target_id) ────────────────
    prep_manifest = args.out_dir / "prep_manifest.tsv"
    with open(prep_manifest, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["seq_idx", "json_file", "md5", "target_id",
                    "priority_tier", "protein_length", "queued_for"])
        w.writerows(prep_rows)

    # ── Summary + SLURM submission lines ────────────────────────────────────
    print()
    print(f"wrote {written:,} json files to {args.out_dir}/")
    print(f"prep manifest: {prep_manifest}")
    print()
    print("─── SLURM submission lines (50 seqs/MSA task, 10 seqs/infer task) ───")
    for tier, label in [(1, "MANE Select"),
                        (2, "Alt GENCODE"),
                        (3, "Novel DB"   )]:
        lo, hi = tier_bounds[tier]
        if lo is None:
            continue
        n_seqs = hi - lo + 1
        msa_tasks   = (n_seqs + 49) // 50
        infer_tasks = (n_seqs + 9)  // 10
        print()
        print(f"# {label}  (seq {lo}..{hi}, {n_seqs:,} sequences)")
        print(f"#   MSA      (~12.5h/task, 24h walltime):")
        print(f"sbatch --export=ALL,SEQ_START={lo},SEQS_PER_TASK=50 \\")
        print(f"       --array=1-{msa_tasks}%50 script.02a.data_pipeline.v1.1.slurm")
        print(f"#   Inference (~2h/task; wait for MSA to substantially complete):")
        print(f"sbatch --export=ALL,SEQ_START={lo},SEQS_PER_TASK=10 \\")
        print(f"       --array=1-{infer_tasks}%50 script.02b.inference.v1.1.slurm")


if __name__ == "__main__":
    main()
