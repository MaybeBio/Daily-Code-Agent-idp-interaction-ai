#!/usr/bin/env python3
"""Concatenate TSV files into one output file.

Usage:
    python concat_tsv.py                  # glob *.tsv in current directory
    python concat_tsv.py a.tsv b.tsv ...  # explicit file list
"""

import sys
from pathlib import Path

output_file = Path("concatenated.tsv")

if len(sys.argv) > 1:
    tsv_files = sorted(Path(f) for f in sys.argv[1:])
else:
    tsv_files = sorted(Path(".").glob("*.tsv"))

if not tsv_files:
    sys.exit("No .tsv files found.")

total_lines = 0
with output_file.open("w") as out_fh:
    for tsv in tsv_files:
        print(f"Reading: {tsv}")
        with tsv.open("r") as in_fh:
            for line in in_fh:
                if line.strip():
                    out_fh.write(line)
                    total_lines += 1

print(f"\nWritten {total_lines:,} lines to {output_file}")

