"""
Demo for the madrigal-predictor package.

Scores the sixteen yeast transcription factors of the held-out classification
test set with the trained MADRIGAL model, and reports the activation domains it
predicts in each. These are the same sixteen proteins (12,981 residues) used to
evaluate the model in the paper, so the demo exercises the model on sequences it
was never trained on.

Requires only the madrigal-predictor package. No GPU and no cluster are needed.

Usage, from the repository root:

    python demo/run_demo.py

Output is written to demo/output/. See demo/README.md for the expected output
and runtime.
"""

from pathlib import Path
import time

t_import = time.perf_counter()
from madrigal_predictor import analyze_fasta_file
import_time = time.perf_counter() - t_import

HERE = Path(__file__).resolve().parent
FASTA = HERE / "demo_seqs.fasta"
OUT = HERE / "output"
OUT.mkdir(exist_ok=True)

THRESHOLD = 0.6     # passed to the model; residues at or above this are called AD
MIN_LENGTH = 10     # applied by this script afterwards, when listing results.
                    # The package has no length filter and returns every run
                    # above the threshold. Set to 1 to list all of them.

print("MADRIGAL demo: sixteen yeast transcription factors (held-out test set)")
print(f"Consensus threshold {THRESHOLD} (package default)")
print(f"Listing domains of at least {MIN_LENGTH} aa; shorter runs are omitted "
      f"by this script\n")

t0 = time.perf_counter()
titles, above_runs, below_runs, above_starts, above_ends = analyze_fasta_file(
    FASTA, threshold=THRESHOLD
)
compute_time = time.perf_counter() - t0

rows = []
no_ad = []
for title, runs, starts, ends in zip(titles, above_runs, above_starts, above_ends):
    kept = [(s, a, b) for s, a, b in zip(runs, starts, ends)
            if len(s) >= MIN_LENGTH]
    if not kept:
        no_ad.append(title)
        print(f"{title:8s}  no activation domain predicted")
        continue
    print(f"{title:8s}  {len(kept)} predicted activation domain"
          f"{'s' if len(kept) > 1 else ''}")
    for seq, start, end in kept:
        # coordinates from the package are 0-based and inclusive;
        # report 1-based residue numbers, as used in the literature
        print(f"          residues {start + 1:4d}-{end + 1:4d}"
              f"  ({len(seq):3d} aa)  {seq}")
        rows.append((title, start + 1, end + 1, len(seq), seq))

outfile = OUT / "predicted_ADs.csv"
with outfile.open("w", encoding="utf-8") as fh:
    fh.write("protein,start,end,length,sequence\n")
    for r in rows:
        fh.write(f"{r[0]},{r[1]},{r[2]},{r[3]},{r[4]}\n")

print(f"\n{len(rows)} predicted activation domains across "
      f"{len(titles) - len(no_ad)} of {len(titles)} proteins.")
if no_ad:
    print(f"No domain predicted in: {', '.join(no_ad)}")
print(f"Written to {outfile}")
print(f"\nModel evaluation: {compute_time:.2f} s "
      f"(plus {import_time:.1f} s to import the model and PyTorch)")
