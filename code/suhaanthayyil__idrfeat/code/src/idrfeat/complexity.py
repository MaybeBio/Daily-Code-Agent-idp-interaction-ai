"""Low-complexity and single-residue-run features.

Low complexity is called with SEG (Wootton and Federhen 1993) through the NCBI ``segmasker``
binary, the field-standard low-complexity tool, rather than a local rule. ``seg_masks`` runs it
once over a set of full sequences and returns a per-residue 0/1 mask per accession, so the
per-segment low-complexity fraction is the mean of the mask over the segment. When ``segmasker``
is not installed, low complexity is left unset (recorded as missing) rather than approximated.
The longest single-residue run is an exact, tool-independent sequence statistic.
"""

from __future__ import annotations

import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np


def segmasker_available() -> bool:
    return shutil.which("segmasker") is not None


def _run_segmasker(fasta_path: str) -> str:
    result = subprocess.run(
        ["segmasker", "-in", fasta_path, "-outfmt", "interval"],
        capture_output=True,
        text=True,
        check=True,
    )
    return result.stdout


def _parse_intervals(output: str, lengths: dict[str, int]) -> dict[str, np.ndarray]:
    masks = {acc: np.zeros(n, dtype=float) for acc, n in lengths.items()}
    current: str | None = None
    for line in output.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            current = line[1:].split()[0]
        elif current is not None and "-" in line:
            a, b = (int(x) for x in line.split("-"))
            mask = masks.get(current)
            if mask is not None:
                mask[a : b + 1] = 1.0
    return masks


def seg_masks(seqs: dict[str, str]) -> dict[str, np.ndarray]:
    """Per-residue low-complexity mask (1 = low complexity) for each accession, via SEG."""
    if not seqs:
        return {}
    with tempfile.TemporaryDirectory() as tmp:
        fasta = Path(tmp) / "input.fasta"
        with fasta.open("w") as fh:
            for acc, seq in seqs.items():
                fh.write(f">{acc}\n{seq}\n")
        output = _run_segmasker(str(fasta))
    return _parse_intervals(output, {acc: len(seq) for acc, seq in seqs.items()})


def seg_lowcomplexity_frac(seq: str) -> float:
    """Fraction of a single sequence masked as low complexity by SEG, or NaN if SEG is absent."""
    if not seq:
        return 0.0
    if not segmasker_available():
        return float("nan")
    mask = seg_masks({"seq": seq}).get("seq")
    return float(mask.mean()) if mask is not None and mask.size else 0.0


def longest_single_run(seq: str) -> int:
    """Length of the longest run of one identical residue."""
    best = run = 0
    prev = None
    for c in seq:
        run = run + 1 if c == prev else 1
        prev = c
        best = max(best, run)
    return best
