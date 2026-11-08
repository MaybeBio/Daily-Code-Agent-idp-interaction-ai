"""Load the joined per-residue pLDDT / DisProt-disorder CSV.

See `docs/research-protocol.md` for the exact sample construction, join
procedure, and exclusion criteria implemented by `scripts/fetch_data.py`,
which produces the CSV this module reads.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Residue:
    """One joined residue: AlphaFold DB pLDDT plus a DisProt disorder label."""

    acc: str
    disprot_id: str
    residue_number: int
    plddt: float
    is_disorder: bool
    is_conditional: bool
    ec_id: str
    ec_name: str


def _row_to_residue(row: dict[str, str]) -> Residue:
    return Residue(
        acc=row["acc"],
        disprot_id=row["disprot_id"],
        residue_number=int(row["residue_number"]),
        plddt=float(row["plddt"]),
        is_disorder=row["is_disorder"].strip().lower() == "true",
        is_conditional=row["is_conditional"].strip().lower() == "true",
        ec_id=row["ec_id"],
        ec_name=row["ec_name"],
    )


def load_residues(path: Path) -> list[Residue]:
    """Load every joined residue row from the CSV produced by scripts/fetch_data.py."""
    with Path(path).open(newline="") as handle:
        reader = csv.DictReader(handle)
        return [_row_to_residue(row) for row in reader]


def split_by_disorder(residues: list[Residue]) -> tuple[list[Residue], list[Residue]]:
    """Return (inside_disorder, outside_disorder)."""
    inside = [r for r in residues if r.is_disorder]
    outside = [r for r in residues if not r.is_disorder]
    return inside, outside


def unique_proteins(residues: list[Residue]) -> set[str]:
    return {r.acc for r in residues}


def evidence_code_counts(residues: list[Residue]) -> dict[str, int]:
    """Count disorder-labeled residues per DisProt evidence code (ec_id)."""
    counts: dict[str, int] = {}
    for r in residues:
        if r.is_disorder and r.ec_id:
            counts[r.ec_id] = counts.get(r.ec_id, 0) + 1
    return counts
