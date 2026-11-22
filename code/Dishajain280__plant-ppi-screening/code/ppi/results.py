"""Parse AlphaFold Server output into per-pair confidence rows.

A server download contains, per job:
  <name>_summary_confidences.json   aggregate scores
  <name>_model.cif                  coordinates
  <name>_full_data_*.json           per-token detail

For screening we only need the summary: ipTM (interface confidence), pTM,
the minimum chain-pair PAE, and the clash/disorder flags. This module reads
those defensively, because the exact key set has changed between releases.
"""

from __future__ import annotations

import csv
import json
from dataclasses import dataclass, asdict
from pathlib import Path

SUMMARY_GLOB = "*summary_confidences*.json"


@dataclass
class Score:
    pair_id: str
    iptm: float
    ptm: float
    min_chain_pair_pae: float
    ranking_score: float
    fraction_disordered: float
    has_clash: bool
    source_file: str

    def as_row(self) -> dict:
        return asdict(self)


def _first_number(payload: dict, *keys, default: float = float("nan")) -> float:
    for key in keys:
        value = payload.get(key)
        if isinstance(value, (int, float)):
            return float(value)
    return default


def _min_pair_pae(payload: dict) -> float:
    """Smallest PAE between chain 0 and chain 1; low values imply an interface."""
    matrix = payload.get("chain_pair_pae_min")
    if isinstance(matrix, list) and matrix:
        row = matrix[0]
        if isinstance(row, list) and len(row) > 1:
            return float(row[1])
    # Newer exports nest it under chain_pair_iptm/pae objects.
    nested = payload.get("chain_pair_pae", {})
    if isinstance(nested, dict):
        for value in nested.values():
            if isinstance(value, (int, float)):
                return float(value)
    return float("nan")


def parse_summary(path: Path) -> Score:
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    pair_id = _pair_id_from_path(Path(path))
    return Score(
        pair_id=pair_id,
        iptm=_first_number(payload, "iptm", "ipTM", "model_iptm"),
        ptm=_first_number(payload, "ptm", "pTM", "model_ptm"),
        min_chain_pair_pae=_min_pair_pae(payload),
        ranking_score=_first_number(payload, "ranking_score"),
        fraction_disordered=_first_number(payload, "fraction_disordered", default=0.0),
        has_clash=bool(payload.get("has_clash", False)),
        source_file=Path(path).name,
    )


def _pair_id_from_path(path: Path) -> str:
    """Recover the pair id from the folder or file name.

    AlphaFold Server writes results as ``fold_<job name>_summary_confidences.json``,
    so the ``fold_`` prefix has to be removed for the id to match the candidate list.
    """
    for part in (path.parent.name, path.name):
        if not part:
            continue
        stem = part
        for marker in ("_summary_confidences", "_model", "_full_data"):
            stem = stem.split(marker)[0]
        if stem.lower().startswith("fold_"):
            stem = stem[5:]
        if stem:
            return stem
    return "unknown"


def collect(results_dir: Path, out_csv: Path) -> list[Score]:
    """Find every summary file under ``results_dir`` and write one CSV row each."""
    results_dir = Path(results_dir)
    scores = [parse_summary(p) for p in sorted(results_dir.rglob(SUMMARY_GLOB))]
    out_csv = Path(out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(Score.__dataclass_fields__))
        writer.writeheader()
        for score in scores:
            writer.writerow(score.as_row())
    return scores


def read_scores(path: Path) -> list[dict]:
    with Path(path).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))
