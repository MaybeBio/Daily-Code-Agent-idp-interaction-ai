"""Filesystem, configuration, and tabular I/O helpers.

Configuration ships as a baked-in default dictionary so the tool runs without a config file;
an optional YAML path deep-merges over the defaults. The bundled example FASTA is packaged data
resolved through importlib.resources, so it is found whether the project is run from a checkout
or installed as a wheel.
"""

from __future__ import annotations

import copy
import gzip
from importlib import resources
from pathlib import Path

import pandas as pd
import yaml

DEFAULT_CONFIG = {
    "seed": 1729,
    "disorder": {
        "primary": "metapredict",  # required disorder predictor that defines IDR segments
        "definitions": ["metapredict", "aiupred"],
        "score_threshold": 0.5,
        "min_segment_length": 10,
        "plddt_threshold": 70.0,
        "binding_threshold": 0.5,
    },
}


def _deep_merge(base: dict, override: dict) -> dict:
    out = copy.deepcopy(base)
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(out.get(key), dict):
            out[key] = _deep_merge(out[key], value)
        else:
            out[key] = value
    return out


def load_config(path: str | Path | None = None) -> dict:
    """Baked-in defaults, deep-merged with an optional YAML override file."""
    if path is None:
        return copy.deepcopy(DEFAULT_CONFIG)
    with open(path) as fh:
        override = yaml.safe_load(fh) or {}
    return _deep_merge(DEFAULT_CONFIG, override)


def _maybe_gzip(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def read_fasta(path: str | Path) -> dict[str, str]:
    """Parse a FASTA into {accession: sequence}.

    Handles `sp|ACC|NAME` and `tr|ACC|NAME` UniProt headers and falls back to the first
    whitespace-delimited token otherwise.
    """
    path = Path(path)
    seqs: dict[str, list[str]] = {}
    acc = None
    with _maybe_gzip(path) as fh:
        for line in fh:
            if line.startswith(">"):
                header = line[1:].strip()
                parts = header.split("|")
                acc = (
                    parts[1] if len(parts) >= 3 and parts[0] in {"sp", "tr"} else header.split()[0]
                )
                seqs[acc] = []
            elif acc is not None:
                seqs[acc].append(line.strip())
    return {a: "".join(chunks) for a, chunks in seqs.items()}


def write_table(df: pd.DataFrame, out: str | Path) -> tuple[Path, Path]:
    """Write the feature table as parquet at ``out`` and CSV alongside it."""
    out = Path(out)
    if out.suffix != ".parquet":
        out = out.with_suffix(".parquet")
    out.parent.mkdir(parents=True, exist_ok=True)
    csv_path = out.with_suffix(".csv")
    df.to_parquet(out, index=False)
    df.to_csv(csv_path, index=False)
    return out, csv_path


def example_fasta_path() -> Path:
    """Path to the bundled synthetic example FASTA."""
    return Path(resources.files("idrfeat").joinpath("data/example.fasta"))
