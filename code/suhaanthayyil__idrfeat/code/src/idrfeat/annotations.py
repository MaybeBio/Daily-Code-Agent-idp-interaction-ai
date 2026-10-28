"""Annotation-database loaders and per-segment annotation aggregation.

Loaders parse each raw database once into accession-keyed indexes; the segment functions then
summarize how an IDR segment overlaps motifs, PTM sites, and curated phase-separation regions.
Every loader is optional. When a database is not provided the caller simply skips it and the
corresponding feature columns are filled with zero counts and false flags.
"""

from __future__ import annotations

import csv
import gzip
import json
import tarfile
from pathlib import Path

import numpy as np

PTM_BUCKETS = ("phospho", "acetyl", "ubiq", "methyl")
_BINDING_ELM_TYPES = frozenset({"LIG", "DOC"})


def _ptm_bucket(ptm_type: str) -> str:
    t = ptm_type.lower()
    if "phospho" in t:
        return "phospho"
    if "acetyl" in t:
        return "acetyl"
    if "ubiquit" in t:
        return "ubiq"
    if "methyl" in t:
        return "methyl"
    return "other"


def load_elm(path: Path) -> dict[str, list[tuple[int, int, str, str]]]:
    """Map accession -> list of (start, end, elmtype, elm_identifier) for true positives."""
    index: dict[str, list[tuple[int, int, str, str]]] = {}
    with open(path, newline="") as fh:
        rows = (line for line in fh if not line.startswith("#"))
        reader = csv.DictReader(rows, delimiter="\t", quotechar='"')
        for row in reader:
            if row.get("InstanceLogic") != "true positive":
                continue
            acc = row["Primary_Acc"]
            index.setdefault(acc, []).append(
                (int(row["Start"]), int(row["End"]), row["ELMType"], row["ELMIdentifier"])
            )
    return index


def _dbptm_lines(path: Path):
    if tarfile.is_tarfile(path):
        with tarfile.open(path, "r:gz") as tf:
            fh = tf.extractfile(tf.getmembers()[0])
            for raw in fh:
                yield raw.decode("utf-8", "replace")
        return
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        yield from fh


def load_dbptm(paths: list[Path]) -> dict[str, list[tuple[int, str]]]:
    """Map accession -> list of (position, ptm_bucket) across the dbPTM files."""
    index: dict[str, list[tuple[int, str]]] = {}
    for path in paths:
        for line in _dbptm_lines(path):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                continue
            acc, pos, ptm_type = fields[1], fields[2], fields[3]
            try:
                position = int(pos)
            except ValueError:
                continue
            index.setdefault(acc, []).append((position, _ptm_bucket(ptm_type)))
    return index


def load_phasepro(path: Path) -> dict[str, list[tuple[int, int]]]:
    """Map accession -> list of (start, end) LLPS regions for human PhaSePro entries."""
    with open(path) as fh:
        data = json.load(fh)
    index: dict[str, list[tuple[int, int]]] = {}
    for entry in data.values():
        if entry.get("organism") != "Homo sapiens":
            continue
        acc = entry["accession"]
        regions = []
        for chunk in entry.get("boundaries", "").split(";"):
            chunk = chunk.strip()
            if "-" not in chunk:
                continue
            start, end = chunk.split("-")
            regions.append((int(start), int(end)))
        if regions:
            index[acc] = regions
    return index


def load_drllps(path: Path) -> set[str]:
    """Human Scaffold and Regulator accessions from the DrLLPS protein table."""
    accs: set[str] = set()
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if row.get("Species") == "Homo sapiens" and row.get("LLPS Type") in {
                "Scaffold",
                "Regulator",
            }:
                accs.add(row["UniProt ID"])
    return accs


def _overlaps(a_start: int, a_end: int, b_start: int, b_end: int) -> bool:
    return a_start <= b_end and b_start <= a_end


def elm_segment_features(
    acc: str,
    start: int,
    end: int,
    elm_index: dict[str, list[tuple[int, int, str, str]]],
    seg_length: int,
) -> dict[str, float | int | bool]:
    """Motif residue coverage, density, distinct classes, and binding-motif overlap."""
    instances = elm_index.get(acc, [])
    covered: set[int] = set()
    classes: set[str] = set()
    binding = False
    for m_start, m_end, elmtype, identifier in instances:
        if not _overlaps(start, end, m_start, m_end):
            continue
        classes.add(identifier)
        if elmtype in _BINDING_ELM_TYPES:
            binding = True
        covered.update(range(max(start, m_start), min(end, m_end) + 1))
    residues = len(covered)
    return {
        "elm_motif_residues": residues,
        "elm_motif_density": residues / seg_length if seg_length else 0.0,
        "elm_classes_distinct": len(classes),
        "overlaps_binding_motif": binding,
    }


def ptm_segment_features(
    acc: str,
    start: int,
    end: int,
    ptm_index: dict[str, list[tuple[int, str]]],
    seg_length: int,
) -> dict[str, float | int]:
    """Total and per-type PTM site counts and densities within the segment."""
    sites = ptm_index.get(acc, [])
    by_bucket = dict.fromkeys(PTM_BUCKETS, 0)
    total = 0
    for pos, bucket in sites:
        if start <= pos <= end:
            total += 1
            if bucket in by_bucket:
                by_bucket[bucket] += 1
    feats: dict[str, float | int] = {
        "ptm_count": total,
        "ptm_density": total / seg_length if seg_length else 0.0,
    }
    for bucket in PTM_BUCKETS:
        feats[f"ptm_{bucket}_count"] = by_bucket[bucket]
        feats[f"ptm_{bucket}_density"] = by_bucket[bucket] / seg_length if seg_length else 0.0
    return feats


def phasepro_overlap(
    acc: str, start: int, end: int, index: dict[str, list[tuple[int, int]]]
) -> bool:
    """Whether the segment overlaps any curated PhaSePro LLPS region for the protein."""
    return any(_overlaps(start, end, r_start, r_end) for r_start, r_end in index.get(acc, []))


def drllps_member(acc: str, accs: set[str]) -> bool:
    """Whether the protein is a DrLLPS scaffold or regulator."""
    return acc in accs


def load_plddt_pdb(path: Path) -> np.ndarray:
    """Per-residue pLDDT from an AlphaFold PDB, read from CA atom B-factors in order."""
    values: list[float] = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                values.append(float(line[60:66]))
    return np.asarray(values, dtype=float)


def load_plddt_dir(directory: Path, accessions: list[str]) -> dict[str, np.ndarray]:
    """Load `{accession}.pdb` AlphaFold models from a directory, skipping any that are absent."""
    out: dict[str, np.ndarray] = {}
    for acc in accessions:
        pdb = Path(directory) / f"{acc}.pdb"
        if pdb.exists():
            out[acc] = load_plddt_pdb(pdb)
    return out
