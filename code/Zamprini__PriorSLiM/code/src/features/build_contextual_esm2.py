"""Encode full-sequence IDR coverage crops with frozen layer-33 ESM2."""
                      

from __future__ import annotations
import os
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))

import argparse
import hashlib
import os
from pathlib import Path
import sys
from typing import Iterable

import numpy as np
import pandas as pd
import torch

from pssm_conservation_idppi import find_pssm_file, load_pssm_with_extras
from chemical_window_idppi_fast import load_regions, parse_fasta

MODEL_NAME = "esm2_t33_650M_UR50D"
REPR_LAYER = 33
EMBEDDING_DIM = 1280
MAX_RESIDUES = 1022
WINDOW_LENGTH = 29
WINDOW_RADIUS = WINDOW_LENGTH // 2
DEFAULT_CONTEXT_MARGIN = 128

SCRIPT_PATH = Path(__file__).resolve()
PROJECT_DIR = Path(os.environ.get("PRIORSLIM_DATA_ROOT", SCRIPT_PATH.parents[2] / "data"))
DEFAULT_METADATA = PROJECT_DIR / "feature_cache/contextual_esm2_sequence_aware_metadata.tsv"
DEFAULT_CHUNK_METADATA = PROJECT_DIR / "feature_cache/contextual_esm2_sequence_aware_chunks.tsv"
DEFAULT_OUTPUT_DIR = PROJECT_DIR / "feature_cache/contextual_esm2"
DEFAULT_PSSM_DIR = PROJECT_DIR / "sample_pssm50"
DEFAULT_REGIONS = PROJECT_DIR / "regions.xlsx"
IDR_COVERAGE_METADATA = PROJECT_DIR / "feature_cache/contextual_esm2_idr_coverage_metadata.tsv"
IDR_COVERAGE_CHUNK_METADATA = PROJECT_DIR / "feature_cache/contextual_esm2_idr_coverage_chunks.tsv"
IDR_COVERAGE_OUTPUT_DIR = PROJECT_DIR / "feature_cache/contextual_esm2_idr_coverage"
DEFAULT_SPLITS = tuple(
    PROJECT_DIR / f"new_data/{split}.tsv" for split in ("train", "val", "test")
)

def base_accession(value: object) -> str:

                                                                        
    import re

    match = re.match(r"^(?P<base>.+)_\d+_\d+_[VnMFC]\d+$", str(value))
    return match.group("base") if match else str(value)

def load_default_windows(split_files: Iterable[Path] = DEFAULT_SPLITS) -> pd.DataFrame:
    frames = []
    for split, path in zip(("train", "val", "test"), split_files):
        frame = pd.read_csv(path, sep="\t")
        required = {"uniprot", "window_start", "window_end", "sequence"}
        missing = required.difference(frame.columns)
        if missing:
            raise ValueError(f"{path} is missing columns: {sorted(missing)}")
        selected = frame.copy()
        selected["split"] = split
        selected["split_row"] = np.arange(len(selected), dtype=np.int64)
        selected["window_id"] = split + ":" + selected["split_row"].astype(str)
        frames.append(selected)
    return pd.concat(frames, ignore_index=True)

def normalize_windows(windows: pd.DataFrame) -> pd.DataFrame:
    required = {"uniprot", "window_start", "window_end"}
    missing = required.difference(windows.columns)
    if missing:
        raise ValueError(f"Window table is missing columns: {sorted(missing)}")
    work = windows.copy().reset_index(drop=True)
    if "window_id" not in work:
        work["window_id"] = [f"window:{index}" for index in range(len(work))]
    if "sequence" not in work:
        work["sequence"] = ""
    work["window_index"] = np.arange(len(work), dtype=np.int64)
    work["source_protein_id"] = work["uniprot"].astype(str)
    work["base_protein_id"] = work["source_protein_id"].map(base_accession)

                                                                     
    work["protein_id"] = work["source_protein_id"]
    work["window_start"] = pd.to_numeric(work["window_start"], errors="raise").astype(int)
    work["window_end"] = pd.to_numeric(work["window_end"], errors="raise").astype(int)
    lengths = work["window_end"] - work["window_start"] + 1
    if not lengths.eq(WINDOW_LENGTH).all():
        bad = work.loc[~lengths.eq(WINDOW_LENGTH), ["window_id", "window_start", "window_end"]]
        raise ValueError(f"All windows must span 29 positions; examples: {bad.head().to_dict('records')}")
    if work["window_id"].astype(str).duplicated().any():
        raise ValueError("window_id must be unique")
    return work

def load_contextual_sequences(
    windows: pd.DataFrame, pssm_dir: Path
) -> tuple[dict[str, str], dict[str, str], dict[str, int]]:
    pass
    base_sequences: dict[str, str] = {}
    for base_id in sorted(set(windows["base_protein_id"].astype(str))):
        pssm_path = find_pssm_file(pssm_dir, base_id)
        if pssm_path is None:
            raise FileNotFoundError(f"No PSSM/natural sequence found for {base_id}")
        _matrix, sequence, _extras = load_pssm_with_extras(pssm_path)
        sequence = str(sequence).strip().upper()
        if not sequence:
            raise ValueError(f"PSSM has no sequence: {pssm_path}")
        base_sequences[base_id] = sequence

    sequences: dict[str, str] = {}
    base_by_sequence: dict[str, str] = {}
    modified_counts: dict[str, int] = {}
    for sequence_id, group in windows.groupby("protein_id", sort=True):
        base_ids = group["base_protein_id"].astype(str).unique()
        if len(base_ids) != 1:
            raise ValueError(f"Sequence entity {sequence_id} maps to multiple base proteins")
        base_id = str(base_ids[0])
        natural = base_sequences[base_id]
        reconstructed = list(natural)
        assigned: dict[int, str] = {}
        for row in group.itertuples(index=False):
            window_sequence = str(row.sequence).strip().upper()
            if len(window_sequence) != WINDOW_LENGTH:
                raise ValueError(
                    f"{row.window_id}: expected a 29-aa sequence, got {len(window_sequence)}"
                )
            for offset, amino_acid in enumerate(window_sequence):
                position = int(row.window_start) + offset
                if not 1 <= position <= len(natural):
                    if amino_acid != "X":
                        raise ValueError(f"{row.window_id}: non-X residue outside protein bounds")
                    continue
                if sequence_id == base_id and amino_acid != natural[position - 1]:
                    raise ValueError(
                        f"Canonical window {row.window_id} disagrees with PSSM at {position}: "
                        f"{amino_acid} != {natural[position - 1]}"
                    )
                previous = assigned.get(position)
                if previous is not None and previous != amino_acid:
                    raise ValueError(
                        f"Conflicting reconstructed residues for {sequence_id} at {position}"
                    )
                assigned[position] = amino_acid
                reconstructed[position - 1] = amino_acid
        sequence = "".join(reconstructed)
        sequences[str(sequence_id)] = sequence
        base_by_sequence[str(sequence_id)] = base_id
        modified_counts[str(sequence_id)] = sum(a != b for a, b in zip(sequence, natural))
        if sequence_id != base_id and modified_counts[str(sequence_id)] == 0:
            raise ValueError(f"Synthetic sequence {sequence_id} contains no reconstructed change")
    return sequences, base_by_sequence, modified_counts

def load_idr_intervals(path: Path) -> dict[str, list[tuple[int, int]]]:
    pass
    idrs, _ordered = load_regions(str(path))
    result: dict[str, list[tuple[int, int]]] = {}
    for region in idrs:
        result.setdefault(str(region.protein_id), []).append((int(region.start), int(region.end)))
    return {protein: sorted(set(intervals)) for protein, intervals in result.items()}

def make_idr_scan_windows(
    sequences: dict[str, str], idrs_by_base: dict[str, list[tuple[int, int]]]
) -> pd.DataFrame:
    pass
    records = []
    for protein_id, sequence in sequences.items():
        intervals = idrs_by_base.get(protein_id)
        if not intervals:
            raise ValueError(f"No IDR annotation for FASTA sequence {protein_id}")
        seen_centers: set[int] = set()
        for idr_start, idr_end in intervals:
            if not 1 <= idr_start <= idr_end <= len(sequence):
                raise ValueError(f"Invalid IDR {protein_id}:{idr_start}-{idr_end} for length {len(sequence)}")
            for center in range(idr_start, idr_end + 1):
                if center in seen_centers:
                    raise ValueError(f"Overlapping IDR intervals for {protein_id} at {center}")
                seen_centers.add(center)
                window_start = center - WINDOW_RADIUS
                window_end = center + WINDOW_RADIUS
                records.append({
                    "uniprot": protein_id,
                    "window_id": f"scan:{protein_id}:{center}",
                    "window_start": window_start,
                    "window_end": window_end,
                    "window_center": center,
                    "sequence": "".join(
                        sequence[position - 1] if 1 <= position <= len(sequence) else "X"
                        for position in range(window_start, window_end + 1)
                    ),
                    "idr_start": idr_start,
                    "idr_end": idr_end,
                })
    if not records:
        raise ValueError("FASTA and region file produced no IDR scoring windows")
    return pd.DataFrame.from_records(records)

def assign_window_idrs(
    windows: pd.DataFrame,
    sequences: dict[str, str],
    base_by_sequence: dict[str, str],
    idrs_by_base: dict[str, list[tuple[int, int]]],
) -> pd.DataFrame:
    pass
    work = windows.copy()
    selected = []
    for row in work.itertuples(index=False):
        sequence_id = str(row.protein_id)
        base_id = base_by_sequence[sequence_id]
        center = (int(row.window_start) + int(row.window_end)) // 2
        length = len(sequences[sequence_id])
        matches = [
            (start, end) for start, end in idrs_by_base.get(base_id, [])
            if start <= center <= end
        ]
        if len(matches) != 1:
            raise ValueError(
                f"{row.window_id}: center {center} has {len(matches)} IDR matches for {base_id}"
            )
        idr_start, idr_end = matches[0]
        if idr_end > length:
            raise ValueError(f"{row.window_id}: IDR end {idr_end} exceeds sequence length {length}")
        if "idr_start" in work and "idr_end" in work:
            stated_start, stated_end = getattr(row, "idr_start"), getattr(row, "idr_end")
            if pd.notna(stated_start) and pd.notna(stated_end):
                stated = (int(stated_start), int(stated_end))
                if stated != matches[0]:
                    allow_value = getattr(row, "allow_idr_region_mismatch", False)
                    allow_mismatch = pd.notna(allow_value) and bool(allow_value)
                    overlaps = max(stated[0], idr_start) <= min(stated[1], idr_end)
                    if not allow_mismatch or not overlaps or not stated[0] <= center <= stated[1]:
                        raise ValueError(
                            f"{row.window_id}: window IDR {stated_start}-{stated_end} "
                            f"disagrees with region file {idr_start}-{idr_end}"
                        )
        selected.append((idr_start, idr_end))
    work["idr_start_used"] = [value[0] for value in selected]
    work["idr_end_used"] = [value[1] for value in selected]
    return work

def centered_crop(group_start: int, group_end: int, protein_length: int) -> tuple[int, int]:
    if protein_length <= MAX_RESIDUES:
        return 1, protein_length

                                       
    crop_start = (group_start + group_end - MAX_RESIDUES + 1) // 2
    crop_start = max(1, min(crop_start, protein_length - MAX_RESIDUES + 1))
    return crop_start, crop_start + MAX_RESIDUES - 1

def _greedy_groups(frame: pd.DataFrame, protein_length: int) -> list[pd.DataFrame]:
    ordered = frame.sort_values(
        ["window_start", "window_end", "window_index"], kind="mergesort"
    )
    if protein_length <= MAX_RESIDUES:
        return [ordered]
    groups: list[list[int]] = []
    current: list[int] = []
    current_left = current_right = 0
    for index, row in ordered.iterrows():
        left, right = int(row.window_start), int(row.window_end)
        if not current:
            current, current_left, current_right = [index], left, right
            continue
        proposed_left = min(current_left, left)
        proposed_right = max(current_right, right)
        if proposed_right - proposed_left + 1 <= MAX_RESIDUES:
            current.append(index)
            current_left, current_right = proposed_left, proposed_right
        else:
            groups.append(current)
            current, current_left, current_right = [index], left, right
    if current:
        groups.append(current)
    return [ordered.loc[indices].copy() for indices in groups]

def _planned_groups(
    windows: pd.DataFrame,
    sequences: dict[str, str],
    crop_mode: str,
    context_margin: int,
):
    pass
    core_length = MAX_RESIDUES - 2 * (WINDOW_RADIUS + context_margin)
    if crop_mode == "idr_coverage" and (context_margin < 0 or core_length <= 0):
        raise ValueError("context_margin must be nonnegative and leave core center positions")
    for protein_id, protein_windows in windows.groupby("protein_id", sort=True):
        protein_length = len(sequences[str(protein_id)])
        if crop_mode == "window":
            for chunk_id, group in enumerate(_greedy_groups(protein_windows, protein_length)):
                group_start = int(group["window_start"].min())
                group_end = int(group["window_end"].max())
                crop_start, crop_end = centered_crop(group_start, group_end, protein_length)
                yield protein_id, chunk_id, group, crop_start, crop_end, {}
            continue
        if protein_length <= MAX_RESIDUES:

            yield protein_id, 0, protein_windows, 1, protein_length, {
                "crop_mode": crop_mode, "context_margin": context_margin,
                "core_start": 1, "core_end": protein_length,
            }
            continue
        chunk_id = 0
        for (idr_start, idr_end), idr_windows in protein_windows.groupby(
            ["idr_start_used", "idr_end_used"], sort=True
        ):
            idr_start, idr_end = int(idr_start), int(idr_end)
            centers = (idr_windows["window_start"] + idr_windows["window_end"]) // 2
            for core_start in range(idr_start, idr_end + 1, core_length):
                core_end = min(idr_end, core_start + core_length - 1)
                group = idr_windows.loc[centers.between(core_start, core_end)].copy()
                if group.empty:
                    continue                                                                
                crop_start, crop_end = centered_crop(core_start, core_end, protein_length)
                yield protein_id, chunk_id, group, crop_start, crop_end, {
                    "crop_mode": crop_mode, "context_margin": context_margin,
                    "idr_start": idr_start, "idr_end": idr_end,
                    "core_start": core_start, "core_end": core_end,
                }
                chunk_id += 1

def build_contextual_plan(
    windows: pd.DataFrame,
    sequences: dict[str, str],
    base_by_sequence: dict[str, str],
    modified_counts: dict[str, int],
    output_dir: Path,
    crop_mode: str = "window",
    idrs_by_base: dict[str, list[tuple[int, int]]] | None = None,
    context_margin: int = DEFAULT_CONTEXT_MARGIN,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    if crop_mode not in {"window", "idr_coverage"}:
        raise ValueError(f"Unknown crop mode: {crop_mode}")
    work = normalize_windows(windows)
    if crop_mode == "idr_coverage":
        if idrs_by_base is None:
            raise ValueError("idr_coverage requires region-file IDR intervals")
        work = assign_window_idrs(work, sequences, base_by_sequence, idrs_by_base)
    window_records: list[dict] = []
    chunk_records: list[dict] = []
    for protein_id, chunk_id, group, crop_start, crop_end, plan_info in _planned_groups(
        work, sequences, crop_mode, context_margin
    ):
        sequence = sequences[str(protein_id)]
        base_protein_id = base_by_sequence[str(protein_id)]
        is_modified = int(modified_counts[str(protein_id)] > 0)
        protein_length = len(sequence)
        if crop_mode == "window" and protein_length > MAX_RESIDUES:
            group_start = int(group["window_start"].min())
            group_end = int(group["window_end"].max())
            if group_end - group_start + 1 > MAX_RESIDUES:
                raise RuntimeError("Greedy group exceeds the ESM crop limit")
        if not 1 <= crop_start <= crop_end <= protein_length or crop_end - crop_start + 1 > MAX_RESIDUES:
            raise RuntimeError("Invalid contextual crop bounds")
        chunk_name = f"{protein_id}__chunk{chunk_id:04d}.pt"
        chunk_path = (output_dir / chunk_name).resolve()
        chunk_sequence = sequence[crop_start - 1 : crop_end]
        sequence_sha256 = hashlib.sha256(chunk_sequence.encode("ascii")).hexdigest()
        needed_crop_rows: set[int] = set()
        for row in group.itertuples(index=False):
            ws, we = int(row.window_start), int(row.window_end)
            valid_start = max(ws, 1)
            valid_end = min(we, protein_length)
            if valid_start > valid_end:
                raise ValueError(f"Window {row.window_id} has no natural-protein residues")
            if valid_start < crop_start or valid_end > crop_end:
                raise RuntimeError(
                    f"Chunk {protein_id}:{chunk_id} does not cover natural part of {row.window_id}"
                )
            left_padding = max(0, 1 - ws)
            right_padding = max(0, we - protein_length)
            valid_length = valid_end - valid_start + 1
            needed_crop_rows.update(range(valid_start - crop_start, valid_end - crop_start + 1))
            if left_padding + valid_length + right_padding != WINDOW_LENGTH:
                raise RuntimeError(f"Window mapping is not 29 residues: {row.window_id}")
            row_info = {}
            if crop_mode == "idr_coverage":
                row_info = {
                    "crop_mode": crop_mode, "context_margin": context_margin,
                    "idr_start": int(row.idr_start_used), "idr_end": int(row.idr_end_used),
                    "core_start": plan_info["core_start"], "core_end": plan_info["core_end"],
                    "window_center": (ws + we) // 2,
                }
            window_records.append({
                "protein_id": protein_id,
                "source_protein_id": row.source_protein_id,
                "base_protein_id": base_protein_id,
                "sequence_is_modified": is_modified,
                "modified_residue_count": modified_counts[str(protein_id)],
                "chunk_id": chunk_id,
                "protein_length": protein_length,
                "crop_start": crop_start,
                "crop_end": crop_end,
                "crop_length": crop_end - crop_start + 1,
                "window_index": int(row.window_index),
                "window_id": str(row.window_id),
                "window_start": ws,
                "window_end": we,

                "window_local_start": ws - crop_start + 1,
                "window_local_end": we - crop_start + 1,
                "valid_protein_start": valid_start,
                "valid_protein_end": valid_end,
                "valid_local_start": valid_start - crop_start + 1,
                "valid_local_end": valid_end - crop_start + 1,
                "window_output_valid_start": left_padding,
                "window_output_valid_end": left_padding + valid_length,
                "left_padding": left_padding,
                "right_padding": right_padding,
                "valid_length": valid_length,
                "embedding_file": chunk_name,
                "embedding_path": str(chunk_path),
                "sequence_sha256": sequence_sha256,
                **row_info,
            })
        embedding_storage = "needed_rows" if is_modified else "full_crop"
        stored_residue_count = (
            len(needed_crop_rows) if embedding_storage == "needed_rows"
            else crop_end - crop_start + 1
        )
        chunk_records.append({
            "protein_id": protein_id,
            "base_protein_id": base_protein_id,
            "sequence_is_modified": is_modified,
            "modified_residue_count": modified_counts[str(protein_id)],
            "chunk_id": chunk_id,
            "protein_length": protein_length,
            "crop_start": crop_start,
            "crop_end": crop_end,
            "crop_length": crop_end - crop_start + 1,
            "window_count": len(group),
            "unique_window_coordinate_count": len(group[["window_start", "window_end"]].drop_duplicates()),
            "embedding_storage": embedding_storage,
            "stored_residue_count": stored_residue_count,
            "embedding_file": chunk_name,
            "embedding_path": str(chunk_path),
            "sequence_sha256": sequence_sha256,
            "esm_model": MODEL_NAME,
            "esm_repr_layer": REPR_LAYER,
            "embedding_dim": EMBEDDING_DIM,
            **plan_info,
        })
    window_metadata = pd.DataFrame(window_records).sort_values("window_index").reset_index(drop=True)
    chunk_metadata = pd.DataFrame(chunk_records).sort_values(["protein_id", "chunk_id"]).reset_index(drop=True)
    if len(window_metadata) != len(work):
        raise RuntimeError("Contextual plan lost windows")
    if not np.array_equal(window_metadata["window_index"].to_numpy(), np.arange(len(work))):
        raise RuntimeError("Contextual metadata is not aligned to input window order")
    return window_metadata, chunk_metadata

def choose_device(requested: str) -> torch.device:
    if requested == "auto":
        if torch.cuda.is_available():
            return torch.device("cuda")
        if getattr(torch.backends, "mps", None) is not None and torch.backends.mps.is_available():
            return torch.device("mps")
        return torch.device("cpu")
    device = torch.device(requested)
    if requested == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("CUDA requested but unavailable")
    if requested == "mps" and not torch.backends.mps.is_available():
        raise RuntimeError("MPS requested but unavailable")
    return device

def _existing_chunk_valid(path: Path, record: pd.Series) -> bool:
    if not path.is_file():
        return False
    try:
        payload = torch.load(path, map_location="cpu", weights_only=True)
        embedding = payload["embedding"]
        storage = str(record.embedding_storage)
        rows = payload.get("embedding_row_indices")
        row_storage_valid = storage == "full_crop" or (
            isinstance(rows, torch.Tensor)
            and rows.dtype == torch.int64
            and rows.ndim == 1
            and len(rows) == int(record.stored_residue_count)
            and (len(rows) == 0 or (
                int(rows.min()) >= 0
                and int(rows.max()) < int(record.crop_length)
                and bool(torch.all(rows[1:] > rows[:-1]))
            ))
        )
        return (
            isinstance(embedding, torch.Tensor)
            and tuple(embedding.shape) == (int(record.stored_residue_count), EMBEDDING_DIM)
            and embedding.dtype == torch.float32
            and payload.get("protein_id") == str(record.protein_id)
            and int(payload.get("crop_start")) == int(record.crop_start)
            and int(payload.get("crop_end")) == int(record.crop_end)
            and payload.get("sequence_sha256") == str(record.sequence_sha256)
            and payload.get("esm_model") == MODEL_NAME
            and int(payload.get("esm_repr_layer")) == REPR_LAYER
            and payload.get("embedding_storage", "full_crop") == storage
            and row_storage_valid
        )
    except (OSError, KeyError, TypeError, ValueError):
        return False

def generate_contextual_embeddings(
    chunk_metadata: pd.DataFrame,
    window_metadata: pd.DataFrame,
    sequences: dict[str, str],
    output_dir: Path,
    device_name: str = "auto",
    batch_size: int = 1,
    force: bool = False,
) -> tuple[int, int, str]:
    if batch_size <= 0:
        raise ValueError("batch_size must be positive")
    output_dir.mkdir(parents=True, exist_ok=True)
    pending = []
    reused = 0
    for index, row in chunk_metadata.iterrows():
        path = output_dir / str(row.embedding_file)
        if not force and _existing_chunk_valid(path, row):
            reused += 1
        else:
            pending.append(index)
    device = choose_device(device_name)
    if not pending:
        return 0, reused, str(device)

    import esm

    model, alphabet = getattr(esm.pretrained, MODEL_NAME)()
    model = model.to(device).eval()
    batch_converter = alphabet.get_batch_converter()
    pending.sort(key=lambda index: (int(chunk_metadata.at[index, "crop_length"]), index))
    generated = 0
    windows_by_chunk = {
        (str(protein_id), int(chunk_id)): group
        for (protein_id, chunk_id), group in window_metadata.groupby(
            ["protein_id", "chunk_id"], sort=False
        )
    }
    with torch.inference_mode():
        for batch_offset in range(0, len(pending), batch_size):
            indices = pending[batch_offset : batch_offset + batch_size]
            batch = []
            for index in indices:
                row = chunk_metadata.loc[index]
                sequence = sequences[str(row.protein_id)][int(row.crop_start) - 1 : int(row.crop_end)]
                if len(sequence) != int(row.crop_length) or len(sequence) > MAX_RESIDUES:
                    raise RuntimeError("Invalid contextual ESM crop")
                batch.append((f"{row.protein_id}|chunk{int(row.chunk_id)}", sequence))
            _labels, _strings, tokens = batch_converter(batch)
            result = model(tokens.to(device), repr_layers=[REPR_LAYER], return_contacts=False)
            representations = result["representations"][REPR_LAYER]
            for batch_index, index in enumerate(indices):
                row = chunk_metadata.loc[index]
                length = int(row.crop_length)
                full_embedding = representations[batch_index, 1 : length + 1]
                storage = str(row.embedding_storage)
                embedding_row_indices = None
                if storage == "needed_rows":
                    group = windows_by_chunk[(str(row.protein_id), int(row.chunk_id))]
                    needed: set[int] = set()
                    for window in group.itertuples(index=False):
                        needed.update(
                            range(int(window.valid_local_start) - 1, int(window.valid_local_end))
                        )
                    embedding_row_indices = torch.tensor(
                        sorted(needed), dtype=torch.int64
                    )
                    embedding = full_embedding.index_select(
                        0, embedding_row_indices.to(full_embedding.device)
                    ).to(device="cpu", dtype=torch.float32).contiguous()
                elif storage == "full_crop":
                    embedding = full_embedding.to(
                        device="cpu", dtype=torch.float32
                    ).contiguous()
                else:
                    raise ValueError(f"Unknown embedding storage mode: {storage}")
                if len(embedding) != int(row.stored_residue_count):
                    raise RuntimeError("Stored contextual row count mismatch")
                payload = {
                    "embedding": embedding,
                    "protein_id": str(row.protein_id),
                    "base_protein_id": str(row.base_protein_id),
                    "sequence_is_modified": int(row.sequence_is_modified),
                    "modified_residue_count": int(row.modified_residue_count),
                    "chunk_id": int(row.chunk_id),
                    "protein_length": int(row.protein_length),
                    "crop_start": int(row.crop_start),
                    "crop_end": int(row.crop_end),
                    "sequence_sha256": str(row.sequence_sha256),
                    "esm_model": MODEL_NAME,
                    "esm_repr_layer": REPR_LAYER,
                    "embedding_dim": EMBEDDING_DIM,
                    "embedding_storage": storage,
                    "embedding_row_indices": embedding_row_indices,
                    "normalized": False,
                }
                path = output_dir / str(row.embedding_file)
                temp = path.with_name(f".{path.name}.{os.getpid()}.tmp")
                torch.save(payload, temp)
                os.replace(temp, path)
                generated += 1
            del result, representations, tokens
            if device.type == "mps" and generated % 250 == 0:
                torch.mps.synchronize()
                torch.mps.empty_cache()
            if generated % 25 == 0 or generated == len(pending):
                print(
                    f"[CONTEXTUAL ESM2] {generated}/{len(pending)} generated; reused={reused}",
                    flush=True,
                )
    return generated, reused, str(device)

def prepare_contextual_esm2(
    windows: pd.DataFrame,
    pssm_dir: Path = DEFAULT_PSSM_DIR,
    output_dir: Path = DEFAULT_OUTPUT_DIR,
    metadata_path: Path = DEFAULT_METADATA,
    chunk_metadata_path: Path = DEFAULT_CHUNK_METADATA,
    device: str = "auto",
    batch_size: int = 1,
    force: bool = False,
    plan_only: bool = False,
    crop_mode: str = "window",
    regions_path: Path = DEFAULT_REGIONS,
    context_margin: int = DEFAULT_CONTEXT_MARGIN,
    sequence_overrides: dict[str, str] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    normalized = normalize_windows(windows)
    if sequence_overrides is None:
        sequences, base_by_sequence, modified_counts = load_contextual_sequences(
            normalized, Path(pssm_dir)
        )
    else:
        sequences = {str(key): str(value).upper() for key, value in sequence_overrides.items()}
        expected_ids = set(normalized["protein_id"].astype(str))
        if expected_ids != set(sequences):
            raise ValueError("FASTA sequence IDs and scoring-window sequence IDs differ")
        base_by_sequence = {key: key for key in sequences}
        modified_counts = {key: 0 for key in sequences}
        for row in normalized.itertuples(index=False):
            actual = "".join(
                sequences[str(row.protein_id)][position - 1]
                if 1 <= position <= len(sequences[str(row.protein_id)]) else "X"
                for position in range(int(row.window_start), int(row.window_end) + 1)
            )
            if str(row.sequence).upper() != actual:
                raise ValueError(f"{row.window_id}: window sequence disagrees with FASTA")
    idrs_by_base = load_idr_intervals(Path(regions_path)) if crop_mode == "idr_coverage" else None
    window_metadata, chunk_metadata = build_contextual_plan(
        normalized, sequences, base_by_sequence, modified_counts, Path(output_dir),
        crop_mode=crop_mode, idrs_by_base=idrs_by_base, context_margin=context_margin,
    )
    Path(metadata_path).parent.mkdir(parents=True, exist_ok=True)
    window_metadata.to_csv(metadata_path, sep="\t", index=False)
    chunk_metadata.to_csv(chunk_metadata_path, sep="\t", index=False)
    if not plan_only:
        generated, reused, selected_device = generate_contextual_embeddings(
            chunk_metadata,
            window_metadata,
            sequences,
            Path(output_dir),
            device_name=device,
            batch_size=batch_size,
            force=force,
        )
        print(
            f"[OK] sequence-aware contextual ESM2 sequences={len(sequences)} "
            f"modified={sum(count > 0 for count in modified_counts.values())} "
            f"chunks={len(chunk_metadata)} "
            f"windows={len(window_metadata)} generated={generated} reused={reused} "
            f"device={selected_device}",
            flush=True,
        )
    else:
        print(
            f"[OK PLAN] sequences={len(sequences)} "
            f"modified={sum(count > 0 for count in modified_counts.values())} "
            f"chunks={len(chunk_metadata)} "
            f"windows={len(window_metadata)} mean_windows_per_chunk="
            f"{len(window_metadata) / max(len(chunk_metadata), 1):.3f}",
            flush=True,
        )
    return window_metadata, chunk_metadata

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--windows", type=Path, default=None)
    parser.add_argument("--fasta", type=Path, default=None,
                        help="Generate every IDR-centered 29-aa window from FASTA (idr_coverage only)")
    parser.add_argument("--crop-mode", choices=("window", "idr_coverage"), default="idr_coverage")
    parser.add_argument("--regions", type=Path, default=DEFAULT_REGIONS)
    parser.add_argument("--context-margin", type=int, default=DEFAULT_CONTEXT_MARGIN,
                        help="Extra ESM residues beyond each 29-aa window at core boundaries")
    parser.add_argument("--pssm-dir", type=Path, default=DEFAULT_PSSM_DIR)
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--metadata-output", type=Path, default=None)
    parser.add_argument("--chunk-metadata-output", type=Path, default=None)
    parser.add_argument("--device", choices=("auto", "cpu", "mps", "cuda"), default="auto")
    parser.add_argument("--batch-size", type=int, default=1)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--plan-only", action="store_true")
    args = parser.parse_args()
    if args.fasta is not None and (args.windows is not None or args.crop_mode != "idr_coverage"):
        parser.error("--fasta requires --crop-mode idr_coverage and cannot be combined with --windows")
    if args.context_margin < 0 or 2 * (WINDOW_RADIUS + args.context_margin) >= MAX_RESIDUES:
        parser.error("--context-margin leaves no valid core centers")
    return args

def main() -> int:
    args = parse_args()
    sequence_overrides = None
    if args.fasta is not None:
        sequence_overrides = parse_fasta(str(args.fasta))
        windows = make_idr_scan_windows(sequence_overrides, load_idr_intervals(args.regions))
    elif args.windows is None:
        windows = load_default_windows()
    else:
        windows = pd.read_csv(args.windows, sep="\t")
    idr_mode = args.crop_mode == "idr_coverage"
    if args.fasta is not None:

        scan_root = PROJECT_DIR / "feature_cache" / "sequence_scans" / args.fasta.stem
        default_output_dir = scan_root / "contextual_esm2"
        default_metadata_path = scan_root / "contextual_esm2_metadata.tsv"
        default_chunk_metadata_path = scan_root / "contextual_esm2_chunks.tsv"
    else:
        default_output_dir = IDR_COVERAGE_OUTPUT_DIR if idr_mode else DEFAULT_OUTPUT_DIR
        default_metadata_path = IDR_COVERAGE_METADATA if idr_mode else DEFAULT_METADATA
        default_chunk_metadata_path = (
            IDR_COVERAGE_CHUNK_METADATA if idr_mode else DEFAULT_CHUNK_METADATA
        )
    prepare_contextual_esm2(
        windows=windows,
        pssm_dir=args.pssm_dir.resolve(),
        output_dir=(args.output_dir or default_output_dir).resolve(),
        metadata_path=(args.metadata_output or default_metadata_path).resolve(),
        chunk_metadata_path=(args.chunk_metadata_output or default_chunk_metadata_path).resolve(),
        device=args.device,
        batch_size=args.batch_size,
        force=args.force,
        plan_only=args.plan_only,
        crop_mode=args.crop_mode,
        regions_path=args.regions.resolve(),
        context_margin=args.context_margin,
        sequence_overrides=sequence_overrides,
    )
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
