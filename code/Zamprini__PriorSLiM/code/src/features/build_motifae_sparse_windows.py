"""Encode contextual ESM2 windows as sparse MotifAE CSR arrays."""
                      

from __future__ import annotations
import os
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))

import argparse
import hashlib
import json
import os
from pathlib import Path
import sys
import time

import numpy as np
import pandas as pd
import torch

from extract_motifae_sparse_f import MotifAEEncoder, choose_device

WINDOW_LENGTH = 29
EXPECTED_INPUT_DIM = 1280
EXPECTED_LATENT_DIM = 40_960

SCRIPT_PATH = Path(__file__).resolve()
PROJECT_DIR = Path(os.environ.get("PRIORSLIM_DATA_ROOT", SCRIPT_PATH.parents[2] / "data"))
DEFAULT_CONTEXT_METADATA = PROJECT_DIR / "feature_cache/contextual_esm2_idr_coverage_metadata.tsv"
DEFAULT_CHUNK_METADATA = PROJECT_DIR / "feature_cache/contextual_esm2_idr_coverage_chunks.tsv"
DEFAULT_OUTPUT_DIR = PROJECT_DIR / "feature_cache/motifae_sparse_windows_idr_coverage"
DEFAULT_CHECKPOINT = Path(os.environ.get(
    "MOTIFAE_CHECKPOINT", PROJECT_DIR / "motifae/MotifAE_k3_step_80000.pt"
)).expanduser()

def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()

def frame_sha256(frame: pd.DataFrame, columns: list[str]) -> str:
    payload = frame[columns].to_csv(index=False, lineterminator="\n").encode("utf-8")
    return hashlib.sha256(payload).hexdigest()

def validate_context_metadata(windows: pd.DataFrame, chunks: pd.DataFrame) -> None:
    window_required = {
        "protein_id", "source_protein_id", "chunk_id", "window_index", "window_id",
        "window_start", "window_end", "valid_local_start", "valid_local_end",
        "window_output_valid_start", "window_output_valid_end", "left_padding",
        "right_padding", "valid_length", "embedding_path",
    }
    chunk_required = {
        "protein_id", "chunk_id", "crop_start", "crop_end", "crop_length",
        "embedding_path", "sequence_sha256", "embedding_storage",
        "stored_residue_count",
    }
    if missing := window_required.difference(windows.columns):
        raise ValueError(f"Context window metadata lacks {sorted(missing)}")
    if missing := chunk_required.difference(chunks.columns):
        raise ValueError(f"Context chunk metadata lacks {sorted(missing)}")
    if windows["window_index"].duplicated().any() or windows["window_id"].duplicated().any():
        raise ValueError("Context window identifiers must be unique")
    lengths = windows["window_end"].astype(int) - windows["window_start"].astype(int) + 1
    if not lengths.eq(WINDOW_LENGTH).all():
        raise ValueError("All contextual windows must span 29 positions")
    mapped = (
        windows["left_padding"].astype(int)
        + windows["valid_length"].astype(int)
        + windows["right_padding"].astype(int)
    )
    if not mapped.eq(WINDOW_LENGTH).all():
        raise ValueError("All contextual windows must map to exactly 29 output residues")
    chunk_keys = set(chunks[["protein_id", "chunk_id"]].itertuples(index=False, name=None))
    window_keys = set(windows[["protein_id", "chunk_id"]].itertuples(index=False, name=None))
    if not window_keys.issubset(chunk_keys):
        raise ValueError("A window references a missing contextual ESM chunk")

def load_chunk(path: Path, expected_row: pd.Series) -> tuple[torch.Tensor, np.ndarray]:
    if not path.is_file():
        raise FileNotFoundError(f"Contextual ESM chunk is missing: {path}")
    payload = torch.load(path, map_location="cpu", weights_only=True)
    embedding = payload.get("embedding")
    if not isinstance(embedding, torch.Tensor):
        raise ValueError(f"Contextual file has no embedding tensor: {path}")
    expected_shape = (int(expected_row.stored_residue_count), EXPECTED_INPUT_DIM)
    if tuple(embedding.shape) != expected_shape or embedding.dtype != torch.float32:
        raise ValueError(f"{path} has {tuple(embedding.shape)}/{embedding.dtype}, expected {expected_shape}/float32")
    if str(payload.get("protein_id")) != str(expected_row.protein_id):
        raise ValueError(f"Protein mismatch in {path}")
    if int(payload.get("chunk_id")) != int(expected_row.chunk_id):
        raise ValueError(f"Chunk mismatch in {path}")
    if payload.get("sequence_sha256") != str(expected_row.sequence_sha256):
        raise ValueError(f"Sequence signature mismatch in {path}")
    storage = str(expected_row.embedding_storage)
    if payload.get("embedding_storage", "full_crop") != storage:
        raise ValueError(f"Embedding storage mismatch in {path}")
    if storage == "full_crop":
        crop_rows = np.arange(int(expected_row.crop_length), dtype=np.int64)
    elif storage == "needed_rows":
        stored_rows = payload.get("embedding_row_indices")
        if not isinstance(stored_rows, torch.Tensor) or stored_rows.dtype != torch.int64:
            raise ValueError(f"Compact contextual file lacks int64 row indices: {path}")
        crop_rows = stored_rows.cpu().numpy()
        if len(crop_rows) != len(embedding) or (
            len(crop_rows) and (
                crop_rows.min() < 0
                or crop_rows.max() >= int(expected_row.crop_length)
                or np.any(crop_rows[1:] <= crop_rows[:-1])
            )
        ):
            raise ValueError(f"Invalid compact contextual row mapping: {path}")
    else:
        raise ValueError(f"Unknown contextual storage mode {storage}: {path}")
    return embedding, crop_rows

def build_unique_windows(windows: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    key_columns = [
        "protein_id", "chunk_id", "window_start", "window_end",
        "valid_local_start", "valid_local_end", "window_output_valid_start",
        "window_output_valid_end", "left_padding", "right_padding", "valid_length",
    ]
    unique = windows.drop_duplicates(key_columns, keep="first").copy()
    unique = unique.sort_values(["protein_id", "chunk_id", "window_start", "window_end"]).reset_index(drop=True)
    unique["sparse_window_index"] = np.arange(len(unique), dtype=np.int64)
    mapping = windows.merge(
        unique[key_columns + ["sparse_window_index"]],
        on=key_columns,
        how="left",
        validate="many_to_one",
        suffixes=("", "_unique"),
    ).sort_values("window_index").reset_index(drop=True)
    if mapping["sparse_window_index"].isna().any():
        raise RuntimeError("Failed to map contextual windows to sparse windows")
    mapping["sparse_window_index"] = mapping["sparse_window_index"].astype(np.int64)
    return unique, mapping

def encode_needed_residues(
    embedding: torch.Tensor,
    needed_rows: np.ndarray,
    encoder: MotifAEEncoder,
    chunk_residues: int,
) -> dict[int, tuple[np.ndarray, np.ndarray]]:
    pass
    result: dict[int, tuple[np.ndarray, np.ndarray]] = {}
    with torch.inference_mode():
        for offset in range(0, len(needed_rows), chunk_residues):
            row_ids = needed_rows[offset : offset + chunk_residues]
            x = embedding[row_ids].to(encoder.device, dtype=torch.float32)
            latent = encoder.encode(x)
            positions = torch.nonzero(latent > 0, as_tuple=False)
            positions_cpu = positions.to("cpu", dtype=torch.int64).numpy()
            values_cpu = latent[positions[:, 0], positions[:, 1]].to(
                "cpu", dtype=torch.float32
            ).numpy()
            counts = np.bincount(
                positions_cpu[:, 0], minlength=len(row_ids)
            ).astype(np.int64, copy=False)
            cursor = 0
            for local_row, crop_row in enumerate(row_ids):
                stop = cursor + int(counts[local_row])
                result[int(crop_row)] = (
                    positions_cpu[cursor:stop, 1].astype(np.int32, copy=False),
                    values_cpu[cursor:stop],
                )
                cursor = stop
            del x, latent, positions
    return result

def build_motifae_sparse_windows(
    context_metadata_path: Path = DEFAULT_CONTEXT_METADATA,
    chunk_metadata_path: Path = DEFAULT_CHUNK_METADATA,
    checkpoint: Path = DEFAULT_CHECKPOINT,
    output_dir: Path = DEFAULT_OUTPUT_DIR,
    device_name: str = "auto",
    chunk_residues: int = 64,
    force: bool = False,
) -> dict:
    if chunk_residues <= 0:
        raise ValueError("chunk_residues must be positive")
    windows = pd.read_csv(context_metadata_path, sep="\t")
    chunks = pd.read_csv(chunk_metadata_path, sep="\t")
    validate_context_metadata(windows, chunks)
    checkpoint = checkpoint.resolve()
    if not checkpoint.is_file():
        raise FileNotFoundError(f"MotifAE checkpoint not found: {checkpoint}")
    output_dir.mkdir(parents=True, exist_ok=True)
    indices_path = output_dir / "motifae_feature_indices.int32.bin"
    values_path = output_dir / "motifae_feature_values.float32.bin"
    indptr_path = output_dir / "motifae_residue_indptr.int64.npy"
    metadata_out = output_dir / "motifae_window_metadata.tsv"
    config_path = output_dir / "motifae_sparse_config.json"

    context_signature = frame_sha256(
        windows,
        [
            "protein_id", "chunk_id", "window_index", "window_start", "window_end",
            "embedding_path", "sequence_sha256",
        ],
    )
    checkpoint_signature = file_sha256(checkpoint)
    if not force and all(path.is_file() for path in (indices_path, values_path, indptr_path, metadata_out, config_path)):
        config = json.loads(config_path.read_text(encoding="utf-8"))
        if (
            config.get("complete") is True
            and config.get("context_metadata_sha256") == context_signature
            and config.get("checkpoint_sha256") == checkpoint_signature
        ):
            print(f"[OK] Reusing complete MotifAE sparse window store: {output_dir}")
            return config

    unique, mapping = build_unique_windows(windows)
    device = choose_device(device_name)
    encoder = MotifAEEncoder(checkpoint, device)
    if encoder.input_dim != EXPECTED_INPUT_DIM or encoder.dict_size != EXPECTED_LATENT_DIM:
        raise ValueError(
            f"Expected MotifAE {EXPECTED_INPUT_DIM}->{EXPECTED_LATENT_DIM}, got "
            f"{encoder.input_dim}->{encoder.dict_size}"
        )

    residue_indptr = np.empty(len(unique) * WINDOW_LENGTH + 1, dtype=np.int64)
    residue_indptr[0] = 0
    unique_nnz = np.zeros(len(unique), dtype=np.int64)
    unique_mean = np.zeros(len(unique), dtype=np.float32)
    unique_q50 = np.zeros(len(unique), dtype=np.float32)
    unique_q90 = np.zeros(len(unique), dtype=np.float32)
    unique_max = np.zeros(len(unique), dtype=np.int32)
    temp_indices = indices_path.with_name(f".{indices_path.name}.{os.getpid()}.tmp")
    temp_values = values_path.with_name(f".{values_path.name}.{os.getpid()}.tmp")
    nnz_total = 0
    residue_pointer = 0
    start_time = time.monotonic()

    chunk_lookup = {
        (str(row.protein_id), int(row.chunk_id)): row
        for row in chunks.itertuples(index=False)
    }
    try:
        with temp_indices.open("wb") as indices_handle, temp_values.open("wb") as values_handle:
            grouped = unique.groupby(["protein_id", "chunk_id"], sort=True)
            for group_number, ((protein_id, chunk_id), group) in enumerate(grouped, start=1):
                chunk_row = pd.Series(chunk_lookup[(str(protein_id), int(chunk_id))]._asdict())
                embedding, crop_rows = load_chunk(
                    Path(str(chunk_row.embedding_path)), chunk_row
                )

                                                                   
                needed_rows = np.arange(len(embedding), dtype=np.int64)
                sparse_by_stored_row = encode_needed_residues(
                    embedding, needed_rows, encoder, chunk_residues
                )
                sparse_by_crop_row = {
                    int(crop_row): sparse_by_stored_row[stored_row]
                    for stored_row, crop_row in enumerate(crop_rows)
                }
                for row in group.sort_values("sparse_window_index").itertuples(index=False):
                    sparse_index = int(row.sparse_window_index)
                    counts = np.zeros(WINDOW_LENGTH, dtype=np.int32)
                    output_start = int(row.window_output_valid_start)
                    crop_start = int(row.valid_local_start) - 1
                    valid_length = int(row.valid_length)
                    for position in range(WINDOW_LENGTH):
                        if output_start <= position < output_start + valid_length:
                            crop_row = crop_start + position - output_start
                            feature_indices, feature_values = sparse_by_crop_row[crop_row]
                            feature_indices.tofile(indices_handle)
                            feature_values.tofile(values_handle)
                            count = len(feature_indices)
                        else:
                            count = 0
                        counts[position] = count
                        nnz_total += count
                        residue_pointer = sparse_index * WINDOW_LENGTH + position + 1
                        residue_indptr[residue_pointer] = nnz_total
                    unique_nnz[sparse_index] = int(counts.sum())
                    unique_mean[sparse_index] = float(counts.mean())
                    unique_q50[sparse_index] = float(np.quantile(counts, 0.50))
                    unique_q90[sparse_index] = float(np.quantile(counts, 0.90))
                    unique_max[sparse_index] = int(counts.max())
                if group_number % 25 == 0 or group_number == len(grouped):
                    elapsed = time.monotonic() - start_time
                    print(
                        f"[MOTIFAE] chunks={group_number}/{len(grouped)} "
                        f"unique_windows={int(group.sparse_window_index.max()) + 1}/{len(unique)} "
                        f"nnz={nnz_total} elapsed_s={elapsed:.1f}",
                        flush=True,
                    )
        if residue_pointer != len(unique) * WINDOW_LENGTH:
            raise RuntimeError("Sparse windows were not written in contiguous index order")
        os.replace(temp_indices, indices_path)
        os.replace(temp_values, values_path)
    finally:
        for temp in (temp_indices, temp_values):
            if temp.exists():
                temp.unlink()

    np.save(indptr_path, residue_indptr)
    stats = pd.DataFrame(
        {
            "sparse_window_index": np.arange(len(unique), dtype=np.int64),
            "stored_nnz": unique_nnz,
            "mean_active_features_per_residue": unique_mean,
            "q50_active_features": unique_q50,
            "q90_active_features": unique_q90,
            "max_active_features": unique_max,
        }
    )
    mapping = mapping.merge(stats, on="sparse_window_index", how="left", validate="many_to_one")
    mapping.to_csv(metadata_out, sep="\t", index=False)

    dataset_weighted_nnz = int(mapping["stored_nnz"].sum())
    config = {
        "complete": True,
        "format": "global_window_residue_csr",
        "window_length": WINDOW_LENGTH,
        "latent_dim": encoder.dict_size,
        "feature_indices": str(indices_path.resolve()),
        "feature_indices_dtype": "int32",
        "feature_values": str(values_path.resolve()),
        "feature_values_dtype": "float32",
        "residue_indptr": str(indptr_path.resolve()),
        "window_metadata": str(metadata_out.resolve()),
        "window_count": len(mapping),
        "unique_sparse_window_count": len(unique),
        "stored_nnz": int(nnz_total),
        "dataset_weighted_nnz": dataset_weighted_nnz,
        "mean_active_features_per_residue": dataset_weighted_nnz / (len(mapping) * WINDOW_LENGTH),
        "q50_window_active_features_per_residue": float(mapping["mean_active_features_per_residue"].quantile(0.50)),
        "q90_window_active_features_per_residue": float(mapping["mean_active_features_per_residue"].quantile(0.90)),
        "max_active_features_in_any_residue": int(mapping["max_active_features"].max()),
        "checkpoint": str(checkpoint),
        "checkpoint_sha256": checkpoint_signature,
        "context_metadata": str(context_metadata_path.resolve()),
        "context_metadata_sha256": context_signature,
        "device_used": str(device),
        "dense_29_by_40960_constructed": False,
        "padding_representation": "empty_sparse_residue",
    }
    config_path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    print(
        f"[OK] MotifAE sparse windows={len(mapping)} unique={len(unique)} "
        f"mean_active={config['mean_active_features_per_residue']:.3f} "
        f"max_active={config['max_active_features_in_any_residue']} output={output_dir}",
        flush=True,
    )
    return config

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--context-metadata", type=Path, default=DEFAULT_CONTEXT_METADATA)
    parser.add_argument("--chunk-metadata", type=Path, default=DEFAULT_CHUNK_METADATA)
    parser.add_argument("--checkpoint", type=Path, default=DEFAULT_CHECKPOINT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--device", choices=("auto", "cpu", "mps", "cuda"), default="auto")
    parser.add_argument("--chunk-residues", type=int, default=64)
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()

def main() -> int:
    args = parse_args()
    build_motifae_sparse_windows(
        context_metadata_path=args.context_metadata.resolve(),
        chunk_metadata_path=args.chunk_metadata.resolve(),
        checkpoint=args.checkpoint.resolve(),
        output_dir=args.output_dir.resolve(),
        device_name=args.device,
        chunk_residues=args.chunk_residues,
        force=args.force,
    )
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
