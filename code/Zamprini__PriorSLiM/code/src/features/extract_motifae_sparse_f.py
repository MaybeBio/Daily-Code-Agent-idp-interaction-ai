"""Frozen MotifAE encoder and regional sparse extraction utilities."""
                      

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
from typing import Any, Mapping

                                                                     
import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F

SCRIPT_PATH = Path(__file__).resolve()
PROJECT_ROOT = Path(os.environ.get("PRIORSLIM_DATA_ROOT", SCRIPT_PATH.parents[2] / "data"))
DEFAULT_METADATA = PROJECT_ROOT / "feature_cache/esm2_region_metadata.tsv"
DEFAULT_EMBEDDING_DIR = PROJECT_ROOT / "feature_cache/esm2_regions"
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "feature_cache/motifae_sparse_f_k3"
EXPECTED_ESM_MODEL = "esm2_t33_650M_UR50D"
EXPECTED_ESM_LAYER = 33
EXPECTED_INPUT_DIM = 1280

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Convert cached raw ESM2-650M layer-33 residue embeddings to "
            "sparse per-residue MotifAE activations."
        )
    )
    parser.add_argument(
        "--checkpoint",
        type=Path,
        required=True,
        help="Released MotifAE checkpoint, e.g. MotifAE_k3_step_80000.pt.",
    )
    parser.add_argument("--metadata", type=Path, default=DEFAULT_METADATA)
    parser.add_argument("--embedding-dir", type=Path, default=DEFAULT_EMBEDDING_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument(
        "--device",
        choices=("auto", "cpu", "mps", "cuda"),
        default="auto",
        help="auto chooses CUDA, then MPS, then CPU.",
    )
    parser.add_argument(
        "--chunk-residues",
        type=int,
        default=64,
        help="Residues encoded at once; 64 uses about 10 MiB for dense f.",
    )
    parser.add_argument(
        "--value-dtype",
        choices=("float32", "float16"),
        default="float32",
        help="Storage dtype for non-zero values; computation remains float32.",
    )
    parser.add_argument(
        "--no-compress",
        action="store_true",
        help="Use np.savez instead of np.savez_compressed.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Recompute valid existing region outputs.",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Process only the first N manifest rows (for a smoke test).",
    )
    parser.add_argument(
        "--progress-every",
        type=int,
        default=25,
        help="Print progress after this many regions.",
    )
    return parser.parse_args()

def choose_device(requested: str) -> torch.device:
    if requested == "auto":
        if torch.cuda.is_available():
            return torch.device("cuda")
        if torch.backends.mps.is_available():
            return torch.device("mps")
        return torch.device("cpu")
    device = torch.device(requested)
    if requested == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("CUDA was requested but torch.cuda.is_available() is false")
    if requested == "mps" and not torch.backends.mps.is_available():
        raise RuntimeError("MPS was requested but torch.backends.mps.is_available() is false")
    return device

def torch_load_weights(path: Path) -> Any:
    try:
        return torch.load(path, map_location="cpu", weights_only=True)
    except TypeError:                                            
        return torch.load(path, map_location="cpu")

def unwrap_state_dict(payload: Any) -> Mapping[str, torch.Tensor]:
    if not isinstance(payload, Mapping):
        raise TypeError(f"Checkpoint must contain a mapping, got {type(payload)!r}")

    state: Mapping[str, Any] = payload
    for wrapper_key in ("state_dict", "model_state_dict", "model"):
        wrapped = state.get(wrapper_key)
        if isinstance(wrapped, Mapping):
            state = wrapped
            break

    clean: dict[str, torch.Tensor] = {}
    for key, value in state.items():
        if not isinstance(value, torch.Tensor):
            continue
        new_key = str(key)
        for prefix in ("module.", "autoencoder.", "model."):
            if new_key.startswith(prefix):
                new_key = new_key[len(prefix) :]
        clean[new_key] = value
    return clean

# MotifAE encoder architecture: Chao.Hou (2025); see MotifAE_LICENSE.
class MotifAEEncoder:
    pass

    def __init__(self, checkpoint: Path, device: torch.device):
        payload = torch_load_weights(checkpoint)
        state = unwrap_state_dict(payload)
        required = ("bias", "encoder.weight", "encoder.bias")
        missing = [key for key in required if key not in state]
        if missing:
            raise KeyError(f"Checkpoint lacks MotifAE encoder tensors: {missing}")

        shared_bias = state["bias"].detach().to(dtype=torch.float32, device="cpu")
        weight = state["encoder.weight"].detach().to(dtype=torch.float32, device="cpu")
        encoder_bias = state["encoder.bias"].detach().to(dtype=torch.float32, device="cpu")
        del payload, state

        if weight.ndim != 2:
            raise ValueError(f"encoder.weight must be 2-D, got {tuple(weight.shape)}")
        self.dict_size, self.input_dim = map(int, weight.shape)
        if tuple(shared_bias.shape) != (self.input_dim,):
            raise ValueError("Shared bias shape does not match encoder input dimension")
        if tuple(encoder_bias.shape) != (self.dict_size,):
            raise ValueError("Encoder bias shape does not match latent dimension")
        if self.input_dim != EXPECTED_INPUT_DIM:
            raise ValueError(
                f"Checkpoint input dimension is {self.input_dim}; expected {EXPECTED_INPUT_DIM}"
            )

        self.device = device
        self.shared_bias = shared_bias.to(device)
        self.weight = weight.to(device)
        self.encoder_bias = encoder_bias.to(device)

    def encode(self, x: torch.Tensor) -> torch.Tensor:

        return F.relu(F.linear(x - self.shared_bias, self.weight, self.encoder_bias))

def checkpoint_fingerprint(path: Path) -> str:
    pass
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(8 * 1024 * 1024):
            digest.update(block)
    return digest.hexdigest()

def validate_metadata(frame: pd.DataFrame) -> None:
    required = {
        "uniprot",
        "region_type",
        "crop_length",
        "embedding_file",
        "embedding_path",
        "esm_model",
        "esm_repr_layer",
        "embedding_dim",
    }
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(f"Metadata is missing columns: {sorted(missing)}")
    if frame["embedding_file"].duplicated().any():
        duplicates = frame.loc[frame["embedding_file"].duplicated(), "embedding_file"].head()
        raise ValueError(f"Metadata has duplicate embedding files: {duplicates.tolist()}")
    if set(frame["esm_model"].astype(str)) != {EXPECTED_ESM_MODEL}:
        raise ValueError(f"Expected only {EXPECTED_ESM_MODEL} embeddings")
    if set(pd.to_numeric(frame["esm_repr_layer"])) != {EXPECTED_ESM_LAYER}:
        raise ValueError(f"Expected only ESM representation layer {EXPECTED_ESM_LAYER}")
    if set(pd.to_numeric(frame["embedding_dim"])) != {EXPECTED_INPUT_DIM}:
        raise ValueError(f"Expected only {EXPECTED_INPUT_DIM}-dimensional embeddings")

def resolve_embedding_path(row: pd.Series, embedding_dir: Path) -> Path:
    recorded = Path(str(row["embedding_path"]))
    if recorded.is_file():
        return recorded
    fallback = embedding_dir / str(row["embedding_file"])
    if fallback.is_file():
        return fallback
    raise FileNotFoundError(
        f"Embedding not found at recorded or fallback path: {recorded}, {fallback}"
    )

def output_path_for(embedding_file: str, output_dir: Path) -> Path:
    stem = Path(embedding_file).stem
    return output_dir / f"{stem}.motifae_f.npz"

def inspect_existing(path: Path, expected_rows: int, dict_size: int) -> dict[str, Any]:
    with np.load(path, allow_pickle=False) as data:
        shape = tuple(int(x) for x in data["shape"])
        if shape != (expected_rows, dict_size):
            raise ValueError(f"Existing output {path} has shape {shape}, expected {(expected_rows, dict_size)}")
        nnz = int(data["indices"].shape[0])
        if int(data["values"].shape[0]) != nnz:
            raise ValueError(f"Existing output {path} has mismatched indices/values")
        indptr = data["indptr"]
        if tuple(indptr.shape) != (expected_rows + 1,) or int(indptr[-1]) != nnz:
            raise ValueError(f"Existing output {path} has invalid indptr")
        value_dtype = str(data["values"].dtype)
        min_value = float(data["values"].min()) if nnz else 0.0
        max_value = float(data["values"].max()) if nnz else 0.0
    return {
        "nnz": nnz,
        "value_dtype": value_dtype,
        "min_nonzero_value": min_value,
        "max_nonzero_value": max_value,
    }

def encode_region_sparse(
    array: np.ndarray,
    encoder: MotifAEEncoder,
    chunk_residues: int,
    value_dtype: np.dtype,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if array.ndim != 2 or array.shape[1] != encoder.input_dim:
        raise ValueError(
            f"Embedding shape must be (L, {encoder.input_dim}), got {array.shape}"
        )
    n_residues = int(array.shape[0])
    counts: list[np.ndarray] = []
    all_indices: list[np.ndarray] = []
    all_values: list[np.ndarray] = []

    with torch.inference_mode():
        for start in range(0, n_residues, chunk_residues):
            stop = min(start + chunk_residues, n_residues)

                                                                                  
            chunk_np = np.asarray(array[start:stop], dtype=np.float32)
            x = torch.tensor(chunk_np, dtype=torch.float32, device=encoder.device)
            latent = encoder.encode(x)
            positive_positions = torch.nonzero(latent > 0, as_tuple=False)

                                                                           
            positions_cpu = positive_positions.to(device="cpu", dtype=torch.int64)
            row_counts = np.bincount(
                positions_cpu[:, 0].numpy(), minlength=stop - start
            ).astype(np.int64, copy=False)
            latent_indices = positions_cpu[:, 1].numpy().astype(np.int32, copy=False)
            latent_values = latent[
                positive_positions[:, 0], positive_positions[:, 1]
            ].to(device="cpu", dtype=torch.float32)

            counts.append(row_counts)
            all_indices.append(latent_indices)
            all_values.append(latent_values.numpy().astype(value_dtype, copy=False))
            del x, latent, positive_positions, positions_cpu, latent_values

    counts_array = np.concatenate(counts) if counts else np.empty(0, dtype=np.int64)
    indptr = np.empty(n_residues + 1, dtype=np.int64)
    indptr[0] = 0
    np.cumsum(counts_array, out=indptr[1:])
    indices = (
        np.concatenate(all_indices).astype(np.int32, copy=False)
        if all_indices
        else np.empty(0, dtype=np.int32)
    )
    values = (
        np.concatenate(all_values).astype(value_dtype, copy=False)
        if all_values
        else np.empty(0, dtype=value_dtype)
    )
    if int(indptr[-1]) != len(indices) or len(indices) != len(values):
        raise RuntimeError("Internal error while constructing sparse arrays")
    return indptr, indices, values

def save_region(
    output_path: Path,
    indptr: np.ndarray,
    indices: np.ndarray,
    values: np.ndarray,
    n_rows: int,
    dict_size: int,
    source_file: str,
    checkpoint_name: str,
    checkpoint_sha256: str,
    compress: bool,
) -> None:
    temp_path = output_path.with_name(f".{output_path.name}.{os.getpid()}.tmp")
    save = np.savez_compressed if compress else np.savez
    try:
        with temp_path.open("wb") as handle:
            save(
                handle,
                indptr=indptr,
                indices=indices,
                values=values,
                shape=np.asarray([n_rows, dict_size], dtype=np.int64),
                source_embedding_file=np.asarray(source_file),
                checkpoint_name=np.asarray(checkpoint_name),
                checkpoint_sha256=np.asarray(checkpoint_sha256),
                representation=np.asarray("per_residue_csr_strict_relu_positive"),
            )
        os.replace(temp_path, output_path)
    finally:
        if temp_path.exists():
            temp_path.unlink()

def write_manifest_atomic(frame: pd.DataFrame, path: Path) -> None:
    temp_path = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    frame.to_csv(temp_path, sep="\t", index=False)
    os.replace(temp_path, path)

def main() -> int:
    args = parse_args()
    if args.chunk_residues <= 0:
        raise ValueError("--chunk-residues must be positive")
    if args.limit is not None and args.limit <= 0:
        raise ValueError("--limit must be positive")
    if args.progress_every <= 0:
        raise ValueError("--progress-every must be positive")

    checkpoint = args.checkpoint.expanduser().resolve()
    metadata_path = args.metadata.expanduser().resolve()
    embedding_dir = args.embedding_dir.expanduser().resolve()
    output_dir = args.output_dir.expanduser().resolve()
    if not checkpoint.is_file():
        raise FileNotFoundError(f"MotifAE checkpoint not found: {checkpoint}")
    if not metadata_path.is_file():
        raise FileNotFoundError(f"ESM metadata not found: {metadata_path}")

    metadata = pd.read_csv(metadata_path, sep="\t")
    validate_metadata(metadata)
    if args.limit is not None:
        metadata = metadata.iloc[: args.limit].copy()

    device = choose_device(args.device)
    output_dir.mkdir(parents=True, exist_ok=True)
    encoder = MotifAEEncoder(checkpoint, device)
    if encoder.dict_size != 40_960:
        print(
            f"[WARN] Checkpoint latent width is {encoder.dict_size}, not released-model 40960",
            file=sys.stderr,
        )
    checkpoint_sha256 = checkpoint_fingerprint(checkpoint)
    storage_dtype = np.dtype(args.value_dtype)
    manifest_path = output_dir / "motifae_sparse_f_manifest.tsv"
    config_path = output_dir / "motifae_sparse_f_config.json"
    config = {
        "checkpoint": str(checkpoint),
        "checkpoint_sha256": checkpoint_sha256,
        "metadata": str(metadata_path),
        "embedding_dir": str(embedding_dir),
        "output_dir": str(output_dir),
        "device": str(device),
        "chunk_residues": args.chunk_residues,
        "value_dtype": args.value_dtype,
        "compressed": not args.no_compress,
        "input_esm_model": EXPECTED_ESM_MODEL,
        "input_esm_layer": EXPECTED_ESM_LAYER,
        "input_dim": encoder.input_dim,
        "latent_dim": encoder.dict_size,
        "sparsity_rule": "strict_relu_positive_f_gt_0",
        "representation": "per_residue_csr",
        "manifest_rows_selected": len(metadata),
    }
    config_path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")

    print(
        f"[INFO] regions={len(metadata)} device={device} input_dim={encoder.input_dim} "
        f"latent_dim={encoder.dict_size} output={output_dir}",
        flush=True,
    )
    start_time = time.monotonic()
    records: list[dict[str, Any]] = []
    generated = 0
    reused = 0
    total_residues = 0
    total_nnz = 0

    for ordinal, (_, row) in enumerate(metadata.iterrows(), start=1):
        embedding_path = resolve_embedding_path(row, embedding_dir)
        embedding = np.load(embedding_path, mmap_mode="r", allow_pickle=False)
        expected_rows = int(row["crop_length"])
        if tuple(embedding.shape) != (expected_rows, encoder.input_dim):
            raise ValueError(
                f"{embedding_path}: shape {embedding.shape}, expected "
                f"{(expected_rows, encoder.input_dim)}"
            )

        output_path = output_path_for(str(row["embedding_file"]), output_dir)
        if output_path.exists() and not args.overwrite:
            sparse_info = inspect_existing(output_path, expected_rows, encoder.dict_size)
            reused += 1
        else:
            indptr, indices, values = encode_region_sparse(
                embedding,
                encoder=encoder,
                chunk_residues=args.chunk_residues,
                value_dtype=storage_dtype,
            )
            save_region(
                output_path,
                indptr,
                indices,
                values,
                n_rows=expected_rows,
                dict_size=encoder.dict_size,
                source_file=str(row["embedding_file"]),
                checkpoint_name=checkpoint.name,
                checkpoint_sha256=checkpoint_sha256,
                compress=not args.no_compress,
            )
            sparse_info = {
                "nnz": len(indices),
                "value_dtype": str(values.dtype),
                "min_nonzero_value": float(values.min()) if len(values) else 0.0,
                "max_nonzero_value": float(values.max()) if len(values) else 0.0,
            }
            generated += 1

        nnz = int(sparse_info["nnz"])
        total_residues += expected_rows
        total_nnz += nnz
        records.append(
            {
                "uniprot": row["uniprot"],
                "region_type": row["region_type"],
                "original_start": row.get("original_start", np.nan),
                "original_end": row.get("original_end", np.nan),
                "crop_start": row.get("crop_start", np.nan),
                "crop_end": row.get("crop_end", np.nan),
                "crop_length": expected_rows,
                "source_embedding_file": row["embedding_file"],
                "source_embedding_path": str(embedding_path),
                "motifae_sparse_file": output_path.name,
                "motifae_sparse_path": str(output_path),
                "latent_dim": encoder.dict_size,
                "nnz": nnz,
                "mean_active_features_per_residue": nnz / max(expected_rows, 1),
                "density": nnz / max(expected_rows * encoder.dict_size, 1),
                "value_dtype": sparse_info["value_dtype"],
                "min_nonzero_value": sparse_info["min_nonzero_value"],
                "max_nonzero_value": sparse_info["max_nonzero_value"],
                "checkpoint_name": checkpoint.name,
                "checkpoint_sha256": checkpoint_sha256,
            }
        )

        if ordinal % args.progress_every == 0 or ordinal == len(metadata):
            elapsed = time.monotonic() - start_time
            rate = total_residues / max(elapsed, 1e-9)
            print(
                f"[PROGRESS] {ordinal}/{len(metadata)} generated={generated} reused={reused} "
                f"residues={total_residues} rate={rate:.1f} residue/s "
                f"mean_active={total_nnz / max(total_residues, 1):.2f}",
                flush=True,
            )

    result_manifest = pd.DataFrame.from_records(records)
    write_manifest_atomic(result_manifest, manifest_path)
    elapsed = time.monotonic() - start_time
    print(
        f"[OK] regions={len(records)} generated={generated} reused={reused} "
        f"residues={total_residues} nnz={total_nnz} "
        f"mean_active={total_nnz / max(total_residues, 1):.2f} "
        f"elapsed_s={elapsed:.1f} manifest={manifest_path}",
        flush=True,
    )
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
