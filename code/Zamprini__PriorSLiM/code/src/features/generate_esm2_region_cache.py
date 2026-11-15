"""Generate regional ESM2 embedding arrays."""
                      

from __future__ import annotations
import os
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import torch

MODEL_NAME = "esm2_t33_650M_UR50D"
REPR_LAYER = 33
EMBEDDING_DIM = 1280
MAX_RESIDUES = 1022

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--region-manifest", required=True, type=Path)
    parser.add_argument("--embedding-dir", required=True, type=Path)
    parser.add_argument("--metadata-output", required=True, type=Path)
    parser.add_argument("--batch-size", type=int, default=1)
    parser.add_argument("--device", choices=("auto", "cpu", "mps", "cuda"), default="auto")
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()

def choose_device(requested: str):
    if requested != "auto":
        return torch.device(requested)
    if torch.cuda.is_available():
        return torch.device("cuda")
    if hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
        return torch.device("mps")
    return torch.device("cpu")

def existing_embedding_is_valid(path: Path, length: int):
    if not path.is_file():
        return False
    try:
        array = np.load(path, mmap_mode="r")
    except (OSError, ValueError):
        return False
    return array.shape == (length, EMBEDDING_DIM) and array.dtype == np.float32

def main():
    args = parse_args()
    manifest = pd.read_csv(args.region_manifest, sep="\t")
    required = {
        "uniprot", "region_type", "original_start", "original_end",
        "crop_start", "crop_end", "original_length", "crop_length",
        "was_truncated", "sequence", "embedding_file",
    }
    missing_columns = required - set(manifest.columns)
    if missing_columns:
        raise ValueError(f"Region manifest is missing columns: {sorted(missing_columns)}")
    if (manifest["crop_length"] > MAX_RESIDUES).any():
        raise ValueError(f"Every ESM input must contain at most {MAX_RESIDUES} residues")
    sequence_lengths = manifest["sequence"].astype(str).str.len()
    if not np.array_equal(sequence_lengths.to_numpy(), manifest["crop_length"].to_numpy()):
        raise ValueError("crop_length does not match the sequence sent to ESM")

    args.embedding_dir.mkdir(parents=True, exist_ok=True)
    pending = []
    for index, row in manifest.iterrows():
        path = args.embedding_dir / str(row["embedding_file"])
        if args.force or not existing_embedding_is_valid(path, int(row["crop_length"])):
            pending.append(index)

                                                                       
    pending.sort(key=lambda index: (int(manifest.at[index, "crop_length"]), index))

    device = choose_device(args.device)
    if pending:
        import esm

        model_loader = getattr(esm.pretrained, MODEL_NAME)
        model, alphabet = model_loader()
        model = model.to(device).eval()
        batch_converter = alphabet.get_batch_converter()
        with torch.no_grad():
            for offset in range(0, len(pending), args.batch_size):
                indices = pending[offset : offset + args.batch_size]
                batch = [
                    (
                        f"{manifest.at[index, 'uniprot']}|{manifest.at[index, 'region_type']}",
                        str(manifest.at[index, "sequence"]),
                    )
                    for index in indices
                ]
                if max(len(sequence) for _label, sequence in batch) > MAX_RESIDUES:
                    raise RuntimeError("An overlength sequence reached the ESM model")
                _labels, _sequences, tokens = batch_converter(batch)
                tokens = tokens.to(device)
                result = model(tokens, repr_layers=[REPR_LAYER], return_contacts=False)
                representations = result["representations"][REPR_LAYER]
                for batch_index, manifest_index in enumerate(indices):
                    length = int(manifest.at[manifest_index, "crop_length"])
                    embedding = (
                        representations[batch_index, 1 : length + 1]
                        .detach().to("cpu", dtype=torch.float32).numpy().copy()
                    )
                    if embedding.shape != (length, EMBEDDING_DIM):
                        raise RuntimeError(
                            f"Unexpected {MODEL_NAME} output shape: {embedding.shape}"
                        )
                    np.save(
                        args.embedding_dir / str(manifest.at[manifest_index, "embedding_file"]),
                        embedding,
                    )

    metadata = manifest.drop(columns=["sequence"]).copy()
    metadata["embedding_path"] = [
        str((args.embedding_dir / filename).resolve())
        for filename in metadata["embedding_file"].astype(str)
    ]
    metadata["esm_model"] = MODEL_NAME
    metadata["esm_repr_layer"] = REPR_LAYER
    metadata["embedding_dim"] = EMBEDDING_DIM
    metadata["esm_max_residues"] = MAX_RESIDUES
    args.metadata_output.parent.mkdir(parents=True, exist_ok=True)
    metadata.to_csv(args.metadata_output, sep="\t", index=False)
    print(
        f"[OK] ESM region cache rows={len(metadata)} newly_generated={len(pending)} "
        f"device={device} output={args.metadata_output}"
    )

if __name__ == "__main__":
    main()
