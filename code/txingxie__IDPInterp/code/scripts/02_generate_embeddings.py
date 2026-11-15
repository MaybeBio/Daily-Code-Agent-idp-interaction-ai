#!/usr/bin/env python3
"""
Generate per-residue PLM embeddings and stream them to memmap files.

Usage:
    python scripts/02_generate_embeddings.py paths=colab model=esm2_8m
    python scripts/02_generate_embeddings.py paths=colab model=esm2_650m
"""

import logging
import sys
from pathlib import Path

import hydra
import numpy as np
import pandas as pd
import torch
from omegaconf import DictConfig
from tqdm import tqdm

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from idpinterp.constants import ESM2_MAX_SEQ_LEN
from idpinterp.models.embedder import get_embedder
from idpinterp.models.memmap_store import (
    MemmapWriter,
    MetadataIndex,
    save_mmap_config,
)

logger = logging.getLogger(__name__)


@hydra.main(version_base=None, config_path="../conf", config_name="config")
def main(cfg: DictConfig) -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(name)s] %(levelname)s: %(message)s",
    )

    logger.info(f"=== Embedding generation ===")
    logger.info(f"Model: {cfg.model.name} ({cfg.model.hf_path})")

    # Set seeds
    torch.manual_seed(cfg.seed)
    np.random.seed(cfg.seed)

    # Load dataset
    dataset_path = Path(cfg.paths.dataset_dir) / f"{cfg.dataset.name}.parquet"
    logger.info(f"Loading dataset: {dataset_path}")
    df = pd.read_parquet(dataset_path)
    logger.info(f"  {len(df)} proteins loaded.")

    # Sort by length descending — process longest first to catch OOM early
    df = df.sort_values("length", ascending=False).reset_index(drop=True)

    # Determine device
    device = "cuda" if torch.cuda.is_available() else "cpu"
    logger.info(f"Device: {device}")

    # Initialize embedder
    embedder = get_embedder(
        model_name=cfg.model.name,
        hf_path=cfg.model.hf_path,
        device=device,
        max_seq_len=cfg.model.max_seq_len,
    )

    max_seq_len = cfg.model.get("max_seq_len", ESM2_MAX_SEQ_LEN)
    batch_size = cfg.model.get("batch_size_inference", 8)

    # Process each layer
    for layer in cfg.model.layers_to_extract:
        logger.info(f"\n--- Processing layer {layer} ---")

        embed_dir = Path(cfg.paths.embedding_dir) / cfg.model.name / cfg.dataset.name
        mmap_path = embed_dir / f"layer{layer}.mmap"
        meta_path = embed_dir / f"layer{layer}_metadata.parquet"

        # Check for existing completed embeddings
        if mmap_path.exists() and meta_path.exists():
            existing_meta = pd.read_parquet(meta_path)
            if len(existing_meta) == len(df):
                logger.info(f"  Layer {layer} already complete. Skipping.")
                continue
            else:
                logger.info(
                    f"  Partial embeddings found ({len(existing_meta)}/{len(df)}). "
                    f"Regenerating from scratch."
                )
        
        mmap_path.unlink(missing_ok=True)
        meta_path.unlink(missing_ok=True)

        # Pre-compute total residues for memmap pre-allocation
        stored_lengths = [min(row["length"], max_seq_len) for _, row in df.iterrows()]
        total_residues = sum(stored_lengths)
        hidden_dim = embedder.get_hidden_dim()

        from idpinterp.models.memmap_store import EMBED_DTYPE
        itemsize = np.dtype(EMBED_DTYPE).itemsize
        logger.info(
            f"  Total residues: {total_residues:,} "
            f"({total_residues * hidden_dim * itemsize / 1e9:.2f} GB, {np.dtype(EMBED_DTYPE).name})"
        )

        # Initialize memmap writer
        writer = MemmapWriter(mmap_path, total_residues, hidden_dim)

        # Track metadata
        protein_ids: list[str] = []
        residue_counts: list[int] = []
        was_truncated: list[bool] = []

        # Process in batches
        sequences = df["sequence"].tolist()
        ids = df["uniprot_id"].tolist()
        original_lengths = df["length"].tolist()

        for batch_start in tqdm(
            range(0, len(sequences), batch_size),
            desc=f"Layer {layer}",
            unit="batch",
            mininterval=30,
        ):
            batch_end = min(batch_start + batch_size, len(sequences))
            batch_seqs = sequences[batch_start:batch_end]
            batch_ids = ids[batch_start:batch_end]
            batch_orig_lens = original_lengths[batch_start:batch_end]

            try:
                embeddings, lengths = embedder.encode_batch(batch_seqs, layer)

                # Write to memmap
                writer.write_batch(embeddings.numpy())

                # Track metadata
                for uid, n_res, orig_len in zip(batch_ids, lengths, batch_orig_lens):
                    protein_ids.append(uid)
                    residue_counts.append(n_res)
                    was_truncated.append(orig_len > max_seq_len)

            except RuntimeError as e:
                if "out of memory" in str(e).lower():
                    logger.error(
                        f"OOM at batch {batch_start}. "
                        f"Reduce batch_size_inference in model config."
                    )
                    torch.cuda.empty_cache()
                    raise
                raise

            # Periodic flush
            if (batch_start // batch_size) % 50 == 0:
                writer.flush()

        writer.finalize()
        save_mmap_config(mmap_path, total_residues, hidden_dim)

        # Save metadata index
        meta_df = MetadataIndex.build(protein_ids, residue_counts, was_truncated)
        MetadataIndex.save(meta_df, meta_path)

        # Clear GPU memory between layers
        torch.cuda.empty_cache()

    logger.info("\n=== Embedding generation complete ===")


if __name__ == "__main__":
    main()
