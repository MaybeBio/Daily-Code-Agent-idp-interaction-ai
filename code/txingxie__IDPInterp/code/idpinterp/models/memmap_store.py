import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

# On-disk dtype for PLM embeddings. float16 halves the embedding footprint and
# read bandwidth vs float32 (ESM-2 itself runs fine in half precision); readers
# convert to float32 in RAM so SAE training math is unchanged.
# Changing this makes previously generated embedding mmaps unreadable; delete and
# regenerate them.
EMBED_DTYPE = np.float16


class MemmapWriter:
    """Writes embeddings sequentially to a memory-mapped file."""

    def __init__(self, file_path: Path, total_rows: int, hidden_dim: int):
        self.file_path = file_path
        self.file_path.parent.mkdir(parents=True, exist_ok=True)
        self.total_rows = total_rows
        self.hidden_dim = hidden_dim

        logger.info(f"Allocating memmap: {file_path} ({total_rows}x{hidden_dim}, {np.dtype(EMBED_DTYPE).name})")
        self.mmap = np.memmap(
            file_path,
            dtype=EMBED_DTYPE,
            mode='w+',
            shape=(total_rows, hidden_dim)
        )
        self.current_idx = 0

    def write_batch(self, batch: np.ndarray) -> None:
        """Write a batch of shape [N, hidden_dim] to the memmap."""
        n_rows = batch.shape[0]
        if self.current_idx + n_rows > self.total_rows:
            raise ValueError(f"Memmap full! Cannot write {n_rows} more rows.")

        self.mmap[self.current_idx : self.current_idx + n_rows] = batch
        self.current_idx += n_rows

    def flush(self) -> None:
        """Flush changes to disk."""
        self.mmap.flush()

    def finalize(self) -> None:
        """Flush and delete the memmap object."""
        self.flush()
        del self.mmap


class MemmapReader:
    """Reads embeddings from a memory-mapped file."""

    def __init__(self, file_path: Path, total_rows: int, hidden_dim: int):
        self.mmap = np.memmap(
            file_path,
            dtype=EMBED_DTYPE,
            mode='r',
            shape=(total_rows, hidden_dim)
        )

    def __getitem__(self, idx):
        return self.mmap[idx]


def save_mmap_config(file_path: Path, total_residues: int, hidden_dim: int) -> None:
    """Save the shape configuration of the memmap alongside it."""
    config_path = file_path.with_suffix(".json")
    with open(config_path, "w") as f:
        json.dump({
            "total_residues": total_residues,
            "hidden_dim": hidden_dim
        }, f)


def load_mmap_config(file_path: Path) -> dict:
    """Load the shape configuration of the memmap."""
    config_path = file_path.with_suffix(".json")
    if not config_path.exists():
        raise FileNotFoundError(f"Missing config for memmap: {config_path}")
    with open(config_path, "r") as f:
        return json.load(f)


class MetadataIndex:
    """Handles mapping between protein IDs and their positions in the memmap."""

    @staticmethod
    def build(
        protein_ids: list[str], 
        residue_counts: list[int], 
        was_truncated: list[bool]
    ) -> pd.DataFrame:
        """Build the metadata index DataFrame."""
        starts = []
        ends = []
        current = 0
        for count in residue_counts:
            starts.append(current)
            current += count
            ends.append(current)

        return pd.DataFrame({
            "uniprot_id": protein_ids,
            "mmap_start_idx": starts,
            "mmap_end_idx": ends,
            "stored_length": residue_counts,
            "was_truncated": was_truncated,
        })

    @staticmethod
    def save(df: pd.DataFrame, path: Path) -> None:
        """Save the metadata index."""
        df.to_parquet(path, engine="pyarrow", index=False)

    @staticmethod
    def load(path: Path) -> pd.DataFrame:
        """Load the metadata index."""
        return pd.read_parquet(path)
