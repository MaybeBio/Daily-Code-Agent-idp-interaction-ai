from pathlib import Path

import numpy as np
import torch
from torch.utils.data import Dataset

from idpinterp.models.memmap_store import EMBED_DTYPE


class MemmapDataset(Dataset):
    """
    A PyTorch Dataset that loads embeddings directly from a memory-mapped file.
    This allows streaming massive datasets (e.g., 50k proteins = ~24GB for ESM2 8M) 
    without loading them entirely into RAM.
    """

    def __init__(self, mmap_path: str | Path, total_rows: int, hidden_dim: int):
        self.mmap_path = Path(mmap_path)
        self.total_rows = total_rows
        self.hidden_dim = hidden_dim
        
        # We don't open the memmap in __init__ because PyTorch DataLoader workers
        # need their own file handles if num_workers > 0.
        self._mmap = None

    def _get_mmap(self) -> np.memmap:
        if self._mmap is None:
            self._mmap = np.memmap(
                self.mmap_path,
                dtype=EMBED_DTYPE,
                mode='r',
                shape=(self.total_rows, self.hidden_dim)
            )
        return self._mmap

    def __len__(self) -> int:
        return self.total_rows

    def __getitem__(self, idx: int) -> torch.Tensor:
        mmap = self._get_mmap()
        # astype(float32) both makes a contiguous copy (memmap must not be
        # mutated) and upcasts the float16 store so SAE math runs in float32.
        return torch.from_numpy(mmap[idx].astype(np.float32))


class BufferedShuffleDataset(torch.utils.data.IterableDataset):
    """
    An IterableDataset that loads embeddings in large contiguous chunks,
    shuffles them in RAM, and yields them one by one. This eliminates the 
    massive I/O bottleneck of purely random disk seeks while preserving 
    the stochasticity required for SAE training.
    """

    def __init__(self, mmap_path: str | Path, total_rows: int, hidden_dim: int, buffer_size: int = 524288, batch_size: int = 2048):
        self.mmap_path = Path(mmap_path)
        self.total_rows = total_rows
        self.hidden_dim = hidden_dim
        self.batch_size = batch_size
        # ~500k rows is roughly 670 MB per worker
        self.buffer_size = min(buffer_size, total_rows)

    def __iter__(self):
        worker_info = torch.utils.data.get_worker_info()
        seed = np.random.randint(0, 2**32 - 1)
        if worker_info is not None:
            seed = (seed + worker_info.id) % (2**32)
        rng = np.random.default_rng(seed)

        mmap = np.memmap(
            self.mmap_path,
            dtype=EMBED_DTYPE,
            mode='r',
            shape=(self.total_rows, self.hidden_dim)
        )

        while True:
            # Pick a random contiguous chunk
            start_idx = rng.integers(0, self.total_rows - self.buffer_size + 1)
            # Read chunk into RAM (fast sequential read; float16 → ~half the bytes)
            chunk = mmap[start_idx : start_idx + self.buffer_size].copy()

            # Shuffle indices instead of the massive 2D array to save CPU time
            indices = np.arange(self.buffer_size)
            rng.shuffle(indices)

            # Yield full batches to completely bypass DataLoader collation overhead
            for i in range(0, self.buffer_size, self.batch_size):
                if i + self.batch_size <= self.buffer_size:
                    batch_indices = indices[i : i + self.batch_size]
                    # upcast float16 store → float32 for SAE math (tiny per-batch cost)
                    yield torch.from_numpy(chunk[batch_indices].astype(np.float32))
