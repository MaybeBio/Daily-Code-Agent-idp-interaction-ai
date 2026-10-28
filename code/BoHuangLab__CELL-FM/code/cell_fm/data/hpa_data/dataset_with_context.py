from typing import List
import os

import h5py
import torch

from .dataset import HPADataset, HPAAllImageDataset
from .collater import collate_fn


def _load_context(data_path):
    h5_path = os.path.join(data_path, 'context', 'interaction.h5')
    with h5py.File(h5_path, 'r') as f:
        query_gene_names = [n.decode() if isinstance(n, bytes) else n for n in f['query_gene_names'][:]]
        embeddings = f['embeddings'][:]  # [N, n_context, emb_dim]
    gene_to_idx = {name: i for i, name in enumerate(query_gene_names)}
    return embeddings, gene_to_idx


class _ContextMixin:
    """Mixin that adds interaction_context lookup to any HPADataset subclass."""

    def _init_context(self):
        self._embeddings, self.gene_to_context_idx = _load_context(self.data_path)
        self.n_context = self._embeddings.shape[1]
        self.context_embedding_dim = self._embeddings.shape[2]

    def _get_context(self, gene_name):
        if gene_name in self.gene_to_context_idx:
            idx = self.gene_to_context_idx[gene_name]
            return torch.from_numpy(self._embeddings[idx].copy())
        return torch.zeros(self.n_context, self.context_embedding_dim, dtype=torch.float32)


class HPADatasetWithContext(_ContextMixin, HPADataset):
    """HPADataset extended with pre-computed protein interaction context."""

    def __init__(self, args, split_key) -> None:
        super().__init__(args, split_key)
        self._init_context()

    def __getitem__(self, index: int) -> dict:
        item = super().__getitem__(index)
        item['interaction_context'] = self._get_context(item['gene_name'])
        return item

    def collate(self, samples: List[dict]) -> dict:
        batch = collate_fn(samples, self.vocab, self.max_protein_sequence_len, 0)
        batch['batched_data']['interaction_context'] = torch.stack(
            [s['interaction_context'] for s in samples]
        )
        return batch


class HPAAllImageDatasetWithContext(_ContextMixin, HPAAllImageDataset):
    """HPAAllImageDataset extended with pre-computed protein interaction context."""

    def __init__(self, args, split_key) -> None:
        super().__init__(args, split_key)
        self._init_context()

    def __getitem__(self, index: int) -> dict:
        item = super().__getitem__(index)
        item['interaction_context'] = self._get_context(item['gene_name'])
        return item
