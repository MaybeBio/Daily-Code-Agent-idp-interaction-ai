"""Data loading (same format as the original repository).

* Embeddings: an .npz with keys `ids` (protein IDs) and `per_residue` (object array of (L, 1024)
  ProtT5 embeddings), optionally `dim`.
* Labels: one TSV per task, each line `<id>\\t<0/1 string of length L>`.

Sequences are always used at full length (no truncation / chunking); mini-batches are dynamically
padded to the longest sequence in the batch and a binary mask marks valid residues.
"""
import os
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import torch
from torch.utils.data import DataLoader, Dataset

TASKS = ["disorder", "protein", "rna", "dna"]


def log(s: str):
    print(f"[LOG] {s}", flush=True)


def load_npz_ids(npz_path: str) -> List[str]:
    with np.load(npz_path, allow_pickle=True) as npz:
        return [str(x) for x in npz["ids"]]


def split_ids(ids: Sequence[str], val_ratio: float = 0.1, seed: int = 42) -> Tuple[List[str], List[str]]:
    rng = np.random.RandomState(seed)
    idx = np.arange(len(ids))
    rng.shuffle(idx)
    n_val = int(len(ids) * val_ratio)
    if n_val == 0 and len(ids) > 1:
        n_val = 1
    return [ids[i] for i in idx[n_val:]], [ids[i] for i in idx[:n_val]]


class IDREmbDataset(Dataset):
    def __init__(self, npz_path: str, label_files: Dict[str, str], id_key: str = "ids",
                 per_res_key: str = "per_residue", include_ids: Optional[set] = None):
        super().__init__()
        npz = np.load(npz_path, allow_pickle=True)
        self.ids: List[str] = list(map(str, npz[id_key]))
        self.per_residue = list(npz[per_res_key])
        self.dim = int(npz["dim"]) if "dim" in npz.files else int(self.per_residue[0].shape[-1])
        self.include_ids = set(include_ids) if include_ids is not None else None
        self.tasks = TASKS

        for t in self.tasks:
            if t not in label_files or label_files[t] is None:
                raise ValueError(f"Missing label file for task '{t}'. Provided keys: {list(label_files.keys())}")
        self.label_maps = {t: self._load_label_txt(label_files[t]) for t in self.tasks}

        self.valid_indices: List[int] = []
        dropped = []
        for i, sid in enumerate(self.ids):
            if self.include_ids is not None and sid not in self.include_ids:
                continue
            L = self.per_residue[i].shape[0]
            mis = {}
            for t in self.tasks:
                arr = self.label_maps[t].get(sid)
                if arr is None:
                    mis[t] = -1
                elif len(arr) != L:
                    mis[t] = len(arr)
            if mis:
                dropped.append((sid, L, mis))
            else:
                self.valid_indices.append(i)

        if dropped:
            log(f"Dropped {len(dropped)} sequences with missing/length-mismatched labels "
                f"(showing up to 30): {dropped[:30]}")
        log(f"Dataset ready ({os.path.basename(npz_path)}). kept={len(self.valid_indices)} / "
            f"total={len(self.ids)}, dim={self.dim}")

    @staticmethod
    def _parse_label_line(line: str) -> Optional[List[int]]:
        line = line.strip()
        if line and set(line) <= set("01"):
            return [int(c) for c in line]
        return None

    def _load_label_txt(self, path: str) -> Dict[str, np.ndarray]:
        if not os.path.isfile(path):
            raise FileNotFoundError(f"label file not found: {path}")
        mp: Dict[str, np.ndarray] = {}
        with open(path, "r", encoding="utf-8") as f:
            for ln in f:
                ln = ln.rstrip("\n")
                if not ln.strip():
                    continue
                sid, payload = ln.split("\t", 1)
                arr = self._parse_label_line(payload)
                if arr is None:
                    raise ValueError(f"Could not parse labels in TSV line: {ln}")
                mp[sid] = np.asarray(arr, dtype=np.int64)
        return mp

    def __len__(self) -> int:
        return len(self.valid_indices)

    def __getitem__(self, idx: int):
        i = self.valid_indices[idx]
        sid = self.ids[i]
        H = np.asarray(self.per_residue[i], dtype=np.float32)
        L = H.shape[0]
        y = np.stack([self.label_maps[t][sid] for t in self.tasks], axis=1).astype(np.float32)
        return {"id": sid, "emb": torch.from_numpy(H), "label": torch.from_numpy(y), "len": L}


def collate_pad(batch: List[dict]):
    Lmax = max(b["len"] for b in batch)
    d = batch[0]["emb"].shape[-1]
    B = len(batch)
    emb = torch.zeros(B, Lmax, d, dtype=torch.float32)
    lab = torch.full((B, Lmax, len(TASKS)), fill_value=-100.0, dtype=torch.float32)
    mask = torch.zeros(B, Lmax, dtype=torch.bool)
    ids = []
    for i, b in enumerate(batch):
        L = b["len"]
        emb[i, :L] = b["emb"]
        lab[i, :L] = b["label"]
        mask[i, :L] = True
        ids.append(b["id"])
    return {"ids": ids, "emb": emb, "label": lab, "mask": mask}


def build_loader(npz_path: str, label_files: Dict[str, str], batch_size: int, shuffle: bool,
                 include_ids: Optional[Sequence[str]] = None, num_workers: int = 2,
                 seed: Optional[int] = None) -> DataLoader:
    ds = IDREmbDataset(npz_path, label_files, include_ids=set(include_ids) if include_ids is not None else None)
    gen = None
    if seed is not None:
        gen = torch.Generator()
        gen.manual_seed(seed)
    return DataLoader(ds, batch_size=batch_size, shuffle=shuffle, num_workers=num_workers,
                      collate_fn=collate_pad, pin_memory=torch.cuda.is_available(), generator=gen)
