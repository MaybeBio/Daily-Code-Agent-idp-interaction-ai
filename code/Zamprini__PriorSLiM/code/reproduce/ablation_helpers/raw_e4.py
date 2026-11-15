from pathlib import Path
import sys
HERE = Path(__file__).resolve().parent
import train_jwarm as t
import model_jwarm as baseline
import numpy as np
import pandas as pd
import torch
from torch import nn
import logging
common = t
digest = t.digest
require = t.require
WINDOW_LENGTH = 29

def read(path):
    return pd.read_csv(path, sep='\t', keep_default_na=False, low_memory=False)
from train_jwarm import digest
import scipy.sparse as sp
import json

def read(path):
    return json.loads(Path(path).read_text())

class RawStore:

    def __init__(self, directory, ids, expected_inputs):
        cfg_path = directory / 'motifae_sparse_config.json'
        self.config = cfg = read(cfg_path)
        require(cfg.get('complete') and cfg.get('dense_29_by_40960_constructed') is False, 'Invalid sparse cache')
        require((cfg['window_length'], cfg['latent_dim']) == (29, 40960), 'Wrong cache dimensions')
        require(digest(cfg_path) == expected_inputs['sparse_config']['sha256'], 'Not the E1 sparse cache')
        paths = {key: directory / Path(cfg[key]).name for key in ['window_metadata', 'residue_indptr', 'feature_indices', 'feature_values']}
        for key, path in paths.items():
            require(digest(path) == expected_inputs[f'sparse_{key}']['sha256'], f'Changed cache: {key}')
        self.metadata = pd.read_csv(paths['window_metadata'], sep='\t', low_memory=False).set_index('window_id', verify_integrity=True)
        self.ptr = np.load(paths['residue_indptr'], mmap_mode='r')
        self.indices = np.memmap(paths['feature_indices'], dtype=np.int32, mode='r')
        self.values = np.memmap(paths['feature_values'], dtype=np.float32, mode='r')
        require(len(self.indices) == len(self.values) == cfg['stored_nnz'], 'Sparse length mismatch')
        self.lookup = np.full(40960, -1, dtype=np.int64)
        self.lookup[ids] = np.arange(len(ids))
        self.aligned = {}

    def align(self, name, ids, meta):
        rows = self.metadata.loc[ids]
        require((rows.crop_mode == 'idr_coverage').all() and (rows.context_margin == 128).all(), 'Requires m128 IDR coverage')
        require(rows.sequence_sha256.notna().all(), 'Missing sequence identity')
        for key in ['window_start', 'window_end']:
            require(np.array_equal(rows[key].to_numpy(), meta[key].to_numpy()), f'Misaligned {key}')
        require(np.array_equal(rows.base_protein_id.astype(str), meta.base_protein_id.astype(str)), 'Protein mismatch')
        self.aligned[name] = rows.reset_index(drop=True)

    def events(self, split, rows, device):
        fi, val, res = ([], [], [])
        valid = np.zeros((len(rows), 29), bool)
        aligned = self.aligned[split]
        for b, row in enumerate(rows):
            m = aligned.iloc[int(row)]
            offset = int(m.sparse_window_index) * 29
            ptr = self.ptr[offset:offset + 30]
            require(len(ptr) == 30 and (np.diff(ptr) >= 0).all(), 'Invalid physical CSR slot')
            start, end = (int(ptr[0]), int(ptr[-1]))
            j = self.lookup[self.indices[start:end]]
            keep = j >= 0
            fi.append(j[keep])
            val.append(np.asarray(self.values[start:end])[keep])
            res.append(np.repeat(np.arange(b * 29, (b + 1) * 29), np.diff(ptr))[keep])
            valid[b, int(m.window_output_valid_start):int(m.window_output_valid_end)] = True
        require(valid.any(1).all(), 'Empty window')
        arrays = [np.concatenate(fi), np.concatenate(val), np.concatenate(res), valid]
        require(np.isfinite(arrays[1]).all() and (arrays[1] >= 0).all(), 'Invalid raw activation')
        return tuple((torch.tensor(a, dtype=t, device=device) for a, t in zip(arrays, [torch.long, torch.float32, torch.long, torch.bool])))

def csr_events(matrix, rows, device):
    x = matrix[rows].tocsr()
    return (torch.tensor(x.indices.copy(), dtype=torch.long, device=device), torch.tensor(x.data.copy(), dtype=torch.float32, device=device), torch.tensor(np.repeat(np.arange(len(rows)), np.diff(x.indptr)), dtype=torch.long, device=device), len(rows))
