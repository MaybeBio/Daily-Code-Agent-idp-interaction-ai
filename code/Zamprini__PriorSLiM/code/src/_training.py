from __future__ import annotations
import argparse
from contextlib import contextmanager
import fcntl
import gc
import hashlib
import json
import logging
import math
import os
from pathlib import Path
import random
import sys
import time
from typing import Iterable
os.environ.setdefault('PYTORCH_ENABLE_MPS_FALLBACK', '1')
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.special import expit
from sklearn.metrics import average_precision_score, roc_auc_score, brier_score_loss, precision_recall_curve, precision_score, recall_score, f1_score
from sklearn.preprocessing import StandardScaler
import torch
from torch.nn import functional as F
PROJECT_DIR = Path(os.environ.get('PRIORSLIM_DATA_ROOT', Path(__file__).resolve().parents[1] / 'data')).expanduser().resolve()
sys.path.insert(0, str(Path(__file__).resolve().parent))
import model as architecture
from model import JointModel, IntrinsicEncoder, AccessibilityResidual, SelectedRefinement, PriorSLiMPredictor
WINDOW_LENGTH = 29
RATIOS = (0.0, 0.001, 0.002, 0.003, 0.005, 0.0075, 0.01, 0.015, 0.02, 0.03, 0.04, 0.05, 0.075, 0.1, 0.15, 0.2, 0.3, 0.5, 0.75, 1.0)
LRS = (0.0001, 0.0003)
BETAS = (0.0, 0.01, 0.1, 1.0)
MAX_EPOCHS, PATIENCE, MIN_IMPROVEMENT = (30, 4, 0.0001)
c = sys.modules[__name__]
engine = c

def jsonable(value):
    if isinstance(value, dict):
        return {str(k): jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [jsonable(v) for v in value]
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, torch.Tensor):
        return jsonable(value.detach().cpu().tolist())
    if isinstance(value, np.ndarray):
        return jsonable(value.tolist())
    if isinstance(value, np.generic):
        return jsonable(value.item())
    if isinstance(value, float) and (not np.isfinite(value)):
        return None
    return value

def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + '.tmp')
    temp.write_text(json.dumps(jsonable(value), indent=2, ensure_ascii=False, allow_nan=False) + '\n')
    temp.replace(path)

def save_pt(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + '.tmp')
    torch.save(value, temp)
    temp.replace(path)

def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()

def identity(value):
    return hashlib.sha256(json.dumps(jsonable(value), sort_keys=True).encode()).hexdigest()

def seed_all(seed):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)

def device_for(name):
    if name == 'auto':
        name = 'cuda' if torch.cuda.is_available() else 'mps' if torch.backends.mps.is_available() else 'cpu'
    if name == 'mps' and (not torch.backends.mps.is_available()):
        raise RuntimeError('MPS unavailable; use --device cpu or another supported device')
    if name == 'cuda' and (not torch.cuda.is_available()):
        raise RuntimeError('CUDA unavailable')
    return torch.device(name)

def batches(n, size, seed=None):
    order = np.arange(n)
    if seed is not None:
        np.random.default_rng(seed).shuffle(order)
    for start in range(0, n, size):
        yield order[start:start + size]

def metrics(y, z, negative_types=None):
    y, z = (np.asarray(y, dtype=int), np.asarray(z, dtype=float))
    if not np.isfinite(z).all() or not set(np.unique(y)) <= {0, 1}:
        raise ValueError('Non-finite logits or non-binary labels')
    p = expit(z)
    result = dict(ap=float(average_precision_score(y, p)), auc=float(roc_auc_score(y, p)) if len(np.unique(y)) == 2 else None, brier=float(brier_score_loss(y, p)), n=len(y), positive_count=int(y.sum()))
    if negative_types is not None:
        kinds = np.asarray(negative_types).astype(str)
        result['ap_by_negative_type'] = {}
        for kind in sorted(set(kinds[y == 0])):
            take = (y == 1) | (y == 0) & (kinds == kind)
            result['ap_by_negative_type'][kind] = float(average_precision_score(y[take], p[take]))
    return result

def fit_preprocessing(data, train):
    arrays = data.arrays
    pssm = np.concatenate((arrays['pssm_residue'][train.tensor_rows], arrays['pssm_valid_mask'][train.tensor_rows, :, None]), axis=-1).reshape(len(train), -1)
    result = {}
    for name in ('pssm', 'pssm_window', 'chemical', 'accessibility'):
        x = pssm if name == 'pssm' else arrays[name][train.tensor_rows]
        if not np.isfinite(x).all():
            raise ValueError(f'Nonfinite values in TRAIN feature {name}; all values must be finite')
        scaler = StandardScaler().fit(x)
        result[name] = {'mean': torch.tensor(scaler.mean_, dtype=torch.float32), 'scale': torch.tensor(scaler.scale_, dtype=torch.float32)}
    return result

class CompactStore:

    def __init__(self, directory):
        d = Path(directory)
        self.ptr = np.load(d / 'indptr.npy', mmap_mode='r')
        self.column = np.memmap(d / 'columns.bin', dtype=np.int32, mode='r') if (d / 'columns.bin').stat().st_size else np.empty(0, np.int32)
        self.value = np.memmap(d / 'values.bin', dtype=np.float32, mode='r') if (d / 'values.bin').stat().st_size else np.empty(0, np.float32)

    def events(self, rows, device='cpu'):
        lengths = self.ptr[np.asarray(rows) + 1] - self.ptr[rows]
        cols = np.concatenate([self.column[self.ptr[r]:self.ptr[r + 1]] for r in rows]).copy() if len(rows) else np.empty(0, np.int32)
        vals = np.concatenate([self.value[self.ptr[r]:self.ptr[r + 1]] for r in rows]).copy() if len(rows) else np.empty(0, np.float32)
        owners = np.repeat(np.arange(len(rows)), lengths)
        return (torch.tensor(cols, dtype=torch.long, device=device), torch.tensor(vals, dtype=torch.float32, device=device), torch.tensor(owners, dtype=torch.long, device=device), len(rows))

def make_compact(store, split, ids, scales, directory, batch_size):
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    if (directory / 'complete.json').exists():
        return
    lookup = np.full(40960, -1, dtype=np.int64)
    lookup[ids.numpy()] = np.arange(len(ids))
    norm = scales.numpy()
    ptr = [0]
    with (directory / 'columns.bin').open('wb') as cf, (directory / 'values.bin').open('wb') as vf:
        for step, take in enumerate(batches(len(split), batch_size)):
            b = store.get_batch(split.sparse_rows[take])
            feature = lookup[b['feature_indices'].numpy()]
            residue = b['residue_indices'].numpy()
            values = b['feature_values'].numpy()
            if not np.isfinite(values).all() or (values < 0).any():
                raise ValueError('Invalid sparse activation')
            keep = (feature >= 0) & (values != 0)
            feature, residue, values = (feature[keep], residue[keep], values[keep])
            column = (feature * 29 + residue % 29).astype(np.int32)
            normalized = (values / (norm[feature] + 1e-08)).astype(np.float32)
            if not np.isfinite(normalized).all():
                raise ValueError('Invalid normalized sparse activation')
            column.tofile(cf)
            normalized.tofile(vf)
            counts = np.bincount(residue // 29, minlength=len(take))
            ptr.extend((ptr[-1] + np.cumsum(counts)).tolist())
            if step % 100 == 0:
                logging.info('Compact %s batch=%d rows=%d/%d', split.name, step, min((step + 1) * batch_size, len(split)), len(split))
    np.save(directory / 'indptr.npy', np.asarray(ptr, dtype=np.int64))
    write_json(directory / 'complete.json', dict(rows=len(split), nnz=ptr[-1], candidate_count=len(ids), normalized=True))

def build_evidence(store, weights, n, directory, batch_size):
    d = Path(directory)
    done = d / 'complete.json'
    weight_hash = hashlib.sha256(np.ascontiguousarray(weights).view(np.uint8)).hexdigest()
    if not done.exists():
        d.mkdir(parents=True, exist_ok=True)
        norm = np.linalg.norm(weights, axis=1)
        unit = np.divide(weights, norm[:, None], out=np.zeros_like(weights), where=norm[:, None] > 0)
        ptr = [0]
        with (d / 'columns.bin').open('wb') as cf, (d / 'values.bin').open('wb') as vf:
            for take in batches(n, batch_size):
                start, end = (int(store.ptr[take[0]]), int(store.ptr[take[-1] + 1]))
                col = np.asarray(store.column[start:end])
                val = np.asarray(store.value[start:end], dtype=np.float64) * unit.reshape(-1)[col]
                rows = np.repeat(np.arange(len(take)), np.diff(store.ptr[take[0]:take[-1] + 2]))
                block = sp.coo_matrix((val, (rows, col // 29)), shape=(len(take), len(weights))).tocsr()
                block.eliminate_zeros()
                block.indices.astype(np.int32).tofile(cf)
                block.data.astype(np.float64).tofile(vf)
                ptr.extend((ptr[-1] + block.indptr[1:]).tolist())
        np.save(d / 'indptr.npy', np.asarray(ptr, dtype=np.int64))
        write_json(done, dict(n=n, features=len(weights), weight_hash=weight_hash))
    elif json.loads(done.read_text())['weight_hash'] != weight_hash:
        raise RuntimeError('Stage-C evidence was generated with different weights')
    values = np.memmap(d / 'values.bin', dtype=np.float64, mode='r') if (d / 'values.bin').stat().st_size else np.empty(0, np.float64)
    cols = np.memmap(d / 'columns.bin', dtype=np.int32, mode='r') if (d / 'columns.bin').stat().st_size else np.empty(0, np.int32)
    ptr = np.load(d / 'indptr.npy', mmap_mode='r')
    return sp.csr_matrix((values, cols, ptr), shape=(n, len(weights)), copy=False)

@contextmanager
def exclusive(path):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('a+') as f:
        try:
            fcntl.flock(f, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise RuntimeError(f'Another process owns {path}; do not start duplicate jobs')
        try:
            yield
        finally:
            fcntl.flock(f, fcntl.LOCK_UN)

def cpu_tree(obj):
    if torch.is_tensor(obj):
        return obj.detach().cpu().clone()
    if isinstance(obj, dict):
        return {k: cpu_tree(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [cpu_tree(v) for v in obj]
    if isinstance(obj, tuple):
        return tuple((cpu_tree(v) for v in obj))
    return obj

def release(device):
    gc.collect()
    if device.type == 'mps':
        torch.mps.synchronize()
        torch.mps.empty_cache()
    elif device.type == 'cuda':
        torch.cuda.empty_cache()

def rng_state(device):
    state = dict(cpu=torch.get_rng_state())
    if device.type == 'mps':
        state['mps'] = torch.mps.get_rng_state()
    elif device.type == 'cuda':
        state['cuda'] = torch.cuda.get_rng_state_all()
    return state

def restore_rng(state, device):
    torch.set_rng_state(state['cpu'])
    if device.type == 'mps' and 'mps' in state:
        torch.mps.set_rng_state(state['mps'])
    elif device.type == 'cuda' and 'cuda' in state:
        torch.cuda.set_rng_state_all(state['cuda'])

def predict(model, forward, split, batch_size):
    model.eval()
    output = np.empty(len(split), np.float64)
    with torch.inference_mode():
        for take in batches(len(split), batch_size):
            z = forward(model, split.name, take)
            if not torch.isfinite(z).all():
                raise RuntimeError('Nonfinite forward')
            output[take] = z.cpu().numpy()
    return output

def fit_stage(model, forward, splits, directory, args, device, *, stage, epochs, patience, lr, wd, stage_signature):
    directory.mkdir(parents=True, exist_ok=True)
    best_path, last_path = (directory / 'best.pt', directory / 'last.pt')
    if last_path.exists() and (not args.resume):
        raise RuntimeError(f'Existing stage {directory}; use --resume or a new output root')
    optimizer = torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=wd)
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer, mode='max', factor=0.5, patience=3)
    start, stale, history, elapsed = (0, 0, [], 0.0)
    if last_path.exists():
        ckpt = load_pt(last_path)
        if ckpt['signature'] != stage_signature:
            raise RuntimeError(f'Stage configuration differs: {directory}')
        if ckpt['complete']:
            model.load_state_dict(ckpt['best_state'])
            if not best_path.exists():
                save_pt(best_path, dict(state=ckpt['best_state'], best_epoch=ckpt['best_epoch'], metrics=ckpt['best_metrics'], signature=stage_signature))
            logging.info('Reuse E1 %s best epoch=%d AP=%.6f', stage, ckpt['best_epoch'], ckpt['best_metrics']['ap'])
            return ckpt
        model.load_state_dict(ckpt['model'])
        optimizer.load_state_dict(ckpt['optimizer'])
        scheduler.load_state_dict(ckpt['scheduler'])
        restore_rng(ckpt['rng'], device)
        start, stale, history, elapsed = (ckpt['epoch'], ckpt['stale'], ckpt['history'], ckpt['seconds'])
        best_state, best_metric, best_epoch = (ckpt['best_state'], ckpt['best_metrics'], ckpt['best_epoch'])
    else:
        z = predict(model, forward, splits['val'], args.batch_size)
        best_metric = metrics(splits['val'].labels, z)
        best_state, best_epoch = (cpu_tree(model.state_dict()), 0)
        history.append(dict(epoch=0, train_bce=None, val_ap=best_metric['ap'], val_auc=best_metric['auc'], seconds=0.0))
    train = splits['train']
    for epoch in range(start + 1, epochs + 1):
        clock = time.monotonic()
        model.train()
        loss_total = 0.0
        for step, take in enumerate(batches(len(train), args.batch_size, args.seed + epoch)):
            optimizer.zero_grad(set_to_none=True)
            z = forward(model, 'train', take)
            y = torch.tensor(train.labels[take], dtype=torch.float32, device=device)
            loss = F.binary_cross_entropy_with_logits(z, y)
            if not torch.isfinite(loss):
                raise RuntimeError(f'{stage} loss is NaN/Inf')
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0, error_if_nonfinite=True)
            optimizer.step()
            loss_total += float(loss.detach().cpu()) * len(take)
            if step and step % 200 == 0:
                logging.info('%s epoch=%d step=%d/%d BCE=%.5f', stage, epoch, step, math.ceil(len(train) / args.batch_size), loss_total / min((step + 1) * args.batch_size, len(train)))
        z = predict(model, forward, splits['val'], args.batch_size)
        val = metrics(splits['val'].labels, z)
        scheduler.step(val['ap'])
        if val['ap'] > best_metric['ap'] + 1e-06:
            best_state, best_metric, best_epoch, stale = (cpu_tree(model.state_dict()), val, epoch, 0)
        else:
            stale += 1
        seconds = time.monotonic() - clock
        elapsed += seconds
        history.append(dict(epoch=epoch, train_bce=loss_total / len(train), val_ap=val['ap'], val_auc=val['auc'], val_brier=val['brier'], seconds=seconds))
        complete = stale >= patience or epoch == epochs
        ckpt = dict(signature=stage_signature, complete=complete, epoch=epoch, stale=stale, model=cpu_tree(model.state_dict()), optimizer=cpu_tree(optimizer.state_dict()), scheduler=scheduler.state_dict(), rng=rng_state(device), best_state=best_state, best_metrics=best_metric, best_epoch=best_epoch, history=history, seconds=elapsed, device=str(device))
        save_pt(last_path, ckpt)
        save_pt(best_path, dict(state=best_state, best_epoch=best_epoch, metrics=best_metric, signature=stage_signature))
        pd.DataFrame(history).to_csv(directory / 'history.tsv', sep='\t', index=False)
        logging.info('%s epoch=%d/%d BCE=%.5f valAP=%.6f valAUC=%.6f best=%d seconds=%.1f', stage, epoch, epochs, loss_total / len(train), val['ap'], val['auc'], best_epoch, seconds)
        release(device)
        if complete:
            break
    model.load_state_dict(best_state)
    return ckpt

def bce(y, z):
    return float(np.mean(np.logaddexp(0.0, z) - y * z))

def solve_proximal(X, y, offset, regularization, initial, *, max_iter, tol, callback=None):
    a = np.asarray(initial, dtype=np.float64).copy()
    if not np.isfinite(regularization) or regularization < 0 or (not np.isfinite(a).all()) or (a < 0).any():
        raise RuntimeError('Invalid Stage C regularization or initial coefficients')
    if not np.isfinite(X.data).all() or not np.isfinite(y).all() or (not np.isfinite(offset).all()):
        raise RuntimeError('Stage C nonfinite input')
    step, stable, last_support = (1.0, 0, a > 1e-10)
    z = offset + X @ a
    obj = bce(y, z) + regularization * a.sum()
    started = time.monotonic()
    rel_obj = rel_a = kkt = float('inf')
    converged = False
    for iteration in range(1, max_iter + 1):
        grad = np.asarray(X.T @ (expit(z) - y)).ravel() / len(y)
        trial = step * 1.2
        old_bce = bce(y, z)
        for _ in range(80):
            candidate = np.maximum(a - trial * (grad + regularization), 0.0)
            delta = candidate - a
            new_z = z + X @ delta
            loss = bce(y, new_z)
            bound = old_bce + grad @ delta + delta @ delta / (2.0 * trial)
            if not np.isfinite(loss):
                raise RuntimeError('Stage C nonfinite objective')
            if loss <= bound + 1e-12:
                break
            trial *= 0.5
        else:
            raise RuntimeError('Stage C line search failed')
        new_obj = loss + regularization * candidate.sum()
        if new_obj > obj + max(1e-08, abs(obj) * 1e-06):
            raise RuntimeError('Stage C objective diverged')
        rel_obj = abs(new_obj - obj) / max(abs(obj), 1e-12)
        rel_a = np.linalg.norm(delta) / max(np.linalg.norm(a), 1e-12)
        a, z, obj, step = (candidate, new_z, new_obj, trial)
        g = np.asarray(X.T @ (expit(z) - y)).ravel() / len(y) + regularization
        support = a > 1e-10
        kkt = float(np.max(np.where(support, np.abs(g), np.maximum(-g, 0.0))))
        stable = stable + 1 if np.array_equal(support, last_support) else 0
        last_support = support
        if kkt < tol and stable >= 10:
            converged = True
            break
        if callback and iteration % 100 == 0:
            callback(iteration, int(support.sum()), kkt)
    if not converged:
        logging.warning('Stage C max_iter reached: K=%d KKT=%.3g; unconverged solution excluded from frontier selection', int((a > 1e-10).sum()), kkt)
    return dict(alpha=a, iterations=iteration, converged=converged, kkt_max=kkt, relative_objective_change=rel_obj, relative_alpha_change=float(rel_a), objective=obj, seconds=time.monotonic() - started)

class SparseStore:

    def __init__(self, store_dir: Path | None=None):
        self.store_dir = (store_dir or PROJECT_DIR / 'feature_cache/motifae_sparse_windows_idr_coverage').resolve()
        config = json.loads((self.store_dir / 'motifae_sparse_config.json').read_text())
        if not config.get('complete') or config.get('dense_29_by_40960_constructed') is not False:
            raise RuntimeError('Sparse MotifAE cache is incomplete or has wrong provenance')
        if config.get('window_length') != 29 or config.get('latent_dim') != 40960:
            raise RuntimeError('Sparse MotifAE cache schema is incompatible')

        def source(key: str) -> Path:
            original = Path(config[key])
            local = self.store_dir / original.name
            return local if local.is_file() else original
        self.paths = {key: source(key) for key in ('window_metadata', 'residue_indptr', 'feature_indices', 'feature_values')}
        self.config = config
        self.metadata = pd.read_csv(source('window_metadata'), sep='\t', low_memory=False)
        self.metadata = self.metadata.sort_values('window_index').reset_index(drop=True)
        if not np.array_equal(self.metadata['window_index'].to_numpy(), np.arange(len(self.metadata))):
            raise RuntimeError('MotifAE window_index must be contiguous')
        if self.metadata['window_id'].duplicated().any():
            raise RuntimeError('MotifAE window_id must be unique')
        self.row_by_id = pd.Series(self.metadata['window_index'].to_numpy(), index=self.metadata['window_id'].astype(str))
        self.indptr = np.load(source('residue_indptr'), mmap_mode='r')
        nnz = int(config['stored_nnz'])
        self.indices = np.memmap(source('feature_indices'), dtype=np.int32, mode='r', shape=(nnz,))
        self.values = np.memmap(source('feature_values'), dtype=np.float32, mode='r', shape=(nnz,))

    def rows_for_ids(self, window_ids: Iterable[str]) -> np.ndarray:
        rows = self.row_by_id.reindex(pd.Index(window_ids, dtype=str))
        if rows.isna().any():
            raise RuntimeError(f'MotifAE cache misses {int(rows.isna().sum())} window IDs')
        return rows.to_numpy(dtype=np.int64)

    def get_batch(self, rows: np.ndarray) -> dict[str, torch.Tensor]:
        rows = np.asarray(rows, dtype=np.int64)
        batch_size = len(rows)
        indices_parts, values_parts, residues_parts = ([], [], [])
        valid_mask = np.zeros((batch_size, WINDOW_LENGTH), dtype=np.bool_)
        for batch_index, row_number in enumerate(rows):
            row = self.metadata.iloc[int(row_number)]
            sparse_index = int(row['sparse_window_index'])
            offset = sparse_index * WINDOW_LENGTH
            pointers = np.asarray(self.indptr[offset:offset + WINDOW_LENGTH + 1], dtype=np.int64)
            counts = np.diff(pointers)
            if np.any(counts < 0):
                raise RuntimeError('MotifAE CSR pointers are non-monotonic')
            start, end = (int(pointers[0]), int(pointers[-1]))
            indices_parts.append(np.asarray(self.indices[start:end], dtype=np.int64))
            values_parts.append(np.asarray(self.values[start:end], dtype=np.float32))
            residues_parts.append(np.repeat(np.arange(batch_index * WINDOW_LENGTH, (batch_index + 1) * WINDOW_LENGTH), counts))
            valid_mask[batch_index, int(row['window_output_valid_start']):int(row['window_output_valid_end'])] = True
        concatenate = lambda parts, dtype: np.concatenate(parts).astype(dtype, copy=False) if parts else np.empty(0, dtype=dtype)
        return {'feature_indices': torch.from_numpy(concatenate(indices_parts, np.int64)), 'feature_values': torch.from_numpy(concatenate(values_parts, np.float32)), 'residue_indices': torch.from_numpy(concatenate(residues_parts, np.int64)), 'valid_mask': torch.from_numpy(valid_mask)}

class SplitData:

    def __init__(self, name: str, metadata: pd.DataFrame, tensor_rows: np.ndarray, sparse_rows: np.ndarray):
        self.name = name
        self.metadata = metadata.reset_index(drop=True)
        self.tensor_rows = np.asarray(tensor_rows, dtype=np.int64)
        self.sparse_rows = np.asarray(sparse_rows, dtype=np.int64)
        self.labels = self.metadata['label'].to_numpy(dtype=np.float32)
        self.negative_types = self.metadata['negative_type'].fillna('').astype(str).to_numpy()

    def __len__(self) -> int:
        return len(self.tensor_rows)

class FeatureData:
    KEYS = ('pssm_residue', 'pssm_valid_mask', 'pssm_window', 'chemical', 'accessibility')

    def __init__(self, store: SparseStore, tensor_path=None, metadata_path=None):
        self.store = store
        archive_path = tensor_path or PROJECT_DIR / 'train/combined_tensor.npz'
        with np.load(archive_path) as archive:
            self.arrays = {name: np.asarray(archive[name]) for name in self.KEYS}
            self.window_ids = np.asarray(archive['window_ids']).astype(str)
        self.combined = pd.read_csv(metadata_path or PROJECT_DIR / 'train/combined_metadata.tsv', sep='\t', low_memory=False)
        if len(self.window_ids) != len(self.combined) or not np.array_equal(self.window_ids, self.combined['window_id'].astype(str).to_numpy()):
            raise RuntimeError('Combined feature archive and metadata are not aligned')

    def split(self, name: str) -> SplitData:
        if name not in ('train', 'val', 'test'):
            raise ValueError('Unknown split')
        frame = self.combined[self.combined['split'] == name].sort_values('split_row')
        provenance = pd.read_csv(PROJECT_DIR / f'new_data/{name}_metadata.tsv', sep='\t', low_memory=False)
        windows = pd.read_csv(PROJECT_DIR / f'new_data/{name}.tsv', sep='\t', low_memory=False)
        if len(frame) != len(provenance) or len(frame) != len(windows):
            raise RuntimeError(f'{name} row counts differ across inputs')
        for feature_column, metadata_column in (('uniprot', 'protein_id'), ('window_start', 'window_start'), ('window_end', 'window_end'), ('sequence', 'window_sequence')):
            if not np.array_equal(frame[feature_column].astype(str).to_numpy(), provenance[metadata_column].astype(str).to_numpy()):
                raise RuntimeError(f'{name}: {feature_column} differs from provenance')
        for column in ('uniprot', 'window_start', 'window_end', 'sequence'):
            if not np.array_equal(frame[column].astype(str).to_numpy(), windows[column].astype(str).to_numpy()):
                raise RuntimeError(f'{name}: window table {column} differs from features')
        if not np.allclose(frame['label'].to_numpy(dtype=float), windows['label'].to_numpy(dtype=float)):
            raise RuntimeError(f'{name}: window table labels differ from features')
        if not np.allclose(frame['label'].to_numpy(dtype=float), provenance['label'].to_numpy(dtype=float)):
            raise RuntimeError(f'{name}: labels differ from provenance')
        rows = frame['tensor_index'].to_numpy(dtype=np.int64)
        return SplitData(name, provenance, rows, self.store.rows_for_ids(frame['window_id'].astype(str)))

    def raw_batch(self, tensor_rows: np.ndarray, device: torch.device):
        return tuple((torch.from_numpy(np.array(self.arrays[name][tensor_rows])).to(device=device, dtype=torch.float32) for name in self.KEYS))

def require(ok, message):
    if not ok:
        raise RuntimeError(message)

def fit_joint(a, inputs, out, signature, device):
    t = inputs.training
    epochs = t['profile_epochs']
    patience = t['residual_patience']
    model = inputs.model(a, device, warm=True)
    require(torch.count_nonzero(model.profile.weights) == 0, 'W must initialize at zero')
    optimizer = torch.optim.AdamW(model.parameter_groups(t))
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer, mode='max', factor=0.5, patience=3)
    directory = out / 'joint'
    directory.mkdir(parents=True, exist_ok=True)
    last = directory / 'last.pt'
    best_path = directory / 'best.pt'
    c.seed_all(a.seed + 200)
    start = stale = 0
    elapsed = 0.0
    history = []
    if last.exists():
        require(a.resume, 'Use --resume to continue')
        ck = c.load_pt(last)
        require(ck['signature'] == signature, 'Training signature changed')
        if ck['complete']:
            model.load_state_dict(ck['best_state'])
            return (model, ck)
        model.load_state_dict(ck['model'])
        optimizer.load_state_dict(ck['optimizer'])
        scheduler.load_state_dict(ck['scheduler'])
        engine.restore_rng(ck['rng'], device)
        start, stale, elapsed, history = (ck['epoch'], ck['stale'], ck['seconds'], ck['history'])
        best_state, best, best_epoch = (ck['best_state'], ck['best_metrics'], ck['best_epoch'])
    else:
        take = np.arange(min(4, len(inputs.splits['val'])))
        model.eval()
        with torch.no_grad():
            z = model(*inputs.batch('val', take, device)).cpu().numpy()
        require(np.allclose(z, inputs.explicit['val'][take], atol=1e-05, rtol=1e-05), 'Warm start differs from pretrained explicit model')
        val = engine.predict(model, lambda m, s, rows: m(*inputs.batch(s, rows, device)), inputs.splits['val'], t['batch_size'])
        best = c.metrics(inputs.splits['val'].labels, val)
        best_state = engine.cpu_tree(model.state_dict())
        best_epoch = 0
        history = [dict(epoch=0, val_ap=best['ap'], val_auc=best['auc'], val_brier=best['brier'])]
    for epoch in range(start + 1, epochs + 1):
        clock = time.monotonic()
        model.train()
        loss_sum = 0.0
        for step, take in enumerate(c.batches(len(inputs.splits['train']), t['batch_size'], a.seed + epoch)):
            optimizer.zero_grad(set_to_none=True)
            z = model(*inputs.batch('train', take, device))
            y = torch.tensor(inputs.splits['train'].labels[take], device=device, dtype=torch.float32)
            loss = F.binary_cross_entropy_with_logits(z, y)
            require(torch.isfinite(loss), 'Nonfinite training loss')
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0, error_if_nonfinite=True)
            optimizer.step()
            loss_sum += float(loss.detach().cpu()) * len(take)
            if step and step % 200 == 0:
                logging.info('%s epoch=%d step=%d', a.arm, epoch, step)
        val = engine.predict(model, lambda m, s, rows: m(*inputs.batch(s, rows, device)), inputs.splits['val'], t['batch_size'])
        m = c.metrics(inputs.splits['val'].labels, val)
        scheduler.step(m['ap'])
        if m['ap'] > best['ap'] + 1e-06:
            best, best_epoch, best_state, stale = (m, epoch, engine.cpu_tree(model.state_dict()), 0)
        else:
            stale += 1
        elapsed += time.monotonic() - clock
        history.append(dict(epoch=epoch, train_bce=loss_sum / len(inputs.splits['train']), val_ap=m['ap'], val_auc=m['auc'], val_brier=m['brier'], seconds=elapsed, lrs=[g['lr'] for g in optimizer.param_groups]))
        ck = dict(signature=signature, model=engine.cpu_tree(model.state_dict()), optimizer=engine.cpu_tree(optimizer.state_dict()), scheduler=scheduler.state_dict(), rng=engine.rng_state(device), epoch=epoch, stale=stale, best_state=best_state, best_metrics=best, best_epoch=best_epoch, seconds=elapsed, history=history, complete=stale >= patience or epoch == epochs)
        c.save_pt(last, ck)
        c.save_pt(best_path, dict(state=best_state, metrics=best, best_epoch=best_epoch, signature=signature))
        pd.DataFrame(history).to_csv(directory / 'history.tsv', sep='\t', index=False)
        logging.info('%s epoch=%d/%d BCE=%.5f valAP=%.6f best=%d', a.arm, epoch, epochs, history[-1]['train_bce'], m['ap'], best_epoch)
        engine.release(device)
        if ck['complete']:
            break
    model.load_state_dict(best_state)
    return (model, ck)

def export_full(a, inputs, out, model, ck, signature, device):
    result = predictions(model, inputs, device, inputs.training['batch_size'])
    W = model.profile.weights.detach().cpu().numpy().astype(np.float64)
    alpha = np.linalg.norm(W, axis=1)
    c.save_pt(out / 'frozen_profiles.pt', dict(weights=torch.from_numpy(W), alpha_init=torch.from_numpy(alpha), feature_ids=inputs.ids, scales=inputs.scales))
    report = {}
    for s, (z, explicit, local) in result.items():
        np.save(out / f'explicit_{s}.npy', explicit)
        replay = explicit.copy()
        store = inputs.compact[s]
        for take in c.batches(len(explicit), inputs.training['batch_size']):
            col, val, owner, n = store.events(take)
            replay[take] += np.bincount(owner.numpy(), weights=W.reshape(-1)[col.numpy()] * val.numpy(), minlength=n)
        err = float(np.max(abs(z - replay)))
        require(np.allclose(z, replay, atol=1e-05, rtol=1e-05), 'Full CPU replay disagrees with forward')
        split = inputs.splits[s]
        m = c.metrics(split.labels, replay, split.negative_types)
        report[s] = dict(full=m, explicit_only=c.metrics(split.labels, explicit, split.negative_types), replay_max_error=err)
        if s == 'val':
            (out / 'predictions').mkdir(exist_ok=True)
            np.savez_compressed(out / 'predictions/full_replay.npz', logits=replay, labels=split.labels, window_ids=split.window_ids.astype(str))
    require(abs(report['val']['full']['ap'] - ck['best_metrics']['ap']) < 0.0001, 'Best AP replay mismatch')
    c.write_json(out / 'full_metrics.json', dict(arm=a.arm, seed=a.seed, candidate_count=5845, K=int((alpha > 1e-10).sum()), best_epoch=ck['best_epoch'], trained_epochs=ck['epoch'], training_seconds=ck['seconds'], metrics=report, inherited_explicit_training_seconds=read(inputs.source / 'explicit_info.json')['seconds'], signature=signature, test_evaluated=False, stage_D_used=False))
    c.write_json(out / 'train_complete.json', dict(signature=signature))

def point(out, pid, alpha, explicit, design, inputs, solver=None, ratio=0.0):
    result = dict(point_id=pid, K=int((alpha > 1e-10).sum()), lambda_ratio=ratio, eligible=True if solver is None else bool(solver['converged']), converged=None if solver is None else bool(solver['converged']), test_evaluated=False)
    if solver:
        result['solver'] = {k: v for k, v in solver.items() if k != 'alpha'}
    for s, split in inputs.splits.items():
        z = explicit[s] + design[s] @ alpha
        result[s] = c.metrics(split.labels, z, split.negative_types)
        if s == 'val':
            (out / 'predictions').mkdir(exist_ok=True)
            np.savez_compressed(out / f'predictions/{pid}.npz', logits=z, labels=split.labels, window_ids=split.window_ids.astype(str))
    c.write_json(out / f'points/{pid}.json', result)
    return result

def compress(a, inputs, out, signature):
    require(read(out / 'train_complete.json')['signature'] == signature, 'Finish joint training first')
    frozen = c.load_pt(out / 'frozen_profiles.pt')
    W = frozen['weights'].numpy()
    alpha = frozen['alpha_init'].numpy()
    explicit = {s: np.load(out / f'explicit_{s}.npy') for s in inputs.splits}
    design = {s: c.build_evidence(inputs.compact[s], W, len(split), out / 'evidence' / s, inputs.training['batch_size']) for s, split in inputs.splits.items()}
    with np.load(out / 'predictions/full_replay.npz') as saved:
        require(np.array_equal(saved['window_ids'], inputs.splits['val'].window_ids), 'Replay IDs differ')
        require(np.allclose(explicit['val'] + design['val'] @ alpha, saved['logits'], atol=1e-08, rtol=0), 'Unit-profile algebra replay failed')
    point(out, 'full_replay', alpha, explicit, design, inputs)
    X = design['train']
    y = inputs.splits['train'].labels.astype(float)
    offset = explicit['train']
    lambda_max = float(np.maximum(-np.asarray(X.T @ (expit(offset) - y)).ravel() / len(y), 0).max())
    for ratio in inputs.training['lambda_ratios']:
        if ratio == 0:
            continue
        pid = 'lambda_' + format(ratio, '.10g').replace('.', 'p')
        path = out / 'solutions' / f'{pid}.pt'
        if path.exists():
            require(a.resume, 'Use --resume')
            saved = c.load_pt(path)
            require(saved['signature'] == signature and saved['ratio'] == ratio, 'Changed proximal solution contract')
            solved = dict(saved['solver'], alpha=saved['alpha'].numpy())
        else:
            if lambda_max <= 0:
                solved = dict(alpha=np.zeros_like(alpha), iterations=0, converged=True, kkt_max=0.0, relative_objective_change=0.0, relative_alpha_change=0.0, objective=engine.bce(y, offset), seconds=0.0)
            else:
                solved = engine.solve_proximal(X, y, offset, ratio * lambda_max, alpha, max_iter=a.solver_max_iter, tol=inputs.training['solver_tol'], callback=lambda i, k, r: logging.info('ratio=%g iter=%d K=%d KKT=%.3g', ratio, i, k, r))
            c.save_pt(path, dict(signature=signature, ratio=ratio, lambda_max=lambda_max, alpha=torch.from_numpy(solved['alpha']), solver=c.jsonable({k: v for k, v in solved.items() if k != 'alpha'})))
        alpha = solved['alpha'].copy()
        result = point(out, pid, alpha, explicit, design, inputs, solved, ratio)
        logging.info('%s ratio=%g converged=%s K=%d valAP=%.6f', a.arm, ratio, result['eligible'], result['K'], result['val']['ap'])
    c.write_json(out / 'compress_complete.json', dict(signature=signature, test_evaluated=False))

def predict_profile(profile, inputs, explicit, split, batch_size, device):
    profile.eval()
    z = np.empty(len(inputs.splits[split]), dtype=np.float64)
    weights = (profile.W0.detach().cpu() + profile.DeltaW.detach().cpu()).numpy().astype(np.float64)
    lookup = profile.lookup.cpu().numpy()
    for take in c.batches(len(z), batch_size):
        col, values, owner, n = inputs.compact[split].events(take, 'cpu')
        col, values, owner = (col.numpy(), values.numpy(), owner.numpy())
        j = lookup[col // 29]
        keep = j >= 0
        delta = np.bincount(owner[keep], weights=weights[j[keep], col[keep] % 29] * values[keep], minlength=n)
        z[take] = explicit[split][take] + delta
    require(np.isfinite(z).all(), 'Nonfinite prediction')
    s = inputs.splits[split]
    return (c.metrics(s.labels, z, s.negative_types), z)

def fit_config(a, inputs, explicit, selected, W0, initial_metrics, signature, out, lr, beta, device):
    directory = out / 'stage_d' / f'lr{lr:g}_beta{beta:g}'
    last = directory / 'last.pt'
    model = SelectedRefinement(selected, W0).to(device)
    require([n for n, p in model.named_parameters() if p.requires_grad] == ['DeltaW'], 'Only DeltaW may train')
    optimizer = torch.optim.AdamW([model.DeltaW], lr=lr, weight_decay=0.0)
    best = dict(metrics=initial_metrics, epoch=0, DeltaW=torch.zeros_like(W0), lr=lr, beta=beta, drift=0.0)
    history = [dict(epoch=0, val_ap=initial_metrics['ap'], drift=0.0)]
    start = stale = 0
    if last.exists():
        require(a.resume, 'Existing refinement; use --resume')
        ck = load_refinement_checkpoint(last)
        require(ck['signature'] == signature, 'Resume signature mismatch')
        if ck['complete']:
            logging.info('Reuse finished lr=%g beta=%g valAP=%.6f', lr, beta, ck['best']['metrics']['ap'])
            return ck['best']
        model.load_state_dict(ck['model'])
        optimizer.load_state_dict(ck['optimizer'])
        start, stale, best, history = (ck['epoch'], ck['stale'], ck['best'], ck['history'])
    c.seed_all(a.seed)
    batch_size = inputs.training['batch_size']
    for epoch in range(start + 1, MAX_EPOCHS + 1):
        model.train()
        loss_sum = 0.0
        for take in c.batches(len(inputs.splits['train']), batch_size, a.seed + epoch):
            events = inputs.compact['train'].events(take, device)
            z0 = torch.tensor(explicit['train'][take], dtype=torch.float32, device=device)
            y = torch.tensor(inputs.splits['train'].labels[take], dtype=torch.float32, device=device)
            optimizer.zero_grad(set_to_none=True)
            loss = F.binary_cross_entropy_with_logits(z0 + model(events), y) + beta * model.drift_loss()
            require(torch.isfinite(loss), 'Nonfinite loss')
            loss.backward()
            require(torch.isfinite(model.DeltaW.grad).all(), 'Nonfinite gradient')
            optimizer.step()
            loss_sum += float(loss.detach().cpu()) * len(take)
        metrics, _ = predict_profile(model, inputs, explicit, 'val', batch_size, device)
        drift = float(model.drift_loss().detach().cpu())
        history.append(dict(epoch=epoch, train_loss=loss_sum / len(inputs.splits['train']), val_ap=metrics['ap'], val_auc=metrics['auc'], drift=drift))
        if metrics['ap'] > best['metrics']['ap'] + MIN_IMPROVEMENT:
            best = dict(metrics=metrics, epoch=epoch, DeltaW=model.DeltaW.detach().cpu().clone(), lr=lr, beta=beta, drift=drift)
            stale = 0
        else:
            stale += 1
        complete = stale >= PATIENCE or epoch == MAX_EPOCHS
        c.save_pt(last, dict(signature=signature, model=engine.cpu_tree(model.state_dict()), optimizer=engine.cpu_tree(optimizer.state_dict()), epoch=epoch, stale=stale, best=best, history=history, complete=complete))
        c.save_pt(directory / 'best.pt', dict(signature=signature, **best))
        pd.DataFrame(history).to_csv(directory / 'history.tsv', sep='\t', index=False)
        logging.info('seed=%d K=%d lr=%g beta=%g epoch=%d BCE+drift=%.6f valAP=%.6f best=%d', a.seed, len(selected), lr, beta, epoch, history[-1]['train_loss'], metrics['ap'], best['epoch'])
        engine.release(device)
        if complete:
            break
    return best

def threshold_metrics(y, z, threshold):
    pred = expit(z) >= threshold
    return dict(threshold=float(threshold), precision=float(precision_score(y, pred, zero_division=0)), recall=float(recall_score(y, pred, zero_division=0)), f1=float(f1_score(y, pred, zero_division=0)))

def val_threshold(y, z):
    precision, recall, thresholds = precision_recall_curve(y, expit(z))
    f1 = 2 * precision[:-1] * recall[:-1] / np.maximum(precision[:-1] + recall[:-1], 1e-12)
    return float(thresholds[np.flatnonzero(f1 == f1.max())[-1]])

def save_predictions(out, name, split, z, explicit):
    np.savez_compressed(out / f'{name}_predictions.npz', window_ids=split.window_ids.astype(str), labels=split.labels, logits=z, probabilities=expit(z), explicit_logits=explicit, motifae_logits=z - explicit, negative_types=split.negative_types.astype(str))
    pd.DataFrame(dict(window_id=split.window_ids, label=split.labels, negative_type=split.negative_types, logit=z, probability=expit(z), explicit_logit=explicit, motifae_logit=z - explicit)).to_csv(out / f'{name}_predictions.tsv', sep='\t', index=False)

def read(path):
    return json.loads(Path(path).read_text())

def load_pt(path):
    allowed = [np.core.multiarray.scalar, np.dtype] + [type(np.dtype(x)) for x in ('float64', 'float32', 'int64', 'int32', 'bool', 'str_')]
    with torch.serialization.safe_globals(allowed):
        return torch.load(path, map_location='cpu', weights_only=True)
load_refinement_checkpoint = load_pt

def validate_alignment(data, store, split):
    split.window_ids = data.window_ids[split.tensor_rows]
    sm = store.metadata.iloc[split.sparse_rows]
    frame = data.combined.iloc[split.tensor_rows]
    for left, right in [('source_protein_id', 'uniprot'), ('window_start', 'window_start'), ('window_end', 'window_end')]:
        require(np.array_equal(sm[left].astype(str), frame[right].astype(str)), f'{split.name}: sparse identity/coordinates differ')
    require(set(sm.crop_mode) == {'idr_coverage'} and set(sm.context_margin) == {128}, 'Requires contextual IDR-coverage m=128')
    require(not sm.sequence_sha256.isna().any(), 'Missing sequence provenance')
    require(len(set(split.window_ids)) == len(split), 'Duplicate window IDs')

def input_fingerprint(a, store):
    paths = {'tensor': a.tensor, 'metadata': a.metadata, 'manifest': a.feature_manifest, 'assets': a.assets, 'sparse_config': a.store_dir / 'motifae_sparse_config.json'}
    paths.update({f'sparse_{k}': v for k, v in store.paths.items()})
    for s in ['train', 'val']:
        for suffix in ['.tsv', '_metadata.tsv']:
            paths[s + suffix] = PROJECT_DIR / 'new_data' / f'{s}{suffix}'
    return {k: dict(sha256=digest(p), bytes=p.stat().st_size) for k, p in paths.items()}

class Inputs:

    def __init__(self, a):
        self.source = a.out_dir / 'stage_a'
        self.raw = SparseStore(a.store_dir)
        self.data = FeatureData(self.raw, a.tensor, a.metadata)
        self.feature_manifest = read(a.feature_manifest)
        self.splits = {s: self.data.split(s) for s in ['train', 'val']}
        for split in self.splits.values():
            validate_alignment(self.data, self.raw, split)
        train_proteins = set(self.raw.metadata.iloc[self.splits['train'].sparse_rows].base_protein_id)
        require(not train_proteins & set(self.raw.metadata.iloc[self.splits['val'].sparse_rows].base_protein_id), 'Train/Validation protein leakage')
        require(self.data.arrays['pssm_residue'].shape[1:] == (29, 6) and self.data.arrays['chemical'].shape[1:] == (58,), 'Feature dimensions differ')
        names = self.feature_manifest['chemical_features']
        require([names[i] for i in architecture.MOTIF_INDICES] == ['motif_min_abs_offset', 'motif_center_enrichment', 'motif_count_log'], 'Motif columns differ')
        self.fingerprint = input_fingerprint(a, self.raw)
        assets = load_pt(a.assets)
        require(assets['format'] == 'jwarm_training_assets_v1', 'Wrong standalone training assets')
        require(assets['sparse_config_sha256'] == self.fingerprint['sparse_config']['sha256'], 'P95 cache provenance does not match the current inputs')
        for s in ['train.tsv', 'train_metadata.tsv']:
            require(assets['train_hashes'][s] == self.fingerprint[s]['sha256'], 'TRAIN-derived candidate/scales provenance differs')
        self.ids, self.scales = (assets['feature_ids'], assets['train_p95_scales'])
        require(len(self.ids) == 5845 and len(torch.unique(self.ids)) == 5845 and (self.scales > 0).all(), 'Wrong candidate universe/scales')
        prepared = a.out_dir / '_prepared'
        signature = identity(self.fingerprint)
        if (prepared / 'input_manifest.json').exists():
            require(read(prepared / 'input_manifest.json')['signature'] == signature, 'Inputs changed; use new output directory')
        else:
            write_json(prepared / 'input_manifest.json', dict(signature=signature, inputs=self.fingerprint))
        if not (prepared / 'preprocessing.pt').exists():
            save_pt(prepared / 'preprocessing.pt', fit_preprocessing(self.data, self.splits['train']))
        self.scalers = load_pt(prepared / 'preprocessing.pt')
        self.compact = {}
        for s, split in self.splits.items():
            make_compact(self.raw, split, self.ids, self.scales, prepared / s, a.batch_size)
            self.compact[s] = CompactStore(prepared / s)
        self.training = {k: getattr(a, k) for k in ['batch_size', 'base_epochs', 'access_epochs', 'profile_epochs', 'base_lr', 'access_lr', 'profile_lr', 'base_wd', 'profile_wd', 'base_patience', 'residual_patience', 'motif_dropout', 'solver_tol']}

    def batch(self, s, take, device):
        return (self.data.raw_batch(self.splits[s].tensor_rows[take], device), self.compact[s].events(take, device))

    def model(self, a, device, warm=True):
        model = JointModel(self.scalers, len(self.ids), a.seed, a.motif_dropout)
        if warm:
            model.load_explicit(load_pt(self.source / 'intrinsic/best.pt')['state'], load_pt(self.source / 'access/best.pt')['state'])
        return model.to(device)

def run_stage_a(a, inputs, signature, device):
    data, splits, scalers = (inputs.data, inputs.splits, inputs.scalers)
    seed_all(a.seed)
    intrinsic = IntrinsicEncoder(scalers, a.motif_dropout).to(device)
    forward = lambda m, s, t: m(*data.raw_batch(splits[s].tensor_rows[t], device)[:4])
    ck = fit_stage(intrinsic, forward, splits, inputs.source / 'intrinsic', a, device, stage='A1 Intrinsic', epochs=a.base_epochs, patience=a.base_patience, lr=a.base_lr, wd=a.base_wd, stage_signature=identity([signature, 'intrinsic']))
    z = {s: predict(intrinsic, forward, split, a.batch_size) for s, split in splits.items()}
    seconds = ck['seconds']
    del intrinsic
    release(device)
    seed_all(a.seed + 100)
    access = AccessibilityResidual().to(device)
    arrays = {s: (data.arrays['accessibility'][split.tensor_rows] - scalers['accessibility']['mean'].numpy()) / scalers['accessibility']['scale'].numpy() for s, split in splits.items()}
    forward_access = lambda m, s, t: torch.tensor(z[s][t], dtype=torch.float32, device=device) + m(torch.tensor(arrays[s][t], dtype=torch.float32, device=device))
    ck = fit_stage(access, forward_access, splits, inputs.source / 'access', a, device, stage='A2 Accessibility', epochs=a.access_epochs, patience=a.residual_patience, lr=a.access_lr, wd=a.base_wd, stage_signature=identity([signature, 'access', digest(inputs.source / 'intrinsic/best.pt')]))
    inputs.explicit = {s: predict(access, forward_access, split, a.batch_size) for s, split in splits.items()}
    for s, zs in inputs.explicit.items():
        np.save(inputs.source / f'explicit_{s}.npy', zs)
    write_json(inputs.source / 'explicit_info.json', dict(seconds=seconds + ck['seconds'], signature=signature))
    del access
    release(device)

def predictions(model, inputs, device, batch_size):
    model.eval()
    result = {}
    with torch.inference_mode():
        for s, split in inputs.splits.items():
            parts = [np.empty(len(split), np.float64) for _ in range(3)]
            for take in batches(len(split), batch_size):
                values = model(*inputs.batch(s, take, device), return_parts=True)
                for dest, z in zip(parts, values):
                    require(torch.isfinite(z).all(), 'Nonfinite prediction')
                    dest[take] = z.cpu().numpy()
            result[s] = parts
    return result

def portable_bundle(locked, state, inputs, seed):
    return dict(format='jwarm_bundle_v1', seed=seed, K=len(locked['feature_ids']), explicit_state={k: v for k, v in state.items() if not k.startswith('profile.')}, preprocessing=inputs.scalers, feature_manifest=inputs.feature_manifest, feature_ids=locked['feature_ids'], train_p95_scales=locked['train_p95_scales'], W_final=locked['W_final'], threshold=locked['threshold'], motif_dropout=inputs.training['motif_dropout'], crop_protocol=dict(mode='idr_coverage', context_margin=128, window_length=29, esm_model='esm2_t33_650M_UR50D', layer=33), selection=dict(K='manual Validation-only operating point', refinement='Validation AP', test_used_for_selection=False))

def run_refine(a, inputs, signature, device):
    out = a.out_dir
    require(read(out / 'compress_complete.json')['signature'] == signature, 'Run frontier first')
    points = [read(p) for p in (out / 'points').glob('lambda_*.json')]
    points = [p for p in points if p['eligible'] and p['converged'] and (p['K'] == a.selected_k)]
    require(bool(points), 'K is not an observed converged support; inspect frontier.tsv')
    point0 = max(points, key=lambda p: (p['val']['ap'], p['lambda_ratio']))
    frozen = load_pt(out / 'frozen_profiles.pt')
    solution = load_pt(out / 'solutions' / f"{point0['point_id']}.pt")
    W = frozen['weights'].double()
    alpha = solution['alpha'].double()
    selected = torch.nonzero(alpha > 1e-10).flatten()
    norm = torch.linalg.vector_norm(W, dim=1)
    W0 = (W[selected] / norm[selected, None] * alpha[selected, None]).float()
    directory = out / f'K{a.selected_k}'
    directory.mkdir(parents=True, exist_ok=True)
    dsig = identity([signature, point0, digest(out / 'solutions' / f"{point0['point_id']}.pt"), LRS, BETAS])
    explicit = {s: np.load(out / f'explicit_{s}.npy') for s in inputs.splits}
    profile = SelectedRefinement(selected, W0).to(device)
    initial, z0 = predict_profile(profile, inputs, explicit, 'val', a.batch_size, device)
    with np.load(out / 'predictions' / f"{point0['point_id']}.npz") as old:
        require(np.array_equal(old['window_ids'], inputs.splits['val'].window_ids), 'Window order changed')
        require(float(np.max(abs(old['logits'] - z0))) < 1e-05, 'Stage C epoch0 replay failed')
    final = directory / 'final_best.pt'
    if final.exists():
        require(a.resume, 'Use --resume')
        locked = load_pt(final)
        require(locked['signature'] == dsig, 'Different selected point')
    else:
        candidates = [dict(metrics=initial, epoch=0, DeltaW=torch.zeros_like(W0), lr=None, beta=None, drift=0.0)]
        for lr in LRS:
            for beta in BETAS:
                candidates.append(fit_config(a, inputs, explicit, selected, W0, initial, dsig, directory, lr, beta, device))
        highest = max((x['metrics']['ap'] for x in candidates))
        tied = [x for x in candidates if highest - x['metrics']['ap'] < 0.0005]
        drift = min((x['drift'] for x in tied))
        tied = [x for x in tied if abs(x['drift'] - drift) <= 1e-12]
        best = max(tied, key=lambda x: (x['beta'] if x['beta'] is not None else -1.0, x['metrics']['ap']))
        profile.DeltaW.data.copy_(best['DeltaW'].to(device))
        met, z = predict_profile(profile, inputs, explicit, 'val', a.batch_size, device)
        require(abs(met['ap'] - best['metrics']['ap']) < 1e-12, 'Deterministic replay failed')
        locked = dict(signature=dsig, seed=a.seed, K=len(selected), **best, selected_candidate_indices=selected, feature_ids=inputs.ids[selected], train_p95_scales=inputs.scales[selected], W0=W0, W_final=W0 + best['DeltaW'], threshold=val_threshold(inputs.splits['val'].labels, z))
        save_pt(final, locked)
        pd.DataFrame([dict(lr=x['lr'], beta=x['beta'], best_epoch=x['epoch'], val_ap=x['metrics']['ap'], drift=x['drift'], selected=x is best) for x in candidates]).to_csv(directory / 'refinement_configs.tsv', sep='\t', index=False)
    state = load_pt(out / 'joint/best.pt')['state']
    bundle = portable_bundle(locked, state, inputs, a.seed)
    bundle_path = directory / 'jwarm_bundle.pt'
    if not bundle_path.exists():
        save_pt(bundle_path, bundle)
    profile.DeltaW.data.copy_(locked['DeltaW'].to(device))
    met, z = predict_profile(profile, inputs, explicit, 'val', a.batch_size, device)
    require(abs(met['ap'] - locked['metrics']['ap']) < 1e-12, 'Locked replay failed')
    met['at_validation_f1_threshold'] = threshold_metrics(inputs.splits['val'].labels, z, locked['threshold'])
    save_predictions(directory, 'val', inputs.splits['val'], z, explicit['val'])
    write_json(directory / 'validation_metrics.json', dict(seed=a.seed, K=len(selected), metrics=met, before_refinement=initial))
    logging.info('Locked K=%d Validation AP=%.6f; bundle=%s', len(selected), met['ap'], bundle_path)
    if a.evaluate_test:
        test_signature = identity([digest(bundle_path), digest(PROJECT_DIR / 'new_data/test.tsv'), digest(PROJECT_DIR / 'new_data/test_metadata.tsv')])
        if (directory / 'test_complete.json').exists():
            done = read(directory / 'test_complete.json')
            require(done['signature'] == test_signature, 'Test inputs changed')
            for name, h in done['outputs'].items():
                require(digest(directory / name) == h, 'Test output changed')
            return
        test = inputs.data.split('test')
        validate_alignment(inputs.data, inputs.raw, test)
        proteins = set(inputs.raw.metadata.iloc[test.sparse_rows].base_protein_id)
        for split in inputs.splits.values():
            require(not proteins & set(inputs.raw.metadata.iloc[split.sparse_rows].base_protein_id), 'Test protein leakage')
        predictor = PriorSLiMPredictor(bundle).to(device)
        z, ze = infer_rows(predictor, inputs.data, inputs.raw, test.tensor_rows, a.batch_size, device)
        met = metrics(test.labels, z, test.negative_types)
        met['at_validation_f1_threshold'] = threshold_metrics(test.labels, z, locked['threshold'])
        save_predictions(directory, 'test', test, z, ze)
        write_json(directory / 'test_metrics.json', dict(metrics=met, seed=a.seed, K=len(selected)))
        names = ['test_predictions.npz', 'test_predictions.tsv', 'test_metrics.json']
        write_json(directory / 'test_complete.json', dict(signature=test_signature, outputs={n: digest(directory / n) for n in names}))

def infer_rows(model, data, store, rows, batch_size, device):
    ids = data.window_ids[rows]
    sparse_rows = store.rows_for_ids(ids)
    total = np.empty(len(rows), np.float64)
    explicit = total.copy()
    model.eval()
    weights = model.network.profile.weights.detach().cpu().numpy().astype(np.float64)
    lookup = model.lookup.cpu().numpy()
    scales = model.scales.cpu().numpy()
    with torch.inference_mode():
        for take in batches(len(rows), batch_size):
            modalities = data.raw_batch(rows[take], device)
            network = model.network
            ze = network.intrinsic(*modalities[:4]) + network.access((modalities[4] - network.access_mean) / network.access_scale)
            explicit[take] = ze.cpu().numpy()
            sparse = store.get_batch(sparse_rows[take])
            f = sparse['feature_indices'].numpy()
            v = sparse['feature_values'].numpy()
            res = sparse['residue_indices'].numpy()
            j = lookup[f]
            keep = j >= 0
            normalized = (v[keep] / (scales[j[keep]] + 1e-08)).astype(np.float32)
            delta = np.bincount(res[keep] // 29, weights=weights[j[keep], res[keep] % 29] * normalized, minlength=len(take))
            total[take] = explicit[take] + delta
    require(np.isfinite(total).all(), 'Prediction contains NaN/Inf')
    return (total, explicit)

def predict_command(a, device):
    bundle = load_pt(a.checkpoint)
    model = PriorSLiMPredictor(bundle).to(device)
    store = SparseStore(a.store_dir)
    data = FeatureData(store, a.tensor, a.metadata)
    manifest = read(a.feature_manifest)
    for key in ['chemical_features', 'pssm_residue_features', 'accessibility_features']:
        if key in bundle['feature_manifest']:
            require(manifest.get(key) == bundle['feature_manifest'][key], f'Feature names differ: {key}')
    rows = np.arange(len(data.window_ids)) if a.split == 'all' else data.combined.index[data.combined['split'] == a.split].to_numpy()
    require(len(rows) > 0, 'No requested windows')
    sm = store.metadata.iloc[store.rows_for_ids(data.window_ids[rows])]
    require(set(sm.crop_mode) == {'idr_coverage'} and set(sm.context_margin) == {128}, 'Inference requires same m128 crop protocol')
    frame = data.combined.iloc[rows]
    for l, r in [('source_protein_id', 'uniprot'), ('window_start', 'window_start'), ('window_end', 'window_end')]:
        require(np.array_equal(sm[l].astype(str), frame[r].astype(str)), 'Prediction sequence/coordinate alignment failed')
    require(not sm.sequence_sha256.isna().any(), 'Missing actual-sequence hash')
    prefix = a.output_prefix
    prefix.parent.mkdir(parents=True, exist_ok=True)
    for suffix in ['.tsv', '.npz', '.json']:
        require(not Path(str(prefix) + suffix).exists(), 'Prediction output exists; choose a new --output-prefix')
    z, ze = infer_rows(model, data, store, rows, a.batch_size, device)
    probabilities = expit(z)
    result = pd.DataFrame(dict(window_id=data.window_ids[rows], logit=z, probability=probabilities, explicit_logit=ze, motifae_logit=z - ze, prediction_at_validation_threshold=probabilities >= model.threshold))
    for key in ['uniprot', 'window_start', 'window_end', 'label']:
        if key in frame:
            result[key] = frame[key].to_numpy()
    result.to_csv(str(prefix) + '.tsv', sep='\t', index=False)
    np.savez_compressed(str(prefix) + '.npz', window_ids=data.window_ids[rows], logits=z, probabilities=probabilities, explicit_logits=ze, motifae_logits=z - ze)
    write_json(Path(str(prefix) + '.json'), dict(checkpoint=str(a.checkpoint), checkpoint_sha256=digest(a.checkpoint), rows=len(z), threshold=model.threshold, threshold_scope='Validation window classifier; not an IDR peak threshold', input_hashes={str(p): digest(p) for p in [a.tensor, a.metadata, a.feature_manifest, a.store_dir / 'motifae_sparse_config.json']}))
    logging.info('Predicted %d windows -> %s', len(z), prefix)

def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('phase', choices=['prepare', 'frontier', 'refine', 'predict'])
    p.add_argument('--seed', type=int, default=2027)
    p.add_argument('--out-dir', type=Path)
    p.add_argument('--assets', type=Path, default=PROJECT_DIR / 'model/jwarm_training_assets.pt')
    p.add_argument('--tensor', type=Path, default=PROJECT_DIR / 'train/combined_tensor.npz')
    p.add_argument('--metadata', type=Path, default=PROJECT_DIR / 'train/combined_metadata.tsv')
    p.add_argument('--feature-manifest', type=Path, default=PROJECT_DIR / 'train/combined_feature_manifest.json')
    p.add_argument('--store-dir', type=Path, default=PROJECT_DIR / 'feature_cache/motifae_sparse_windows_idr_coverage')
    p.add_argument('--device', choices=['auto', 'cpu', 'mps', 'cuda'], default='auto')
    p.add_argument('--threads', type=int, default=4)
    p.add_argument('--batch-size', type=int, default=128)
    p.add_argument('--solver-max-iter', type=int, default=10000)
    p.add_argument('--selected-k', type=int)
    p.add_argument('--evaluate-test', action='store_true')
    p.add_argument('--resume', action='store_true')
    p.add_argument('--checkpoint', type=Path, default=PROJECT_DIR / 'model/jwarm_seed2027_K743.pt')
    p.add_argument('--split', choices=['all', 'train', 'val', 'test'], default='all')
    p.add_argument('--output-prefix', type=Path)
    a = p.parse_args(argv)
    a.out_dir = (a.out_dir or PROJECT_DIR / 'train/outputs_jwarm' / f'seed{a.seed}').resolve()
    for k in ['assets', 'tensor', 'metadata', 'feature_manifest', 'store_dir', 'checkpoint']:
        setattr(a, k, getattr(a, k).expanduser().resolve())
    if a.output_prefix:
        a.output_prefix = a.output_prefix.expanduser().resolve()
    if a.phase == 'predict' and a.output_prefix is None:
        p.error('predict needs --output-prefix')
    if a.phase == 'refine' and (not a.selected_k):
        p.error('refine needs an observed --selected-k')
    if a.evaluate_test and a.phase != 'refine':
        p.error('--evaluate-test only applies after refinement')
    require(min(a.threads, a.batch_size, a.solver_max_iter) > 0, 'Positive sizes required')
    for k, v in dict(base_epochs=60, access_epochs=30, profile_epochs=20, base_lr=0.0001, access_lr=0.0003, profile_lr=0.0003, base_wd=0.001, profile_wd=0.0001, base_patience=8, residual_patience=4, motif_dropout=0.4, solver_tol=1e-05, arm='J_WARM').items():
        setattr(a, k, v)
    return a

def main(argv=None):
    a = parse_args(argv)
    torch.set_num_threads(a.threads)
    device = device_for(a.device)
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s', handlers=[logging.StreamHandler(sys.stdout)], force=True)
    if a.phase == 'predict':
        predict_command(a, device)
        return 0
    protected = [PROJECT_DIR / 'model', PROJECT_DIR / 'feature_cache', PROJECT_DIR / 'new_data', PROJECT_DIR / 'prior_compression_experiments', PROJECT_DIR / 'train/robust_slim_variants']
    require(not any((a.out_dir == p or p in a.out_dir.parents or a.out_dir in p.parents for p in protected)), 'Training output overlaps protected inputs or results')
    a.out_dir.mkdir(parents=True, exist_ok=True)
    logging.getLogger().addHandler(logging.FileHandler(a.out_dir / 'run.log'))
    with exclusive(a.out_dir / 'process.lock'):
        inputs = Inputs(a)
        inputs.training['lambda_ratios'] = list(RATIOS)
        config = dict(architecture='J_WARM', seed=a.seed, training=inputs.training, solver_max_iter=a.solver_max_iter, inputs=inputs.fingerprint, objective='classification_BCE_only', device=str(device), code={p.name: digest(p) for p in [Path(__file__), Path(architecture.__file__)]})
        signature = identity(config)
        if (a.out_dir / 'config.json').exists():
            require(read(a.out_dir / 'config.json')['signature'] == signature, 'Training inputs/code/settings changed')
        else:
            write_json(a.out_dir / 'config.json', dict(signature=signature, **config))
        if a.phase == 'prepare':
            return 0
        if a.phase == 'refine':
            run_refine(a, inputs, signature, device)
            return 0
        if not (a.out_dir / 'train_complete.json').exists():
            run_stage_a(a, inputs, signature, device)
            model, ck = fit_joint(a, inputs, a.out_dir, signature, device)
            export_full(a, inputs, a.out_dir, model, ck, signature, device)
            del model
            release(device)
        else:
            require(read(a.out_dir / 'train_complete.json')['signature'] == signature, 'Joint checkpoint contract differs')
        if not (a.out_dir / 'compress_complete.json').exists():
            compress(a, inputs, a.out_dir, signature)
        else:
            require(read(a.out_dir / 'compress_complete.json')['signature'] == signature, 'Compression contract differs')
        rows = []
        for path in (a.out_dir / 'points').glob('*.json'):
            q = read(path)
            rows.append(dict(point=q['point_id'], K=q['K'], lambda_ratio=q['lambda_ratio'], converged=q['converged'], eligible=q['eligible'], val_ap=q['val']['ap'], val_auc=q['val']['auc']))
        pd.DataFrame(rows).sort_values('lambda_ratio').to_csv(a.out_dir / 'frontier.tsv', sep='\t', index=False)
        logging.info('Frontier complete. Manually select K from frontier.tsv; no Test evaluated.')
    return 0
if __name__ == '__main__':
    raise SystemExit(main())
