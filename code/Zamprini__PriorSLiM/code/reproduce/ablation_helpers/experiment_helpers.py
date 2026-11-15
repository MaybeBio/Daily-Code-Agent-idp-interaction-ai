from pathlib import Path
import sys
import json
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

class IntrinsicAblation(nn.Module):

    def __init__(self, baseline, flags, scalers, motif_dropout=0.4):
        super().__init__()
        self.flags = flags.copy()
        self.motif_dropout = motif_dropout
        self.motif_indices = baseline.MOTIF_INDICES
        self.chemistry_indices = baseline.CHEMISTRY_INDICES
        for key, values in scalers.items():
            for field in ('mean', 'scale'):
                self.register_buffer(f'{key}_{field}', values[field].clone())
        self.conservation = baseline.ConservationEncoder() if flags['conservation'] else None
        self.motif = baseline.MotifEncoder() if flags['motif'] else None
        self.chemistry = baseline.ChemistryEncoder() if flags['chemistry'] else None
        self.cm = baseline._interaction_fusion() if self.conservation is not None and self.motif is not None else None
        self.qm = baseline._interaction_fusion() if self.chemistry is not None and self.motif is not None else None
        count = sum((m is not None for m in (self.conservation, self.motif, self.chemistry, self.cm, self.qm)))
        if not count:
            raise ValueError('A uses an intercept, not IntrinsicAblation')
        self.fusion = nn.Sequential(nn.Linear(count * 32, 96), nn.LayerNorm(96), nn.GELU(), nn.Dropout(0.2), nn.Linear(96, 64), nn.LayerNorm(64), nn.GELU(), nn.Dropout(0.15))
        self.head = nn.Linear(64, 1)

    def forward(self, residue, mask, peak, chemical):
        pieces = []
        cons = motif = chem = None
        if self.conservation is not None:
            p = torch.cat((residue, mask.unsqueeze(-1)), dim=-1).reshape(len(residue), -1)
            p = ((p - self.pssm_mean) / self.pssm_scale).reshape(-1, 29, 7)
            cons = self.conservation(p, (peak - self.pssm_window_mean) / self.pssm_window_scale)
            pieces.append(cons)
        if self.motif is not None:
            idx = list(self.motif_indices)
            m = (chemical[:, idx] - self.chemical_mean[idx]) / self.chemical_scale[idx]
            if self.training and self.motif_dropout:
                m = m * (torch.rand((len(m), 1), device=m.device) >= self.motif_dropout)
            motif = self.motif(m)
            pieces.append(motif)
        if self.chemistry is not None:
            idx = list(self.chemistry_indices)
            chem = self.chemistry((chemical[:, idx] - self.chemical_mean[idx]) / self.chemical_scale[idx])
            pieces.append(chem)
        if self.cm is not None:
            pieces.append(self.cm(torch.cat((cons, motif, cons * motif), dim=1)))
        if self.qm is not None:
            pieces.append(self.qm(torch.cat((chem, motif, chem * motif), dim=1)))
        return self.head(self.fusion(torch.cat(pieces, dim=1))).squeeze(-1)

class EnsembleResidual(nn.Module):

    def __init__(self):
        super().__init__()
        self.groups = ((0, 12), (12, 24), (24, 36), (36, 38))
        self.branches = nn.ModuleList()
        for start, end in self.groups:
            hidden = 4 if end - start == 2 else 8
            branch = nn.Sequential(nn.Linear(end - start, hidden), nn.LayerNorm(hidden), nn.GELU(), nn.Linear(hidden, hidden), nn.GELU(), nn.Linear(hidden, 1))
            nn.init.zeros_(branch[-1].weight)
            nn.init.zeros_(branch[-1].bias)
            self.branches.append(branch)

    def forward(self, x):
        if x.ndim != 2 or x.shape[1] != 38:
            raise ValueError('Ensemble input must be [B,38]')
        return sum((branch(x[:, start:end]).squeeze(-1) for branch, (start, end) in zip(self.branches, self.groups)))

def ensemble_inputs(project, data, splits):
    from sklearn.preprocessing import StandardScaler
    path = project / 'feature_cache/gussidr_ensemble_features.tsv'
    keys = ['uniprot', 'idr_start', 'idr_end']
    names = [f'feat_{i}' for i in range(1, 39)]
    cache = pd.read_csv(path, sep='\t', low_memory=False)
    if cache.columns.tolist() != keys + names or cache.duplicated(keys).any():
        raise ValueError('Expected unique exact protein/IDR keys and full38 GussIDR schema')
    raw = {}
    for name, split in splits.items():
        left = data.combined.iloc[split.tensor_rows][keys].copy()
        if left.isna().any().any():
            raise ValueError('Missing ensemble join keys')
        left['_order'] = np.arange(len(left))
        joined = left.merge(cache, on=keys, how='left', sort=False, validate='many_to_one', indicator=True)
        if not joined._merge.eq('both').all() or not np.array_equal(joined._order, np.arange(len(left))):
            raise ValueError(f'{name}: missing exact sequence-derived ensemble rows; no fallback allowed')
        raw[name] = joined[names].to_numpy(dtype=np.float32)
        if not np.isfinite(raw[name]).all():
            raise ValueError(f'{name}: nonfinite ensemble input')
    scaler = StandardScaler().fit(raw['train'])
    values = {s: scaler.transform(x).astype(np.float32) for s, x in raw.items()}
    return (values, dict(source_sha256=digest(path), feature_names=names, groups=dict(asphericity=names[:12], end_to_end=names[12:24], rg=names[24:36], global_polymer=names[36:]), join_keys=keys, canonical_fallback=False, normalization='TRAIN_only', mean=scaler.mean_, scale=scaler.scale_, rows={s: len(x) for s, x in raw.items()}))

def _observed_candidate_tracks(local: dict[str, torch.Tensor], lookup: np.ndarray, candidate_count: int) -> tuple[np.ndarray, np.ndarray]:
    feature_ids = local['feature_indices'].numpy().astype(np.int64, copy=False)
    values = local['feature_values'].numpy().astype(np.float32, copy=False)
    residue_indices = local['residue_indices'].numpy().astype(np.int64, copy=False)
    compact = lookup[feature_ids]
    keep = (values != 0) & (compact >= 0)
    if not np.any(keep):
        return (np.empty(0, dtype=np.int64), np.empty(0, dtype=np.float32))
    values = values[keep]
    residue_indices = residue_indices[keep]
    compact = compact[keep]
    windows = residue_indices // WINDOW_LENGTH
    positions = residue_indices % WINDOW_LENGTH
    pair_keys = windows * int(candidate_count) + compact
    unique_keys, inverse = np.unique(pair_keys, return_inverse=True)
    tracks = np.zeros((len(unique_keys), WINDOW_LENGTH), dtype=np.float32)
    np.add.at(tracks, (inverse, positions), values)
    summaries = np.max((tracks[:, :-2] + tracks[:, 1:-1] + tracks[:, 2:]) / 3.0, axis=1)
    return (unique_keys % int(candidate_count), summaries.astype(np.float32, copy=False))

def load_subset(args):
    legacy = t
    store = t.SparseStore(args.store_dir)
    data = t.FeatureData(store, args.tensor, args.metadata)
    assets = common.load_pt(args.project_root / 'model/jwarm_training_assets.pt')
    ids = assets['feature_ids'].long().clone()
    del assets
    require(ids.shape == (5845,) and len(torch.unique(ids)) == 5845, 'Expected fixed5845 candidates')
    feature_manifest = json.loads((args.project_root / 'train/combined_feature_manifest.json').read_text())
    require([feature_manifest['chemical_features'][i] for i in baseline.MOTIF_INDICES] == ['motif_min_abs_offset', 'motif_center_enrichment', 'motif_count_log'], 'Motif schema changed')
    require(data.arrays['pssm_residue'].shape[1:] == (29, 6) and data.arrays['chemical'].shape[1:] == (58,), 'Unexpected feature shape')
    manifest = json.loads((args.subset_dir / 'subset_manifest.json').read_text())
    stats = manifest['statistics']
    require(stats['seed'] == args.seed and stats['requested_percent'] == args.percent, 'Wrong subset')
    for name, field in (('train.tsv', 'train_sha256'), ('train_metadata.tsv', 'metadata_sha256'), ('train_row_index.tsv', 'row_index_sha256'), ('val.tsv', 'validation_sha256')):
        require(common.digest(args.subset_dir / name) == manifest[field], f'Subset file changed: {name}')
    full_index = read(args.datasets_dir / 'train_full_row_index.tsv')
    splits, indices = ({}, {})
    for split in ('train', 'val'):
        ix = read(args.subset_dir / f'{split}_row_index.tsv')
        rows = ix.tensor_index.to_numpy(dtype=np.int64)
        sparse_rows = store.rows_for_ids(ix.window_id)
        frame = data.combined.iloc[rows].reset_index(drop=True)
        sm = store.metadata.iloc[sparse_rows].reset_index(drop=True)
        table = read(args.subset_dir / f'{split}.tsv')
        meta = read(args.subset_dir / f'{split}_metadata.tsv')
        require(len(ix) == len(table) == len(meta) and len(np.unique(rows)) == len(rows), 'Invalid subset lengths/rows')
        require(frame.split.eq(split).all() and np.array_equal(frame.split_row, ix.source_split_row), 'Wrong original split rows')
        for source in (frame, sm):
            require(np.array_equal(source.window_id, ix.window_id), 'Feature window-ID mismatch')
        require(np.array_equal(sm.sparse_window_index, ix.sparse_window_index), 'CSR physical slot mapping differs')
        require(np.array_equal(sm.base_protein_id, ix.base_protein_id), 'Biological protein mismatch')
        for col, mcol in (('uniprot', 'protein_id'), ('window_start', 'window_start'), ('window_end', 'window_end'), ('sequence', 'window_sequence')):
            require(np.array_equal(frame[col], table[col]) and np.array_equal(frame[col], meta[mcol]), 'Sequence/coordinate mismatch: ' + col)
        require(np.array_equal(sm.source_protein_id, frame.uniprot) and np.array_equal(sm.window_start, frame.window_start) and np.array_equal(sm.window_end, frame.window_end), 'Sparse sequence/crop mapping differs')
        require(np.array_equal(meta.label.astype(int), frame.label.astype(int)) and np.array_equal(table.label.astype(int), frame.label.astype(int)), 'Labels mismatch')
        require(sm.crop_mode.eq('idr_coverage').all() and sm.context_margin.eq(128).all() and sm.sequence_sha256.notna().all(), 'Expected sequence-aware m128 cache')
        splitdata = legacy.SplitData(split, meta, rows, sparse_rows)
        splitdata.window_ids = ix.window_id.to_numpy(dtype=str)
        splits[split], indices[split] = (splitdata, ix)
    expected = full_index[full_index.homology_component.isin(indices['train'].homology_component)]
    require(np.array_equal(expected.source_split_row, indices['train'].source_split_row), 'Partial homology component')
    require(not set(indices['train'].homology_component) & set(indices['val'].homology_component), 'Homology leakage')
    require(not set(indices['train'].base_protein_id) & set(indices['val'].base_protein_id), 'Protein leakage')
    original_val = read(args.project_root / 'new_data/val.tsv')
    require(np.array_equal(indices['val'].source_split_row, np.arange(len(original_val))), 'Validation subset changed')
    require(common.digest(args.subset_dir / 'val.tsv') == common.digest(args.project_root / 'new_data/val.tsv') and common.digest(args.subset_dir / 'val_metadata.tsv') == common.digest(args.project_root / 'new_data/val_metadata.tsv'), 'Validation must remain fixed')
    require(len(splits['train']) == stats['windows'], 'TRAIN size mismatch')
    return (baseline, legacy, store, data, splits, ids, feature_manifest, stats)

def fit_subset_p95(legacy, store, train, ids, output, batch_size):
    count = len(ids)
    lookup = np.full(40960, -1, dtype=np.int64)
    lookup[ids.numpy()] = np.arange(count)
    capacities = np.zeros(count, dtype=np.int64)

    def summaries():
        for step, take in enumerate(common.batches(len(train), batch_size)):
            j, q = legacy._observed_candidate_tracks(store.get_batch(train.sparse_rows[take]), lookup, count)
            require(np.isfinite(q).all() and (q >= 0).all(), 'Nonfinite/negative MotifAE summary')
            keep = q > 0
            if step % 100 == 0:
                logging.info('Subset P95 pass batch=%d/%d', step + 1, (len(train) + batch_size - 1) // batch_size)
            yield (j[keep], q[keep])
    for j, q in summaries():
        capacities += np.bincount(j, minlength=count)
    offsets = np.r_[0, capacities.cumsum()]
    output.parent.mkdir(parents=True, exist_ok=True)
    temp = output.with_suffix('.nonzero.tmp')
    values = np.memmap(temp, mode='w+', dtype=np.float32, shape=(max(1, int(offsets[-1])),))
    cursors = np.zeros(count, dtype=np.int64)
    try:
        for j, q in summaries():
            order = np.argsort(j, kind='stable')
            j, q = (j[order], q[order])
            unique, first, sizes = np.unique(j, return_index=True, return_counts=True)
            within = np.arange(len(j)) - np.repeat(first, sizes)
            values[offsets[j] + cursors[j] + within] = q
            cursors[unique] += sizes
        require(np.array_equal(cursors, capacities), 'P95 two-pass mismatch')
        scales = np.ones(count, dtype=np.float32)
        for j in np.flatnonzero(capacities):
            scales[j] = np.quantile(values[offsets[j]:offsets[j + 1]], 0.95)
    finally:
        del values
        temp.unlink(missing_ok=True)
    require(np.isfinite(scales).all() and (scales > 0).all(), 'Invalid P95')
    frame = pd.DataFrame(dict(feature_id=ids.numpy(), scale_p95=scales, nonzero_summary_count=capacities, unseen_in_subset=capacities == 0, summary='max_3_consecutive_mean', scale_split='selected_train_only'))
    frame.to_csv(output, sep='\t', index=False)
    logging.info('Subset P95 ready: seen=%d unseen=%d', int((capacities > 0).sum()), int((capacities == 0).sum()))
    return frame
