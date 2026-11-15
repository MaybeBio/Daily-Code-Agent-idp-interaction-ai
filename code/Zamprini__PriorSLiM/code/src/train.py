from __future__ import annotations
import argparse
import gzip
import hashlib
from contextlib import contextmanager
from functools import partial
import logging
from pathlib import Path
import sys
import tempfile
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _training as t
import numpy as np
import pandas as pd
import torch
ROOT = t.PROJECT_DIR
MIN_ACTIVE_WINDOWS = 10
TRAIN_SEED = 42
EFFICIENCY_AP_TOLERANCE = 0.005
TRAIN_KEYS = ['batch_size', 'base_epochs', 'access_epochs', 'profile_epochs', 'base_lr', 'access_lr', 'profile_lr', 'base_wd', 'profile_wd', 'base_patience', 'residual_patience', 'motif_dropout', 'solver_tol']

def table(path):
    return pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False)

def verify_file(path, expected):
    t.require(path.is_file(), f'Missing input: {path}')
    t.require(t.digest(path) == expected, f'Input changed: {path}')

def verify_sparse_metadata(path, expected, metadata_path):
    if t.digest(path) == expected:
        return
    reference = metadata_path.parent.parent / 'feature_inputs/jobs/benchmark/reference_sparse_metadata.tsv.gz'
    t.require(reference.is_file(), f'Sparse mapping changed without a reconstruction reference: {path}')
    with gzip.open(reference, 'rb') as handle:
        original = handle.read()
    t.require(hashlib.sha256(original).hexdigest() == expected, 'Reference sparse mapping checksum differs')
    columns = ['window_id', 'window_index', 'protein_id', 'base_protein_id',
               'sparse_window_index', 'window_start', 'window_end', 'chunk_id',
               'crop_start', 'crop_end', 'sequence_sha256', 'valid_local_start',
               'valid_local_end', 'left_padding', 'right_padding']
    old = pd.read_csv(reference, sep='\t', dtype=str, keep_default_na=False, usecols=columns)
    new = pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False, usecols=columns)
    t.require(old[columns].equals(new[columns]), 'Reconstructed sparse window identities or coordinates differ')

def verify_fold_contract(a):
    manifest = t.read(a.cv_dir / 'manifest.json')
    t.require(manifest['n_splits'] == 10, 'Expected ten grouped development folds')
    t.require(manifest['protocol'] == 'development_grouped_10fold_9train_1val_fixed_test', 'Wrong fold protocol')
    folder = a.cv_dir / f'fold_{a.fold:02d}'
    for role in ('train', 'val'):
        for suffix in ('.tsv', '_metadata.tsv', '_indices.tsv'):
            p = folder / f'{role}{suffix}'
            verify_file(p, manifest['output_sha256'][str(p.relative_to(a.cv_dir))])
    verify_file(a.metadata, manifest['sources']['combined']['sha256'])
    return manifest

def load_fold_splits(a, data, raw):
    folder = a.cv_dir / f'fold_{a.fold:02d}'
    splits = {}
    for role in ('train', 'val'):
        meta, windows = (table(folder / f'{role}_metadata.tsv'), table(folder / f'{role}.tsv'))
        rows = meta.tensor_index.to_numpy(dtype=np.int64)
        t.require(len(windows) == len(meta) and len(set(rows)) == len(rows), f'{role}: rows')
        t.require(meta['split'].eq(role).all() and meta.cv_fold.astype(int).eq(a.fold).all(), f'{role}: role')
        t.require(np.all((rows >= 0) & (rows < len(data.combined))), f'{role}: tensor bounds')
        cache = data.combined.iloc[rows]
        for x, y in [('uniprot', 'protein_id'), ('sequence', 'window_sequence')]:
            t.require(np.array_equal(cache[x].astype(str), meta[y]), f'{role}: cache {x}')
            t.require(np.array_equal(windows[x], meta[y]), f'{role}: table {x}')
        for key in ('window_start', 'window_end', 'label'):
            t.require(np.array_equal(cache[key].astype(float), meta[key].astype(float)), f'{role}: cache {key}')
            t.require(np.array_equal(windows[key].astype(float), meta[key].astype(float)), f'{role}: table {key}')
        t.require(np.array_equal(cache.window_id.astype(str), meta.window_id), f'{role}: window IDs')
        t.require(not cache['split'].eq('test').any(), 'Test row entered CV')
        sparse_rows = raw.rows_for_ids(meta.window_id)
        t.require(np.array_equal(sparse_rows, meta.sparse_row_index.astype(int)), 'Sparse row mapping changed')
        split = t.SplitData(role, meta, rows, sparse_rows)
        t.validate_alignment(data, raw, split)
        t.require(np.array_equal(raw.metadata.iloc[sparse_rows].base_protein_id.astype(str), meta.base_protein_id), 'Sparse parent protein mismatch')
        t.require(meta.heldout_fold.astype(int).eq(a.fold).all() if role == 'val' else meta.heldout_fold.astype(int).ne(a.fold).all(), 'Wrong held-out fold')
        splits[role] = split
    for col in ('window_id', 'component_id', 'base_protein_id'):
        t.require(not set(splits['train'].metadata[col]) & set(splits['val'].metadata[col]), f'Cross-fold {col}')
    t.require(len(splits['train']) + len(splits['val']) == 87801, 'Development pool size changed')
    return splits

def observed_summaries(local):
    f = local['feature_indices'].numpy().astype(np.int64, copy=False)
    r = local['residue_indices'].numpy().astype(np.int64, copy=False)
    v = local['feature_values'].numpy().astype(np.float32, copy=False)
    t.require(np.isfinite(v).all() and (v >= 0).all(), 'Invalid sparse activations')
    t.require(((f >= 0) & (f < 40960)).all(), 'Feature ID out of bounds')
    keep = v > 0
    f, r, v = (f[keep], r[keep], v[keep])
    keys, inverse = np.unique(r // 29 * 40960 + f, return_inverse=True)
    tracks = np.zeros((len(keys), 29), np.float32)
    np.add.at(tracks, (inverse, r % 29), v)
    q = np.max((tracks[:, :-2] + tracks[:, 1:-1] + tracks[:, 2:]) / 3.0, axis=1)
    return (keys % 40960, q)

def fit_train_assets(raw, train, directory, batch_size, min_windows=MIN_ACTIVE_WINDOWS):
    directory.mkdir(parents=True, exist_ok=True)
    counts = np.zeros(40960, np.int64)

    def iterate():
        for step, take in enumerate(t.batches(len(train), batch_size)):
            j, q = observed_summaries(raw.get_batch(train.sparse_rows[take]))
            t.require(np.isfinite(q).all() and (q > 0).all(), 'Invalid positive summary')
            if step % 100 == 0:
                logging.info('Fold TRAIN feature statistics batch=%d rows=%d/%d', step, min((step + 1) * batch_size, len(train)), len(train))
            yield (j, q)
    for j, q in iterate():
        counts += np.bincount(j, minlength=40960)
    ids = np.flatnonzero(counts >= min_windows)
    t.require(len(ids) > 0, 'No TRAIN candidates passed active-window filter')
    lookup = np.full(40960, -1, np.int64)
    lookup[ids] = np.arange(len(ids))
    capacities = counts[ids]
    offsets = np.r_[0, capacities.cumsum()]
    cursors = np.zeros(len(ids), np.int64)
    with tempfile.TemporaryDirectory(prefix='p95_', dir=directory) as temp:
        values = np.memmap(Path(temp) / 'nonzero.bin', mode='w+', dtype=np.float32, shape=(int(offsets[-1]),))
        try:
            for f, q in iterate():
                j = lookup[f]
                keep = j >= 0
                j, q = (j[keep], q[keep])
                order = np.argsort(j, kind='stable')
                j, q = (j[order], q[order])
                unique, first, sizes = np.unique(j, return_index=True, return_counts=True)
                within = np.arange(len(j)) - np.repeat(first, sizes)
                values[offsets[j] + cursors[j] + within] = q
                cursors[unique] += sizes
            t.require(np.array_equal(cursors, capacities), 'P95 pass counts differ')
            scales = np.array([np.quantile(values[offsets[k]:offsets[k + 1]], 0.95) for k in range(len(ids))], np.float32)
        finally:
            del values
    t.require(np.isfinite(scales).all() and (scales > 0).all(), 'Invalid P95')
    pd.DataFrame(dict(feature_id=np.arange(40960), active_train_windows=counts, selected=counts >= min_windows)).to_csv(directory / 'train_feature_frequency.tsv', sep='\t', index=False)
    pd.DataFrame(dict(feature_id=ids, train_p95_scale=scales, active_train_windows=counts[ids])).to_csv(directory / 'train_candidate_scales.tsv', sep='\t', index=False)
    logging.info('Fold-local candidates=%d (not forced to 5845)', len(ids))
    return (torch.tensor(ids, dtype=torch.long), torch.tensor(scales))

class FoldInputs(t.Inputs):

    def __init__(self, a):
        self.source = a.out_dir / 'stage_a'
        split_manifest = verify_fold_contract(a)
        self.raw = t.SparseStore(a.store_dir)
        verify_sparse_metadata(self.raw.paths['window_metadata'], split_manifest['sources']['sparse']['sha256'], a.metadata)
        self.data = t.FeatureData(self.raw, a.tensor, a.metadata)
        self.splits = load_fold_splits(a, self.data, self.raw)
        self.feature_manifest = t.read(a.feature_manifest)
        names = self.feature_manifest['chemical_features']
        t.require([names[i] for i in t.architecture.MOTIF_INDICES] == ['motif_min_abs_offset', 'motif_center_enrichment', 'motif_count_log'], 'Motif columns differ')
        t.require(self.data.arrays['pssm_residue'].shape[1:] == (29, 6) and self.data.arrays['chemical'].shape[1:] == (58,), 'Wrong biological feature dimensions')
        paths = {'tensor': a.tensor, 'metadata': a.metadata, 'feature_manifest': a.feature_manifest, 'split_manifest': a.cv_dir / 'manifest.json', 'sparse_config': a.store_dir / 'motifae_sparse_config.json'}
        paths.update({f'sparse_{k}': p for k, p in self.raw.paths.items()})
        for s in ('train', 'val'):
            for suffix in ('.tsv', '_metadata.tsv'):
                paths[s + suffix] = a.cv_dir / f'fold_{a.fold:02d}' / f'{s}{suffix}'
        self.fingerprint = {k: dict(sha256=t.digest(p), bytes=p.stat().st_size) for k, p in paths.items()}
        prep = a.out_dir / '_prepared'
        prep_signature = t.identity(dict(inputs=self.fingerprint, fold=a.fold, min_active_windows=MIN_ACTIVE_WINDOWS, adapter=t.digest(__file__), engine=t.digest(t.__file__)))
        manifest = prep / 'input_manifest.json'
        if manifest.exists():
            t.require(t.read(manifest)['signature'] == prep_signature, 'Prepared inputs changed; use new output directory')
        else:
            t.write_json(manifest, dict(signature=prep_signature, inputs=self.fingerprint, fold=a.fold))
        assets = prep / 'fold_training_assets.pt'
        if not assets.exists():
            ids, scales = fit_train_assets(self.raw, self.splits['train'], prep, a.batch_size)
            t.save_pt(assets, dict(feature_ids=ids, train_p95_scales=scales, preprocessing=t.fit_preprocessing(self.data, self.splits['train']), signature=prep_signature, train_window_ids=self.splits['train'].window_ids.tolist(), min_active_windows=MIN_ACTIVE_WINDOWS))
        saved = t.load_pt(assets)
        t.require(saved['signature'] == prep_signature, 'Wrong fold-local assets')
        t.require(saved['train_window_ids'] == self.splits['train'].window_ids.tolist(), 'Assets trained on other rows')
        self.ids, self.scales, self.scalers = (saved['feature_ids'], saved['train_p95_scales'], saved['preprocessing'])
        t.require(len(self.ids) == len(torch.unique(self.ids)) and (self.scales > 0).all(), 'Bad candidate IDs/scales')
        self.compact = {}
        for name, split in self.splits.items():
            t.make_compact(self.raw, split, self.ids, self.scales, prep / name, a.batch_size)
            self.compact[name] = t.CompactStore(prep / name)
        self.training = {k: getattr(a, k) for k in TRAIN_KEYS}
        self.training['lambda_ratios'] = list(t.RATIOS)

@contextmanager
def dynamic_refinement(candidate_count):
    original = t.SelectedRefinement
    t.SelectedRefinement = partial(original, candidate_count=candidate_count)
    try:
        yield
    finally:
        t.SelectedRefinement = original

def choose_automatic_points(points):
    by_k = {}
    for p in points:
        if not (p.get('eligible') is True and p.get('converged') is True and (int(p['K']) > 0) and (float(p['lambda_ratio']) > 0)):
            continue
        t.require(np.isfinite(p['val']['ap']), 'Nonfinite Validation AP in eligible point')
        k = int(p['K'])
        if k not in by_k or (p['val']['ap'], p['lambda_ratio']) > (by_k[k]['val']['ap'], by_k[k]['lambda_ratio']):
            by_k[k] = p
    t.require(bool(by_k), 'No converged positive-K Stage C point; cannot auto-select or refine')
    performance = min(by_k.values(), key=lambda p: (-p['val']['ap'], int(p['K']), -p['lambda_ratio']))
    cutoff = performance['val']['ap'] - EFFICIENCY_AP_TOLERANCE
    efficient = [p for p in by_k.values() if p['val']['ap'] >= cutoff - 1e-12]
    efficiency = min(efficient, key=lambda p: (int(p['K']), -p['val']['ap'], -p['lambda_ratio']))
    return dict(efficiency=efficiency, performance=performance)

def lock_automatic_selection(a, signature):
    out = a.out_dir
    t.require(t.read(out / 'compress_complete.json')['signature'] == signature, 'Finish this fold Stage C first')
    paths = sorted((out / 'points').glob('lambda_*.json'))
    selections = choose_automatic_points([t.read(p) for p in paths])
    selected = {name: dict(K=int(p['K']), point_id=p['point_id'], lambda_ratio=p['lambda_ratio'], val_ap=p['val']['ap'], solution_sha256=t.digest(out / 'solutions' / f"{p['point_id']}.pt")) for name, p in selections.items()}
    reference = t.read(out / 'points/full_replay.json')
    record = dict(training_signature=signature, fold=a.fold, seed=TRAIN_SEED, policy=dict(eligible='converged positive-lambda positive-K compression points only', performance='maximum Validation AP; ties prefer smaller K', efficiency='minimum K with Validation AP >= performance AP - tolerance', efficiency_ap_tolerance=EFFICIENCY_AP_TOLERANCE, selection_stage='C_before_D', test_used=False, full_replay_eligible=False), full_replay_reference=dict(K=reference['K'], val_ap=reference['val']['ap']), points_sha256={p.name: t.digest(p) for p in paths}, selected=selected)
    path = out / 'automatic_k_selection.json'
    if path.exists():
        t.require(t.identity(t.read(path)) == t.identity(record), 'Locked automatic K selection changed')
    else:
        t.write_json(path, record)
    pd.DataFrame([dict(operating_point=name, **row) for name, row in selected.items()]).to_csv(out / 'automatic_k_selection.tsv', sep='\t', index=False)
    for name, p in selected.items():
        logging.info('Fold %02d %s K=%d Stage-C valAP=%.6f', a.fold, name, p['K'], p['val_ap'])
    return record

def refine_automatic(a, inputs, signature, device):
    record = lock_automatic_selection(a, signature)
    by_k = {}
    for name, p in record['selected'].items():
        by_k.setdefault(p['K'], []).append(name)
    for k, roles in by_k.items():
        a.selected_k = k
        logging.info('Fold %02d Stage D K=%d roles=%s', a.fold, k, ','.join(roles))
        with dynamic_refinement(len(inputs.ids)):
            t.run_refine(a, inputs, signature, device)
        path = a.out_dir / f'K{k}' / 'jwarm_bundle.pt'
        bundle = t.load_pt(path)
        selection = dict(K='automatic Validation-only Stage-C operating point', operating_points=roles, refinement='Validation AP', test_used_for_selection=False, policy=record['policy'], stage_c=record['selected'][roles[0]])
        if bundle.get('cv_fold') != a.fold or t.identity(bundle.get('selection', {})) != t.identity(selection):
            bundle['selection'] = selection
            bundle['cv_fold'] = a.fold
            t.save_pt(path, bundle)
    t.write_json(a.out_dir / 'automatic_refinement_complete.json', dict(training_signature=signature, selected=record['selected'], unique_refined_K=sorted(by_k), test_evaluated=False))

def summarize(a):
    rows = []
    frontier = []
    selected_rows = []
    for fold in a.folds:
        out = a.output_root / f'seed{a.seed}' / f'fold_{fold:02d}'
        path = out / 'full_metrics.json'
        if path.exists():
            report = t.read(path)
            v = report['metrics']['val']['full']
            rows.append(dict(fold=fold, stage='B_full', K=report['K'], candidate_count=report['candidate_count'], **v))
        path = out / 'frontier.tsv'
        if path.exists():
            f = pd.read_csv(path, sep='\t')
            f.insert(0, 'fold', fold)
            frontier.append(f)
        selection = out / 'automatic_k_selection.json'
        if selection.exists():
            record = t.read(selection)
            for name, p in record['selected'].items():
                selected_rows.append(dict(fold=fold, operating_point=name, **p))
                path = out / f"K{p['K']}" / 'validation_metrics.json'
                if path.exists():
                    r = t.read(path)
                    rows.append(dict(fold=fold, stage='D_refined', operating_point=name, K=r['K'], **r['metrics']))
    directory = a.output_root / f'seed{a.seed}'
    directory.mkdir(parents=True, exist_ok=True)
    if rows:
        flat = [{k: v for k, v in r.items() if not isinstance(v, dict)} for r in rows]
        pd.DataFrame(flat).to_csv(directory / 'cv_validation_summary.tsv', sep='\t', index=False)
        t.write_json(directory / 'cv_validation_metrics.json', dict(results=rows, test_evaluated=False, caveat='Validation also used for early stopping/K selection; not an unbiased outer-CV estimate'))
    if frontier:
        pd.concat(frontier, ignore_index=True).to_csv(directory / 'all_folds_frontier.tsv', sep='\t', index=False)
    if selected_rows:
        pd.DataFrame(selected_rows).to_csv(directory / 'automatic_k_summary.tsv', sep='\t', index=False)
    print(f'Summary written to: {directory}', flush=True)

def args_from_cli(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('phase', choices=['check', 'prepare', 'frontier', 'refine', 'summary'])
    p.add_argument('--fold', type=int)
    p.add_argument('--folds', type=int, nargs='+', default=list(range(1, 11)))
    p.add_argument('--seed', type=int, default=TRAIN_SEED, choices=[TRAIN_SEED])
    p.add_argument('--cv-dir', type=Path, default=ROOT / 'new_data/cv10_homology_seed2026')
    p.add_argument('--output-root', type=Path, default=ROOT / 'train/outputs_jwarm_cv10')
    p.add_argument('--device', choices=['cpu', 'mps', 'cuda', 'auto'], default='auto')
    p.add_argument('--threads', type=int, default=4)
    p.add_argument('--batch-size', type=int, default=128)
    p.add_argument('--solver-max-iter', type=int, default=10000)
    p.add_argument('--resume', action='store_true')
    p.add_argument('--tensor', type=Path, default=ROOT / 'train/combined_tensor.npz')
    p.add_argument('--metadata', type=Path, default=ROOT / 'train/combined_metadata.tsv')
    p.add_argument('--feature-manifest', type=Path, default=ROOT / 'train/combined_feature_manifest.json')
    p.add_argument('--store-dir', type=Path, default=ROOT / 'feature_cache/motifae_sparse_windows_idr_coverage')
    opts = p.parse_args(argv)
    t.require(all((1 <= f <= 10 for f in opts.folds)), 'Folds must be 1..10')
    if opts.phase not in ('check', 'summary'):
        t.require(opts.fold in range(1, 11), 'Training needs --fold 1..10')
    a = t.parse_args(['prepare', '--seed', str(opts.seed), '--device', opts.device, '--batch-size', str(opts.batch_size), '--threads', str(opts.threads), '--solver-max-iter', str(opts.solver_max_iter)])
    for k, v in vars(opts).items():
        setattr(a, k, v.resolve() if isinstance(v, Path) else v)
    a.evaluate_test = False
    a.out_dir = a.output_root / f'seed{a.seed}' / f'fold_{a.fold or 1:02d}'
    protected = [ROOT / 'model', ROOT / 'feature_cache', ROOT / 'new_data', ROOT / 'prior_compression_experiments', ROOT / 'train/robust_slim_variants']
    t.require(not any((a.output_root == q or q in a.output_root.parents or a.output_root in q.parents for q in protected)), 'Output overlaps protected files')
    return a

def main(argv=None):
    a = args_from_cli(argv)
    torch.set_num_threads(a.threads)
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s', stream=sys.stdout, force=True)
    if a.phase == 'summary':
        summarize(a)
        return 0
    if a.phase == 'check':
        raw = t.SparseStore(a.store_dir)
        data = t.FeatureData(raw, a.tensor, a.metadata)
        for fold in [a.fold] if a.fold else a.folds:
            a.fold = fold
            manifest = verify_fold_contract(a)
            verify_sparse_metadata(raw.paths['window_metadata'], manifest['sources']['sparse']['sha256'], a.metadata)
            splits = load_fold_splits(a, data, raw)
            logging.info('Fold %02d check passed: TRAIN=%d VAL=%d; no fitting/training', fold, len(splits['train']), len(splits['val']))
        audit = a.cv_dir / 'sequence_overlap_audit.json'
        if audit.exists():
            logging.warning('Split overlap audit: %s', t.read(audit))
        return 0
    device = t.device_for(a.device)
    a.out_dir.mkdir(parents=True, exist_ok=True)
    logging.getLogger().addHandler(logging.FileHandler(a.out_dir / 'run.log'))
    with t.exclusive(a.out_dir / 'process.lock'):
        inputs = FoldInputs(a)
        config = dict(architecture='J_WARM_CV10', fold=a.fold, seed=a.seed, training=inputs.training, solver_max_iter=a.solver_max_iter, inputs=inputs.fingerprint, device=str(device), candidate_rule=dict(min_active_train_windows=MIN_ACTIVE_WINDOWS, count=len(inputs.ids)), objective='classification_BCE_only', test_evaluated=False, automatic_k=dict(efficiency_ap_tolerance=EFFICIENCY_AP_TOLERANCE, scope='converged_nonzero_compression_points', selected_before_refinement=True), code={str(p.name): t.digest(p) for p in [Path(__file__), Path(t.__file__), Path(t.architecture.__file__)]})
        signature = t.identity(config)
        if (a.out_dir / 'config.json').exists():
            t.require(t.read(a.out_dir / 'config.json')['signature'] == signature, 'Training contract changed')
        else:
            t.write_json(a.out_dir / 'config.json', dict(signature=signature, **config))
        if a.phase == 'prepare':
            return 0
        if a.phase == 'refine':
            refine_automatic(a, inputs, signature, device)
            return 0
        if not (a.out_dir / 'train_complete.json').exists():
            t.run_stage_a(a, inputs, signature, device)
            model, ck = t.fit_joint(a, inputs, a.out_dir, signature, device)
            t.export_full(a, inputs, a.out_dir, model, ck, signature, device)
            report = t.read(a.out_dir / 'full_metrics.json')
            report['candidate_count'] = len(inputs.ids)
            t.write_json(a.out_dir / 'full_metrics.json', report)
            del model
            t.release(device)
        else:
            t.require(t.read(a.out_dir / 'train_complete.json')['signature'] == signature, 'Wrong joint checkpoint')
        report = t.read(a.out_dir / 'full_metrics.json')
        if report['candidate_count'] != len(inputs.ids):
            report['candidate_count'] = len(inputs.ids)
            t.write_json(a.out_dir / 'full_metrics.json', report)
        if not (a.out_dir / 'compress_complete.json').exists():
            t.compress(a, inputs, a.out_dir, signature)
        else:
            t.require(t.read(a.out_dir / 'compress_complete.json')['signature'] == signature, 'Wrong Stage C checkpoint')
        rows = []
        for path in (a.out_dir / 'points').glob('*.json'):
            q = t.read(path)
            rows.append(dict(point=q['point_id'], K=q['K'], lambda_ratio=q['lambda_ratio'], converged=q['converged'], eligible=q['eligible'], val_ap=q['val']['ap'], val_auc=q['val']['auc']))
        pd.DataFrame(rows).sort_values('lambda_ratio').to_csv(a.out_dir / 'frontier.tsv', sep='\t', index=False)
        lock_automatic_selection(a, signature)
        logging.info('Fold %02d A-C and automatic K selection complete; no Test evaluated', a.fold)
    return 0
if __name__ == '__main__':
    raise SystemExit(main())
