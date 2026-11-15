from __future__ import annotations
import argparse
import logging
import sys
from pathlib import Path
from types import SimpleNamespace
import experiment_helpers as h
from experiment_helpers import t, baseline as m, np, pd, torch, nn
import training_steps as steps
HERE = Path(__file__).resolve().parent
FLAGS = {'A': (False, False, False, False, False), 'B': (True, True, True, False, False), 'C': (True, True, True, True, False), 'D': (True, True, True, True, False), 'B1': (False, True, True, False, False), 'B2': (True, False, True, False, False), 'B3': (True, True, False, False, False), 'CE': (True, True, True, True, True)}
NAMES = ('conservation', 'motif', 'chemistry', 'accessibility', 'ensemble')
TRAINING_KEYS = ('batch_size', 'base_epochs', 'access_epochs', 'profile_epochs', 'base_lr', 'access_lr', 'profile_lr', 'base_wd', 'profile_wd', 'base_patience', 'residual_patience', 'motif_dropout', 'solver_tol')

class JointAblation(nn.Module):

    def __init__(self, scalers, variant, features, seed, prior, motif_dropout=0.4):
        super().__init__()
        self.flags = dict(zip(NAMES, FLAGS[variant]))
        t.seed_all(seed)
        self.intrinsic = h.IntrinsicAblation(m, self.flags, scalers, motif_dropout) if any(FLAGS[variant][:3]) else None
        self.register_buffer('intercept', torch.tensor(float(prior), dtype=torch.float32))
        t.seed_all(seed + 100)
        self.access = m.AccessibilityResidual() if self.flags['accessibility'] else None
        self.register_buffer('access_mean', scalers['accessibility']['mean'].clone())
        self.register_buffer('access_scale', scalers['accessibility']['scale'].clone())
        t.seed_all(seed + 150)
        self.ensemble = h.EnsembleResidual() if self.flags['ensemble'] else None
        t.seed_all(seed + 200)
        self.profile = m.PositionalReadout(features)

    def forward(self, modalities, events, return_parts=False):
        explicit = self.intrinsic(*modalities[:4]) if self.intrinsic is not None else self.intercept.expand(len(modalities[0]))
        if self.access is not None:
            explicit = explicit + self.access((modalities[4] - self.access_mean) / self.access_scale)
        if self.ensemble is not None:
            explicit = explicit + self.ensemble(modalities[5])
        residual = self.profile(events)
        return (explicit + residual, explicit, residual) if return_parts else explicit + residual

    def parameter_groups(self, training):
        result = []
        for module, lr, wd in [(self.intrinsic, 'base_lr', 'base_wd'), (self.access, 'access_lr', 'base_wd'), (self.ensemble, 'access_lr', 'base_wd'), (self.profile, 'profile_lr', 'profile_wd')]:
            if module is not None:
                result.append(dict(params=module.parameters(), lr=training[lr], weight_decay=training[wd]))
        return result

def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('phase', choices=['check', 'prepare', 'run'])
    p.add_argument('--experiment', choices=['E1', 'E3'], required=True)
    p.add_argument('--variant', choices=list(FLAGS), default='C')
    p.add_argument('--seed', type=int, choices=range(2026, 2031), default=2026)
    p.add_argument('--percent', type=int, choices=[10, 25, 50, 75], default=10)
    p.add_argument('--project-root', type=Path, default=HERE.parents[1])
    p.add_argument('--output-root', type=Path)
    p.add_argument('--e1-results', type=Path, default=HERE / 'E1/results')
    p.add_argument('--datasets-dir', type=Path, default=HERE.parent / 'E3/datasets')
    p.add_argument('--store-dir', type=Path)
    p.add_argument('--device', choices=['cpu', 'mps', 'cuda', 'auto'], default='auto')
    p.add_argument('--threads', type=int, default=4)
    p.add_argument('--batch-size', type=int, default=128)
    p.add_argument('--solver-max-iter', type=int, default=10000)
    p.add_argument('--resume', action='store_true')
    a = p.parse_args(argv)
    a.project_root = a.project_root.expanduser().resolve()
    t.PROJECT_DIR = a.project_root
    a.output_root = (a.output_root or HERE / a.experiment / 'results').expanduser().resolve()
    for k in ['e1_results', 'datasets_dir']:
        setattr(a, k, getattr(a, k).expanduser().resolve())
    a.store_dir = (a.store_dir or a.project_root / 'feature_cache/motifae_sparse_windows_idr_coverage').resolve()
    a.tensor = a.project_root / 'train/combined_tensor.npz'
    a.metadata = a.project_root / 'train/combined_metadata.tsv'
    a.feature_manifest = a.project_root / 'train/combined_feature_manifest.json'
    a.assets = a.project_root / 'model/jwarm_training_assets.pt'
    a.scope_root = a.output_root if a.experiment == 'E1' else a.output_root / f'seed{a.seed}' / f'pct{a.percent}'
    a.out_dir = a.scope_root / f'{a.variant}_seed{a.seed}'
    a.subset_dir = a.datasets_dir / f'seed{a.seed}' / f'pct{a.percent}'
    defaults = t.parse_args(['frontier'])
    for key in TRAINING_KEYS:
        if not hasattr(a, key):
            setattr(a, key, getattr(defaults, key))
    a.arm = a.variant
    t.require(min(a.threads, a.batch_size, a.solver_max_iter) > 0, 'Positive sizes required')
    protected = [a.project_root / x for x in ['model', 'train', 'feature_cache', 'new_data']]
    protected += [HERE.parent / x for x in ['E1', 'E2', 'E3', 'E4']]
    for path in protected:
        t.require(not (a.output_root == path or path in a.output_root.parents or a.output_root in path.parents), 'Output overlaps protected inputs or experiment results')
    return a

def code_hashes():
    paths = [Path(__file__), HERE / 'experiment_helpers.py', HERE / 'training_steps.py', Path(t.__file__), Path(m.__file__)]
    return {p.name: t.digest(p) for p in paths}

def ensure_signature(path, payload):
    signature = t.identity(payload)
    if path.exists():
        t.require(t.read(path)['signature'] == signature, f'Input/config changed: {path}; use a new output directory')
    else:
        t.write_json(path, dict(signature=signature, **payload))
    return signature

class Inputs:

    def __init__(self, a, prepare=True):
        self.a = a
        self.source = a.out_dir / 'pretraining'
        self.training = {k: getattr(a, k) for k in TRAINING_KEYS}
        self.training['lambda_ratios'] = list(t.RATIOS)
        self.feature_manifest = t.read(a.feature_manifest)
        t.require([self.feature_manifest['chemical_features'][i] for i in m.MOTIF_INDICES] == ['motif_min_abs_offset', 'motif_center_enrichment', 'motif_count_log'], 'Chemical motif schema changed')
        assets = t.load_pt(a.assets)
        t.require(assets['format'] == 'jwarm_training_assets_v1', 'Wrong candidate asset format')
        if a.experiment == 'E1':
            self.raw = t.SparseStore(a.store_dir)
            self.data = t.FeatureData(self.raw, a.tensor, a.metadata)
            self.splits = {s: self.data.split(s) for s in ['train', 'val']}
            self.ids = assets['feature_ids']
            self.statistics = {'requested_percent': 100}
        else:
            _, _, self.raw, self.data, self.splits, self.ids, _, self.statistics = h.load_subset(a)
        t.require(len(self.ids) == 5845 and len(torch.unique(self.ids)) == 5845, 'Wrong candidate universe')
        for split in self.splits.values():
            t.validate_alignment(self.data, self.raw, split)
        proteins = {s: set(self.raw.metadata.iloc[d.sparse_rows].base_protein_id) for s, d in self.splits.items()}
        t.require(not proteins['train'] & proteins['val'], 'TRAIN/Validation protein overlap')
        self.fingerprint = t.input_fingerprint(a, self.raw)
        t.require(assets['sparse_config_sha256'] == self.fingerprint['sparse_config']['sha256'], 'Candidate/cache provenance changed')
        for key, expected in assets['train_hashes'].items():
            t.require(self.fingerprint[key]['sha256'] == expected, 'Full TRAIN provenance changed')
        if a.experiment == 'E3':
            for s in ['train', 'val']:
                for suffix in ['.tsv', '_metadata.tsv', '_row_index.tsv']:
                    path = a.subset_dir / f'{s}{suffix}'
                    self.fingerprint['subset_' + s + suffix] = {'sha256': t.digest(path), 'bytes': path.stat().st_size}
        self.reference = None
        if a.experiment == 'E3':
            self.reference = self.check_reference()
        positive_fraction = float(np.clip(self.splits['train'].labels.mean(), 1e-06, 1 - 1e-06))
        self.prior = float(np.log(positive_fraction / (1 - positive_fraction)))
        self.ensemble, self.ensemble_mapping = (None, None)
        if FLAGS[a.variant][-1]:
            path = a.project_root / 'feature_cache/gussidr_ensemble_features.tsv'
            t.require(path.is_file(), 'CE requires full38 sequence-aware ensemble cache')
        if not prepare:
            return
        directory = a.scope_root / '_prepared'
        with t.exclusive(directory / 'process.lock'):
            prep_signature = ensure_signature(directory / 'input_manifest.json', dict(inputs=self.fingerprint, code=code_hashes(), statistics=self.statistics, normalization='selected_TRAIN_only', candidate_universe='fixed_full_TRAIN5845', test_loaded=False))
            norm = directory / 'normalization.pt'
            if norm.exists():
                saved = t.load_pt(norm)
                t.require(saved['signature'] == prep_signature, 'Normalization signature mismatch')
            else:
                if a.experiment == 'E1':
                    scales = assets['train_p95_scales'].clone()
                    unseen = []
                else:
                    frame = h.fit_subset_p95(h, self.raw, self.splits['train'], self.ids, directory / 'subset_p95.tsv', a.batch_size)
                    scales = torch.tensor(frame.scale_p95.to_numpy(), dtype=torch.float32)
                    unseen = np.flatnonzero(frame.unseen_in_subset.to_numpy()).tolist()
                saved = dict(signature=prep_signature, scales=scales, scalers=t.fit_preprocessing(self.data, self.splits['train']), unseen=unseen)
                t.save_pt(norm, saved)
            self.scales, self.scalers, self.unseen = (saved['scales'], saved['scalers'], saved['unseen'])
            t.require(self.scales.shape == (5845,) and torch.isfinite(self.scales).all() and (self.scales > 0).all(), 'Invalid candidate P95 scales')
            self.compact = {}
            for s, split in self.splits.items():
                t.make_compact(self.raw, split, self.ids, self.scales, directory / s, a.batch_size)
                self.compact[s] = t.CompactStore(directory / s)
                meta = split.metadata.copy()
                meta['window_id'] = split.window_ids
                meta['base_protein_id'] = self.raw.metadata.iloc[split.sparse_rows].base_protein_id.to_numpy()
                meta.to_csv(directory / f'{s}_metadata.tsv', sep='\t', index=False)
        if FLAGS[a.variant][-1]:
            self.ensemble, self.ensemble_mapping = h.ensemble_inputs(a.project_root, self.data, self.splits)

    def check_reference(self):
        a = self.a
        source = a.e1_results / f'{a.variant}_seed{a.seed}'
        t.require((source / 'complete.json').exists(), f'Completed full-training-set J_WARM E1 outputs are required: {source}')
        cfg = t.read(source / 'config.json')
        done = t.read(source / 'complete.json')
        t.require(cfg['architecture'] == 'J_WARM' and cfg['experiment'] == 'E1' and (cfg['signature'] == done['signature']), 'Expected completed J_WARM E1 outputs with matching configuration signatures')
        t.require(cfg['training'] == self.training and cfg['code'] == code_hashes(), '100% protocol/code mismatch')
        prior = t.read(a.e1_results / '_prepared/input_manifest.json')['inputs']
        for key, value in self.fingerprint.items():
            if not key.startswith('subset_'):
                t.require(prior[key]['sha256'] == value['sha256'], 'E1/E3 representation mismatch: ' + key)
        point = source / 'points' / ('explicit.json' if a.variant == 'D' else 'full_replay.json')
        return dict(source=str(source), point_sha256=t.digest(point), config_sha256=t.digest(source / 'config.json'), metric=t.read(point), endpoint='full_uncompressed_not_best_K')

    def batch(self, s, take, device):
        modalities = self.data.raw_batch(self.splits[s].tensor_rows[take], device)
        if self.ensemble is not None:
            modalities = (*modalities, torch.tensor(self.ensemble[s][take], device=device))
        return (modalities, self.compact[s].events(take, device))

    def model(self, a, device, warm=True):
        model = JointAblation(self.scalers, a.variant, len(self.ids), a.seed, self.prior, a.motif_dropout)
        if warm:
            for name, path in self.warm_paths.items():
                getattr(model, name).load_state_dict(t.load_pt(path)['state'], strict=True)
        return model.to(device)

def pretrain(a, inputs, device):
    flags = dict(zip(NAMES, FLAGS[a.variant]))
    splits = inputs.splits
    shared = a.scope_root / '_pretrain' / f'seed{a.seed}'
    warm = {}
    seconds = 0.0
    sig = t.identity(dict(inputs=inputs.fingerprint, training=inputs.training, code=code_hashes(), seed=a.seed))
    z = {s: np.full(len(split), inputs.prior, dtype=np.float64) for s, split in splits.items()}
    if any(FLAGS[a.variant][:3]):
        name = 'intrinsic_' + ''.join((str(int(x)) for x in FLAGS[a.variant][:3]))
        directory = shared / name
        t.seed_all(a.seed)
        model = h.IntrinsicAblation(m, flags, inputs.scalers, a.motif_dropout).to(device)
        forward = lambda mod, s, rows: mod(*inputs.data.raw_batch(splits[s].tensor_rows[rows], device)[:4])
        ck = steps.fit_stage(model, forward, splits, directory, a, device, stage=name, epochs=a.base_epochs, patience=a.base_patience, lr=a.base_lr, wd=a.base_wd, stage_signature=t.identity([sig, name]))
        z = {s: t.predict(model, forward, split, a.batch_size) for s, split in splits.items()}
        warm['intrinsic'] = directory / 'best.pt'
        seconds += ck['seconds']
        del model
        t.release(device)
    if flags['accessibility']:
        directory = shared / 'access'
        t.seed_all(a.seed + 100)
        model = m.AccessibilityResidual().to(device)
        arrays = {s: (inputs.data.arrays['accessibility'][split.tensor_rows] - inputs.scalers['accessibility']['mean'].numpy()) / inputs.scalers['accessibility']['scale'].numpy() for s, split in splits.items()}
        forward = lambda mod, s, rows: torch.tensor(z[s][rows], device=device, dtype=torch.float32) + mod(torch.tensor(arrays[s][rows], device=device))
        ck = steps.fit_stage(model, forward, splits, directory, a, device, stage='access', epochs=a.access_epochs, patience=a.residual_patience, lr=a.access_lr, wd=a.base_wd, stage_signature=t.identity([sig, 'access', t.digest(warm['intrinsic'])]))
        z = {s: t.predict(model, forward, split, a.batch_size) for s, split in splits.items()}
        warm['access'] = directory / 'best.pt'
        seconds += ck['seconds']
        del model
        t.release(device)
    if flags['ensemble']:
        directory = shared / 'ensemble'
        t.seed_all(a.seed + 150)
        model = h.EnsembleResidual().to(device)
        forward = lambda mod, s, rows: torch.tensor(z[s][rows], device=device, dtype=torch.float32) + mod(torch.tensor(inputs.ensemble[s][rows], device=device))
        ck = steps.fit_stage(model, forward, splits, directory, a, device, stage='ensemble', epochs=a.access_epochs, patience=a.residual_patience, lr=a.access_lr, wd=a.base_wd, stage_signature=t.identity([sig, 'ensemble', inputs.ensemble_mapping, t.digest(warm['access'])]))
        z = {s: t.predict(model, forward, split, a.batch_size) for s, split in splits.items()}
        warm['ensemble'] = directory / 'best.pt'
        seconds += ck['seconds']
        del model
        t.release(device)
    inputs.explicit, inputs.warm_paths = (z, warm)
    t.write_json(inputs.source / 'explicit_info.json', dict(seconds=seconds, paths=warm, signature=sig))

def run(a, inputs, device):
    out = a.out_dir
    payload = dict(architecture='J_WARM', experiment=a.experiment, variant=a.variant, seed=a.seed, percent=100 if a.experiment == 'E1' else a.percent, training=inputs.training, solver_max_iter=a.solver_max_iter, code=code_hashes(), inputs=inputs.fingerprint, statistics=inputs.statistics, flags=dict(zip(NAMES, FLAGS[a.variant])), sparse=a.variant != 'D', stage_C=a.experiment == 'E1' and a.variant != 'D', stage_D=False, test_evaluated=False, objective='classification_BCE_only', ensemble=inputs.ensemble_mapping, reference_100pct=inputs.reference)
    signature = ensure_signature(out / 'config.json', payload)
    if (out / 'complete.json').exists():
        t.require(t.read(out / 'complete.json')['signature'] == signature, 'Stale completion')
        return
    if not (out / 'train_complete.json').exists():
        pretrain(a, inputs, device)
        if a.variant == 'D':
            point = {'point_id': 'explicit', 'K': 0, 'eligible': True, 'converged': None, 'lambda_ratio': 0.0}
            for s, split in inputs.splits.items():
                point[s] = t.metrics(split.labels, inputs.explicit[s], split.negative_types)
                np.save(out / f'explicit_{s}.npy', inputs.explicit[s])
            t.write_json(out / 'points/explicit.json', point)
            (out / 'predictions').mkdir(exist_ok=True)
            np.savez_compressed(out / 'predictions/explicit.npz', window_ids=inputs.splits['val'].window_ids, labels=inputs.splits['val'].labels, logits=inputs.explicit['val'])
            t.write_json(out / 'train_complete.json', dict(signature=signature))
        else:
            model, ck = steps.fit_joint(a, inputs, out, signature, device)
            t.require(not inputs.unseen or bool((model.profile.weights[inputs.unseen] == 0).all()), 'Unseen subset features acquired weights')
            t.export_full(a, inputs, out, model, ck, signature, device)
            full = t.read(out / 'full_metrics.json')
            t.write_json(out / 'points/full_replay.json', dict(point_id='full_replay', K=full['K'], eligible=True, converged=None, lambda_ratio=0.0, **{s: full['metrics'][s]['full'] for s in inputs.splits}))
            del model
            t.release(device)
    else:
        t.require(t.read(out / 'train_complete.json')['signature'] == signature, 'Training signature changed')
    if a.experiment == 'E1' and a.variant != 'D':
        if not (out / 'compress_complete.json').exists():
            t.compress(a, inputs, out, signature)
        else:
            t.require(t.read(out / 'compress_complete.json')['signature'] == signature, 'Compression signature changed')
    rows = []
    for path in (out / 'points').glob('*.json'):
        q = t.read(path)
        rows.append(dict(point=q['point_id'], K=q['K'], eligible=q['eligible'], converged=q['converged'], val_ap=q['val']['ap'], val_auc=q['val']['auc'], train_ap=q['train']['ap']))
    pd.DataFrame(rows).to_csv(out / 'frontier.tsv', sep='\t', index=False)
    t.write_json(out / 'complete.json', dict(signature=signature, architecture='J_WARM', validation_only=True))

def main(argv=None):
    a = parse_args(argv)
    torch.set_num_threads(a.threads)
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s', stream=sys.stdout)
    if a.phase == 'check':
        inputs = Inputs(a, prepare=False)
        logging.info('Checked %s %s seed%d rows=%s; no training/statistics fitting', a.experiment, a.variant, a.seed, {s: len(d) for s, d in inputs.splits.items()})
        return
    a.out_dir.mkdir(parents=True, exist_ok=True)
    logging.getLogger().addHandler(logging.FileHandler(a.out_dir / 'run.log'))
    with t.exclusive(a.out_dir / 'process.lock'):
        inputs = Inputs(a)
        if a.phase == 'run':
            run(a, inputs, t.device_for(a.device))
if __name__ == '__main__':
    main()
