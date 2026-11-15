from __future__ import annotations
import argparse
import copy
import logging
from pathlib import Path
import sys
from functools import partial
import numpy as np
import scipy
import sklearn
import pandas as pd
import torch
HERE = Path(__file__).resolve().parent
ROOT = Path(__import__('os').environ['PRIORSLIM_DATA_ROOT'])
import run_e1_e3 as e1
import training_steps as steps
import train_cv10 as cv
from experiment_helpers import t, baseline as m, ensemble_inputs
ARMS = ['A', 'B', 'C', 'D', 'B1', 'B2', 'B3', 'CE', 'S', 'J_SCRATCH', 'H3', 'AROMATIC_ADD', 'AROMATIC_SPLIT']
SEEDS = [2026, 2027, 2028, 2029, 2030]

def arguments(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('phase', choices=['prepare', 'smoke', 'train', 'compress', 'refine'])
    p.add_argument('--arm', choices=ARMS, default='C')
    p.add_argument('--seed', type=int, choices=SEEDS, default=2026)
    p.add_argument('--device', choices=['cpu', 'mps', 'cuda', 'auto'], default='mps')
    p.add_argument('--output-root', type=Path, default=HERE / 'results')
    p.add_argument('--resume', action='store_true')
    p.add_argument('--solver-max-iter', type=int, default=10000)
    opts = p.parse_args(argv)
    a = t.parse_args(['frontier', '--seed', str(opts.seed), '--device', opts.device])
    for k, v in vars(opts).items():
        setattr(a, k, v)
    a.output_root = a.output_root.resolve()
    t.require(not a.output_root.is_relative_to(ROOT), 'Output must be outside input data')
    a.project_root = ROOT
    a.cv_dir = ROOT / 'new_data/cv10_homology_seed2026'
    a.fold = 6
    a.scope_root = a.output_root
    a.out_dir = a.scope_root / f'{a.arm}_seed{a.seed}'
    a.variant = a.arm if a.arm in e1.FLAGS else 'C'
    a.experiment = 'E1'
    a.evaluate_test = False
    if a.arm == 'J_SCRATCH':
        a.profile_epochs = a.base_epochs
        a.residual_patience = a.base_patience
    return a

def code_hashes():
    return {**e1.code_hashes(), Path(__file__).name: t.digest(__file__), 'train_cv10.py': t.digest(cv.__file__)}

class Sequential(e1.JointAblation):

    def freeze_explicit(self):
        for name, parameter in self.named_parameters():
            parameter.requires_grad_(name.startswith('profile.'))
        return self

    def train(self, mode=True):
        super().train(mode)
        for module in (self.intrinsic, self.access):
            if module is not None:
                module.eval()
        return self

    def parameter_groups(self, training):
        return [dict(params=self.profile.parameters(), lr=training['profile_lr'], weight_decay=training['profile_wd'])]

class FoldInputs(e1.Inputs):

    def __init__(self, a):
        self.a = a
        self.source = a.out_dir / 'pretraining'
        self.training = {k: getattr(a, k) for k in e1.TRAINING_KEYS}
        self.training['lambda_ratios'] = list(t.RATIOS)
        self.feature_manifest = t.read(a.feature_manifest)
        cv.verify_fold_contract(a)
        self.raw = t.SparseStore(a.store_dir)
        self.data = t.FeatureData(self.raw, a.tensor, a.metadata)
        self.splits = cv.load_fold_splits(a, self.data, self.raw)
        self.reference = None
        self.statistics = dict(requested_percent=100, fold=6, test_used=False)
        source = ROOT / 'train/outputs_jwarm_cv10/seed42/fold_06/_prepared'
        old = t.read(source / 'input_manifest.json')
        paths = dict(tensor=a.tensor, metadata=a.metadata, feature_manifest=a.feature_manifest, split_manifest=a.cv_dir / 'manifest.json', sparse_config=a.store_dir / 'motifae_sparse_config.json')
        paths.update({f'sparse_{k}': v for k, v in self.raw.paths.items()})
        for s in ('train', 'val'):
            for suffix in ('.tsv', '_metadata.tsv'):
                paths[s + suffix] = a.cv_dir / 'fold_06' / (s + suffix)
        self.fingerprint = {k: dict(sha256=t.digest(p), bytes=p.stat().st_size) for k, p in paths.items()}
        t.require(self.fingerprint == old['inputs'], 'Fold6 inputs changed: cannot reuse TRAIN assets')
        assets = t.load_pt(source / 'fold_training_assets.pt')
        t.require(assets['signature'] == old['signature'], 'Asset signature mismatch')
        t.require(assets['train_window_ids'] == self.splits['train'].window_ids.tolist(), 'Assets fitted on other TRAIN rows')
        self.ids, self.scales = (assets['feature_ids'], assets['train_p95_scales'])
        self.scalers = copy.deepcopy(assets['preprocessing'])
        self.unseen = []
        t.require(len(self.ids) == len(torch.unique(self.ids)) and (self.scales > 0).all(), 'Invalid candidates/scales')
        self.ensemble = self.ensemble_mapping = None
        self.compact = {}
        prep = a.scope_root / '_prepared'
        prep.mkdir(parents=True, exist_ok=True)
        record = dict(inputs=self.fingerprint, fold=6, candidate_count=len(self.ids), normalization='fold06_TRAIN_only', source_asset_sha256=t.digest(source / 'fold_training_assets.pt'), split_paths={s: str(a.cv_dir / 'fold_06' / f'{s}_metadata.tsv') for s in self.splits}, test_used=False)
        e1.ensure_signature(prep / 'input_manifest.json', record)
        for s, split in self.splits.items():
            directory = source / s
            info = t.read(directory / 'complete.json')
            t.require(info['rows'] == len(split) and info['candidate_count'] == len(self.ids), 'Compact schema mismatch')
            self.compact[s] = store = t.CompactStore(directory)
            sample = np.unique([0, len(split) // 2, len(split) - 1])
            local = self.raw.get_batch(split.sparse_rows[sample])
            lookup = np.full(40960, -1, np.int64)
            lookup[self.ids.numpy()] = np.arange(len(self.ids))
            feature = lookup[local['feature_indices'].numpy()]
            valid = (feature >= 0) & (local['feature_values'].numpy() != 0)
            residue = local['residue_indices'].numpy()[valid]
            expected_col = feature[valid] * 29 + residue % 29
            expected_val = (local['feature_values'].numpy()[valid] / (self.scales.numpy()[feature[valid]] + 1e-08)).astype(np.float32)
            col, val, owner, n = store.events(sample)
            t.require(np.array_equal(col.numpy(), expected_col) and np.array_equal(val.numpy(), expected_val) and np.array_equal(owner.numpy(), residue // 29), 'Compact event replay mismatch')
            link = prep / s
            if link.exists():
                t.require(link.resolve() == directory.resolve(), 'Wrong prepared symlink')
            else:
                link.symlink_to(directory, target_is_directory=True)
            metapath = prep / f'{s}_metadata.tsv'
            if not metapath.exists():
                split.metadata.to_csv(metapath, sep='\t', index=False)
            else:
                frame = pd.read_csv(metapath, sep='\t', low_memory=False)
                t.require(np.array_equal(frame.window_id.astype(str), split.window_ids), 'Prepared metadata differs')
        self.prior = float(np.log(self.splits['train'].labels.mean() / (1 - self.splits['train'].labels.mean())))
        self.custom_intrinsic = self.custom_joint = None
        if a.arm == 'H3':
            self.add_h3()
        if a.arm.startswith('AROMATIC_'):
            self.add_aromatic(a.arm.removeprefix('AROMATIC_'))
        if a.arm == 'CE':
            self.ensemble, self.ensemble_mapping = ensemble_inputs(ROOT, self.data, self.splits)
        logging.info('Fold6 %s seed=%d TRAIN=%d Validation=%d candidates=%d; no Test split', a.arm, a.seed, len(self.splits['train']), len(self.splits['val']), len(self.ids))

    def add_h3(self):
        from model_h3 import H3Intrinsic, H3Joint, FEATURE_NAMES
        extra = np.array([[seq.count(aa) / 29.0 for aa in 'SYF'] for seq in self.data.combined.sequence.astype(str)], np.float32)
        self.extend_chemical(extra, FEATURE_NAMES)
        self.custom_intrinsic, self.custom_joint = (H3Intrinsic, H3Joint)

    def extend_chemical(self, extra, names):
        from sklearn.preprocessing import StandardScaler
        fit = StandardScaler().fit(extra[self.splits['train'].tensor_rows])
        self.data.arrays['chemical'] = np.column_stack([self.data.arrays['chemical'], extra]).astype(np.float32)
        for field, values in [('mean', fit.mean_), ('scale', fit.scale_)]:
            self.scalers['chemical'][field] = torch.cat([self.scalers['chemical'][field], torch.tensor(values, dtype=torch.float32)])
        self.feature_manifest = copy.deepcopy(self.feature_manifest)
        self.feature_manifest['chemical_features'] += names

    def add_aromatic(self, variant):
        from model_aromatic import AromaticIntrinsic, AromaticJoint, EXTRA_NAMES
        directory = ROOT / 'prior_compression_experiments/AROMATIC_SPLIT/feature_cache'
        manifest = t.read(directory / 'manifest.json')
        cv.verify_file(directory / 'features.npz', manifest['values_sha256'])
        for path, digest in manifest['contract']['files'].items():
            cv.verify_file(Path(path), digest)
        with np.load(directory / 'features.npz') as z:
            ids = np.concatenate([s.window_ids for s in self.splits.values()])
            t.require(np.array_equal(z['window_ids'], ids), 'Aromatic fold window order differs')
            rows = np.concatenate([s.tensor_rows for s in self.splits.values()])
            extra = np.zeros((len(self.data.combined), 18), np.float32)
            extra[rows] = z['values']
            names = self.feature_manifest['chemical_features']
            pooled_names = [f'{metric}_aromatic{suffix}' for metric in ('density', 'X', 'center_enrich') for suffix in ('', '_percentile_IDR')]
            removed = [names.index(name) for name in pooled_names]
            t.require(np.max(np.abs(z['pooled'] - self.data.arrays['chemical'][rows][:, removed])) < 1e-06, 'Pooled aromatic descriptors differ from the input feature arrays')
        self.extend_chemical(extra, EXTRA_NAMES)
        layout = [i for i in range(76) if variant == 'ADD' or i not in removed]
        self.scalers['aromatic_layout'] = layout
        self.data.arrays['chemical'] = self.data.arrays['chemical'][:, layout].copy()
        self.scalers['chemical'] = {k: v[layout].clone() for k, v in self.scalers['chemical'].items()}
        self.feature_manifest['chemical_features'] = [self.feature_manifest['chemical_features'][i] for i in layout]
        self.custom_intrinsic, self.custom_joint = (AromaticIntrinsic, AromaticJoint)

    def model(self, a, device, warm=True):
        if self.custom_joint:
            model = self.custom_joint(self.scalers, len(self.ids), a.seed, a.motif_dropout)
        else:
            cls = Sequential if a.arm == 'S' else e1.JointAblation
            model = cls(self.scalers, a.variant, len(self.ids), a.seed, self.prior, a.motif_dropout)
        if warm and a.arm != 'J_SCRATCH':
            for name, path in self.warm_paths.items():
                getattr(model, name).load_state_dict(t.load_pt(path)['state'], strict=True)
        if a.arm == 'S':
            model.freeze_explicit()
        return model.to(device)

def pretrain(a, inputs, device):
    if a.arm == 'J_SCRATCH':
        model = inputs.model(a, device, warm=False)
        inputs.explicit = {s: t.predict(model, lambda net, role, rows: net(*inputs.batch(role, rows, device)), split, a.batch_size) for s, split in inputs.splits.items()}
        inputs.warm_paths = {}
        t.write_json(inputs.source / 'explicit_info.json', dict(seconds=0.0, paths={}, initialization='scratch'))
        del model
        t.release(device)
    elif inputs.custom_intrinsic:
        original = t.IntrinsicEncoder
        t.IntrinsicEncoder = inputs.custom_intrinsic
        try:
            t.run_stage_a(a, inputs, t.identity([code_hashes(), inputs.fingerprint, a.arm, a.seed]), device)
        finally:
            t.IntrinsicEncoder = original
        inputs.warm_paths = {name: inputs.source / sub / 'best.pt' for name, sub in [('intrinsic', 'intrinsic'), ('access', 'access')]}
    else:
        original_training, original_a = (inputs.training, copy.copy(a))
        try:
            e1.pretrain(a, inputs, device)
        finally:
            inputs.training = original_training
    for s, z in inputs.explicit.items():
        inputs.source.mkdir(parents=True, exist_ok=True)
        np.save(inputs.source / f'explicit_{s}.npy', z)

def contract(a, inputs):
    payload = dict(architecture='J_WARM', display_name='PriorSLiM', experiment='E1', fold=6, variant=a.variant, arm=a.arm, seed=a.seed, training=inputs.training, inputs=inputs.fingerprint, code=code_hashes(), objective='classification_BCE_only', candidate_count=len(inputs.ids), test_evaluated=False, split_paths={s: str(a.cv_dir / 'fold_06' / f'{s}_metadata.tsv') for s in inputs.splits}, strategy='frozen_warm_explicit' if a.arm == 'S' else 'joint_scratch' if a.arm == 'J_SCRATCH' else 'joint_warm', chemical_names=inputs.feature_manifest['chemical_features'], ensemble=inputs.ensemble_mapping)
    return e1.ensure_signature(a.out_dir / 'config.json', payload)

def smoke(a, inputs, device):
    model = inputs.model(a, device, warm=False)
    rows = np.arange(4)
    modalities, events = inputs.batch('val', rows, device)
    model.eval()
    z, explicit, residual = model(modalities, events, return_parts=True)
    t.require(z.shape == (4,) and torch.isfinite(z).all(), 'Forward shape/finiteness')
    t.require(torch.equal(z, explicit) and torch.count_nonzero(residual) == 0, 'Zero W replay')
    y = torch.tensor(inputs.splits['val'].labels[rows], device=device)
    torch.nn.functional.binary_cross_entropy_with_logits(z, y).backward()
    t.require(torch.isfinite(model.profile.weights.grad).all() and model.profile.weights.grad.abs().sum() > 0, 'Profile gradient')
    if a.arm == 'C':
        reference = m.JointModel(inputs.scalers, len(inputs.ids), a.seed, a.motif_dropout).to(device).eval()
        reference_z = reference(modalities, events)
        t.require(torch.allclose(z, reference_z, atol=1e-06, rtol=0), 'C architecture differs from baseline')
    if a.arm == 'S':
        model.train()
        t.require(not model.intrinsic.training and (not model.access.training), 'Frozen S dropout must be disabled')
        t.require(all((p.grad is None for n, p in model.named_parameters() if not n.startswith('profile.'))), 'S explicit gradients')
    t.write_json(a.scope_root / 'smoke' / f'{a.arm}.json', dict(passed=True, candidate_count=len(inputs.ids), rows={s: len(d) for s, d in inputs.splits.items()}, test_used=False, device=str(device)))

def execute(a, inputs, device):
    out = a.out_dir
    signature = contract(a, inputs)
    if a.phase == 'train':
        if (out / 'train_complete.json').exists():
            t.require(t.read(out / 'train_complete.json')['signature'] == signature, 'Training signature differs')
            logging.info('Reuse completed training %s', out.name)
            return
        pretrain(a, inputs, device)
        if a.arm == 'D':
            q = dict(point_id='explicit', K=0, eligible=True, converged=None, lambda_ratio=0.0)
            for s, split in inputs.splits.items():
                q[s] = t.metrics(split.labels, inputs.explicit[s], split.negative_types)
                np.save(out / f'explicit_{s}.npy', inputs.explicit[s])
            t.write_json(out / 'points/explicit.json', q)
            t.write_json(out / 'train_complete.json', dict(signature=signature))
            t.write_json(out / 'complete.json', dict(signature=signature, validation_only=True))
        else:
            model, ck = steps.fit_joint(a, inputs, out, signature, device)
            t.export_full(a, inputs, out, model, ck, signature, device)
            full = t.read(out / 'full_metrics.json')
            full['candidate_count'] = len(inputs.ids)
            t.write_json(out / 'full_metrics.json', full)
            t.write_json(out / 'points/full_replay.json', dict(point_id='full_replay', K=full['K'], eligible=True, converged=None, lambda_ratio=0.0, **{s: full['metrics'][s]['full'] for s in inputs.splits}))
    elif a.phase == 'compress':
        if a.arm == 'D':
            return
        if not (out / 'compress_complete.json').exists():
            t.compress(a, inputs, out, signature)
        else:
            t.require(t.read(out / 'compress_complete.json')['signature'] == signature, 'Compression signature differs')
        rows = []
        for path in (out / 'points').glob('*.json'):
            q = t.read(path)
            rows.append(dict(point=q['point_id'], K=q['K'], eligible=q['eligible'], converged=q['converged'], val_ap=q['val']['ap'], val_auc=q['val']['auc'], train_ap=q['train']['ap']))
        pd.DataFrame(rows).to_csv(out / 'frontier.tsv', sep='\t', index=False)
        t.write_json(out / 'complete.json', dict(signature=signature, validation_only=True))
    elif a.phase == 'refine':
        t.require(a.arm in ['C', 'S', 'J_SCRATCH'], 'Stage D scoped to E2 strategies')
        points = [t.read(p) for p in (out / 'points').glob('lambda_*.json')]
        selected = cv.choose_automatic_points(points)
        t.write_json(out / 'automatic_k_selection.json', dict(seed=a.seed, fold=6, test_used=False, policy='performance=max converged Validation AP; efficiency=min K within 0.005', selected=selected))
        with cv.dynamic_refinement(len(inputs.ids)):
            for k in sorted({int(p['K']) for p in selected.values()}):
                a.selected_k = k
                t.run_refine(a, inputs, signature, device)

def main(argv=None):
    a = arguments(argv)
    a.out_dir.mkdir(parents=True, exist_ok=True)
    torch.set_num_threads(a.threads)
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s', handlers=[logging.StreamHandler(sys.stdout), logging.FileHandler(a.out_dir / 'run.log')], force=True)
    with t.exclusive(a.out_dir / 'process.lock'):
        inputs = FoldInputs(a)
        if a.phase == 'prepare':
            return
        device = t.device_for(a.device)
        if a.phase == 'smoke':
            smoke(a, inputs, device)
        else:
            execute(a, inputs, device)
if __name__ == '__main__':
    main()
