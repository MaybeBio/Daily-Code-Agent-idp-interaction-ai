from pathlib import Path
from experiment_helpers import t, np, pd, torch
from raw_e4 import RawStore, csr_events
import scipy.sparse as sp
HERE = Path(__file__).resolve().parent
read = t.read
require = t.require

class Split:

    def __init__(self, name, frame):
        require(name in ('train', 'val'), 'E4 is Validation-only: expected train or val split')
        self.name = name
        self.labels = frame.label.to_numpy(dtype=np.float32)
        self.negative_types = frame.negative_type.fillna('').to_numpy(dtype=str)
        self.window_ids = frame.window_id.to_numpy(dtype=str)

    def __len__(self):
        return len(self.labels)

def load_design(directory, n, weights):
    cfg = read(directory / 'complete.json')
    import hashlib
    require(cfg['weight_hash'] == hashlib.sha256(np.ascontiguousarray(weights).view(np.uint8)).hexdigest(), 'Evidence/profile hash mismatch')
    ptr = np.load(directory / 'indptr.npy', mmap_mode='r')
    data = np.memmap(directory / 'values.bin', dtype=np.float64, mode='r') if (directory / 'values.bin').stat().st_size else np.empty(0, np.float64)
    col = np.memmap(directory / 'columns.bin', dtype=np.int32, mode='r') if (directory / 'columns.bin').stat().st_size else np.empty(0, np.int32)
    require(len(ptr) == n + 1 and ptr[-1] == len(data) == len(col), 'Invalid evidence CSR')
    return sp.csr_matrix((data, col, ptr), shape=(n, len(weights)), copy=False)

class ExperimentData:

    def __init__(self, args, with_raw=True):
        self.args = args
        self.source = source = args.e1_results / f'{args.variant}_seed{args.seed}'
        cfg = read(source / 'config.json')
        done = read(source / 'complete.json')
        require(cfg['architecture'] == 'J_WARM' and cfg['experiment'] == 'E1', 'Expected NEW J_WARM E1, not old S')
        require(cfg['signature'] == done['signature'] and done['validation_only'] and (not cfg['test_evaluated']), 'Incomplete/invalid J_WARM source')
        require(cfg['seed'] == args.seed and cfg['variant'] == args.variant, 'Wrong source identity')
        require(cfg['objective'] == 'classification_BCE_only', 'Different objective')
        self.manifest = manifest = read(args.e1_results / '_prepared/input_manifest.json')
        require(manifest['inputs'] == cfg['inputs'], 'Prepared/source inputs differ')
        self.meta = {s: pd.read_csv(args.e1_results / f'_prepared/{s}_metadata.tsv', sep='\t', low_memory=False) for s in ['train', 'val']}
        self.splits = {s: Split(s, frame) for s, frame in self.meta.items()}
        for s, frame in self.meta.items():
            original = args.project_root / f'new_data/{s}_metadata.tsv'
            require(t.digest(original) == cfg['inputs'][f'{s}_metadata.tsv']['sha256'], 'Split metadata changed')
            df = pd.read_csv(original, sep='\t', low_memory=False)
            require(np.array_equal(df.sample_id, frame.sample_id) and np.array_equal(df.label, frame.label), 'Split row/label mismatch')
        require(not set(self.meta['train'].base_protein_id) & set(self.meta['val'].base_protein_id), 'Protein leakage')
        self.explicit = {s: np.load(source / f'explicit_{s}.npy') for s in self.splits}
        frozen = t.load_pt(source / 'frozen_profiles.pt')
        self.ids = frozen['feature_ids'].numpy().astype(np.int64)
        self.scales = frozen['scales'].numpy()
        self.weights = frozen['weights'].numpy().astype(np.float64)
        self.alpha = frozen['alpha_init'].numpy().astype(np.float64)
        require(self.weights.shape == (5845, 29) and len(set(self.ids)) == 5845, 'Wrong candidate universe')
        require(np.allclose(self.alpha, np.linalg.norm(self.weights, axis=1), atol=1e-12), 'Profile amplitudes differ')
        self.design = {s: load_design(source / 'evidence' / s, len(split), self.weights) for s, split in self.splits.items()}
        self.evidence = {s: x.multiply(self.alpha).tocsr() for s, x in self.design.items()}
        full = read(source / 'points/full_replay.json')
        for s, split in self.splits.items():
            z = self.explicit[s] + np.asarray(self.evidence[s].sum(1)).ravel()
            require(np.isfinite(z).all() and abs(t.metrics(split.labels, z)['ap'] - full[s]['ap']) < 1e-10, 'J_WARM full replay failed')
        with np.load(source / 'predictions/full_replay.npz') as stored:
            require(np.array_equal(stored['window_ids'], self.splits['val'].window_ids), 'Validation IDs mismatch')
            require(np.allclose(stored['logits'], self.explicit['val'] + np.asarray(self.evidence['val'].sum(1)).ravel(), atol=1e-08, rtol=0), 'Full logits mismatch')
        self.raw = RawStore(args.store_dir, self.ids, cfg['inputs']) if with_raw else None
        self.vectors = self.bias = None
        if with_raw:
            for s, split in self.splits.items():
                self.raw.align(s, split.window_ids, self.meta[s])
            require(t.digest(args.decoder_checkpoint) == self.raw.config['checkpoint_sha256'], 'Decoder != cached encoder checkpoint')
            state = t.load_pt(args.decoder_checkpoint)
            require(state['decoder.weight'].shape == (1280, 40960) and state['bias'].shape == (1280,), 'Wrong official decoder orientation')
            self.vectors = state['decoder.weight'].T[self.ids].contiguous().float()
            self.bias = state['bias'].clone().float()
        paths = [source / x for x in ['config.json', 'complete.json', 'explicit_train.npy', 'explicit_val.npy', 'frozen_profiles.pt', 'points/full_replay.json', 'joint/best.pt']]
        paths += [args.e1_results / f'_prepared/{s}_metadata.tsv' for s in self.splits]
        for s in self.splits:
            paths += [source / 'evidence' / s / x for x in ['complete.json', 'indptr.npy', 'columns.bin', 'values.bin']]
        prox = args.proximal_results / source.name
        require(read(prox / 'complete.json')['signature'] == cfg['signature'], 'Proximal source not same joint model')
        paths += [f for p in sorted((prox / 'points').glob('lambda_*.json')) for f in [p, prox / 'solutions' / (p.stem + '.pt')]]
        self.fingerprint = dict(source_assets={str(p): t.digest(p) for p in paths}, recorded_inputs=cfg['inputs'], decoder_sha256=t.digest(args.decoder_checkpoint), candidate_ids=t.identity(self.ids.tolist()), test_loaded=False, explicit_source='post_joint_J_WARM_not_pretraining_A')

    def proximal_points(self):
        source = self.args.proximal_results / self.source.name
        points = []
        for path in sorted((source / 'points').glob('lambda_*.json')):
            q = read(path)
            if not (q['converged'] is True and q['eligible'] is True):
                continue
            solution = source / 'solutions' / (path.stem + '.pt')
            alpha = t.load_pt(solution)['alpha'].numpy().astype(float)
            require(alpha.shape == self.alpha.shape and np.isfinite(alpha).all() and (alpha >= 0).all(), 'Invalid alpha')
            require(int((alpha > 1e-10).sum()) == q['K'], 'Support differs')
            for s, split in self.splits.items():
                require(abs(t.metrics(split.labels, self.explicit[s] + self.design[s] @ alpha)['ap'] - q[s]['ap']) < 1e-09, 'Proximal replay failed')
            points.append((path.stem, alpha, q, dict(point_sha256=t.digest(path), solution_sha256=t.digest(solution))))
        require(points, 'No converged proximal points; finish J_WARM E1 first')
        return points
