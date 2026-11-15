import sys
from pathlib import Path
import numpy as np
import scipy
import sklearn
import torch
import pandas as pd
HERE = Path(__file__).resolve().parent
import train_e4 as engine
from e4_data import ExperimentData, Split, load_design, RawStore, t, read, require

class FoldE4Data(ExperimentData):

    def __init__(self, args, with_raw=True):
        self.args = args
        self.source = source = args.e1_results / f'C_seed{args.seed}'
        cfg = read(source / 'config.json')
        done = read(source / 'complete.json')
        require(cfg['fold'] == 6 and cfg['arm'] == 'C' and (cfg['seed'] == args.seed), 'Not fold6 PriorSLiM')
        require(cfg['signature'] == done['signature'] and done['validation_only'] and (not cfg['test_evaluated']), 'Incomplete source')
        require(cfg['objective'] == 'classification_BCE_only', 'Wrong objective')
        self.manifest = read(args.e1_results / '_prepared/input_manifest.json')
        require(self.manifest['inputs'] == cfg['inputs'], 'Prepared/source inputs differ')
        self.meta = {s: pd.read_csv(args.e1_results / f'_prepared/{s}_metadata.tsv', sep='\t', low_memory=False) for s in ['train', 'val']}
        self.splits = {s: Split(s, frame) for s, frame in self.meta.items()}
        for s, frame in self.meta.items():
            path = Path(cfg['split_paths'][s])
            require(t.digest(path) == cfg['inputs'][f'{s}_metadata.tsv']['sha256'], 'Fold split changed')
            original = pd.read_csv(path, sep='\t', low_memory=False)
            for key in ['window_id', 'label', 'component_id', 'window_sequence']:
                require(np.array_equal(original[key].astype(str), frame[key].astype(str)), f'Fold mismatch {key}')
        for key in ['window_id', 'component_id', 'base_protein_id']:
            require(not set(self.meta['train'][key]) & set(self.meta['val'][key]), 'Split leakage')
        frozen = t.load_pt(source / 'frozen_profiles.pt')
        self.ids = frozen['feature_ids'].numpy().astype(np.int64)
        self.scales = frozen['scales'].numpy()
        self.weights = frozen['weights'].numpy().astype(np.float64)
        self.alpha = frozen['alpha_init'].numpy().astype(np.float64)
        require(self.weights.shape == (cfg['candidate_count'], 29) and len(set(self.ids)) == cfg['candidate_count'], 'Candidate schema')
        require(np.allclose(self.alpha, np.linalg.norm(self.weights, axis=1), atol=1e-12), 'Profile amplitudes')
        self.explicit = {s: np.load(source / f'explicit_{s}.npy') for s in self.splits}
        self.design = {s: load_design(source / 'evidence' / s, len(split), self.weights) for s, split in self.splits.items()}
        self.evidence = {s: x.multiply(self.alpha).tocsr() for s, x in self.design.items()}
        full = read(source / 'points/full_replay.json')
        for s, split in self.splits.items():
            z = self.explicit[s] + np.asarray(self.evidence[s].sum(1)).ravel()
            require(np.isfinite(z).all() and abs(t.metrics(split.labels, z)['ap'] - full[s]['ap']) < 1e-10, 'Full replay failed')
        with np.load(source / 'predictions/full_replay.npz') as z:
            require(np.array_equal(z['window_ids'], self.splits['val'].window_ids), 'Replay IDs')
            require(np.allclose(z['logits'], self.explicit['val'] + np.asarray(self.evidence['val'].sum(1)).ravel(), atol=1e-08, rtol=0), 'Replay logits')
        self.raw = RawStore(args.store_dir, self.ids, cfg['inputs']) if with_raw else None
        self.vectors = self.bias = None
        if with_raw:
            for s, split in self.splits.items():
                self.raw.align(s, split.window_ids, self.meta[s])
            require(t.digest(args.decoder_checkpoint) == self.raw.config['checkpoint_sha256'], 'Decoder/cache checkpoint mismatch')
            state = t.load_pt(args.decoder_checkpoint)
            require(state['decoder.weight'].shape == (1280, 40960) and state['bias'].shape == (1280,), 'Decoder orientation')
            self.vectors = state['decoder.weight'].T[self.ids].contiguous().float()
            self.bias = state['bias'].clone().float()
        paths = [source / name for name in ['config.json', 'complete.json', 'explicit_train.npy', 'explicit_val.npy', 'frozen_profiles.pt', 'points/full_replay.json', 'joint/best.pt']]
        paths += [args.e1_results / f'_prepared/{s}_metadata.tsv' for s in self.splits]
        for s in self.splits:
            paths += [source / 'evidence' / s / n for n in ['complete.json', 'indptr.npy', 'columns.bin', 'values.bin']]
        paths += [p for q in sorted((source / 'points').glob('lambda_*.json')) for p in [q, source / 'solutions' / (q.stem + '.pt')]]
        paths.append(Path(__file__))
        self.fingerprint = dict(source_assets={str(p): t.digest(p) for p in paths}, recorded_inputs=cfg['inputs'], decoder_sha256=t.digest(args.decoder_checkpoint), candidate_ids=t.identity(self.ids.tolist()), test_loaded=False, explicit_source='post_joint_PriorSLiM_fold06', fold=6)

def main():
    require(len(sys.argv) == 3, 'Usage: run_e4_fold06.py PHASE SEED')
    phase, seed = sys.argv[1:]
    require(int(seed) in range(2026, 2031), 'Expected five registered seeds')
    engine.ExperimentData = FoldE4Data
    engine.main([phase, '--seed', seed, '--variant', 'C', '--device', 'mps', '--resume', '--e1-results', str(HERE / 'results'), '--proximal-results', str(HERE / 'results'), '--output-root', str(HERE / 'E4'), '--budgets', '64', '128', '256', '512', '1024', '2048', '4096', '5829'])
if __name__ == '__main__':
    main()
