"""Predict from aligned feature files or run the small numerical example."""
import argparse
from pathlib import Path
import numpy as np
import torch
import _training as engine
from model import PriorSLiMPredictor

DEFAULT_CHECKPOINT = Path(__file__).resolve().parents[1] / 'data/models/model.pt'


def predict_example(checkpoint, fixture):
    model = PriorSLiMPredictor(engine.load_pt(checkpoint))
    with np.load(fixture, allow_pickle=False) as arrays:
        modalities = tuple(torch.tensor(arrays[k], dtype=torch.float32) for k in
                           ('pssm_residue', 'pssm_valid_mask', 'pssm_window', 'chemical', 'accessibility'))
        sparse = {k[7:]: torch.tensor(arrays[k]) for k in arrays.files if k.startswith('sparse_')}
        with torch.inference_mode():
            logits = model(modalities, sparse).numpy()
        np.testing.assert_allclose(logits, arrays['expected_logit'], rtol=0, atol=1e-5)
        return logits


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--checkpoint', type=Path, default=DEFAULT_CHECKPOINT,
                   help='Override the default data/models/model.pt for experiment reproduction')
    p.add_argument('--example', action='store_true')
    p.add_argument('--fixture', type=Path, default=Path(__file__).resolve().parents[1] / 'examples/synthetic_windows.npz')
    p.add_argument('--tensor', type=Path)
    p.add_argument('--metadata', type=Path)
    p.add_argument('--feature-manifest', type=Path)
    p.add_argument('--store-dir', type=Path)
    p.add_argument('--output-prefix', type=Path)
    p.add_argument('--device', choices=['cpu', 'cuda', 'mps', 'auto'], default='cpu')
    p.add_argument('--batch-size', type=int, default=128)
    p.add_argument('--threads', type=int, default=2)
    args = p.parse_args(argv)
    if not args.example:
        missing = [k for k in ('tensor', 'metadata', 'feature_manifest', 'store_dir', 'output_prefix') if getattr(args, k) is None]
        if missing:
            p.error('Required for feature-file prediction: ' + ', '.join('--' + k.replace('_', '-') for k in missing))
    return args


def main():
    args = parse_args()
    torch.set_num_threads(args.threads)
    if args.example:
        logits = predict_example(args.checkpoint, args.fixture)
        print(f'PASS: {len(logits)} synthetic windows; {args.checkpoint.name}; absolute tolerance 1e-5')
        return
    args.split = 'all'
    engine.predict_command(args, engine.device_for(args.device))


if __name__ == '__main__':
    main()
