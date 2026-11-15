"""Build the five explicit PriorSLiM arrays using the training accessibility reference."""
import os
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))


import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
import build_all_features as features

def build(windows, pssm_dir, regions, cache_dir, reference):
    features.PSSM_DIR = Path(pssm_dir)
    features.FRAGMENT_FILE = Path(regions)
    features.GUSSIDR_CACHE_FILE = Path(cache_dir) / 'gussidr_ensemble_features.tsv'
    features.NARDINI_CACHE_FILE = Path(cache_dir) / 'nardini_features.tsv'
    features.SKIP_BAD_WINDOWS = False
    if windows.window_id.duplicated().any():
        raise ValueError('window_id must be unique')
    residue, peak, mask, pmeta, _, perrors = features.build_pssm_tensor(windows)
    chemistry, cmeta, cerrors = features.build_chemical_tensors(windows)
    access, ameta, aerrors = features.build_accessibility_tensors(windows, rg_norm_ref=reference)
    for metadata, errors in [(pmeta, perrors), (cmeta, cerrors), (ameta, aerrors)]:
        if len(errors) or not np.array_equal(metadata.window_id.astype(str), windows.window_id.astype(str)):
            raise ValueError('Feature extraction failed or changed window ordering')
    arrays = dict(pssm_residue=residue, pssm_valid_mask=mask, pssm_window=peak, chemical=chemistry, accessibility=access)
    for key, value in arrays.items():
        if not np.isfinite(value).all():
            raise ValueError(f'Non-finite {key}')
    return arrays

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--windows', type=Path, required=True, help='TSV with window_id, sequence, coordinates and IDR bounds')
    parser.add_argument('--pssm-dir', type=Path, required=True)
    parser.add_argument('--regions', type=Path, required=True)
    parser.add_argument('--cache-dir', type=Path, required=True)
    parser.add_argument('--training-manifest', type=Path, required=True)
    parser.add_argument('--output-prefix', type=Path, required=True)
    args = parser.parse_args()
    output = args.output_prefix
    if any(Path(str(output) + suffix).exists() for suffix in ['_tensor.npz','_metadata.tsv','_feature_manifest.json']):
        parser.error('Use a new output prefix')
    windows = pd.read_csv(args.windows, sep='\t', low_memory=False)
    manifest = json.loads(args.training_manifest.read_text())
    arrays = build(windows, args.pssm_dir, args.regions, args.cache_dir, manifest['accessibility_rg_norm_ref'])
    output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(str(output) + '_tensor.npz', **arrays, window_ids=windows.window_id.to_numpy(dtype=str))
    windows = windows.copy()
    windows['tensor_index'] = np.arange(len(windows))
    windows.to_csv(str(output) + '_metadata.tsv', sep='\t', index=False)
    manifest['release_explicit_arrays'] = {k:list(v.shape) for k,v in arrays.items()}
    manifest['explicit_array_names'] = list(arrays)
    Path(str(output) + '_feature_manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(f'Built {len(windows)} windows. Sparse MotifAE inputs must be prepared separately.')

if __name__ == '__main__':
    main()
