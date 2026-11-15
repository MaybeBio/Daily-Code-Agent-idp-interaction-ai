"""Rebuild experiment features from the released PSSMs and sequence-level inputs."""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import shutil
import sys
from pathlib import Path

FEATURE_CODE = Path(__file__).resolve().parent / 'features'
sys.path.insert(0, str(FEATURE_CODE))


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024**2), b''):
            h.update(block)
    return h.hexdigest()


def under(root, relative):
    path = (root / relative).resolve()
    if not path.is_relative_to(root.resolve()):
        raise ValueError(f'Path outside data root: {relative}')
    return path


def read_sequences(path):
    sequences, name = {}, None
    with gzip.open(path, 'rt') as handle:
        for line in handle:
            line = line.strip()
            if line.startswith('>'):
                name = line[1:]
                if name in sequences:
                    raise ValueError(f'Duplicate sequence: {name}')
                sequences[name] = ''
            elif line:
                if name is None:
                    raise ValueError('Sequence before FASTA header')
                sequences[name] += line
    return sequences


def verify_inputs(root, jobs):
    import csv
    checked = set()
    cohorts = {job['experiment'] for job in jobs}
    for job in jobs:
        inputs = under(root, job['inputs'])
        for name, expected in job['files'].items():
            path = under(inputs, name)
            if path not in checked:
                if digest(path) != expected:
                    raise ValueError(f'Input checksum mismatch: {path}')
                checked.add(path)
    with (root / 'pssm_resources/files.tsv').open() as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            if Path(row['path']).parts[0] in cohorts:
                path = under(root / 'pssm_resources', row['path'])
                if digest(path) != row['sha256']:
                    raise ValueError(f'PSSM checksum mismatch: {path}')
                checked.add(path)
    print(f'Verified {len(checked)} input files for {len(jobs)} feature stores', flush=True)


def validate_sequences(windows, context, chunks, sequences):
    if windows.window_id.duplicated().any():
        raise ValueError('Duplicate scoring-window IDs')
    if list(windows.window_id) != list(context.window_id):
        raise ValueError('Scoring windows and contextual mapping differ')
    if set(windows.uniprot) != set(sequences):
        raise ValueError('Sequence entities differ from scoring windows')
    for row in windows.itertuples(index=False):
        seq = sequences[row.uniprot]
        expected = ''.join(seq[i-1] if 1 <= i <= len(seq) else 'X'
                           for i in range(row.window_start, row.window_end + 1))
        if expected != row.sequence:
            raise ValueError(f'Window/sequence mismatch: {row.window_id}')
    for row in chunks.itertuples(index=False):
        seq = sequences[row.protein_id]
        crop = seq[int(row.crop_start)-1:int(row.crop_end)]
        if len(seq) != row.protein_length or hashlib.sha256(crop.encode()).hexdigest() != row.sequence_sha256:
            raise ValueError(f'Crop/sequence mismatch: {row.protein_id}, chunk {row.chunk_id}')


def code_signature():
    paths = [Path(__file__), *FEATURE_CODE.glob('*.py'), *FEATURE_CODE.glob('*.json')]
    return {p.name: digest(p) for p in sorted(paths)}


def build_job(root, job, args, code_hashes):
    import numpy as np
    import pandas as pd
    from build_priorslim_features import build
    from build_contextual_esm2 import generate_contextual_embeddings
    from build_motifae_sparse_windows import build_motifae_sparse_windows

    inputs = under(root, job['inputs'])
    output = under(root, job['output'])
    sparse = under(root, job['sparse_output'])
    marker = output / 'feature_rebuild.json'
    signature = dict(inputs=job['files'], pssm_dir=job['pssm_dir'], code=code_hashes)
    state = json.loads(marker.read_text()) if marker.exists() else {}
    if state and state['signature'] != signature:
        raise ValueError(f'{job["id"]}: inputs or code changed; use a new data directory')
    if not state and (output / 'combined_tensor.npz').exists():
        raise FileExistsError(f'{output}: existing features have no rebuild record; use a new data directory')
    if state.get('explicit_complete'):
        for name, expected in state['outputs'].items():
            if digest(under(root, name)) != expected:
                raise ValueError(f'Completed output changed: {name}')
    if state.get('explicit_complete') and (args.explicit_only or state.get('sparse_complete')):
        print(f'Reusing {job["id"]}', flush=True)
        return
    windows = pd.read_csv(inputs / 'windows.tsv.gz', sep='\t', low_memory=False)
    context = pd.read_csv(inputs / 'context_windows.tsv.gz', sep='\t', low_memory=False)
    chunks = pd.read_csv(inputs / 'context_chunks.tsv', sep='\t', low_memory=False)
    sequences = read_sequences(inputs / 'sequences.fasta.gz')
    validate_sequences(windows, context, chunks, sequences)
    if len(windows) != job['windows']:
        raise ValueError(f'{job["id"]}: window count differs')
    output.mkdir(parents=True, exist_ok=True)
    state = state or dict(signature=signature, outputs={})

    def save_state():
        temp = marker.with_suffix('.tmp')
        temp.write_text(json.dumps(state, indent=2) + '\n')
        os.replace(temp, marker)

    def record(paths):
        for path in paths:
            state['outputs'][str(path.relative_to(root))] = digest(path)

    save_state()
    if not state.get('explicit_complete'):
        manifest = json.loads((inputs / 'feature_manifest.json').read_text())
        work = windows.copy()
        work['protein_sequence'] = work.uniprot.map(sequences)
        arrays = build(work, under(root, job['pssm_dir']), inputs / 'regions.xlsx', inputs,
                       manifest['accessibility_rg_norm_ref'])
        tensor = output / 'combined_tensor.npz'
        temporary = output / 'combined_tensor.partial.npz'
        np.savez_compressed(temporary, **arrays, window_ids=windows.window_id.to_numpy(dtype=str))
        os.replace(temporary, tensor)
        metadata = output / 'combined_metadata.tsv'
        template = inputs / 'combined_metadata.tsv.gz'
        if template.exists():
            with gzip.open(template, 'rb') as source, metadata.open('wb') as target:
                shutil.copyfileobj(source, target)
        else:
            windows.to_csv(metadata, sep='\t', index=False)
        manifest['release_explicit_arrays'] = {key: list(value.shape) for key, value in arrays.items()}
        manifest_path = output / 'combined_feature_manifest.json'
        manifest_path.write_text(json.dumps(manifest, indent=2) + '\n')
        record([tensor, metadata, manifest_path])
        state['explicit_complete'] = True
        save_state()
    if args.explicit_only:
        return

    # Use the archived crop plan so shard boundaries cannot change sequence context.
    cache = sparse.parent
    embeddings = cache / 'contextual_esm2_idr_coverage'
    for frame in (context, chunks):
        frame['embedding_path'] = frame.embedding_file.map(lambda name: str(under(embeddings, name)))
    cache.mkdir(parents=True, exist_ok=True)
    context_path = cache / 'contextual_esm2_idr_coverage_metadata.tsv'
    chunks_path = cache / 'contextual_esm2_idr_coverage_chunks.tsv'
    context.to_csv(context_path, sep='\t', index=False)
    chunks.to_csv(chunks_path, sep='\t', index=False)
    generate_contextual_embeddings(chunks, context, sequences, embeddings,
                                   device_name=args.device, batch_size=args.batch_size)
    build_motifae_sparse_windows(context_path, chunks_path, args.checkpoint, sparse,
                                device_name=args.device)
    mapping = pd.read_csv(sparse / 'motifae_window_metadata.tsv', sep='\t', usecols=['window_id'])
    if set(mapping.window_id) != set(windows.window_id) or len(mapping) != len(windows):
        raise ValueError('Sparse store and explicit windows differ')
    record([context_path, chunks_path, *sorted(sparse.glob('motifae_*'))])
    state['sparse_complete'] = True
    state['device'] = args.device
    state['motifae_checkpoint_sha256'] = digest(args.checkpoint)
    save_state()
    if not args.keep_esm:
        shutil.rmtree(embeddings)
    print(f'Completed {job["id"]}: {len(windows)} windows', flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-root', type=Path, required=True, help='Directory containing feature_inputs/ and pssm_resources/')
    parser.add_argument('--experiment', nargs='+', default=['all'], help='Names shown by --list, or all')
    parser.add_argument('--job', help='One store ID shown by --list; useful for a small first run')
    parser.add_argument('--list', action='store_true', help='List experiments and store IDs without running encoders')
    parser.add_argument('--check-inputs', action='store_true', help='Verify inputs and stop')
    parser.add_argument('--explicit-only', action='store_true', help='Build only the five explicit arrays; no upstream weights required')
    parser.add_argument('--checkpoint', type=Path, default=os.environ.get('MOTIFAE_CHECKPOINT'), help='MotifAE k3 checkpoint, or set MOTIFAE_CHECKPOINT')
    parser.add_argument('--device', choices=['cpu', 'mps', 'cuda'], default='cpu')
    parser.add_argument('--batch-size', type=int, default=1, help='ESM2 crop batch size')
    parser.add_argument('--keep-esm', action='store_true', help='Keep intermediate ESM2 embeddings after sparse encoding')
    args = parser.parse_args()
    root = args.data_root.resolve()
    manifest = json.loads((root / 'feature_inputs/manifest.json').read_text())
    all_jobs = manifest['jobs']
    available = {job['experiment'] for job in all_jobs}
    unknown = set(args.experiment) - available - {'all'}
    if unknown:
        parser.error('Unknown experiments: ' + ', '.join(sorted(unknown)))
    jobs = [j for j in all_jobs if ('all' in args.experiment or j['experiment'] in args.experiment)
            and (args.job is None or j['id'] == args.job)]
    if not jobs:
        parser.error('No matching feature stores')
    if args.list:
        for j in jobs:
            print(f'{j["id"]}\t{j["windows"]} windows\t{j["output"]}')
        return
    verify_inputs(root, jobs)
    if args.check_inputs:
        return
    if not args.explicit_only:
        if not args.checkpoint or not args.checkpoint.is_file():
            parser.error('Set --checkpoint or MOTIFAE_CHECKPOINT to MotifAE_k3_step_80000.pt')
        resources = json.loads((Path(__file__).parent / 'configs/resources.json').read_text())
        expected = resources['external_resources']['motifae']['files'][0]['sha256']
        if digest(args.checkpoint) != expected:
            parser.error('MotifAE checkpoint checksum differs from the released k3 model')
    hashes = code_signature()
    for index, job in enumerate(jobs, 1):
        print(f'[{index}/{len(jobs)}] {job["id"]}', flush=True)
        build_job(root, job, args, hashes)


if __name__ == '__main__':
    main()
