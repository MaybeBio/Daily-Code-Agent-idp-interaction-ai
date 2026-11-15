"""Restore losslessly packed MotifAE activations from the feature archive using NumPy."""
import argparse
from contextlib import ExitStack
import hashlib
import json
from pathlib import Path
import tempfile

import numpy as np

FORMAT = 'motifae_residue_dictionary_v1'
OUTPUT_NAMES = ('motifae_feature_indices.int32.bin', 'motifae_feature_values.float32.bin')


def file_sha256(path):
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024**2), b''):
            digest.update(block)
    return digest.hexdigest()


def local_path(root, relative):
    path = (root / relative).resolve()
    if Path(relative).is_absolute() or not path.is_relative_to(root.resolve()):
        raise ValueError(f'Path is outside the data directory: {relative}')
    return path


def checked_file(root, record):
    path = local_path(root, record['path'])
    if path.stat().st_size != record['bytes'] or file_sha256(path) != record['sha256']:
        raise ValueError(f'Checksum mismatch: {path}')
    return path


def restore_store(root, record, verify_only=False, block_rows=4096):
    root = Path(root).resolve()
    if record['format'] != FORMAT or block_rows < 1:
        raise ValueError('Unsupported feature format or invalid block size')
    archive = checked_file(root, record)
    ptr_path = checked_file(root, record['residue_indptr'])
    outputs = record['restores']
    paths = [local_path(root, r['path']) for r in outputs]
    if tuple(p.name for p in paths) != OUTPUT_NAMES or any(p.parent != archive.parent for p in paths):
        raise ValueError('Unexpected sparse output paths')
    pending = []
    for index, (path, output) in enumerate(zip(paths, outputs)):
        if path.exists():
            checked_file(root, output)
        else:
            pending.append(index)
    if not pending:
        return 'verified'
    if verify_only:
        raise FileNotFoundError(f'Sparse activations have not been restored: {archive.parent}')

    with np.load(archive, allow_pickle=False) as data:
        ids = data['feature_indices']
        values = data['feature_values']
        offsets = data['unique_indptr']
        references = data['residue_rows']
    if ids.dtype != np.dtype('<i4') or values.dtype != np.dtype('<f4'):
        raise ValueError('Expected int32 feature IDs and float32 activations')
    if offsets.dtype != np.dtype('<i8') or references.dtype != np.dtype('<i4'):
        raise ValueError('Invalid sparse dictionary index types')
    if any(a.ndim != 1 for a in [ids, values, offsets, references]):
        raise ValueError('Sparse dictionary arrays must be one-dimensional')
    if len(ids) != len(values) or len(offsets) < 1 or offsets[0] != 0 or offsets[-1] != len(ids):
        raise ValueError('Invalid sparse dictionary lengths')
    if np.any(offsets[1:] < offsets[:-1]):
        raise ValueError('Sparse dictionary offsets are not monotonic')
    if len(references) and (references.min() < 0 or references.max() >= len(offsets) - 1):
        raise ValueError('Residue reference is outside the dictionary')
    ptr = np.load(ptr_path, mmap_mode='r', allow_pickle=False)
    if ptr.dtype != np.dtype('<i8') or ptr.shape != (len(references) + 1,) or ptr[0] != 0:
        raise ValueError('Invalid original residue offsets')
    if any(r['bytes'] != int(ptr[-1]) * 4 for r in outputs):
        raise ValueError('Sparse output sizes do not match residue offsets')

    temporary = {}
    try:
        with ExitStack() as stack:
            handles, digests = {}, {}
            for index in pending:
                handle = stack.enter_context(tempfile.NamedTemporaryFile(
                    prefix='.' + paths[index].name + '.', suffix='.part',
                    dir=archive.parent, delete=False))
                temporary[index] = Path(handle.name)
                handles[index] = handle
                digests[index] = hashlib.sha256()
            arrays = (ids, values)
            for start in range(0, len(references), block_rows):
                stop = min(start + block_rows, len(references))
                rows = references[start:stop]
                counts = offsets[rows + 1] - offsets[rows]
                if not np.array_equal(counts, ptr[start + 1:stop + 1] - ptr[start:stop]):
                    raise ValueError('Dictionary row sizes disagree with original residue offsets')
                for index in pending:
                    source = arrays[index]
                    block = np.concatenate([source[offsets[j]:offsets[j + 1]] for j in rows])
                    content = block.tobytes()
                    handles[index].write(content)
                    digests[index].update(content)
            for index in pending:
                if digests[index].hexdigest() != outputs[index]['sha256']:
                    raise ValueError(f'Restored checksum mismatch: {paths[index]}')
        for index in pending:
            if temporary[index].stat().st_size != outputs[index]['bytes']:
                raise ValueError(f'Restored file size mismatch: {paths[index]}')
        for index in pending:
            temporary[index].replace(paths[index])
    finally:
        for path in temporary.values():
            path.unlink(missing_ok=True)
    return 'restored'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-root', type=Path, required=True)
    parser.add_argument('--store', help='Restore only this store directory, relative to --data-root')
    parser.add_argument('--verify-only', action='store_true', help='Verify existing restored files without writing')
    args = parser.parse_args()
    root = args.data_root.expanduser().resolve()
    try:
        manifest = json.loads((root / 'features_manifest.json').read_text())
        entries = [r for r in manifest['files'] if r.get('format') == FORMAT]
        if args.store:
            entries = [r for r in entries if Path(r['path']).parent == Path(args.store)]
        if not entries:
            raise ValueError('No matching packed feature stores in features_manifest.json')
        restored = 0
        for number, record in enumerate(entries, 1):
            restored += restore_store(root, record, args.verify_only) == 'restored'
            if number % 25 == 0 or number == len(entries):
                print(f'Checked {number}/{len(entries)} stores; restored {restored}.', flush=True)
        print('All sparse activations match the reference SHA-256 checksums.')
    except (OSError, ValueError, KeyError) as error:
        parser.exit(1, str(error) + '\n')


if __name__ == '__main__':
    main()
