#!/usr/bin/env python3
"""Accumulate saved soft-count contributions on absolute slab-time blocks."""
from pathlib import Path
import hashlib
import json
import re
import numpy as np
import analyze as a


def fingerprint(path):
    size = path.stat().st_size
    digest = hashlib.sha256()
    digest.update(str(size).encode())
    with path.open('rb') as handle:
        for offset in (0, max(0, size // 2 - 32768), max(0, size - 65536)):
            handle.seek(offset)
            digest.update(handle.read(65536))
    return digest.hexdigest()


def accumulate(q, start, origin, counts, masses, lag=200, unit=1000):
    for absolute_block in range(start // unit, (start + len(q) - 1) // unit + 1):
        lo = max(0, absolute_block * unit - start)
        hi = min(len(q), (absolute_block + 1) * unit - start)
        bi = absolute_block - origin // unit
        masses[bi] += q[lo:hi].sum(axis=0)
        stop = min(hi, len(q) - lag)
        if stop > lo:
            counts[bi] += q[lo:stop].T @ q[lo + lag:stop + lag]


def main():
    provenance, exclusions = [], []
    for system in a.s.SYSTEMS:
        for role in a.s.ROLES:
            directory = a.s.model_dir(system, 'condensate', role)
            local = directory / 'extracted_500ns_bulk_chains'
            source = a.s.density_paths(system, role)[0].parent / 'extracted_500ns_bulk_chains'
            meta_path = local / 'segments_meta.json'
            if not meta_path.exists():
                meta_path = source / 'segments_meta.json'
            metadata = json.loads(meta_path.read_text())
            by_chain = {int(v['chain_id']): v for v in metadata}
            assert len(by_chain) == len(metadata)
            files = sorted(local.glob('chain*_bulk_500ns.dcd'))
            inter = directory / a.s.intermediate_name(system)
            with np.load(inter / 'combined_meta.npz') as z:
                lengths = z['lengths']
            assert len(files) == len(lengths)
            windows = []
            for ti, (file, length) in enumerate(zip(files, lengths)):
                chain = int(re.search(r'chain(\d+)', file.name)[1])
                prefix = dict(**a.condition(system, 'condensate', role), trajectory=ti,
                              chain_id=chain, qcache_frames=int(length), dcd=str(file.relative_to(a.ROOT)),
                              metadata_path=str(meta_path.relative_to(a.ROOT)))
                if chain not in by_chain:
                    exclusions.append(dict(**prefix, reason='chain absent from segments_meta.json'))
                    continue
                entry = by_chain[chain]
                assert int(entry['end']) - int(entry['start']) + 1 == length, prefix
                sampled_hash = fingerprint(file)
                check = 'local metadata'
                if meta_path.parent != local:
                    assert sampled_hash == fingerprint(source / file.name), prefix
                    check = 'source metadata; same size and 3x64KiB sampled DCD SHA256'
                windows.append((ti, int(entry['start']), int(length)))
                provenance.append(dict(**prefix, start_frame=int(entry['start']), end_frame=int(entry['end']),
                                       sampled_dcd_sha256=sampled_hash, provenance_check=check,
                                       metadata_sha256=hashlib.sha256(meta_path.read_bytes()).hexdigest()))
            origin = min(w[1] for w in windows) // 10000 * 10000
            end = (max(w[1] + w[2] for w in windows) + 9999) // 10000 * 10000
            tag = f'{system}_condensate_{a.s.TEMPS[system][role]}'
            with np.load(a.CACHE / f'{tag}_n5_lag50.npz') as z:
                attempts = z['attempts']
            counts = np.zeros((len(attempts), (end - origin) // a.UNIT, 5, 5))
            masses = np.zeros((len(attempts), (end - origin) // a.UNIT, 5))
            for ai, att in enumerate(attempts):
                perm = np.argsort(np.load(directory / f'state_map_att{att}_raw_to_final.npy'))
                for ti, start, length in windows:
                    file = inter / f'qcache_combined_n5_att{att}_lag50/traj{ti}.npy'
                    raw = np.load(file, mmap_mode='r')
                    assert raw.shape == (length, 5)
                    q = np.asarray(raw[:, perm], dtype=np.float64)
                    accumulate(q, start, origin, counts[ai], masses[ai])
            np.savez_compressed(a.OUT / f'{tag}_absolute_statistics.npz', counts=counts, masses=masses,
                                attempts=attempts, absolute_origin=origin, absolute_end=end,
                                unit_frames=a.UNIT, excluded_trajectories=np.asarray([
                                    ti for ti in range(len(files)) if ti not in {w[0] for w in windows}], dtype=int))
            print(tag, 'absolute blocks', counts.shape[1], 'windows', len(windows), flush=True)
            a.write_csv('absolute_time_provenance.csv', provenance)
            a.write_csv('absolute_time_exclusions.csv', exclusions)
    (a.OUT / 'absolute_statistics_metadata.json').write_text(json.dumps(dict(
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), unit_frames=a.UNIT,
        dt_ns=.01, block_rule='absolute source-frame index; lagged pairs grouped by start frame',
        limitation='conditional on selected windows and recorded provenance; sampled DCD hashes are not full-file hashes'), indent=2) + '\n')


if __name__ == '__main__':
    main()
