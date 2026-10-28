#!/usr/bin/env python3
"""Rebuild selected-window intermolecular observables from raw kinetic DCDs."""
import argparse
from concurrent.futures import ProcessPoolExecutor
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import sys
import time

os.environ.setdefault('MPLCONFIGDIR', '/tmp/mpl-interchain-recompute')
import numpy as np
from MDAnalysis.coordinates.DCD import DCDReader
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PILOT_PATH = HERE.parent / 'selection_density/corrected_time_pilot/run_corrected_contacts.py'
spec = importlib.util.spec_from_file_location('corrected_pilot', PILOT_PATH)
pilot = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pilot)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def dump(path, obj):
    tmp = path.with_suffix(path.suffix + '.tmp')
    tmp.write_text(json.dumps(obj, indent=2, allow_nan=False) + '\n')
    tmp.replace(path)


def partner_contact_counts(xyz, box, chain_len, chains, cutoff=.66):
    """Count bead pairs, retaining their partner-chain identity."""
    if not np.isfinite(xyz).all() or not np.isfinite(box).all() or np.any(box <= 0):
        raise ValueError('Nonfinite coordinates or invalid periodic box')
    pos = np.asarray(xyz, dtype=float) % box
    tree = cKDTree(pos, boxsize=box)
    result = {}
    for chain in chains:
        neighbors = tree.query_ball_point(pos[chain*chain_len:(chain+1)*chain_len], cutoff)
        ids = np.concatenate(neighbors).astype(np.int64) // chain_len
        counts = np.bincount(ids, minlength=200).astype(np.int32)
        counts[chain] = 0
        result[chain] = counts
    return result


def process_case(case):
    legacy, legacy_path = pilot.load_legacy()
    condition = next(c for c in legacy.CONDITIONS if f'{c.system}_{c.temp[:-1]}' == case)
    out = HERE / 'raw' / case
    out.mkdir(parents=True, exist_ok=True)
    parts = pilot.inspect_parts(condition)
    segments = json.loads(condition.metadata.read_text())
    if len(set(s['id'] for s in segments)) != len(segments):
        raise ValueError('Duplicate segment IDs')
    schedule, windows = pilot.build_schedule(segments, parts, 'full', 70)
    # Cache reuse requires identical source versions and analysis configuration.
    identity = dict(case=case, script_sha256=sha(__file__), pilot_sha256=sha(PILOT_PATH),
                    operator_sha256=sha(legacy_path), metadata_sha256=sha(condition.metadata),
                    topology_sha256=sha(condition.topology), parts=parts,
                    stride_raw=50, raw_dt_ns=.01, cutoff_nm=.66, last_parts=70,
                    chain_id_base=0, schedule_frames=len(schedule))
    identity_hash = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()
    dump(out / 'inputs.json', identity)
    dump(out / 'windows.json', windows)
    by_part = {}
    for frame in sorted(schedule):
        by_part.setdefault(frame//10000 + 1, []).append(frame)
    segment_index = {s['id']: i for i, s in enumerate(segments)}
    dump(out / 'segments.json', segments)
    chain_len = parts[0]['atoms']//200
    atom_chain = np.repeat(np.arange(200), chain_len)
    start = time.monotonic()
    filled = 0
    checks = []
    for part, frames in by_part.items():
        dest = out / f'part{part:03d}.npz'
        if dest.exists():
            with np.load(dest, allow_pickle=False) as z:
                if str(z['identity_hash']) != identity_hash:
                    raise ValueError(f'Incompatible checkpoint: {dest}')
                filled += len(z['frame'])
            continue
        rows, counts = [], []
        info = parts[part-1]
        with DCDReader(str(ROOT / info['path'])) as reader:
            for j, frame in enumerate(frames):
                local = frame % 10000
                ts = reader[local]
                if ts.dimensions is None or not np.allclose(ts.dimensions[3:], 90):
                    raise ValueError(f'Invalid box at {case}:{frame}')
                header = float(ts.time/1000)
                expected = info['first_header_time_ns'] + local*info['raw_dt_ns']
                if not np.isclose(header, expected, atol=.001, rtol=0):
                    raise ValueError(f'Unexpected header time at {case}:{frame}')
                targets = schedule[frame]
                chains = sorted({s['chain_id'] for s, _ in targets})
                xyz, box = ts.positions.astype(float)*.1, ts.dimensions[:3].astype(float)*.1
                result = partner_contact_counts(xyz, box, chain_len, chains)
                # Compare the retained partner-resolved counts with the original operator.
                if j == 0:
                    reference = legacy.count_contacts_for_chains(xyz, box, chain_len, chains, .66, atom_chain)
                    for c in chains:
                        if int(result[c].sum()) != reference[c]:
                            raise AssertionError(f'Contact operator mismatch {case}:{frame}:{c}')
                    checks.append(dict(part=part, frame=frame, chains=len(chains), exact_match=True))
                for s, _ in targets:
                    rows.append((frame, segment_index[s['id']], s['chain_id'], header))
                    counts.append(result[s['chain_id']])
        arr = np.asarray(rows)
        tmp = dest.with_suffix('.tmp.npz')
        np.savez_compressed(tmp, identity_hash=np.asarray(identity_hash),
                            frame=arr[:, 0].astype(np.int64), segment=arr[:, 1].astype(np.int32),
                            chain=arr[:, 2].astype(np.int32), time_ns=arr[:, 3],
                            partner_contacts=np.asarray(counts, dtype=np.int32))
        tmp.replace(dest)
        filled += len(rows)
        dump(out/'progress.json', dict(case=case, last_completed_part=part,
             completed_chain_frames=filled, elapsed_seconds=time.monotonic()-start,
             planned_source_frames=len(schedule), status='running'))
        print(f'{case}: part {part}/{max(by_part)}; {filled} chain-frames; '
              f'{time.monotonic()-start:.1f}s', flush=True)
    all_arrays = {key: [] for key in ['frame', 'segment', 'chain', 'time_ns', 'partner_contacts']}
    for part in by_part:
        with np.load(out/f'part{part:03d}.npz', allow_pickle=False) as z:
            for key in all_arrays:
                all_arrays[key].append(z[key])
    arrays = {k: np.concatenate(v) for k, v in all_arrays.items()}
    if len(arrays['frame']) != sum(len(v) for v in schedule.values()):
        raise AssertionError('Incomplete raw series')
    for i, s in enumerate(segments):
        got = arrays['frame'][arrays['segment'] == i]
        expected = [f for f in sorted(schedule) if any(t['id']==s['id'] for t,_ in schedule[f])]
        if not np.array_equal(got, expected):
            raise AssertionError(f'Incomplete or reordered segment {s["id"]}')
    np.savez_compressed(out/'timeseries.npz', **arrays, identity_hash=np.asarray(identity_hash),
                        chain_len=np.asarray(chain_len), dt_ns=np.asarray(.5))
    summary = dict(case=case, status='complete', identity_hash=identity_hash,
                   source_frames=len(schedule), chain_frames=len(arrays['frame']),
                   sampled_segments=len(np.unique(arrays['segment'])),
                   selected_part_start=min(by_part), selected_part_end=max(by_part),
                   time_start_ns=float(arrays['time_ns'].min()), time_end_ns=float(arrays['time_ns'].max()),
                   header_jump_parts=[p['part'] for p in parts if p['header_boundary_discontinuous']],
                   contact_operator_checks=checks, failed_rows=0,
                   elapsed_seconds=time.monotonic()-start,
                   limitation='Historical selected windows, not new bulk selection or independent replicas.')
    dump(out/'completed.json', summary)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cases', default='FUS_280,FUS_325,A1_230,A1_260,TDP43_260,TDP43_300')
    parser.add_argument('--workers', type=int, default=1)
    args = parser.parse_args()
    cases = args.cases.split(',')
    if len(cases) != len(set(cases)):
        parser.error('Duplicate case')
    if args.workers == 1:
        results = [process_case(c) for c in cases]
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            results = list(pool.map(process_case, cases))
    print(json.dumps(results, indent=2), flush=True)


if __name__ == '__main__':
    main()
