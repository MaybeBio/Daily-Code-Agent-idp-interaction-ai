#!/usr/bin/env python3
"""Correct raw-frame contact reader. Defaults to six diagnostic frames per case.

No ACF or full trajectory export is exposed by this pilot entry point.
"""
import argparse
import csv
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import sys

import numpy as np
from MDAnalysis.coordinates.DCD import DCDReader

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
RAW_DT_NS = 0.01
RAW_STRIDE = 50
FRAMES_PER_PART = 10000
FIELDS = ['case', 'mode', 'segment_id', 'chain_id', 'source_part', 'source_path',
          'raw_local_frame', 'raw_global_frame', 'nominal_index_time_ns', 'physical_time_ns', 'header_time_ns',
          'source_window_start_raw', 'source_window_end_raw', 'within_source_window',
          'pilot_role', 'delta_from_previous_row_ns', 'continuity_run',
          'N_inter', 'status', 'error']


def load_legacy():
    # This module defines the original operator but has a guarded main().
    os.environ['MPLCONFIGDIR'] = str(HERE / '.mplconfig')
    path = ROOT / 'ppt_5_11_figures/compute_interchain_contact_relaxation_parts.py'
    spec = importlib.util.spec_from_file_location('corrected_contact_legacy', path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module, path


def global_frame(part, local):
    return (part - 1) * FRAMES_PER_PART + local


def time_ns(frame):
    # OpenMM first report is after the first report interval, not at t=0.
    return (frame + 1) * RAW_DT_NS


def validate_parts(parts):
    ids = [p['part'] for p in parts]
    if not ids or ids != list(range(1, ids[-1] + 1)):
        raise ValueError(f'Missing/noncontiguous/duplicate part IDs: {ids}')
    atoms = parts[0]['atoms']
    for i, p in enumerate(parts):
        if p['atoms'] != atoms or atoms % 200:
            raise ValueError(f'Inconsistent atom count: {p}')
        if not np.isclose(p['raw_dt_ns'], RAW_DT_NS, rtol=1e-6, atol=1e-9):
            raise ValueError(f'Wrong raw dt: {p}')
        if not (0 < p['frames'] <= FRAMES_PER_PART):
            raise ValueError(f'Invalid part length: {p}')
        if i < len(parts) - 1 and p['frames'] != FRAMES_PER_PART:
            raise ValueError(f'Short internal part would create a time gap: {p}')
        span=p['last_header_time_ns']-p['first_header_time_ns']
        if not np.isclose(span,(p['frames']-1)*p['raw_dt_ns'],rtol=0,atol=.001):
            raise ValueError(f'Inconsistent within-part header time span: {p}')
        p['header_offset_from_index_ns']=p['first_header_time_ns']-time_ns(global_frame(p['part'],0))
        delta=None if i==0 else p['first_header_time_ns']-parts[i-1]['last_header_time_ns']
        p['header_boundary_delta_ns']=delta
        p['header_boundary_discontinuous']=bool(delta is not None and not np.isclose(delta,RAW_DT_NS,rtol=0,atol=.001))


def inspect_parts(condition):
    parts = []
    paths = sorted(condition.directory.glob(condition.part_glob),
                   key=lambda p: int(re.search(r'part(\d+)', p.name)[1]))
    for path in paths:
        part = int(re.search(r'part(\d+)', path.name)[1])
        with DCDReader(str(path)) as reader:
            first_time = float(reader[0].time / 1000)
            last_time = float(reader[-1].time / 1000)
            parts.append(dict(part=part, path=str(path.relative_to(ROOT)),
                 frames=len(reader), atoms=reader.n_atoms, bytes=path.stat().st_size,
                 mtime_ns=path.stat().st_mtime_ns, raw_dt_ns=float(reader.dt/1000),
                 first_header_time_ns=first_time, last_header_time_ns=last_time))
    validate_parts(parts)
    return parts


def sampled_frames(start, end):
    first = start + (-start) % RAW_STRIDE
    return range(first, end+1, RAW_STRIDE)


def build_schedule(segments, parts, mode, last_parts):
    first_part = max(1, parts[-1]['part'] - last_parts + 1) if last_parts else 1
    lo = global_frame(first_part, 0)
    hi = global_frame(parts[-1]['part'], parts[-1]['frames']-1)
    schedule, stage = {}, []
    for s in segments:
        if not (0 <= s['chain_id'] < 200 and 0 <= s['start'] <= s['end']):
            raise ValueError(f'Invalid zero-based metadata: {s}')
        a, b = max(lo, s['start']), min(hi, s['end'])
        indices = sampled_frames(a, b)
        stage.append(dict(segment_id=s['id'], chain_id=s['chain_id'],
            original_start_raw=s['start'], original_end_raw=s['end'],
            clipped_start_raw=a, clipped_end_raw=b, n_planned_samples=len(indices),
            clipped=bool(a != s['start'] or b != s['end'])))
        if mode == 'full':
            for frame in indices:
                schedule.setdefault(frame, []).append((s, 'metadata_window'))
    if mode == 'pilot':
        # Exercise two within-part steps and a boundary, inside one metadata window.
        choice = None
        for s in segments:
            a, b = max(lo, s['start']), min(hi, s['end'])
            boundary = ((a + 2*RAW_STRIDE + FRAMES_PER_PART-1)//FRAMES_PER_PART)*FRAMES_PER_PART
            if boundary + RAW_STRIDE <= b:
                choice = s, boundary
                break
        if choice is None:
            raise ValueError('No source window spans a complete diagnostic part boundary')
        s, boundary = choice
        for frame in [boundary-100, boundary-50, boundary, boundary+50]:
            schedule[frame] = [(s, 'boundary_window')]
        tail = hi-hi%RAW_STRIDE
        for frame in [tail-RAW_STRIDE, tail]:
            # Tail may not be inside this chain's historical selected window.
            schedule.setdefault(frame, []).append((s, 'tail_reader_only'))
    return schedule, stage


def execute_frames(case, mode, parts, schedule, operator, writer):
    index = {p['part']: p for p in parts}
    chain_len = parts[0]['atoms']//200
    atom_chain_index = np.repeat(np.arange(200), chain_len)
    previous, previous_time, run_ids, broken = {}, {}, {}, set()
    reader, current = None, None
    stats = dict(source_frames=len(schedule), contact_rows=0, failed_rows=0)
    try:
        for frame, targets in sorted(schedule.items()):
            part, local = frame//FRAMES_PER_PART+1, frame%FRAMES_PER_PART
            info = index[part]
            if local >= info['frames']:
                raise ValueError(f'Scheduled frame past actual EOF: {frame}')
            status, error, header_time, counts = 'ok', '', '', {}
            try:
                if current != part:
                    if reader is not None:
                        reader.close()
                        reader = None
                    reader = DCDReader(str(ROOT/info['path']))
                    current = part
                ts = reader[local]
                if ts.dimensions is None or not np.allclose(ts.dimensions[3:], 90):
                    raise ValueError('Requires orthorhombic periodic box')
                xyz = ts.positions.astype(np.float64)*.1
                if not np.isfinite(xyz).all():
                    raise ValueError('Nonfinite coordinates')
                header_time = float(ts.time/1000)
                expected=info['first_header_time_ns']+local*info['raw_dt_ns']
                if not np.isclose(header_time, expected, rtol=0, atol=.001):
                    raise ValueError('Sample header time mismatch')
                counts = operator(xyz, ts.dimensions[:3].astype(float)*.1,
                    chain_len, sorted({s['chain_id'] for s,_ in targets}), .66, atom_chain_index)
            except Exception as exc:
                status, error = 'error', f'{type(exc).__name__}: {exc}'
            for s, role in targets:
                key = (s['id'], s['chain_id'])
                prev = previous.get(key)
                physical_delta=header_time-previous_time[key] if header_time!='' and key in previous_time else ''
                gap = prev is not None and (frame-prev != RAW_STRIDE or physical_delta=='' or not np.isclose(physical_delta,.5,rtol=0,atol=.001))
                if key not in run_ids:
                    run_ids[key] = 0
                elif gap or key in broken:
                    run_ids[key] += 1
                writer.writerow(dict(case=case, mode=mode, segment_id=s['id'], chain_id=s['chain_id'],
                    source_part=part, source_path=info['path'], raw_local_frame=local,
                    raw_global_frame=frame, nominal_index_time_ns=f'{time_ns(frame):.8f}', physical_time_ns=header_time,
                    header_time_ns=header_time, source_window_start_raw=s['start'],
                    source_window_end_raw=s['end'], within_source_window=s['start']<=frame<=s['end'],
                    pilot_role=role, delta_from_previous_row_ns=physical_delta,
                    continuity_run=run_ids[key], N_inter=counts.get(s['chain_id'], ''),status=status,error=error))
                previous[key] = frame
                if header_time!='':
                    previous_time[key]=header_time
                else:
                    previous_time.pop(key,None)
                if status == 'error':
                    broken.add(key)
                    stats['failed_rows'] += 1
                else:
                    broken.discard(key)
                stats['contact_rows'] += 1
    finally:
        if reader is not None:
            reader.close()
    return stats


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode', choices=['pilot'], default='pilot')
    parser.add_argument('--cases', help='Comma-separated identifiers, e.g. FUS_280,A1_260; default all six')
    parser.add_argument('--last-parts', type=int, default=70,
                        help='Clip metadata to last N parts, like legacy scope; 0 means all parts')
    parser.add_argument('--output-name', help='New child directory name; default results_pilot or results_full')
    args = parser.parse_args()
    if args.last_parts < 0:
        parser.error('--last-parts must be nonnegative')
    name = args.output_name or f'results_{args.mode}'
    if name in ('.','..') or Path(name).name != name:
        parser.error('--output-name must be one child directory name')
    out = HERE/name
    if out.exists():
        parser.error(f'{out} exists; choose a new --output-name to preserve evidence')
    legacy, source = load_legacy()
    requested = set(args.cases.split(',')) if args.cases else None
    conditions = [c for c in legacy.CONDITIONS if requested is None or f'{c.system}_{c.temp[:-1]}' in requested]
    if not conditions or (requested and len(conditions)!=len(requested)):
        parser.error('Unknown case identifier')
    out.mkdir()
    report = dict(mode=args.mode, raw_dt_ns=RAW_DT_NS, raw_stride=RAW_STRIDE,
        sample_dt_ns=.5, metadata_frame_divisor=1, chain_id_base=0,
        nominal_first_frame_time_ns=.01, physical_time_source='DCD header; retained separately from numbered index time', last_parts=args.last_parts, n_chains=200, cutoff_nm=.66,
        original_contact_operator=str(source), operator_file_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        runner_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        limitation='Reader/contact-operator validation only. No ACF, no relaxation inference, no new MD.',
        gap_policy='Keep error rows and raw timestamps; never drop NaNs or compress gaps; no ACF implemented.',
        pilot_tail_policy='Tail rows are reader diagnostics and may be outside selected residence window.', cases=[])
    failures = False
    for c in conditions:
        case = f'{c.system}_{c.temp[:-1]}'
        print(f'[{case}] checking part continuity, headers and actual EOF', flush=True)
        record = dict(case=case, metadata=str(c.metadata), topology=str(c.topology))
        try:
            parts = inspect_parts(c)
            segments = json.loads(c.metadata.read_text())
            if len({s['id'] for s in segments}) != len(segments):
                raise ValueError('Duplicate metadata segment IDs')
            schedule, stage = build_schedule(segments, parts, args.mode, args.last_parts)
            (out/f'{case}_parts.json').write_text(json.dumps(parts, indent=2)+'\n')
            (out/f'{case}_windows.json').write_text(json.dumps(stage, indent=2)+'\n')
            with (out/f'{case}_contacts.csv').open('w',newline='') as f:
                writer = csv.DictWriter(f,fieldnames=FIELDS)
                writer.writeheader()
                stats = execute_frames(case,args.mode,parts,schedule,legacy.count_contacts_for_chains,writer)
            record.update(status='ok' if stats['failed_rows']==0 else 'sample_errors',
                 n_parts=len(parts),last_part_frames=parts[-1]['frames'],
                 header_discontinuity_parts=[p['part'] for p in parts if p['header_boundary_discontinuous']],
                 last_raw_global_frame=global_frame(parts[-1]['part'],parts[-1]['frames']-1),
                 metadata_sha256=hashlib.sha256(c.metadata.read_bytes()).hexdigest(), **stats)
            failures |= stats['failed_rows'] > 0
        except Exception as exc:
            record.update(status='failed',error=f'{type(exc).__name__}: {exc}')
            failures = True
        report['cases'].append(record)
        (out/'run_manifest.json').write_text(json.dumps(report,indent=2)+'\n')
        print(f'[{case}] {record}',flush=True)
    if failures:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
