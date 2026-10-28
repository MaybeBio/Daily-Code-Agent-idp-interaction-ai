#!/usr/bin/env python3
"""Versioned middle-60% C5 check; preserve the earlier archived-rule analysis."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil

os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np
from MDAnalysis.coordinates.DCD import DCDReader
import extract

BASE = Path(__file__).resolve().parent
OUT = BASE / 'central60_2026_10_10'


def prepare():
    (OUT / 'raw').mkdir(parents=True, exist_ok=True)
    original = json.loads((BASE / 'plan.json').read_text())
    revised = []
    for item in original:
        item = dict(item)
        stem = f"{item['case']}_block{item['block']}"
        src = BASE / 'raw' / (stem + '.npz')
        dst = OUT / 'raw' / src.name
        meta = json.loads(src.with_suffix('.json').read_text())
        if item['case'] != 'FUS_325':
            if not dst.exists():
                os.link(src, dst)
            shutil.copy2(src.with_suffix('.json'), dst.with_suffix('.json'))
        else:
            item['selection_params'] = [5.0, 20, 80]
            item['selection_revision'] = 'central60_2026_10_10'
            ident = hashlib.sha256(json.dumps(item, sort_keys=True).encode()).hexdigest()
            if not dst.exists():
                flags, frames = [], []
                for part in item['parts']:
                    path = extract.ROOT / part['path']
                    assert path.stat().st_size == part['bytes']
                    index = int(path.stem.split('part')[-1])
                    with DCDReader(str(path)) as reader:
                        for j in range(0, 10000, 20):
                            ts = reader[j]
                            com = extract.selector.get_com_nm(ts.positions, 200)
                            box = ts.dimensions[:3].astype(float) * .1
                            flags.append(extract.bulk_flags(com, box, item['selection_params']))
                            frames.append((index-1)*10000+j)
                    print('Reselected', stem, path.name, flush=True)
                with np.load(src) as old:
                    arrays = {key: old[key] for key in old.files}
                np.testing.assert_array_equal(arrays['frames'], frames)
                arrays['bulk'] = np.asarray(flags)
                arrays['identity'] = np.asarray(ident)
                np.savez_compressed(dst, **arrays)
            else:
                with np.load(dst) as cached:
                    assert str(cached['identity']) == ident
            meta.update(identity=ident, selection_params=item['selection_params'],
                        source_cache_sha256=extract.digest(src),
                        change='Only spatial eligibility recomputed; contacts and mass COM unchanged.')
            extract.dump(dst.with_suffix('.json'), meta)
        revised.append(item)
    extract.dump(OUT / 'plan.json', revised)
    for name in ('training_coordinate_provenance.csv', 'training_provenance.json'):
        shutil.copy2(BASE / name, OUT / name)
    extract.dump(OUT / 'selection_source_parameters.json', {
        'analysis_spatial_parameters': {
            key: ([5.0, 20, 80] if key == 'FUS_325' else list(value))
            for key, value in extract.CONFIG.items()},
        'historical_FUS_protocol': 'Author reports collaborator used central 60%, 2026-10-10.',
        'scope': 'Author clarification is not an independent run-log verification. '
                 'The new C5 calculation explicitly uses central 60% for every condition. '
                 'FUS high-temperature cluster cutoff remains 5 nm in this sensitivity analysis.',
        'previous_audit_retained_at': str(BASE / 'selection_source_parameters.json')})
    ids = {}
    for path in sorted((OUT / 'raw').glob('*.npz')):
        with np.load(path) as z:
            bulk = z['bulk']
            ids[path.stem] = {'500ns': np.flatnonzero(bulk.all(axis=0)).tolist(),
                             '100ns': [np.flatnonzero(a.all(axis=0)).tolist()
                                       for a in bulk.reshape(5, 500, 200)]}
    extract.dump(OUT / 'new_residence_chain_ids.json', ids)


def summarize():
    def read(path):
        with path.open() as f:
            return list(csv.DictReader(f))
    summary = read(BASE / 'selection_summary.csv')
    mobility = read(OUT / 'mobility_summary.csv')
    counts = read(OUT / 'residence_counts.csv')
    occupancy = read(OUT / 'vamp_sampled_residence.csv')
    for r in summary:
        case = r['case']
        spatial = next(m for m in mobility if m['case'] == case and m['lag_ns'] == '10'
                       and m['coordinates'] == 'xy' and m['selected'] == 'original_ids_bulk_origins')
        r['bulk_matched_mobility_ratio'] = spatial['ratio']
        c = [x for x in counts if x['case'] == case]
        o = [x for x in occupancy if x['case'] == case]
        for duration, label in [('500', 'long'), ('100', 'short')]:
            r[f'min_{label}_chains'] = min(int(x[f'n_{duration}']) for x in c)
            r[f'max_{label}_chains'] = max(int(x[f'n_{duration}']) for x in c)
        for rule in ('source', 'uniform'):
            r[f'sampled_VAMP_{rule}_rule_bulk_fraction'] = (
                sum(int(x['observed_frames'])*float(x[f'{rule}_rule_bulk_fraction']) for x in o)
                / sum(int(x['observed_frames']) for x in o))
    with (OUT / 'selection_summary.csv').open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(summary[0]))
        writer.writeheader()
        writer.writerows(summary)


def run():
    import analyze
    import plot
    import verify
    analyze.HERE = plot.HERE = verify.HERE = OUT
    analyze.CONFIG = verify.CONFIG = {
        key: ((5.0, 20, 80) if key == 'FUS_325' else value)
        for key, value in extract.CONFIG.items()}
    analyze.main()
    summarize()
    plot.main()
    verify.main()
    # Assert all observations except FUS spatial eligibility are unchanged.
    for path in (OUT / 'raw').glob('FUS_325*.npz'):
        with np.load(path) as new, np.load(BASE / 'raw' / path.name) as old:
            for key in old.files:
                if key not in ('bulk', 'identity'):
                    np.testing.assert_array_equal(new[key], old[key])
            assert np.all(~new['bulk'] | old['bulk'])
    print('PASS: FUS central-60 flags are a subset; all cached observables unchanged.', flush=True)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--prepare-only', action='store_true')
    args = ap.parse_args()
    prepare()
    if not args.prepare_only:
        run()
