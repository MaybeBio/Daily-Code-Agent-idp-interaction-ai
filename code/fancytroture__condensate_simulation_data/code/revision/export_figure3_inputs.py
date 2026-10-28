#!/usr/bin/env python3
"""Export compact Figure 3 inputs from the archived analysis workspace."""
import argparse
import csv
import hashlib
import importlib.util
import json
import os
from pathlib import Path

os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('OMP_NUM_THREADS', '1')
os.environ.setdefault('MPLCONFIGDIR', '/tmp/jcp_revision_export_mpl')
import numpy as np

HERE = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workspace', type=Path, required=True)
    parser.add_argument('--output', type=Path, default=HERE / 'figure3_inputs')
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    source = HERE / 'project_figure/scripts/make_figure2_dynamics_story.py'
    spec = importlib.util.spec_from_file_location('figure3_original', source)
    style = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(style)
    style.RG_DIRS = {s: args.workspace / path.relative_to(HERE) for s, path in style.RG_DIRS.items()}
    style.INTRA_DIR = args.workspace / style.INTRA_DIR.relative_to(HERE)
    arrays, response, hashes = {}, {}, {}

    def record(path):
        hashes[str(path.relative_to(args.workspace))] = hashlib.sha256(path.read_bytes()).hexdigest()

    for system in style.SYSTEM_ORDER:
        for family in style.FAMILY_ORDER:
            for role in style.ROLE_ORDER:
                prefix = f'{system}__{family}__{role}__'
                rg = style.rg_cache(system, family, role)
                keep = rg['time_axis'] <= 60
                for key in ('time_axis', 'acf_rg'):
                    arrays[prefix + key] = rg[key][keep]
                intra = style.load_intra_cache(system, family, role)
                for key in ('time_ns', 'acf_mean', 'acf_sem'):
                    arrays[prefix + key] = intra[key]
                temp = style.TEMPS[system][role]
                record(style.RG_DIRS[system] / f'cache_{family}_{temp}.npz')
                record(style.INTRA_DIR / f'intrachain_{system}_{family}_{temp}.npz')
            response[f'{system}__{family}'] = style.rg_heating_response_interval(system, family)
            print(f'Exported {system} {family} Rg response', flush=True)
        for role in style.ROLE_ORDER:
            temp = style.TEMPS[system][role].removesuffix('K')
            for metric in ('inter', 'partner'):
                path = args.workspace / 'revision_2026_09_30/interchain_recompute/analysis' / f'{system}_{temp}_{metric}_curves.npz'
                record(path)
                with np.load(path, allow_pickle=False) as data:
                    assert data['block_250_boot'].shape[0] == 5000
                    np.testing.assert_allclose(np.diff(data['time_ns']), .5, atol=1e-12)
                    np.testing.assert_allclose(data['mean'], data['unit_curves'].mean(axis=0), atol=1e-12)
                    prefix = f'{system}__{role}__{metric}__'
                    arrays[prefix + 'time_ns'] = data['time_ns']
                    arrays[prefix + 'mean'] = data['mean']
                    arrays[prefix + 'interval'] = np.quantile(data['block_250_boot'], [.025, .975], axis=0, method='inverted_cdf')
    np.savez_compressed(args.output / 'curves.npz', **arrays)
    (args.output / 'rg_heating_response_percent.json').write_text(json.dumps(response, indent=2) + '\n')
    path = args.workspace / 'storyline_figures_for_ppt_corrected_v3/intra_inter_contact_relaxation/intrachain_contact_relaxation_ratio.csv'
    record(path)
    (args.output / 'intrachain_contact_relaxation_ratio.csv').write_bytes(path.read_bytes())
    (args.output / 'provenance.json').write_text(json.dumps({'input_sha256': hashes,
        'rg_bootstrap': 5000, 'intermolecular_bootstrap': 5000, 'origin_block_ns': 250,
        'intermolecular_dt_ns': .5, 'intermolecular_interval_quantile_method': 'inverted_cdf'}, indent=2) + '\n')


if __name__ == '__main__':
    main()
