#!/usr/bin/env python3
"""Render revised Figure 3 using compact, validated analysis summaries."""
import argparse
import csv
import importlib.util
import json
import os
from pathlib import Path

os.environ.setdefault('MPLCONFIGDIR', '/tmp/jcp_revision_figure3_mpl')
import numpy as np

HERE = Path(__file__).resolve().parent


def read(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, default=HERE.parent / 'outputs/revision_20261018')
    args = parser.parse_args()
    path = HERE / 'project_figure/scripts/make_figure2_dynamics_story.py'
    spec = importlib.util.spec_from_file_location('figure3_style', path)
    plot = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(plot)
    with np.load(HERE / 'figure3_inputs/curves.npz', allow_pickle=False) as archive:
        data = dict(archive)
    rg_ratios = json.loads((HERE / 'figure3_inputs/rg_heating_response_percent.json').read_text())
    ratios = read(HERE / 'revision_2026_10_08/JCP_revision/final_dynamics_20261013/Fig3_updated_ratios.csv')
    plot.rg_cache = lambda s, f, r: {k: data[f'{s}__{f}__{r}__{k}'] for k in ('time_axis', 'acf_rg')}
    plot.load_intra_cache = lambda s, f, r: {k: data[f'{s}__{f}__{r}__{k}'] for k in ('time_ns', 'acf_mean', 'acf_sem')}
    plot.rg_heating_response_interval = lambda s, f: rg_ratios[f'{s}__{f}']

    def read_ratio(path):
        if path == plot.INTER_RATIO_BOOTSTRAP:
            metric = 'inter'
        elif path == plot.RESULTS / 'interchain_partner_survival_ratio.csv':
            metric = 'partner'
        elif path.name == 'intrachain_contact_relaxation_ratio.csv':
            return read(HERE / 'figure3_inputs/intrachain_contact_relaxation_ratio.csv')
        else:
            raise ValueError(f'Unexpected Figure 3 input: {path}')
        return [dict(system=r['system'], ratio_th_over_tl=r['ratio'],
                     ratio_ci_low=r['ci_low'], ratio_ci_high=r['ci_high'])
                for r in ratios if r['metric'] == metric]

    plot.read_csv = read_ratio

    def draw_curve(ax, system, metric, **kwargs):
        threshold = .2 if metric == 'inter' else .5
        for role in plot.ROLE_ORDER:
            prefix = f'{system}__{role}__{metric}__'
            t, y = data[prefix + 'time_ns'], data[prefix + 'mean']
            lo, hi = data[prefix + 'interval']
            color = plot.TEMP_COLORS[role]
            ax.plot(t, y, color=color, lw=1.8, zorder=3, label=plot.temp_label(role))
            ax.fill_between(t, lo, hi, color=color, alpha=.14, zorder=2, lw=0)
            plot.draw_threshold_time(ax, plot.first_crossing_time(t, y, threshold), threshold,
                                     color=color, linestyle='-', lw=1.15)
        ax.axhline(threshold, color='#555555', ls=':', lw=.9)
        ax.set_xlabel('Lag time (ns)')
        ax.set_xlim(0, 30 if metric == 'inter' else 50)
        ax.set_ylim((-.05, 1.03) if metric == 'inter' else (0, 1.03))
        plot.clean_axis(ax)

    plot.plot_interchain_acf = lambda ax, s, **kw: draw_curve(ax, s, 'inter', **kw)
    plot.plot_partner_survival = lambda ax, s, **kw: draw_curve(ax, s, 'partner', **kw)
    plot.setup_style()
    plot.build_figure(args.output_dir)
    for suffix in ('pdf', 'png'):
        (args.output_dir / f'Figure2_dynamics_story_columnwise_layout_v22.{suffix}').replace(args.output_dir / f'Fig3.{suffix}')


if __name__ == '__main__':
    main()
