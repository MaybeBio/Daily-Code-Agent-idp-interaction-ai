#!/usr/bin/env python3
"""Reproduce the revised SI figures from the accompanying compact input tables."""
import argparse
import csv
import importlib.util
import json
import os
from pathlib import Path
import runpy
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
REV = ROOT / 'revision_2026_10_08'
JCP = REV / 'JCP_revision'
OUTPUT = ROOT.parent / 'outputs/revision_20261018'
GROUPS = ('main', 'phase', 'states', 'models', 'sensitivity', 'mobility', 'tj', 'events')


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def read(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def collect(source, label, response=None):
    for suffix in ('pdf', 'png'):
        path = source.with_suffix('.' + suffix)
        if not path.is_file():
            raise FileNotFoundError(path)
        shutil.copy2(path, OUTPUT / f'{label}.{suffix}')
        if response:
            shutil.copy2(path, OUTPUT / f'{response}.{suffix}')


def run_group(group):
    OUTPUT.mkdir(parents=True, exist_ok=True)
    if group == 'main':
        subprocess.run([sys.executable, str(ROOT / 'make_figure3.py')], check=True)
    elif group == 'phase':
        path = JCP / 'phase_diagram_20261013/plot_phase_diagram.py'
        runpy.run_path(str(path), run_name='__main__')
        collect(path.parent / 'Fig_SI_phase_diagram', 'FigS2')
    elif group == 'states':
        path = REV / 'R2C2_state_remodeling/response_package/build_package.py'
        module = load('common_state_figure', path)
        module.plt.rcParams.update({'font.size': 13, 'axes.labelsize': 14,
                                   'axes.titlesize': 15, 'xtick.labelsize': 12,
                                   'ytick.labelsize': 12, 'pdf.fonttype': 42,
                                   'ps.fonttype': 42})
        rows = read(path.parent / 'Table_R2C2_2_metrics.csv')
        for row in rows:
            for key in ('estimate', 'ci_low', 'ci_high'):
                row[key] = float(row[key])
        module.validate_metrics()
        module.combined_plot(rows)
        collect(path.parent / 'Fig_R2C2_common_transition_and_population', 'FigS8', 'FigR2')
    elif group == 'models':
        path = REV / 'R2C3_model_consistency/analyze.py'
        module = load('model_consistency_figure', path)
        module.self_tests()
        module.plot(read(path.parent / 'Table_C3_1_all_populations.csv'),
                    read(path.parent / 'Table_C3_2_state_agreement.csv'),
                    read(path.parent / 'condition_summary.csv'))
        collect(path.parent / 'Fig_C3_state_population_and_agreement', 'FigS9', 'FigR3')
    elif group == 'sensitivity':
        path = JCP / 'response_revision_20261016/plot_figures.py'
        module = load('final_ratio_figures', path)
        module.sensitivity()
        module.residence()
        for name, si, response in (
            ('Fig_R4_all_observable_thresholds', 'FigS11', 'FigR4'),
            ('Fig_R5_all_observable_integrals', 'FigS12', 'FigR5'),
            ('Fig_R7_all_observable_residence', 'FigS14', 'FigR7'),
        ):
            collect(module.OUT / name, si, response)
        for figure in ('R4', 'R5', 'R7'):
            assert {r['metric'] for r in module.VALUES if r['figure'] == figure} == set(module.METRICS)
        (OUTPUT / 'sensitivity_plotted_values.json').write_text(json.dumps(module.VALUES, indent=2) + '\n')
    elif group == 'mobility':
        path = JCP / 'remove_R6B_20261016/plot_mobility.py'
        runpy.run_path(str(path), run_name='__main__')
        collect(path.parent / 'figures/Fig_R6_mobility', 'FigS13', 'FigR6')
    elif group == 'tj':
        path = JCP / 'tj_notation_20261017/make_figures.py'
        runpy.run_path(str(path), run_name='__main__')
        collect(path.parent / 'figures/Fig_control_density_only', 'FigS16')
        collect(path.parent / 'figures/Fig_R1_all_observable_controls', 'FigS17', 'FigR1')
    elif group == 'events':
        # This is the plotting block from update_dynamics_figures.py, evaluated
        # on its saved 2 ns partner-event summary instead of reopening raw series.
        import numpy as np
        style = load('event_style', ROOT / 'project_figure/scripts/make_figure2_dynamics_story.py')
        style.setup_style()
        plt = style.plt
        events = read(JCP / 'final_dynamics_20261013/FigS12_partner_events.csv')
        fig, axes = plt.subplots(1, 3, figsize=(11.4, 3.8), layout='constrained')
        for ax, key, ylabel, label in zip(axes, ['turnover', 'net_change', 'compensation'],
                [r'Fractional turnover, $\mathcal{T}$', r'Fractional change, $\mathcal{D}$',
                 r'Compensation index, $\mathcal{C}$'], ['A', 'B', 'C']):
            values = []
            for role, shift in [('Tl', -.12), ('Th', .12)]:
                y = [float(next(r[key] for r in events if (r['system'], r['role']) == (s, role))) for s in style.SYSTEM_ORDER]
                ax.plot(np.arange(3) + shift, y, 'o', ms=8, color=style.TEMP_COLORS[role], label=style.temp_label(role))
                values += y
            delta = max(values) - min(values)
            ax.set_ylim(min(values) - .2 * delta, max(values) + .35 * delta)
            ax.set_xticks(np.arange(3), style.display_systems(), fontsize=11)
            ax.set_ylabel(ylabel, fontsize=14)
            ax.tick_params(axis='y', labelsize=12)
            style.box_axis(ax)
            ax.text(-.17, 1.04, label, transform=ax.transAxes, fontsize=20, fontweight='bold')
        axes[0].legend(ncol=2, loc='upper left', fontsize=11)
        for ext in ('pdf', 'png'):
            fig.savefig(OUTPUT / f'FigS15.{ext}', dpi=220)
        plt.close(fig)
    print(f'Completed {group}', flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--only', choices=GROUPS)
    args = parser.parse_args()
    os.environ.setdefault('MPLBACKEND', 'Agg')
    os.environ.setdefault('MPLCONFIGDIR', '/tmp/jcp_revision_figures_mpl')
    os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
    os.environ.setdefault('OMP_NUM_THREADS', '1')
    if args.only:
        run_group(args.only)
    else:
        # Each figure family gets the original independent Matplotlib settings.
        for group in GROUPS:
            subprocess.run([sys.executable, str(Path(__file__).resolve()), '--only', group], check=True)
        print(f'Figures written to {OUTPUT}')


if __name__ == '__main__':
    main()
