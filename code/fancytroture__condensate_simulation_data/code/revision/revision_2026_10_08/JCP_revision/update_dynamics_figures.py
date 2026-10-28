#!/usr/bin/env python3
"""Build final dynamics figures from the validated physical-time contact series."""
import csv
import hashlib
import importlib.util
import json
import os
from pathlib import Path

os.environ.setdefault('MPLCONFIGDIR', '/tmp/jcp_final_figures_mpl')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT = HERE / 'final_dynamics_20261013'
CORRECTED = ROOT / 'revision_2026_09_30/interchain_recompute'
OUT.mkdir(exist_ok=True)


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def rows(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def write_csv(path, data):
    with path.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(data[0]))
        w.writeheader()
        w.writerows(data)


def crossing(t, curves, threshold):
    result = []
    for y in np.atleast_2d(curves):
        indices = np.flatnonzero(y <= threshold)
        if not len(indices):
            result.append(np.inf)
            continue
        i = indices[0]
        result.append(float(t[0]) if i == 0 else
                      float(np.interp(threshold, y[i-1:i+1][::-1], t[i-1:i+1][::-1])))
    return np.asarray(result)


def main():
    plot = load_module('original_dynamics', ROOT / 'project_figure/scripts/make_figure2_dynamics_story.py')
    time_methods = load_module('physical_time_analysis', CORRECTED / 'analyze.py')
    sources = {}
    curve_data = {}
    comparison = []
    source_ratios = rows(CORRECTED / 'analysis/temperature_ratios.csv')
    ratios = {}
    for system in plot.SYSTEM_ORDER:
        for metric, threshold in [('inter', .2), ('partner', .5)]:
            values = {}
            for role in plot.ROLE_ORDER:
                temp = plot.TEMPS[system][role].removesuffix('K')
                path = CORRECTED / 'analysis' / f'{system}_{temp}_{metric}_curves.npz'
                sources[str(path.relative_to(ROOT))] = hashlib.sha256(path.read_bytes()).hexdigest()
                with np.load(path, allow_pickle=False) as z:
                    data = {k: z[k] for k in ['time_ns', 'mean', 'unit_curves', 'block_250_boot']}
                t, y, boot = data['time_ns'], data['mean'], data['block_250_boot']
                np.testing.assert_allclose(np.diff(t), .5, atol=1e-12)
                np.testing.assert_allclose(y, data['unit_curves'].mean(axis=0), atol=1e-12)
                curve_data[system, role, metric] = data
                values[role] = crossing(t, y, threshold)[0], crossing(t, boot, threshold)
            selected = next(r for r in source_ratios if
                            (r['system'], r['metric'], r['method'], r['quantity']) ==
                            (system, metric, 'origin_block_250ns', f'crossing_{threshold:g}'))
            estimate = values['Th'][0] / values['Tl'][0]
            assert np.isfinite(values['Th'][1]).all() and np.isfinite(values['Tl'][1]).all()
            lo, hi = np.quantile(values['Th'][1] / values['Tl'][1], [.025, .975], method='inverted_cdf')
            np.testing.assert_allclose([estimate, lo, hi],
                                       [float(selected[k]) for k in ['ratio', 'ci_low', 'ci_high']], atol=1e-12)
            ratios[system, metric] = dict(system=system, ratio_th_over_tl=estimate,
                                          ratio_ci_low=lo, ratio_ci_high=hi)
            comparison.append(dict(system=system, metric=metric, threshold=threshold,
                                   ratio=estimate, ci_low=lo, ci_high=hi))

    def draw_curve(ax, system, metric, show_ylabel=False, show_legend=False):
        threshold = .2 if metric == 'inter' else .5
        for role in plot.ROLE_ORDER:
            d = curve_data[system, role, metric]
            t, y = d['time_ns'], d['mean']
            lo, hi = np.quantile(d['block_250_boot'], [.025, .975], axis=0, method='inverted_cdf')
            color = plot.TEMP_COLORS[role]
            ax.plot(t, y, color=color, lw=1.8, zorder=3, label=plot.temp_label(role))
            ax.fill_between(t, lo, hi, color=color, alpha=.14, zorder=2, lw=0)
            plot.draw_threshold_time(ax, crossing(t, y, threshold)[0], threshold,
                                     color=color, linestyle='-', lw=1.15)
        ax.axhline(threshold, color='#555555', ls=':', lw=.9)
        ax.set_xlabel('Lag time (ns)')
        ax.set_xlim(0, 30 if metric == 'inter' else 50)
        ax.set_ylim((-.05, 1.03) if metric == 'inter' else (0, 1.03))
        plot.clean_axis(ax)

    plot.plot_interchain_acf = lambda ax, system, **kw: draw_curve(ax, system, 'inter', **kw)
    plot.plot_partner_survival = lambda ax, system, **kw: draw_curve(ax, system, 'partner', **kw)
    original_read = plot.read_csv

    def revised_read(path):
        if path == plot.INTER_RATIO_BOOTSTRAP:
            return [ratios[s, 'inter'] for s in plot.SYSTEM_ORDER]
        if path == plot.RESULTS / 'interchain_partner_survival_ratio.csv':
            return [ratios[s, 'partner'] for s in plot.SYSTEM_ORDER]
        return original_read(path)

    plot.read_csv = revised_read
    plot.setup_style()
    plot.build_figure(OUT)
    for ext in ['pdf', 'png']:
        (OUT / f'Figure2_dynamics_story_columnwise_layout_v22.{ext}').replace(OUT / f'Fig3.{ext}')
    write_csv(OUT / 'Fig3_updated_ratios.csv', comparison)

    # The fixed-lag partner-event figure must use the same physical-time series.
    events = []
    for system in plot.SYSTEM_ORDER:
        for role in plot.ROLE_ORDER:
            temp = plot.TEMPS[system][role].removesuffix('K')
            path = CORRECTED / 'raw' / f'{system}_{temp}' / 'timeseries.npz'
            sources[str(path.relative_to(ROOT))] = hashlib.sha256(path.read_bytes()).hexdigest()
            with np.load(path, allow_pickle=False) as z:
                frame, seg, time, contacts = [z[k] for k in ['frame', 'segment', 'time_ns', 'partner_contacts']]
            loss, gain, n0 = [], [], []
            for s in np.unique(seg):
                ids = np.flatnonzero(seg == s)
                ids = ids[np.argsort(frame[ids])]
                for run in time_methods.runs(frame[ids], time[ids]):
                    ix = ids[run]
                    if len(ix) < 8:
                        continue
                    np.testing.assert_allclose(time[ix[4:]] - time[ix[:-4]], 2, atol=.004, rtol=0)
                    present = contacts[ix] > 0
                    a, b = present[:-4], present[4:]
                    loss.extend(np.sum(a & ~b, axis=1))
                    gain.extend(np.sum(b & ~a, axis=1))
                    n0.extend(np.sum(a, axis=1))
            loss, gain, n0 = np.asarray(loss), np.asarray(gain), np.asarray(n0)
            events.append(dict(system=system, role=role, lag_ns=2, pairs=len(n0),
                               turnover=float(np.mean((loss + gain) / np.maximum(n0, 1))),
                               net_change=float(np.mean(np.abs(gain - loss) / np.maximum(n0, 1))),
                               compensation=float(1 - np.abs(gain-loss).sum() / (gain+loss).sum())))
    write_csv(OUT / 'FigS12_partner_events.csv', events)
    fig, axes = plt.subplots(1, 3, figsize=(11.4, 3.8), layout='constrained')
    for ax, key, ylabel, label in zip(axes, ['turnover', 'net_change', 'compensation'],
                                    [r'Fractional turnover, $\mathcal{T}$', r'Fractional change, $\mathcal{D}$',
                                     r'Compensation index, $\mathcal{C}$'], ['A', 'B', 'C']):
        values = []
        for role, shift in [('Tl', -.12), ('Th', .12)]:
            y = [next(r[key] for r in events if (r['system'], r['role']) == (s, role)) for s in plot.SYSTEM_ORDER]
            ax.plot(np.arange(3)+shift, y, 'o', ms=8, color=plot.TEMP_COLORS[role], label=plot.temp_label(role))
            values += y
        delta = max(values) - min(values)
        ax.set_ylim(min(values)-.2*delta, max(values)+.35*delta)
        ax.set_xticks(np.arange(3), plot.display_systems(), fontsize=11)
        ax.set_ylabel(ylabel, fontsize=14)
        ax.tick_params(axis='y', labelsize=12)
        plot.box_axis(ax)
        ax.text(-.17, 1.04, label, transform=ax.transAxes, fontsize=20, fontweight='bold')
    axes[0].legend(ncol=2, loc='upper left', fontsize=11)
    for ext in ['pdf', 'png']:
        fig.savefig(OUT / f'FigS12.{ext}', dpi=220)
    plt.close(fig)
    (OUT / 'figure_validation.json').write_text(json.dumps(dict(
        figure3_ratios=comparison, partner_events=events, input_sha256=sources,
        curve_means_and_bootstrap_intervals_reproduced=True,
        physical_spacing_ns=.5, partner_event_lag_ns=2,
        main_figure_intramolecular_data_unchanged=True), indent=2)+'\n')
    print(json.dumps(dict(ratios=comparison, partner_events=events), indent=2))


if __name__ == '__main__':
    main()
