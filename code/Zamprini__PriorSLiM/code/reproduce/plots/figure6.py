from pathlib import Path
import os, csv, json, hashlib
from support import DATA, OUT, save as export
ROOT = OUT
OLD = DATA / 'fig06'
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
source = OLD / 'main_results.tsv'
with (OLD / 'plotted_values.tsv').open() as f:
    rows = list(csv.DictReader(f, delimiter='\t'))
with source.open() as f:
    raw = {(r['evaluation_method'], r['subclass']): r for r in csv.DictReader(f, delimiter='\t')}
for r in rows:
    original = raw['site' if r['panel'] == 'a' else 'region', r['subclass']]
    for k in ['delta', 'ci_low', 'ci_high', 'q_joint']:
        assert abs(float(r[k]) - float(original[k])) < 1e-14
plt.rcParams.update({'font.family': 'Arial', 'font.size': 7, 'axes.labelsize': 7, 'xtick.labelsize': 6.8, 'axes.linewidth': 0.55, 'axes.edgecolor': '#26333B', 'text.color': '#26333B', 'xtick.color': '#26333B', 'pdf.fonttype': 42, 'svg.fonttype': 'none', 'savefig.facecolor': 'white'})
audit = {}
for panel in ['a', 'b']:
    data = [r for r in rows if r['panel'] == panel]
    assert len(data) == 10
    fig = plt.figure(figsize=(88 / 25.4, 68 / 25.4), dpi=180)
    ax = fig.add_axes([0.375, 0.16, 0.6, 0.72])
    ax.set_xlim(-32, 85)
    ax.set_ylim(9.55, -0.8)
    ax.set_yticks([])
    ax.set_xticks([-20, 0, 20, 40, 60, 80])
    ax.tick_params(axis='x', length=2, pad=2)
    ax.spines[['top', 'left', 'right']].set_visible(False)
    ax.axvline(0, c='#82929C', lw=0.65, zorder=0)
    ax.grid(axis='x', color='#EDF0F2', lw=0.5)
    ax.set_axisbelow(True)
    ax.set_xlabel('Hit-rate difference (percentage points)', labelpad=4)
    ax.text(-0.04, 1.012, 'Annotation', transform=ax.transAxes, ha='right', fontsize=7, color='#64747F')
    labels = []
    highlighted = []
    for i, r in enumerate(data):
        delta, lo, hi, q = [float(r[k]) for k in ['delta', 'ci_low', 'ci_high', 'q_joint']]
        significant = q < 0.05
        orange = significant and delta > 0 if panel == 'a' else significant and delta < 0
        if orange:
            highlighted.append(r['label'])
        labels.append(ax.text(-0.04, i, r['label'], transform=ax.get_yaxis_transform(), ha='right', va='center', fontsize=7.2, color='#B44D12' if orange else '#26333B', fontweight='bold' if orange else 'normal', clip_on=False))
        ax.plot([100 * lo, 100 * hi], [i, i], color='#276F91', lw=0.65, alpha=0.8)
        ax.plot(100 * delta, i, 'o', ms=3.6, mfc='#276F91' if significant else 'white', mec='#276F91', mew=0.8)
    for y in [1.5, 5.5, 6.5] if panel == 'a' else [4.5]:
        ax.axhline(y, color='#E4E9ED', lw=0.45)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    bounds = [t.get_window_extent(renderer) for t in labels]
    assert min((b.x0 for b in bounds)) > 0, 'Annotation clipped at left edge'
    gap = ax.get_window_extent(renderer).x0 - max((b.x1 for b in bounds))
    assert gap > 4
    export(fig, 6, panel)
    audit[panel] = {'orange_labels': highlighted, 'label_alignment': 'right', 'axes_left_fraction': 0.375, 'plot_width_fraction': 0.6, 'label_to_plot_gap_points': gap * 72 / fig.dpi, 'min_label_left_points': min((b.x0 for b in bounds)) * 72 / fig.dpi}
    plt.close(fig)
assert len(audit['a']['orange_labels']) == 5 and len(audit['b']['orange_labels']) == 4
audit['source_results_sha256'] = sha(source)
audit['source_values_sha256'] = sha(OLD / 'plotted_values.tsv')
audit['display_rule'] = 'a: positive delta and joint q<0.05. b: negative delta and joint q<0.05. Marker fill retains statistical significance irrespective of direction.'
from pathlib import Path
import csv, json, hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
SRC = DATA / 'fig06/primary360_tests.tsv'
with SRC.open() as f:
    rows = list(csv.DictReader(f, delimiter='\t'))
plt.rcParams.update({'font.family': 'Arial', 'font.size': 15, 'text.color': '#263540', 'axes.labelcolor': '#263540', 'ytick.color': '#263540', 'axes.edgecolor': '#61717C', 'axes.linewidth': 0.7, 'svg.fonttype': 'path', 'pdf.fonttype': 42})
audit = []
specs = [('c', 'file1_union_vs_file1_outside', 'Peak-density ratio', (0.95, 1.5), [1, 1.1, 1.2, 1.3, 1.4, 1.5]), ('d', 'file2_multi_vs_one', 'Adjusted odds ratio', (0.8, 2.3), [1, 1.25, 1.5, 1.75, 2, 2.25])]
for letter, contrast, label, lim, ticks in specs:
    r = sorted([r for r in rows if r['min_run'] == '3' and r['stratum'] == 'Domain' and (r['contrast'] == contrast)], key=lambda r: int(r['fold']))
    assert len(r) == 10 and all((x['status'] == 'estimated' for x in r))
    values = np.array([float(x['effect_ratio']) for x in r])
    sig = [float(x['q_bh_primary360']) < 0.05 for x in r]
    fig, ax = plt.subplots(figsize=(7.1, 3.2))
    fig.subplots_adjust(left=0.17, right=0.97, bottom=0.23, top=0.95)
    ax.spines[['top', 'right', 'bottom']].set_visible(False)
    ax.set_axisbelow(True)
    ax.grid(axis='y', lw=0.6, color='#E5EAED')
    ax.axhline(1, ls=(0, (4, 4)), lw=1, color='#8A99A3')
    offsets = np.array([-0.14, -0.105, -0.07, -0.035, 0, 0.035, 0.07, 0.105, 0.14, 0.175])
    color = '#287C9D' if letter == 'c' else '#8B6FA8'
    ax.scatter(1 + offsets, values, s=68, facecolors=[color if s else 'white' for s in sig], edgecolors=color, linewidths=1.1, zorder=3)
    ax.plot([0.77, 1.23], [np.median(values)] * 2, color='#263540', lw=1.8, zorder=2)
    ax.set_xlim(0.6, 1.4)
    ax.set_ylim(*lim)
    ax.set_xticks([])
    ax.set_yticks(ticks)
    ax.tick_params(axis='y', length=0, pad=8, labelsize=13)
    ax.set_ylabel(label, labelpad=12, fontsize=15)
    if letter == 'd':
        fig.legend(handles=[Line2D([], [], ls='', marker='o', ms=6, color='#263540', markerfacecolor='#263540', label='BH q < 0.05'), Line2D([], [], ls='', marker='o', ms=6, color='#263540', markerfacecolor='white', markeredgewidth=1.1, label='BH q ≥ 0.05'), Line2D([], [], color='#263540', lw=1.8, label='Median')], loc='lower center', bbox_to_anchor=(0.55, 0.005), ncol=3, frameon=False, fontsize=11, columnspacing=1.3, handlelength=1.4)
    export(fig, 6, letter)
    plt.close(fig)
    audit.append({'panel': letter, 'contrast': contrast, 'folds': [int(x['fold']) for x in r], 'estimates': values.tolist(), 'q': [float(x['q_bh_primary360']) for x in r], 'median': float(np.median(values)), 'significant': sum(sig)})
