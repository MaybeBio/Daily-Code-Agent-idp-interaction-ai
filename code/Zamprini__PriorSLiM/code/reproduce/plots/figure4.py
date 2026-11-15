from pathlib import Path
import os, csv, hashlib, json
from support import DATA, save as export
ROOT = Path(__file__).resolve().parent
FIG = DATA / 'fig04'
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
COLORS = ['#3B7097', '#59A092', '#D18066', '#9783AD', '#939CA6']
INK = '#27323A'
DS = ['ASHI', 'Phospho_WT', 'Nonhuman', 'MOTIF']
HEAD = ['Human IDRs (ASHI)', 'Phospho-screen WT', 'Viral IDRs', 'Swiss-Prot MOTIF']
METHODS = ['PriorSLiM', 'SLiMMine', 'SHARK-Capture', 'ELM']
MARKERS = ['+', 's', 'D', '^', 'o']
HEIGHTS = {'a': 340, 'b': 660, 'c': 300, 'd': 520}
plt.rcParams.update({'font.family': 'Arial', 'font.size': 7.5, 'axes.labelsize': 8, 'axes.titlesize': 8, 'xtick.labelsize': 7, 'ytick.labelsize': 7, 'text.color': INK, 'axes.labelcolor': INK, 'xtick.color': INK, 'ytick.color': INK, 'axes.edgecolor': INK, 'axes.linewidth': 0.55, 'pdf.fonttype': 42, 'svg.fonttype': 'none', 'savefig.facecolor': 'white', 'lines.solid_capstyle': 'round'})

def table(relative):
    p = FIG / Path(relative).name
    with p.open() as f:
        rows = list(csv.DictReader(f, delimiter='\t'))
    for row in rows:
        for k, v in row.items():
            try:
                row[k] = float(v) if v else float('nan')
            except ValueError:
                pass
    return rows

def figure(letter):
    fig = plt.figure(figsize=(7.2, HEIGHTS[letter] / 200))
    fig.text(0.017, 0.975, {'b': 'c', 'c': 'd', 'd': 'e'}[letter], fontsize=11, fontweight='bold', va='top')
    return fig

def style(ax, ylim, ticks, ylabel=None):
    ax.set_xlim(0, 4)
    ax.set_ylim(*ylim)
    ax.set_yticks(ticks)
    for side in ['top', 'right']:
        ax.spines[side].set_visible(False)
    ax.tick_params(axis='both', length=2, width=0.5, pad=2)
    if ylabel:
        ax.set_ylabel(ylabel, labelpad=7)
    for x in [1, 2, 3]:
        ax.axvline(x, color='#AEB6BD', ls=(0, (3, 3)), lw=0.55, zorder=0)

def headings(ax):
    for i, h in enumerate(HEAD):
        ax.text(i + 0.5, 1.075, h, transform=ax.get_xaxis_transform(), ha='center', va='bottom', fontsize=8)

def tx(i, value, count):
    return i + 0.07 + 0.86 * value / count

def method_handles():
    return [Line2D([], [], ls='none', marker=m, color=c, markersize=5, label=n) for n, c, m in zip(METHODS, COLORS, MARKERS)]

def save(fig, letter, audit):
    export(fig, 4, {'b': 'c', 'c': 'd', 'd': 'e'}[letter])

def panel_b():
    peak = [r for r in table('four_cohorts/phospho_update_20261010/peak_summary.tsv') if r['method'] in METHODS]
    top = [r for r in table('four_cohorts/phospho_update_20261010/topk_summary.tsv') if r['geometry'] == 'native']
    assert len(peak) == 16 and len(top) == 24
    assert all((np.isfinite(r['enrichment_mean']) for r in peak))
    fig = figure('b')
    upper = fig.add_axes([0.09, 0.57, 0.88, 0.31])
    lower = fig.add_axes([0.09, 0.18, 0.88, 0.28])
    style(upper, (0, 4), np.arange(0, 4.1, 0.5), 'Observed / expected')
    headings(upper)
    upper.axhline(1, color='#939CA6', ls=(0, (3, 3)), lw=0.7)
    xt = []
    labels = []
    audit = []
    for i, ds in enumerate(DS):
        for j, m in enumerate(METHODS):
            r = next((r for r in peak if r['dataset'] == ds and r['method'] == m))
            x = tx(i, j + 0.5, len(METHODS))
            y = r['enrichment_mean']
            c = COLORS[j]
            xt.append(x)
            labels.append(['Prior-\nSLiM', 'SLiM-\nMine', 'SHARK\nCapture', 'ELM'][j])
            upper.errorbar(x, y, yerr=[[y - r['enrichment_ci_low']], [r['enrichment_ci_high'] - y]], fmt=MARKERS[j], ms=6 if j == 0 else 4, mew=0.8, color=c, lw=0.6, capsize=0, zorder=3)
            dx, dy = (5, 0) if j < 3 else (0, -8 if j == 3 else 7)
            if j == 1:
                dx, dy = (0, -8)
            if y == 0:
                dx, dy = (0, 7)
            upper.annotate(f'{y:.2f}×', (x, y), xytext=(dx, dy), textcoords='offset points', fontsize=7, weight='bold', color=c, ha='left' if dx else 'center', va='center', zorder=5)
            audit.append({'dataset': ds, 'method': m, 'metric': 'enrichment', 'value': y, 'label': f'{y:.2f}×'})
    upper.set_xticks(xt)
    upper.set_xticklabels(labels, fontsize=5.5)
    upper.tick_params(axis='x', length=0, pad=6)
    style(lower, (0, 100), [0, 25, 50, 75, 100], 'Recovered loci (%)')
    xt = []
    labels = []
    for i, ds in enumerate(DS):
        for j, m in enumerate(METHODS[:2]):
            r = sorted([r for r in top if r['dataset'] == ds and r['method'] == m], key=lambda r: r['k'])
            x = np.array([tx(i, k + 0.5, 3) for k in range(3)])
            y = 100 * np.array([r['recall_mean'] for r in r])
            sd = 100 * np.array([r['recall_fold_sd'] for r in r]) if j == 0 else None
            lower.errorbar(x, y, yerr=sd, fmt=MARKERS[j] + '-', color=COLORS[j], lw=1.1, ms=6 if j == 0 else 4, mew=0.8, capsize=0, elinewidth=0.6, zorder=3)
            for xx, yy, rr in zip(x, y, r):
                lower.annotate(f'{yy:.1f}%', (xx, yy), xytext=(0, 8 if j == 0 else -6), textcoords='offset points', ha='center', va='center', fontsize=7, color=COLORS[j], weight='bold')
                audit.append({'dataset': ds, 'method': m, 'metric': f"Top-{int(rr['k'])}", 'value': yy, 'label': f'{yy:.1f}%'})
        xt.extend([tx(i, k + 0.5, 3) for k in range(3)])
        labels.extend(['1', '3', '5'])
    lower.set_xticks(xt)
    lower.set_xticklabels(labels)
    lower.tick_params(axis='x', length=0, pad=6)
    lower.set_xlabel('Maximum candidates per IDR (Top-K)', labelpad=4)
    fig.legend(handles=method_handles(), loc='center', bbox_to_anchor=(0.53, 0.04), ncol=4, frameon=False, fontsize=7.7, columnspacing=2.5, handletextpad=0.5)
    assert len(audit) == 40
    save(fig, 'b', {'point_labels': audit, 'count': len(audit), 'undefined_enrichment_points': 0, 'methods': METHODS, 'topk_protocol': 'native_thresholds', 'prior_filter': 'probability >= fold-specific Validation single-window threshold', 'slimmine_filter': 'final_score > 0.9 AND flank_ratio < 0.5', 'geometry': 'native', 'topk': [1, 3, 5]})

def panel_c():
    rows = table('four_cohorts/source_data/pattern_summary.tsv')
    folds = table('four_cohorts/source_data/pattern_per_fold.tsv')
    fig = figure('c')
    ax = fig.add_axes([0.09, 0.24, 0.88, 0.55])
    style(ax, (30, 90), np.arange(30, 91, 10), 'Target–background\nwin rate (%)')
    headings(ax)
    ax.axhline(50, color=COLORS[4], ls=(0, (3, 3)), lw=0.7)
    ax.set_xticks([])
    audit = []
    for i, ds in enumerate(DS):
        r = next((r for r in rows if r['dataset'] == ds))
        f = sorted([r for r in folds if r['dataset'] == ds], key=lambda r: r['fold'])
        assert len(f) == 10
        y = 100 * r['win_mean']
        lo = 100 * r['win_ci_low']
        hi = 100 * r['win_ci_high']
        ax.scatter(i + 0.34 + np.arange(10) * 0.028, [100 * r['pairwise_win'] for r in f], s=13, color='#A4C0CF', zorder=3)
        ax.errorbar(i + 0.5, y, yerr=[[y - lo], [hi - y]], fmt='D', ms=5, color=COLORS[0], lw=0.8, capsize=0, zorder=4)
        ax.annotate(f'{y:.2f}%', (i + 0.5, max(hi, y)), xytext=(0, 5), textcoords='offset points', ha='center', fontsize=8, color=COLORS[0], weight='bold')
        audit.append({'dataset': ds, 'mean': y, 'ci': [lo, hi], 'folds': [100 * r['pairwise_win'] for r in f]})
    handles = [Line2D([], [], ls='none', marker='o', color='#A4C0CF', ms=4, label='Individual fold'), Line2D([], [], ls='none', marker='D', color=COLORS[0], ms=4, label='Ten-fold mean'), Line2D([], [], ls='--', color=COLORS[4], lw=0.7, label='Chance (50%)')]
    fig.legend(handles=handles, loc='center', bbox_to_anchor=(0.53, 0.09), ncol=3, frameon=False, fontsize=7.5, columnspacing=3)
    save(fig, 'c', audit)

def panel_d():
    rows = table('combined_center_profiles_20261010/combined_profiles.tsv')
    fig = figure('d')
    ax = fig.add_axes([0.09, 0.36, 0.88, 0.47])
    support = fig.add_axes([0.09, 0.205, 0.88, 0.095], sharex=ax)
    style(ax, (25, 80), np.arange(25, 81, 5), 'Mean within-IDR\npercentile (%)')
    style(support, (0, 100), [0, 100])
    fig.text(0.027, 0.252, 'Available\nloci (%)', ha='center', va='center', fontsize=6.5)
    ax.tick_params(axis='x', bottom=False, labelbottom=False)
    ax.axhline(50, color='#9A8B6F', ls=(0, (4, 3)), lw=0.75, zorder=1)
    xt = []
    xl = []
    audit = []
    for i, ds in enumerate(DS):
        q = sorted([r for r in rows if r['dataset'] == ds], key=lambda r: r['offset_aa'])
        assert len(q) == 101
        x = tx(i, np.array([r['offset_aa'] for r in q]) + 50, 100)
        v = lambda k: np.array([r[k] for r in q])
        for prefix, c in [('all', COLORS[0]), ('unhit', COLORS[2])]:
            ax.fill_between(x, v(prefix + '_ci_low'), v(prefix + '_ci_high'), color=c, alpha=0.17, lw=0, zorder=1)
            ax.plot(x, v(prefix + '_mean'), color=c, lw=1.5, zorder=4)
        ax.plot(x, v('unhit_conditional_null'), color=COLORS[4], ls=(0, (3, 2)), lw=0.8, zorder=2)
        a, b = (100 * v('all_available_fraction'), 100 * v('unhit_available_fraction'))
        support.fill_between(x, 0, np.minimum(a, b), color=COLORS[0], alpha=0.25, lw=0, zorder=1)
        support.fill_between(x, np.minimum(a, b), np.maximum(a, b), color=COLORS[2], alpha=0.3, lw=0, zorder=2)
        support.plot(x, a, color=COLORS[0], lw=1.4, zorder=3)
        support.plot(x, b, color=COLORS[2], lw=1.4, zorder=4)
        ax.text(i + 0.5, 1.18, HEAD[i], transform=ax.get_xaxis_transform(), ha='center', fontsize=8)
        center = next((r for r in q if r['offset_aa'] == 0))
        ax.text(i + 0.27, 1.095, f"All {center['all_mean']:.1f}%", transform=ax.get_xaxis_transform(), ha='center', color=COLORS[0], fontsize=7)
        ax.text(i + 0.74, 1.095, f"Unhit {center['unhit_mean']:.1f}%", transform=ax.get_xaxis_transform(), ha='center', color=COLORS[2], fontsize=7)
        xt.extend(tx(i, np.array([0, 25, 50, 75, 100]), 100))
        xl.extend(['−50', '−25', '0', '25', '50'])
        audit.append({'dataset': ds, 'all_curve': v('all_mean').tolist(), 'unhit_curve': v('unhit_mean').tolist(), 'availability_all': a.tolist(), 'availability_unhit': b.tolist(), 'ci_clipped_all': int(np.sum((v('all_ci_low') < 25) | (v('all_ci_high') > 80))), 'ci_clipped_unhit': int(np.sum((v('unhit_ci_low') < 25) | (v('unhit_ci_high') > 80)))})
    support.set_xticks(xt)
    support.set_xticklabels(xl, fontsize=6.7)
    support.tick_params(axis='x', length=0, pad=4)
    support.set_xlabel('Position relative to the experimental peptide / annotation center (aa)', labelpad=6, fontsize=7.5)
    handles = [Line2D([], [], color=COLORS[0], lw=1.5, label='All loci (including hits)'), Line2D([], [], color=COLORS[2], lw=1.5, label='Unhit loci'), Patch(fc='#AAB5BD', alpha=0.45, label='95% CI'), Line2D([], [], color='#9A8B6F', lw=0.8, ls='--', label='50% reference'), Line2D([], [], color=COLORS[4], lw=0.8, ls='--', label='Unhit shift null')]
    fig.legend(handles=handles, loc='center', bbox_to_anchor=(0.53, 0.05), ncol=5, frameon=False, fontsize=7, columnspacing=1.3, handlelength=2)
    save(fig, 'd', {'axis': [25, 80], 'step': 5, 'curves': audit, 'availability_fill': '0 to min = overlap; min to max = nonoverlap'})

def panel_topk():
    import pandas as pd
    import statistics
    rows = pd.read_csv(FIG / 'pooled_summary.tsv', sep='\t')
    folds = pd.read_csv(FIG / 'pooled_per_fold.tsv', sep='\t')
    fig, ax = plt.subplots(figsize=(6, 3.6), layout='constrained')
    for method, label, color in zip(['PriorSLiM', 'Explicit_only', 'MotifAE_only'], ['PriorSLiM', 'Explicit-only', 'MotifAE-only'], COLORS):
        q = rows[(rows.method == method) & (rows['mode'] == 'single_window')].sort_values('k')
        assert q.k.tolist() == list(range(1, 11))
        for r in q.itertuples():
            f = folds[(folds.method == method) & (folds['mode'] == 'single_window') & (folds.k == r.k)]
            assert len(f) == 10 and r.loci == 4286
            assert abs(f.recall.mean() - r.recall_mean) < 1e-12
            assert abs(f.recall.std(ddof=1) - r.recall_sd) < 1e-12
        ax.errorbar(q.k, 100 * q.recall_mean, yerr=100 * q.recall_sd, color=color, marker='o', label=label, capsize=2)
    ax.set(xlabel='Maximum candidates per IDR (K)', ylabel='Recovered sites (%)', xticks=range(1, 11))
    ax.legend(frameon=False)
    export(fig, 4, 'b')
from matplotlib.patches import Patch
panel_topk()
panel_b()
panel_c()
panel_d()
