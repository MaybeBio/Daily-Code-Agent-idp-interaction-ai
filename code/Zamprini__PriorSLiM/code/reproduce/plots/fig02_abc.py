from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.text import Text
from matplotlib.colors import LinearSegmentedColormap, Normalize, ListedColormap, to_rgb, to_rgba
from sklearn.metrics import precision_recall_curve
METHODS = ['PriorSLiM', 'SLiMMine', 'SHARK', 'ELM']
COLORS = {'PriorSLiM': '#276D91', 'SLiMMine': '#BE8345', 'SHARK': '#679995', 'ELM': '#737B82'}
INK = '#27323A'
GRAY = '#A6ADB1'
TYPES = ['natural_background', 'motif_matched', 'high_conservation', 'chemistry_context_matched', 'critical_mutation', 'offset_core_excluded']
TYPE_LABELS = ['Natural', 'Motif', 'Conserved', 'Chemistry', 'Critical', 'Offset', 'Recall']

def style(ax):
    ax.spines[['top', 'right']].set_visible(False)
    ax.tick_params(length=2, width=0.8, pad=1.5)
RC = {'font.family': 'Arial', 'font.size': 7, 'axes.labelsize': 7, 'axes.titlesize': 8, 'xtick.labelsize': 7, 'ytick.labelsize': 7, 'axes.labelpad': 2, 'pdf.fonttype': 42, 'svg.fonttype': 'none', 'axes.linewidth': 1, 'text.color': INK, 'axes.labelcolor': INK, 'xtick.color': INK, 'ytick.color': INK}

def plot_a(y, models, metrics, strata, within):
    w, h = (65, 54)
    fig = plt.figure(figsize=(w / 25.4, h / 25.4))
    ax = fig.add_axes([10 / w, 13 / h, 53 / w, 31 / h])
    grid = np.linspace(0, 1, 501)
    for name in ['PriorSLiM', 'SLiMMine', 'SHARK']:
        group = [m for m in models if m['method'] == name]
        aps = metrics[metrics.method.eq(name)].ap.mean()
        if len(group) > 1:
            curves = []
            for m in group:
                p, r, _ = precision_recall_curve(y, m['score'])
                curves.append(np.interp(grid, r[::-1], p[::-1]))
            curves = np.array(curves)
            mean = curves.mean(0)
            sd = curves.std(0, ddof=1)
            ax.fill_between(grid, np.clip(mean - sd, 0, 1), np.clip(mean + sd, 0, 1), color=COLORS[name], alpha=0.12, lw=0)
            ax.plot(grid, mean, color=COLORS[name], lw=1.2, label=f'{name} ({aps:.3f})', linestyle='-' if name == 'PriorSLiM' else '--')
        else:
            p, r, _ = precision_recall_curve(y, group[0]['score'])
            ax.plot(r, p, color=COLORS[name], lw=1.1, label=f'{name} ({aps:.3f})')
    elm = metrics[metrics.method.eq('ELM')].iloc[0]
    ax.scatter(elm.recall, elm.precision, marker='D', s=15, color=COLORS['ELM'], edgecolors='white', linewidths=0.45, zorder=5, label='ELM (binary)')
    ax.axhline(y.mean(), color=GRAY, ls=':', lw=0.8)
    ax.set(xlim=(0, 1.025), ylim=(0, 1.03), xticks=[0, 0.5, 1], yticks=[0, 0.5, 1], xlabel='Recall', ylabel='Precision')
    ax.set_title('Test-set discrimination', pad=6)
    ax.legend(frameon=False, fontsize=7, loc='upper right', handlelength=1.2, labelspacing=0.25, borderaxespad=0.2, handletextpad=0.4)
    style(ax)
    return fig

def plot_b(y, models, metrics, strata, within):
    w, h = (62, 54)
    fig = plt.figure(figsize=(w / 25.4, h / 25.4))
    ax = fig.add_axes([16 / w, 14 / h, 44 / w, 30 / h])
    values = []
    for name in METHODS:
        s = strata[strata.method.eq(name)]
        values.append([s[s.negative_type.eq(k)].fpr.mean() * 100 for k in TYPES] + [s.recall.mean() * 100])
    vals = np.array(values).T
    cmap = LinearSegmentedColormap.from_list('rates', ['#F4F7F8', '#96BECE', COLORS['PriorSLiM']])
    ax.imshow(vals, aspect='auto', cmap=cmap, vmin=0, vmax=100)
    for i in range(7):
        for j in range(4):
            ax.text(j, i, f'{vals[i, j]:.1f}', ha='center', va='center', fontsize=7, color='white' if vals[i, j] > 65 else INK)
    ax.axhline(5.5, color='white', lw=2)
    ax.set(yticks=range(7), yticklabels=TYPE_LABELS, xticks=range(4), xticklabels=METHODS)
    ax.tick_params(axis='x', rotation=60, length=0, pad=1)
    ax.tick_params(axis='y', length=0, pad=2)
    ax.set_title('False positives and recall (%)', pad=6)
    for spine in ax.spines.values():
        spine.set_visible(False)
    return fig

def plot_c(y, models, metrics, strata, within):
    w, h = (41, 54)
    fig = plt.figure(figsize=(w / 25.4, h / 25.4))
    ax = fig.add_axes([9 / w, 14 / h, 30 / w, 30 / h])
    groups = within.groupby(['method', 'family']).auc.mean()
    rng = np.random.default_rng(20261013)
    values = [groups.loc[name].to_numpy() for name in METHODS]
    bp = ax.boxplot(values, positions=range(4), widths=0.5, patch_artist=True, showfliers=False, medianprops=dict(color=INK, lw=0.9), whiskerprops=dict(lw=0.8), capprops=dict(lw=0.8))
    for i, name in enumerate(METHODS):
        bp['boxes'][i].set(facecolor=COLORS[name], alpha=0.24)
        ax.scatter(i + rng.uniform(-0.2, 0.2, 62), values[i], s=3.1, color=COLORS[name], alpha=0.55, lw=0)
        ax.scatter(i, np.mean(values[i]), marker='D', s=12, color=COLORS[name], edgecolor='white', lw=0.4, zorder=4)
    ax.axhline(0.5, color=GRAY, ls='--', lw=0.8)
    ax.set(xticks=range(4), xticklabels=METHODS, xlim=(-0.55, 3.55), ylim=(0, 1.04), yticks=[0, 0.5, 1], ylabel='AUROC across 62 classes')
    ax.tick_params(axis='x', rotation=60, length=0, pad=1)
    ax.set_title('Within-motif ranking', pad=6)
    style(ax)
    return fig
