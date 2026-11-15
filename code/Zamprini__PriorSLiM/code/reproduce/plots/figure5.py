from pathlib import Path
import os, sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
import pandas as pd
import numpy as np
from support import DATA, save as export
OUT = DATA / 'fig05'
BLUE, ORANGE, INK, ACCESS = ('#3B7097', '#D99F68', '#27323A', '#59A092')
PURPLE, ROSE = ('#9783AD', '#BB6576')
plt.rcParams.update({'font.family': 'Arial', 'font.size': 8, 'axes.labelsize': 7.6, 'axes.linewidth': 0.65, 'xtick.labelsize': 6.4, 'ytick.labelsize': 7, 'text.color': INK, 'axes.labelcolor': INK, 'xtick.color': INK, 'ytick.color': INK, 'svg.fonttype': 'path', 'pdf.fonttype': 42})

def clean(ax):
    ax.spines[['top', 'right']].set_visible(False)
    ax.tick_params(length=2.3, width=0.6, pad=2)

def direction(letter):
    group = '正向对照' if letter == 'a' else '主验证'
    d = pd.read_csv(DATA / 'experiments/mutation.tsv', sep='\t')
    d = d[d.group.eq(group)].sort_values(['paper', 'gene', 'mutation']).reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(3.42, 2.35))
    fig.subplots_adjust(left=0.15, right=0.985, top=0.84, bottom=0.37)
    handles = [Line2D([], [], ls='', marker='^', ms=4, color=INK, label='Gain'), Line2D([], [], ls='', marker='v', ms=4, color=INK, label='Loss'), Line2D([], [], ls='', marker='o', ms=4, color=INK, label='Concordant'), Line2D([], [], ls='', marker='o', ms=4, color=INK, mfc='white', label='Opposite')]
    if letter == 'a':
        fig.legend(handles=handles, loc='upper left', bbox_to_anchor=(0.15, 0.98), ncol=4, frameon=False, fontsize=5.9, columnspacing=0.65, handletextpad=0.25, handlelength=0.8, borderaxespad=0)
    nA = d.paper.eq('A').sum()
    x = np.arange(len(d), dtype=float)
    x[nA:] += 0.9
    ax.axhline(0, color='#9DA6AB', lw=0.65, zorder=0)
    ax.axvline(nA - 0.05, color='#CDD3D6', lw=0.6, zorder=0)
    for i, r in d.iterrows():
        color = BLUE if r.paper == 'A' else PURPLE
        ax.vlines(x[i], 0, r.delta_final_logit, color='#BFC9CE', lw=0.65, zorder=1)
        ax.scatter(x[i], r.delta_final_logit, marker='^' if r.expected_sign > 0 else 'v', s=23, facecolor=color if r.correct_final_logit else 'white', edgecolor=color, lw=0.85, zorder=3)
    for paper in ['A', 'B']:
        ix = np.flatnonzero(d.paper.eq(paper))
        g = d.iloc[ix]
        ax.text(x[ix].mean(), 1.025, paper, ha='center', transform=ax.get_xaxis_transform(), fontsize=7.1, color=BLUE if paper == 'A' else PURPLE)
    ax.set_xticks(x, [r.gene + ' ' + r.mutation.split('.')[1] for r in d.itertuples()], rotation=62, ha='right', rotation_mode='anchor', fontsize=5.6)
    ax.set_xlim(-0.65, x[-1] + 0.65)
    ax.set_ylim((-5.35, 1.4) if letter == 'a' else (-2.55, 1.0))
    ax.set_yticks([-4, -2, 0] if letter == 'a' else [-2, -1, 0])
    ax.set_ylabel('Centre Δlogit', labelpad=2)
    clean(ax)
    save(fig, letter)

def branches():
    d = pd.read_csv(DATA / 'experiments/mutation.tsv', sep='\t')
    d = d[d.group.eq('主验证')].sort_values(['paper', 'gene', 'mutation']).reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(6.94, 2.05))
    fig.subplots_adjust(left=0.075, right=0.992, top=0.86, bottom=0.34)
    x = np.arange(len(d), dtype=float)
    nA = d.paper.eq('A').sum()
    x[nA:] += 0.7
    pos = np.zeros(len(d))
    neg = pos.copy()
    for field, color, label in [('intrinsic_logit', BLUE, 'Intrinsic'), ('accessibility_logit_change', ACCESS, 'Accessibility'), ('motifae_logit', ORANGE, 'MotifAE')]:
        y = d['delta_' + field].to_numpy()
        ax.bar(x, y, bottom=np.where(y >= 0, pos, neg), width=0.61, color=color, label=label, linewidth=0)
        pos += np.maximum(y, 0)
        neg += np.minimum(y, 0)
    ax.axhline(0, color='#9DA6AB', lw=0.6)
    ax.axvline(nA - 0.15, color='#CDD3D6', lw=0.6)
    for i, r in d.iterrows():
        ax.scatter(x[i], r.delta_final_logit, marker='^' if r.expected_sign > 0 else 'v', s=22, facecolor=INK if r.correct_final_logit else 'white', edgecolor=INK, lw=0.75, zorder=4)
    ax.set_xticks(x, [r.gene + ' ' + r.mutation.split('.')[1] for r in d.itertuples()], rotation=50, ha='right', rotation_mode='anchor', fontsize=6.4)
    ax.set_xlim(-0.7, x[-1] + 0.7)
    ax.set_ylim(-2.65, 1.2)
    ax.set_yticks([-2, -1, 0, 1])
    ax.set_ylabel('Branch Δlogit', labelpad=3)
    fig.legend(*ax.get_legend_handles_labels(), loc='upper right', bbox_to_anchor=(0.99, 0.99), ncol=3, frameon=False, fontsize=7, handlelength=1.1, columnspacing=1)
    clean(ax)
    save(fig, 'c')
CASES = {'d': ('P52952.A61V', 'Main validation: concordant gain'), 'e': ('Q14444.F372A', 'Main validation: concordant loss'), 'f': ('Q9Y2J4.R71I', 'Main validation: opposite direction')}

def curve(letter):
    mutation, subtitle = CASES[letter]
    d = pd.read_csv(OUT / 'curve_plot_data.tsv', sep='\t')
    d = d[d.mutation.eq(mutation)].sort_values('window_center')
    r = d.iloc[0]
    x = d.window_center.to_numpy() - r.position
    wt = d.final_logit_WT.to_numpy()
    mut = d.final_logit_MUT.to_numpy()
    dz = mut - wt
    i = np.flatnonzero(x == 0).item()
    fig = plt.figure(figsize=(2.23, 2.24))
    grid = fig.add_gridspec(2, 1, height_ratios=[2, 1], hspace=0.13, left=0.22, right=0.98, bottom=0.22, top=0.87)
    ax = fig.add_subplot(grid[0])
    da = fig.add_subplot(grid[1], sharex=ax)
    ax.plot(x, wt, color=BLUE, lw=1.2, marker='o', ms=2.0, label='WT')
    ax.plot(x, mut, color=ROSE, lw=1.2, marker='s', ms=1.9, label='MUT')
    ax.fill_between(x, wt, mut, color=ROSE, alpha=0.1, lw=0)
    for yy, col in [(wt[i], BLUE), (mut[i], ROSE)]:
        ax.scatter(0, yy, s=18, facecolors='white', edgecolors=col, lw=1, zorder=5)
    lo = min(wt.min(), mut.min())
    hi = max(wt.max(), mut.max())
    pad = max((hi - lo) * 0.15, 0.25)
    ax.set_ylim(lo - pad, hi + pad)
    ax.yaxis.set_major_locator(MaxNLocator(3))
    ax.set_ylabel('Final logit', labelpad=3)
    ax.tick_params(labelbottom=False)
    if letter == 'd':
        fig.legend(*ax.get_legend_handles_labels(), loc='upper left', bbox_to_anchor=(0.22, 0.99), ncol=2, frameon=False, fontsize=6.8, handlelength=1.2, handletextpad=0.4, columnspacing=0.9, borderaxespad=0)
    da.axhline(0, color='#A6ADB1', lw=0.6)
    da.plot(x, dz, color=PURPLE, lw=1.1)
    da.fill_between(x, 0, dz, color=PURPLE, alpha=0.16, lw=0)
    da.scatter(0, dz[i], s=21, facecolor=PURPLE if dz[i] * r.expected_sign > 0 else 'white', edgecolor=PURPLE, lw=0.85, zorder=5)
    lim = max(abs(dz).max() * 1.22, 0.15)
    da.set_ylim(-lim, lim)
    da.yaxis.set_major_locator(MaxNLocator(3))
    da.set_ylabel('Δlogit', labelpad=3)
    da.set_xlabel('Window centre relative\nto mutation (aa)', fontsize=6.9, labelpad=2)
    da.set_xticks(sorted(set([int(x.min()), 0, int(x.max())])))
    for a in [ax, da]:
        a.axvline(0, color='#9DA6AB', lw=0.65, ls=(0, (3, 2)), zorder=0)
        a.set_xlim(x.min() - 0.4, x.max() + 0.4)
        clean(a)
    save(fig, letter)

def save(fig, letter):
    export(fig, 5, letter)
for letter in 'ab':
    direction(letter)
branches()
for letter in 'def':
    curve(letter)
