from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.text import Text
from matplotlib.colors import LinearSegmentedColormap, Normalize, ListedColormap, to_rgb, to_rgba
import fig02_evidence as old
BLUE = '#276D91'
TEAL = '#679995'
ORANGE = '#BE8345'
GRAY = '#A6ADB1'
INK = '#27323A'
DIV = LinearSegmentedColormap.from_list('evidence', [BLUE, '#F8F8F5', ORANGE])
SEQ = LinearSegmentedColormap.from_list('intensity', ['#F4F7F8', '#96BECE', BLUE])
RET = LinearSegmentedColormap.from_list('retention', [ORANGE, '#E7DED1', '#A3BDCA', BLUE])
PLACEMENT = {'a': [3, 3, 65, 54], 'b': [71, 3, 62, 54], 'c': [136, 3, 41, 54], 'd': [3, 60, 91, 53], 'e': [97, 60, 80, 53], 'f': [3, 118, 174, 62], 'g': [3, 184, 174, 62]}

def style(ax):
    ax.spines[['top', 'right']].set_visible(False)
    ax.tick_params(length=2, pad=2)

def f_selection(ax, d):
    colors = {'FP_corrected': BLUE, 'FN_rescued': ORANGE}
    names = {'FP_corrected': 'Corrected negatives', 'FN_rescued': 'Rescued positives'}
    for kind, g in d['discovery'].groupby('kind'):
        ax.scatter(g.mean_bce, g.proteins, c=colors[kind], s=7, alpha=0.22, lw=0)
    offsets = {5719: (-12, 7), 27995: (-9, 9), 35819: (-34, -7), 27176: (-37, -2), 10467: (18, 0), 36920: (-55, 13)}
    for r in d['locked']['f']:
        fid = int(r['feature_id'])
        co = colors[r['kind']]
        ax.scatter(r['mean_bce'], r['proteins'], s=23 if fid in [5719, 27995] else 15, c=co, edgecolor=INK, lw=0.6, zorder=5)
        ax.annotate(str(fid), (r['mean_bce'], r['proteins']), xytext=offsets[fid], textcoords='offset points', fontsize=7, fontweight='bold' if fid in [5719, 27995] else 'normal', color=co, arrowprops=dict(arrowstyle='-', color=co, lw=0.6))
    ax.set_xscale('symlog', linthresh=0.0001)
    ax.set_xticks([-0.01, 0, 0.001, 0.1])
    ax.set_xticklabels(['-0.01', '0', '0.001', '0.1'])
    ax.set(xlabel='Deletion ΔBCE', ylabel='Proteins with |Δp| ≥ 0.05', ylim=(-4, 158))
    ax.set_title('Influence across discovery windows', loc='left', fontsize=8, pad=8)
    ax.axvline(0, c=GRAY, lw=1)
    ax.legend(handles=[Line2D([], [], marker='o', ls='', color=colors[k], markersize=3, label=names[k]) for k in ['FP_corrected', 'FN_rescued']], fontsize=7, frameon=False, loc='upper left', handlelength=0.7, handletextpad=0.4, borderpad=0, labelspacing=0.2)
    style(ax)

def f_cases(fig, hm, score, cb, d):
    cases = d['fcases']
    tr = d['ftracks']
    arr = np.array([tr[(tr.window_id == r.window_id) & (tr.feature_id == r.feature_id) & (tr.kind == r.kind)].sort_values('position').contribution.to_numpy() for r in cases.itertuples()])
    assert arr.shape == (8, 29) and np.isfinite(arr).all()
    bound = max(0.01, float(abs(arr).max()))
    im = hm.imshow(arr, cmap=DIV, vmin=-bound, vmax=bound, aspect='auto')
    hm.set(xticks=[0, 14, 28], xticklabels=['-14', '0', '+14'], xlabel='Window position', yticks=range(8), yticklabels=[r.base_protein_id for r in cases.itertuples()])
    hm.tick_params(axis='y', length=0, pad=2)
    hm.axvline(14, c=GRAY, ls=':', lw=1)
    hm.axhline(3.5, c='white', lw=2)
    hm.set_title('Follow-up position contributions', fontsize=8, pad=8)
    hm.text(-0.3, 0.76, '5719', transform=hm.transAxes, rotation=90, va='center', ha='center', color=BLUE, fontweight='bold')
    hm.text(-0.3, 0.25, '27995', transform=hm.transAxes, rotation=90, va='center', ha='center', color=ORANGE, fontweight='bold')
    for n, r in enumerate(cases.itertuples()):
        co = BLUE if r.kind == 'FP_corrected' else ORANGE
        score.plot([r.probability, r.without_probability], [n, n], c=co, lw=1)
        score.scatter(r.probability, n, c=co, s=14, zorder=4)
        score.scatter(r.without_probability, n, c=INK, marker='x', s=14, lw=1)
        score.scatter(r.threshold, n, c=GRAY, marker='|', s=45, lw=1)
    score.set(ylim=(7.5, -0.5), xlim=(0, 1), xticks=[0, 0.5, 1], xticklabels=['0', '0.5', '1'], yticks=[], xlabel='Probability')
    score.set_title('Latent deletion', fontsize=8, pad=8)
    style(score)
    colorbar = fig.colorbar(im, cax=cb, orientation='horizontal', ticks=[-bound, 0, bound])
    colorbar.ax.set_xticklabels([f'{-bound:.2f}', '0', f'{bound:.2f}'])
    colorbar.set_label('Signed logit contribution', labelpad=1, fontsize=7)
    colorbar.ax.tick_params(length=1, pad=1)
    fig.legend(handles=[Line2D([], [], marker='o', ls='', color=BLUE, label='Full', ms=3), Line2D([], [], marker='x', ls='', color=INK, label='Deleted', ms=3), Line2D([], [], marker='|', ls='', color=GRAY, label='Threshold', ms=5)], loc='lower right', bbox_to_anchor=(0.995, 0.02), ncol=1, frameon=False, fontsize=7, borderpad=0, handlelength=0.7, labelspacing=0.1, columnspacing=1)

def g_coverage(fig, ax, cb, d):
    broad = d['broad']
    g = broad[broad.M_center_bce >= 1e-05]
    sc = ax.scatter(100 * g.train_active_fraction, 100 * g.relative_dependency_drop, c=g.retained_D_count, cmap=RET, vmin=0, vmax=5, s=8, alpha=0.68, lw=0)
    ax.set_yscale('symlog', linthresh=20)
    ax.axhline(80, c=ORANGE, ls='--', lw=1)
    ax.axhline(0, c=GRAY, lw=1)
    label_rows = broad[broad.near_zero_stable & (broad.retained_D_count <= 1)].sort_values('train_active_fraction', ascending=False).head(4)
    offsets = {37849: (5, 15), 27779: (4, 4), 13621: (-10, 18), 40054: (-3, -20)}
    for r in label_rows.itertuples():
        x = 100 * r.train_active_fraction
        y = 100 * r.relative_dependency_drop
        ax.scatter(x, y, facecolor='none', edgecolor=INK, s=30, lw=0.8, zorder=5)
        ax.annotate(str(r.feature_id), (x, y), xytext=offsets.get(r.feature_id, (4, 8)), textcoords='offset points', color=ORANGE if r.feature_id == 27779 else INK, fontweight='bold' if r.feature_id == 27779 else 'normal', fontsize=7, arrowprops=dict(arrowstyle='-', color=GRAY, lw=0.6))
    ax.set(xlabel='TRAIN activation coverage (%)', ylabel='Deletion ΔBCE decrease (%)')
    ax.set_title('Coverage, dependence\nand retention', fontsize=8, pad=7)
    ax.set_xticks([0, 25, 50, 75])
    style(ax)
    bar = fig.colorbar(sc, cax=cb, orientation='horizontal', ticks=range(6))
    bar.set_label('Stage-D retention (folds out of 5)', fontsize=7, labelpad=1)
    bar.ax.tick_params(length=1, pad=1)

def g_example(ax, d):
    g = d['gfolds'].query('feature_id == 27779 and fold >= 6').sort_values('fold')
    assert len(g) == 5 and int(g.retained_D.sum()) == 0
    for r in g.itertuples():
        ax.plot([0, 1], [r.M_center_bce * 10000.0, r.B_center_bce * 10000.0], c=GRAY, alpha=0.75, lw=1)
    for i, (col, co) in enumerate([('M_center_bce', ORANGE), ('B_center_bce', BLUE)]):
        ax.scatter(np.full(len(g), i), g[col] * 10000.0, c=co, s=14, zorder=3)
        ax.plot([i - 0.12, i + 0.12], [g[col].mean() * 10000.0] * 2, c=INK, lw=1.5, zorder=4)
    ax.axhline(0, c=GRAY, ls='--', lw=1)
    ax.set(xticks=[0, 1], xticklabels=['MotifAE\nonly', 'Joint'], ylabel='Deletion ΔBCE (× 1e-4)', xlim=(-0.35, 1.35), ylim=(-0.55, 2.8), yticks=[0, 1, 2])
    ax.set_title('ID27779\nfive follow-up folds', fontsize=8, pad=7)
    ax.text(0.5, 0.96, '52% TRAIN coverage\nStage-D retained: 0/5', transform=ax.transAxes, fontsize=7, ha='center', va='top')
    style(ax)

def build(letter, d):
    _, _, w, h = PLACEMENT[letter]
    fig = plt.figure(figsize=(w / 25.4, h / 25.4), facecolor='white')

    def ax(x, y, ww, hh):
        return fig.add_axes([x / w, 1 - (y + hh) / h, ww / w, hh / h])
    if letter == 'd':
        axs = [ax(15, 12, 19, 29), ax(43, 12, 20, 29), ax(71, 12, 19, 29)]
        old.draw_d(axs, d)
        for a, title in zip(axs, ['AP gain', 'Positives', 'Negatives']):
            a.set_title(title, fontsize=8, pad=6)
        axs[0].set_xlabel('AP gain (pp)')
        axs[0].set_xticks([0, 5, 10, 15])
        axs[1].set(xlabel='% positives', xlim=(-15, 25), xticks=[-10, 0, 10, 20])
        axs[2].set(xlabel='% negatives', xlim=(-5, 15), xticks=[0, 10])
        fig.legend(handles=[Line2D([], [], color=BLUE, lw=4, label='Corrected'), Line2D([], [], color=ORANGE, lw=4, label='New errors')], loc='upper right', bbox_to_anchor=(1, 1.02), ncol=2, frameon=False, fontsize=7, handlelength=0.7, columnspacing=0.8)
    elif letter == 'e':
        a = ax(34, 9, 33, 33)
        c = ax(75, 9, 1.5, 33)
        old.draw_e(a, c, d)
        a.images[0].set_cmap(DIV)
        a.set_title('Corrected versus still wrong', fontsize=8, pad=6)
        a.set_xticklabels(['E to F\nPos.', 'E to F\nNeg.', 'M to F\nPos.', 'M to F\nNeg.'])
        c.set_ylabel('')
        c.yaxis.set_ticks_position('left')
        c.set_title('SD', fontsize=7, pad=3)
    elif letter == 'f':
        s = ax(15, 12, 52, 31)
        f_selection(s, d)
        hm = ax(95, 12, 47, 31)
        score = ax(150, 12, 22, 31)
        cb = ax(104, 53, 29, 1.4)
        f_cases(fig, hm, score, cb, d)
    elif letter == 'g':
        a = ax(15, 12, 40, 31)
        old.draw_g(a, d)
        a.set_title('Dependence after fusion', fontsize=8, pad=7)
        a.set_xticklabels(['0', '1e-4', '1e-2'])
        a.set_yticklabels(['0', '1e-4', '1e-2'])
        a.set_xlabel('MotifAE-only deletion ΔBCE')
        for t in a.texts:
            if 'near-zero' in t.get_text():
                t.set_text('5880 IDs tested\n1459 broad; 26 near-zero')
        b = ax(79, 12, 44, 31)
        cb = ax(85, 53, 32, 1.4)
        g_coverage(fig, b, cb, d)
        c = ax(145, 12, 28, 31)
        g_example(c, d)
    for t in fig.findobj(Text):
        t.set_fontfamily('Arial')
        t.set_fontsize(max(7, t.get_fontsize()))
    for a in fig.axes:
        for c in a.collections:
            c.set_rasterized(False)
        for s in a.spines.values():
            s.set_linewidth(1)
        for line in a.lines:
            line.set_linewidth(max(1, line.get_linewidth()))
        a.tick_params(width=1)
    return fig
