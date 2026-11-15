from pathlib import Path
import os, sys, json, hashlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.text import Text, Annotation
from matplotlib.lines import Line2D
from matplotlib.colors import ListedColormap, to_rgb
from matplotlib.patches import Rectangle
ROOT = Path(__file__).resolve().parent
PACKAGE = ROOT.parent
sys.path.insert(0, str(PACKAGE.parent))
import fig02_data as src
import fig02_detail as detail
from support import OUT, save
USED = {}

class Context:

    def table(self, name):
        p = PACKAGE / name
        USED[name] = hashlib.sha256(p.read_bytes()).hexdigest()
        return pd.read_csv(p, sep='\t', float_precision='round_trip')

    def json(self, name):
        p = PACKAGE / name
        USED[name] = hashlib.sha256(p.read_bytes()).hexdigest()
        return json.loads(p.read_text())
CTX = Context()
LAYOUT = {'a': [0.18, 0.25, 3.8, 5.0], 'b': [4.16, 0.25, 4.4, 5.0], 'c': [8.74, 0.25, 6.08, 5.0], 'd': [0.18, 5.55, 5.6, 6.1], 'e': [5.96, 5.55, 8.86, 6.1], 'f': [0.18, 11.95, 5.6, 6.4], 'g': [5.96, 11.95, 8.86, 6.4]}
plt.rcParams.update(src.RC)
plt.rcParams.update({'font.size': 16, 'axes.labelsize': 16, 'axes.titlesize': 18, 'xtick.labelsize': 16, 'ytick.labelsize': 16, 'svg.fonttype': 'none', 'axes.labelpad': 7, 'savefig.facecolor': 'white'})
detail.DIV = detail.LinearSegmentedColormap.from_list('evidence', [detail.BLUE, '#F8F8F5', detail.ORANGE])

def canvas(k):
    w, h = LAYOUT[k][2:]
    return plt.figure(figsize=(w, h), facecolor='white')

def ax(fig, x, y, w, h):
    fw, fh = fig.get_size_inches()
    return fig.add_axes([x / fw, 1 - (y + h) / fh, w / fw, h / fh])

def title(fig, text, x, y=0.38):
    _, h = fig.get_size_inches()
    fig.text(x / fig.get_figwidth(), 1 - y / h, text, ha='center', va='top', fontsize=18)

def common(fig):
    for a in fig.axes:
        a.tick_params(width=1.4, pad=4)
        for s in a.spines.values():
            s.set_linewidth(1.4)
        for collection in a.collections:
            collection.set_rasterized(False)
        for text in a.texts:
            if isinstance(text, Annotation):
                text.set_position(tuple((v * 1.7 for v in text.get_position())))
    for t in fig.findobj(Text):
        t.set_fontfamily('Arial')
        t.set_fontsize(max(16, t.get_fontsize()))
    fig.canvas.draw()

def panel_a():
    fig = src.abc.plot_a(*src.load_abc(CTX))
    fig.set_size_inches(*LAYOUT['a'][2:])
    a = fig.axes[0]
    a.set_position([0.8 / 3.8, 1 - (0.95 + 3.0) / 5, 2.82 / 3.8, 3.0 / 5])
    a.set_title('')
    title(fig, 'Test-set discrimination', 2.0)
    a.legend(frameon=False, fontsize=16, loc='upper right', handlelength=1.1, labelspacing=0.22, borderaxespad=0.08, handletextpad=0.3)
    for line in a.lines:
        line.set_linewidth(2)
    return fig

def panel_b():
    _, _, _, strata, _ = src.load_abc(CTX)
    values = np.array([[*[strata[(strata.method == m) & (strata.negative_type == k)].fpr.mean() * 100 for k in src.abc.TYPES], strata[strata.method == m].recall.mean() * 100] for m in src.abc.METHODS]).T
    fig = canvas('b')
    a = ax(fig, 1.12, 0.95, 3.13, 3.0)
    cells = np.zeros((7, 4, 3))
    for r in range(7):
        best = np.argmin(values[r]) if r < 6 else np.argmax(values[r])
        for c in range(4):
            intensity = 1 - values[r, c] / 100 if r < 6 else values[r, c] / 100
            start = np.array(to_rgb('#f5f7f8'))
            end = np.array(to_rgb('#9abbd0' if r < 6 else '#e0c7a3'))
            cells[r, c] = start + intensity * (end - start)
            a.text(c, r, f'{values[r, c]:.1f}', ha='center', va='center', fontweight='bold' if c == best else 'normal')
        a.add_patch(Rectangle((best - 0.48, r - 0.48), 0.96, 0.96, fill=False, lw=1.4, ec='#3f7093' if r < 6 else '#9b793e'))
    a.imshow(cells, aspect='auto')
    a.axhline(5.5, c='white', lw=2)
    a.set(yticks=range(7), yticklabels=src.abc.TYPE_LABELS, xticks=range(4), xticklabels=src.abc.METHODS)
    a.tick_params(axis='x', rotation=60, length=0)
    a.tick_params(axis='y', length=0)
    for s in a.spines.values():
        s.set_visible(False)
    title(fig, 'False positives ↓ / recall ↑ (%)', 2.3)
    return fig

def panel_c():
    fig = canvas('c')
    a = ax(fig, 2.82, 0.95, 2.65, 3.0)
    cb = ax(fig, 5.86, 0.95, 0.12, 3.0)
    detail.old.draw_e(a, cb, src.load_detail(CTX, 'e'))
    a.images[0].set_cmap(detail.DIV)
    a.set_title('')
    a.set_xticklabels(['E→F\nPos.', 'E→F\nNeg.', 'M→F\nPos.', 'M→F\nNeg.'])
    cb.set_ylabel('')
    cb.yaxis.set_ticks_position('left')
    cb.set_title('SD', fontsize=16, pad=8)
    title(fig, 'Corrected versus still wrong', 3.65)
    for l in a.lines:
        if l.get_marker() == '.':
            l.set_markersize(6)
    return fig

def panel_d():
    fig = canvas('d')
    a = ax(fig, 0.95, 0.95, 4.4, 3.45)
    detail.f_selection(a, src.load_detail(CTX, 'f'))
    a.set_title('', loc='left')
    title(fig, 'Influence across discovery windows', 2.95)
    for c in a.collections:
        c.set_sizes(c.get_sizes() * 5)
    a.legend(handles=[Line2D([], [], marker='o', ls='', color=c, markersize=6, label=l) for c, l in [(detail.BLUE, 'Corrected negatives'), (detail.ORANGE, 'Rescued positives')]], frameon=False, loc='upper left', fontsize=16, handlelength=0.6, handletextpad=0.4, borderpad=0)
    for t in a.texts:
        if isinstance(t, Annotation) and t.get_text() in ['27995', '10467']:
            t.set_ha('right')
            t.set_position((-8, 12 if t.get_text() == '27995' else -20))
    return fig

def panel_e():
    fig = canvas('e')
    hm = ax(fig, 1.45, 0.95, 4.55, 3.45)
    score = ax(fig, 6.65, 0.95, 2.0, 3.45)
    cb = ax(fig, 2.0, 5.12, 3.2, 0.12)
    x = np.linspace(-1, 1, 1025)
    white = np.array(to_rgb('#FAFAF8'))
    end = np.where((x < 0)[:, None], to_rgb('#276D91'), to_rgb('#C4513E'))
    detail.DIV = ListedColormap(white + np.abs(x)[:, None] ** 0.5 * (end - white))
    detail.f_cases(fig, hm, score, cb, src.load_detail(CTX, 'f'))
    hm.set_title('')
    score.set_title('')
    title(fig, 'Follow-up position contributions', 3.6)
    title(fig, 'Latent deletion', 7.65)
    for t in hm.texts:
        t.set_x(-0.27)
    for c in score.collections:
        c.set_sizes(c.get_sizes() * 3.6)
    for l in score.lines:
        l.set_linewidth(2)
    for legend in fig.legends:
        legend.remove()
    fig.legend(handles=[Line2D([], [], marker='o', ls='', color=detail.BLUE, label='Full', ms=7), Line2D([], [], marker='x', ls='', color=detail.INK, label='Deleted', ms=7), Line2D([], [], marker='|', ls='', color=detail.GRAY, label='Threshold', ms=11)], loc='lower right', bbox_to_anchor=(0.995, 0.005), frameon=False, fontsize=16, labelspacing=0.15, handlelength=0.7)
    return fig

def panel_f():
    fig = canvas('f')
    a = ax(fig, 0.95, 0.95, 4.4, 4.0)
    d = src.load_detail(CTX, 'g')
    allg = d['allg']
    a.set_xscale('symlog', linthresh=1e-05)
    a.set_yscale('symlog', linthresh=1e-05)
    limit = max(abs(allg[['M_center_bce', 'B_center_bce']].to_numpy(float)).max(), 0.001) * 1.1
    a.add_patch(Rectangle((0.0001, -1e-05), limit - 0.0001, 2e-05, facecolor='#F4E6CF', edgecolor='#B17C3E', lw=1.1, zorder=0))
    a.scatter(allg.M_center_bce, allg.B_center_bce, s=12, c=detail.BLUE, alpha=0.32, linewidth=0, zorder=2)
    a.plot([-limit, limit], [-limit, limit], c='#777777', ls='--', lw=1.2, zorder=1)
    a.axhline(0, c='#C8CDD1', lw=1, zorder=1)
    a.axvline(0, c='#C8CDD1', lw=1, zorder=1)
    a.set(xlabel='MotifAE-only deletion ΔBCE', ylabel='Joint deletion ΔBCE')
    a.set_xticks([0, 0.0001, 0.01])
    a.set_yticks([0, 0.0001, 0.01])
    a.set_xticklabels(['0', '1e-4', '1e-2'])
    a.set_yticklabels(['0', '1e-4', '1e-2'])
    a.text(0.02, 0.97, '5,880 candidates', ha='left', va='top', transform=a.transAxes)
    row = d['broad'][d['broad'].feature_id == 27779].iloc[0]
    a.scatter([row.M_center_bce], [row.B_center_bce], s=62, facecolors='none', edgecolors='#333B43', lw=1.2, zorder=4)
    a.annotate('27779', (row.M_center_bce, row.B_center_bce), xytext=(10, -30), textcoords='offset points', color='#333B43', arrowprops=dict(arrowstyle='-', lw=0.8, color='#69747E'))
    detail.style(a)
    fig.text(0.95 / 5.6, 1 - 5.85 / 6.4, 'Shaded: MotifAE-only ΔBCE ≥ 1e-4\n|Joint ΔBCE| ≤ 1e-5', ha='left', va='top', fontsize=16, color='#805B2C')
    title(fig, 'Dependence after fusion', 2.95)
    region = (allg.M_center_bce >= 0.0001) & (allg.B_center_bce.abs() <= 1e-05)
    return fig

def panel_g():
    fig = canvas('g')
    a = ax(fig, 1.15, 0.95, 3.92, 4.0)
    b = ax(fig, 6.6, 0.95, 2.05, 4.0)
    d = src.load_detail(CTX, 'g')
    broad = d['broad']
    g = broad[broad.M_center_bce >= 1e-05]
    colors = ['#B77A35', '#C3CDD3', '#A9BDC9', '#7C9FB3', '#507F9C', '#185D82']
    for k in [1, 2, 3, 4, 5, 0]:
        rows = g[g.retained_D_count == k]
        kwargs = dict(facecolors='none', edgecolors=colors[k], s=42, linewidth=1.25, alpha=1, zorder=4) if k == 0 else dict(c=colors[k], s=25 if k == 5 else 19, linewidth=0.45, edgecolors='white', alpha=1 if k == 5 else 0.5, zorder=3 if k == 5 else 2)
        a.scatter(100 * rows.train_active_fraction, 100 * rows.relative_dependency_drop, **kwargs)
    a.set_yscale('symlog', linthresh=20)
    a.axhline(80, c='#B89561', ls='--', lw=1, zorder=0)
    a.axhline(0, c=detail.GRAY, lw=1, zorder=0)
    labels = broad[broad.near_zero_stable & (broad.retained_D_count <= 1)].sort_values('train_active_fraction', ascending=False).head(4)
    offsets = {37849: (5, 15), 27779: (4, 4), 13621: (-10, 18), 40054: (-3, -20)}
    for r in labels.itertuples():
        xy = (100 * r.train_active_fraction, 100 * r.relative_dependency_drop)
        a.scatter(*xy, facecolor='none', edgecolor=detail.INK, s=100, lw=0.9, zorder=5)
        a.annotate(str(r.feature_id), xy, xytext=offsets[r.feature_id], textcoords='offset points', color=colors[0] if r.feature_id == 27779 else detail.INK, fontweight='bold' if r.feature_id == 27779 else 'normal', arrowprops=dict(arrowstyle='-', color=detail.GRAY, lw=0.7))
    a.set(xlabel='TRAIN activation coverage (%)', ylabel='Deletion ΔBCE decrease (%)', xticks=[0, 25, 50, 75])
    detail.style(a)
    handles = [Line2D([], [], ls='', marker='o', ms=7, markerfacecolor='none' if k == 0 else colors[k], markeredgecolor=colors[k], alpha=1 if k in [0, 5] else 0.5, label=str(k)) for k in range(6)]
    fig.legend(handles=handles, ncol=6, loc='upper center', bbox_to_anchor=(3.11 / 8.86, 1 - 5.62 / 6.4), frameon=False, title='Stage-D retention (folds out of 5)', title_fontsize=16, fontsize=16, columnspacing=0.65, handletextpad=0.3, handlelength=0.8, borderpad=0, labelspacing=0.3)
    detail.g_example(b, d)
    a.set_title('')
    b.set_title('')
    title(fig, 'Coverage, dependence\nand retention', 3.04, 0.1)
    title(fig, 'ID27779\nfive follow-up folds', 7.48, 0.1)
    for c in b.collections:
        c.set_sizes(c.get_sizes() * 4)
    for l in b.lines:
        l.set_linewidth(2)
    for t in b.texts:
        if 'coverage' in t.get_text():
            t.set_text('Retained: 0/5')
    return fig
for letter in 'abcdefg':
    fig = globals()['panel_' + letter]()
    common(fig)
    save(fig, 2, letter)
