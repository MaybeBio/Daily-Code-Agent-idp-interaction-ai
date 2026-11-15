from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, to_rgba
from support import DATA as INPUT, save as export
DATA = INPUT / 'fig03'
SEEDS = list(range(2026, 2031))
plt.rcParams.update({'font.family': 'Arial', 'font.size': 7.2, 'axes.spines.top': False, 'axes.spines.right': False, 'svg.fonttype': 'none'})

def table(path):
    return pd.read_csv(path, sep='\t')

def read(name):
    return table(DATA / name)

def title(fig, text):
    fig.text(0.045, 0.965, text, ha='left', va='top', fontsize=8)

def save(fig, key):
    w, h = LAYOUT['panels'][key][2:]
    fig.set_size_inches(w / 25.4, h / 25.4)
    export(fig, 3, key)
COL = dict(A='#7864A4', B='#769D99', C='#167F83', CE='#335D73', B1='#B68C5A', B2='#879CA8', B3='#A0A0A0')
ARMS = ['A', 'B1', 'B2', 'B3', 'B', 'C', 'CE']
NAMES = dict(A='MotifAE-only', B='Intrinsic', C='PriorSLiM', CE='+ Ensemble', B1='− Conservation', B2='− Motif', B3='− Chemistry')
INVERSE = dict(zip(['MotifAE-only', 'Intrinsic without conservation', 'Intrinsic without motif', 'Intrinsic without chemistry', 'Intrinsic + MotifAE', 'PriorSLiM', 'PriorSLiM + ensemble'], ARMS))
CM = LinearSegmentedColormap.from_list('muted_teal', ['#F3F5F4', '#BCD6D2', '#649F9A', '#167F83'])
LAYOUT = {'canvas_mm': [180, 234], 'panels': {'a': [3, 7, 92, 58], 'b': [99, 7, 78, 58], 'c': [3, 72, 68, 58], 'd': [77, 72, 100, 58], 'e': [3, 137, 174, 32], 'f': [3, 176, 174, 54]}}

def plot_a():
    d = table(DATA / 'joint_deletion_endpoint.tsv')
    d['arm'] = d.model.map(INVERSE)
    controls = table(DATA / 'joint_deletion_curves.tsv')
    controls = controls[(controls['set'] == 'matched_control') & (controls['mode'] == 'mean_preserving') & (controls.deleted_count == 398)]
    controls = controls.groupby(['seed', 'model']).ap_loss.mean().reset_index()
    controls['arm'] = controls.model.map(INVERSE)
    fig, ax = plt.subplots(figsize=(85 / 25.4, 55 / 25.4), layout='constrained')
    for i, arm in enumerate(ARMS):
        q = d[d.arm == arm].sort_values('seed')
        assert q.seed.tolist() == list(range(2026, 2031))
        v = q.ap_loss.to_numpy()
        offsets = np.linspace(-0.15, 0.15, 5)
        ax.scatter(v, i + offsets, s=10, c=COL[arm], zorder=3)
        ax.errorbar(v.mean(), i, xerr=v.std(ddof=1), fmt='|', ms=10, capsize=2, lw=1.1, color=COL[arm])
        c = controls[controls.arm == arm].sort_values('seed').ap_loss.to_numpy()
        assert len(c) == 5
        ax.scatter(c, i + offsets, s=11, marker='s', facecolors='white', edgecolors=COL[arm], linewidths=0.75, zorder=2)
    ax.set_yticks(range(7), [NAMES[a] for a in ARMS])
    ax.invert_yaxis()
    ax.set_xlim(0, 0.15)
    ax.set_xticks([0, 0.05, 0.1, 0.15])
    ax.set_ylim(6.6, -1.4)
    ax.set_xlabel('Validation AP loss after deletion')
    ax.set_title('Dependence on 398 target features', loc='left')
    ax.plot([], [], marker='o', color='#666666', ls='', ms=3, label='Targets')
    ax.plot([], [], marker='s', mfc='white', mec='#666666', ls='', ms=3, label='Matched controls')
    ax.legend(frameon=False, loc='upper right', ncol=2, handletextpad=0.3, columnspacing=0.7)
    save(fig, 'a')

def b_retention():
    d = read('retention_near600_display.tsv')
    arms = ['B1', 'B2', 'B3', 'B', 'C', 'CE']
    d['retained_count_difference'] = d.retained - d.A_retained
    q = d.pivot(index='arm', columns='seed', values='retained_count_difference').loc[arms, SEEDS]
    check = d.set_index(['arm', 'seed'])
    for arm in arms:
        for seed in SEEDS:
            r = check.loc[arm, seed]
            assert np.isclose(r.retention_difference * r.set_size, r.retained - r.A_retained, atol=1e-10)
    assert not q.isna().any().any()
    fig = plt.figure(figsize=(85 / 25.4, 60 / 25.4))
    ax = fig.add_axes([0.31, 0.24, 0.64, 0.58])
    assert q.to_numpy().max() <= 0
    im = ax.imshow(-q.to_numpy(), cmap=CM, vmin=0, vmax=55, aspect='auto')
    for y in range(6):
        for x in range(5):
            v = int(q.iloc[y, x])
            ax.text(x, y, str(v).replace('-', '−'), ha='center', va='center', fontsize=7, color='white' if v <= -34 else '#26383C')
    ax.set_xticks(range(5), [str(s) for s in SEEDS])
    ax.set_xlabel('Seed', labelpad=2)
    ax.set_yticks(range(6), ['− Cons.', '− Motif', '− Chem.', 'Intrinsic', 'PriorSLiM', '+ Ensemble'])
    ax.set_xticks(np.arange(6) - 0.5, minor=True)
    ax.set_yticks(np.arange(7) - 0.5, minor=True)
    ax.grid(which='minor', color='white', lw=1.3)
    ax.tick_params(which='both', length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    cax = fig.add_axes([0.31, 0.095, 0.64, 0.025])
    cb = fig.colorbar(im, cax=cax, orientation='horizontal', ticks=[0, 25, 50])
    cb.outline.set_visible(False)
    cb.ax.tick_params(labelsize=6, length=2, pad=1)
    cb.set_label('Fewer retained targets vs MotifAE-only', fontsize=6.3, labelpad=1)
    title(fig, 'Retention at matched observed $K$')
    save(fig, 'b')

def c_association():
    q = pd.read_csv(DATA / 'c_signed_contribution_composition.tsv', sep='\t')
    pairs = [(39522, 'Y'), (39522, 'F'), (8196, 'S'), (14927, 'S')]
    fig = plt.figure()
    ax = fig.add_axes([0.35, 0.29, 0.61, 0.56])
    for i, (fid, aa) in enumerate(pairs):
        part = q[(q.feature_id == fid) & (q.amino_acid == aa)]
        assert len(part) == 2
        vals = [part[part.phase.str.startswith(t)].partial_r.item() for t in ['Discovery', 'Follow-up']]
        ax.plot(vals, [i - 0.07, i + 0.07], color='#B6C0C5', lw=0.8)
        ax.scatter(vals[0], i - 0.07, color='#879CA8', s=17, marker='o', label='Discovery, folds 1–5' if i == 0 else None)
        ax.scatter(vals[1], i + 0.07, color='#167F83', s=16, marker='D', label='Follow-up, folds 6–10' if i == 0 else None)
    ax.set_xlim(-0.38, 0.3)
    ax.set_ylim(3.6, -0.6)
    ax.axvline(0, color='#B0B8BC', lw=0.65)
    ax.set_yticks(range(4), [f'ID{fid} · {aa}' for fid, aa in pairs])
    ax.set_xticks([-0.3, 0, 0.2])
    ax.set_xlabel('Partial correlation', labelpad=3)
    fig.text(0.045, 0.965, 'Composition-linked evidence', va='top', fontsize=8)
    hs, ls = ax.get_legend_handles_labels()
    fig.legend(hs, ls, loc='lower center', bbox_to_anchor=(0.55, 0.03), frameon=False, ncol=1, fontsize=6.2, handletextpad=0.4, labelspacing=0.3)
    save(fig, 'c')

def d_case():
    q = pd.read_csv(DATA / 'd_original_absolute_values.tsv', sep='\t')
    arms = ['C', 'H3', 'AROMATIC_ADD', 'AROMATIC_SPLIT']
    labels = ['PriorSLiM', '+ SYF density', '+ FYW detail', 'Split aromatic']
    colors = ['#167F83', '#BE8345', '#879CA8', '#7864A4']
    fig = plt.figure()
    axes = [fig.add_axes([0.235, 0.25, 0.3, 0.48]), fig.add_axes([0.645, 0.25, 0.3, 0.48])]
    metrics = [('mean_abs_contribution', 1.0, 'Contribution', 'Mean |logit contribution|', 0.02, [0, 0.01, 0.02]), ('deletion_bce_increase', 10000.0, 'Deletion effect', 'BCE increase ($\\times10^{-4}$)', 8.0, [0, 4, 8])]
    records = []
    for ax, (col, scale, head, xlabel, xmax, ticks) in zip(axes, metrics):
        for i, (arm, color) in enumerate(zip(arms, colors)):
            t = q[q.arm == arm].sort_values('seed')
            assert t.seed.tolist() == SEEDS
            vals = t[col].to_numpy() * scale
            ax.barh(i, vals.mean(), height=0.52, color=to_rgba(color, 0.48), edgecolor=color, lw=0.6, zorder=1)
            ax.errorbar(vals.mean(), i, xerr=vals.std(ddof=1), fmt='none', ecolor='#354047', lw=0.8, capsize=2, zorder=3)
            ax.scatter(vals, i + np.linspace(-0.14, 0.14, 5), s=8, facecolors='white', edgecolors='#455258', linewidths=0.5, zorder=4)
            records.append({'arm': arm, 'metric': col, 'mean': float(vals.mean() / scale), 'sd': float(vals.std(ddof=1) / scale), 'n': 5})
        ax.set_ylim(3.6, -0.6)
        ax.set_xlim(0, xmax)
        ax.set_xticks(ticks)
        if col == 'mean_abs_contribution':
            ax.set_xticklabels(['0', '.01', '.02'])
        ax.set_title(head, fontsize=7.2, pad=5)
        ax.set_xlabel(xlabel, fontsize=6.2, labelpad=3)
        ax.set_yticks(range(4), labels if ax is axes[0] else [''] * 4)
        ax.tick_params(axis='y', length=0, pad=3)
        ax.spines['left'].set_visible(False)
    fig.text(0.04, 0.965, 'ID39522 after explicit chemistry changes', va='top', fontsize=8)
    fig.text(0.235, 0.075, 'Bars: mean ± s.d.   Circles: five seeds', fontsize=6.4, color='#4D5B60')
    save(fig, 'd')

def e_training(cmap):
    q = pd.read_csv(DATA / 'e_training_strategy_seed_values.tsv', sep='\t')
    arms = ['J_SCRATCH', 'S', 'C']
    roles = ['full', 'performance', 'efficiency']
    stats = q.groupby(['arm', 'role']).ap.agg(['mean', 'std', 'size'])
    assert (stats['size'] == 5).all()
    vals = np.array([[stats.loc[(a, r), 'mean'] for r in roles] for a in arms])
    fig = plt.figure()
    ax = fig.add_axes([0.15, 0.13, 0.73, 0.53])
    im = ax.imshow(vals, cmap=cmap, vmin=0.51, vmax=0.537, aspect='auto')
    for i, a in enumerate(arms):
        for j, r in enumerate(roles):
            m, s = stats.loc[(a, r), ['mean', 'std']]
            ax.text(j, i, f'{m:.4f} ± {s:.4f}', ha='center', va='center', fontsize=7.3, color='white' if m > 0.528 else '#26383C')
    ax.set_xticks(range(3), ['Full readout', 'Performance refined', 'Efficiency refined'])
    ax.xaxis.tick_top()
    ax.tick_params(axis='x', pad=3)
    ax.set_yticks(range(3), ['Joint scratch', 'Sequential', 'PriorSLiM'])
    ax.set_xticks(np.arange(4) - 0.5, minor=True)
    ax.set_yticks(np.arange(4) - 0.5, minor=True)
    ax.tick_params(which='both', length=0)
    ax.grid(which='minor', color='white', lw=1.3)
    for s in ax.spines.values():
        s.set_visible(False)
    ca = fig.add_axes([0.91, 0.13, 0.012, 0.53])
    cb = fig.colorbar(im, cax=ca, ticks=[0.51, 0.52, 0.53])
    cb.outline.set_visible(False)
    cb.ax.tick_params(length=2, labelsize=6.3, pad=2)
    fig.text(0.025, 0.965, 'Training strategy', fontsize=8, va='top')
    fig.text(0.99, 0.965, 'Validation AP, five-seed mean ± s.d.', fontsize=7, va='top', ha='right', color='#48575B')
    save(fig, 'e')

def f_curves():
    q = pd.read_csv(DATA / 'f_budget_envelope.tsv', sep='\t')
    styles = {'proximal': ('#007F73', 'o', '-', 'Proximal'), 'binary': ('#7B3294', 's', (0, (5, 2)), 'Positional binary gate'), 'topk': ('#D55E00', '^', (0, (3, 1, 1, 1)), 'Top-K'), 'random': ('#41464C', 'x', (0, (1, 2)), 'Random'), 'decoder': ('#0072B2', 'D', (0, (6, 2, 1, 2)), 'Binary gate + decoder')}
    fig = plt.figure()
    ax = fig.add_axes([0.09, 0.24, 0.65, 0.64])
    summaries = []
    for key, (color, marker, line, label) in styles.items():
        m = q[q.method == key].groupby('budget').val_ap.agg(['mean', 'std', 'size'])
        m = m[m['size'] == 5]
        x = m.index.to_numpy(float)
        y = m['mean'].to_numpy()
        sd = m['std'].to_numpy()
        ax.fill_between(x, y - sd, y + sd, color=color, alpha=0.055, lw=0, zorder=1)
        ax.plot(x, y, color=color, marker=marker, ls=line, lw=1.65, ms=3.7, markerfacecolor='white' if key != 'proximal' else color, markeredgewidth=0.85, label=label, zorder=3 if key == 'proximal' else 2)
        summaries.extend((dict(method=key, K=int(k), mean_ap=float(r['mean']), sd=float(r['std']), n=int(r['size'])) for k, r in m.iterrows()))
    ax.set_xscale('log', base=2)
    ax.set_xticks([64, 128, 256, 512, 1024, 2048, 4096, 5829], ['64', '128', '256', '512', '1,024', '2,048', '4,096', '5,829'])
    ax.set_xlim(55, 6600)
    ax.set_ylim(0.39, 0.56)
    ax.set_yticks([0.4, 0.45, 0.5, 0.55])
    ax.set_xlabel('Maximum retained features, $K$', labelpad=4)
    ax.set_ylabel('Validation AP', labelpad=4)
    ax.legend(loc='center left', bbox_to_anchor=(1.015, 0.5), frameon=False, fontsize=6.8, handlelength=2.8, handletextpad=0.65, labelspacing=1.25, borderaxespad=0)
    fig.text(0.025, 0.965, 'Selection and readout controls', va='top', fontsize=8)
    save(fig, 'f')
plot_a()
b_retention()
c_association()
d_case()
e_training(CM)
f_curves()
