"""
fig1_split.py — 2026-10-13 reframe: Figure 1 keeps only the pipeline (a-d, saved as fig1_pipeline) and the three analysis
choices move to the Supplement (former e-g, relabelled a-c, saved as figS_pitfalls). Panels are drawn exactly as in
fig1_pitfalls.py (2026-10-11), which stays unchanged. Original description below.
  a-c  One CALVADOS-2 snapshot of Q9DUC6 (87 residues) as a C-alpha trace, its contacts, and its largest contact cluster as
       a network with the hubs ringed (pipeline/figS_pipeline_1-3.png, drawn for the deck).
  d    Hub occupancy of every residue over 1,000 snapshots (hub_frames_v3_cal/per_residue.csv), canonical stickers F, W,
       Y, R orange and G, S light blue.
  e    Positional profile: within-protein z of mean betweenness vs relative chain position, |i-j| >= 2
       (200 frames) vs |i-j| >= 3 (1,000 frames); the same 1,198 CALVADOS-2 proteins.
  f    Residue-type "taxonomy": mean z per amino acid, |i-j| >= 2 raw, |i-j| >= 2 after
       removing the within-protein positional trend, and |i-j| >= 3 after the same adjustment.
  g    Sampling: correlation between per-residue mean betweenness from two disjoint random sets of n snapshots
       (100 proteins; frames_reliability.py), with the Spearman-Brown projection to 1,000 snapshots.
Inputs: ../../redo_evidence/hub_frames_sep2/, hub_frames/, hub_frames_v3_cal/, frames_reliability_summary.txt.
  ~/.venvs/bender/bin/python fig1_split.py  ->  fig1_pipeline.pdf/.png, figS_pitfalls.pdf/.png
"""
import csv
import os
import re
from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np

import figstyle as fs

HERE = os.path.dirname(os.path.abspath(__file__))
EV = fs.EV
AAS = 'ACDEFGHIKLMNPQRSTVWY'
C2, C3 = fs.RULE_COLORS['sep2'], fs.RULE_COLORS['sep3']


def load(d):
    by = defaultdict(list)
    for r in csv.DictReader(open(os.path.join(EV, d, 'per_residue.csv'))):
        by[r['UniProt_ID']].append((int(r['pos']), r['aa'], float(r['mean_btw']), float(r['z_posadj'])))
    out = {}
    for u, rows in by.items():
        rows.sort()
        b = np.array([x[2] for x in rows])
        sd = b.std() or 1.0
        out[u] = ([x[1] for x in rows], (b - b.mean()) / sd, np.array([x[3] for x in rows]))
    return out


def main():
    d2, d3 = load('hub_frames_sep2'), load('hub_frames')
    ids = sorted(set(d2) & set(d3))
    fig = plt.figure(figsize=(7.2, 1.3))
    top = fig.add_gridspec(1, 4, width_ratios=[1, 1, 1, 1.2], wspace=0.32, left=0.02, right=0.99, top=0.86, bottom=0.25)
    heads = ('one snapshot', 'contacts < 0.8 nm, |i−j| ≥ 3', 'network, hubs ringed')
    for k in range(3):
        a_ = fig.add_subplot(top[0, k])
        img = plt.imread(os.path.join(HERE, 'pipeline', f'figS_pipeline_{k + 1}.png'))[..., :3]
        ys_, xs_ = np.nonzero(img.min(-1) < 0.97)            # crop the white margin so the drawing fills the panel
        img = img[max(0, ys_.min() - 15):ys_.max() + 15, max(0, xs_.min() - 15):xs_.max() + 15]
        a_.imshow(img)
        a_.set_axis_off()
        fs.panel(a_, 'abc'[k], x=0.0, y=1.02, title=heads[k])
        if k < 3:
            a_.annotate('', xy=(1.13, 0.5), xytext=(1.0, 0.5), xycoords='axes fraction',
                        arrowprops=dict(arrowstyle='->', color='0.45', lw=0.8))
    ad = fig.add_subplot(top[0, 3])
    occ_rows = [r for r in csv.DictReader(open(os.path.join(EV, 'hub_frames_v3_cal', 'per_residue.csv')))
                if r['UniProt_ID'] == 'Q9DUC6']
    occ_rows.sort(key=lambda r: int(r['pos']))
    occ = np.array([float(r['occupancy']) for r in occ_rows])
    cls = ['#d95f02' if r['aa'] in 'FWYR' else '#9ecae1' if r['aa'] in 'GS' else '#aeaeb2' for r in occ_rows]
    ad.bar(np.arange(1, len(occ) + 1), occ, width=1.0, color=cls, lw=0)
    ad.axhline(occ.mean(), color='0.25', lw=0.6, ls=(0, (3, 2)))
    ad.text(2, occ.mean() + 0.01, f'mean {occ.mean():.2f}', fontsize=5.5, color='0.25', va='bottom')
    ad.set_xlim(0.5, len(occ) + 0.5)
    ad.set_xticks([1, 40, len(occ)])
    ad.set_xlabel('residue', labelpad=1)
    ad.set_ylabel('hub occupancy', labelpad=1)
    ad.tick_params(labelsize=6)
    fs.panel(ad, 'd', x=-0.3, y=1.02, title='occupancy, 1,000 snapshots')
    fs.save(fig, 'fig1_pipeline')

    fig = plt.figure(figsize=(7.2, 1.9))
    bot = fig.add_gridspec(1, 3, width_ratios=[1, 1.45, 0.9], wspace=0.38, left=0.075, right=0.99, top=0.9, bottom=0.2)
    ax = [fig.add_subplot(bot[0, k]) for k in range(3)]

    # A: positional profile
    bins = np.linspace(0, 1, 21)
    mids = (bins[:-1] + bins[1:]) / 2
    for d, col, lab in ((d2, C2, '|i−j| ≥ 2'), (d3, C3, '|i−j| ≥ 3')):
        prof = []
        for u in ids:
            z = d[u][1]
            x = np.arange(len(z)) / max(1, len(z) - 1)
            k = np.clip(np.digitize(x, bins) - 1, 0, 19)
            prof.append([z[k == j].mean() if (k == j).any() else np.nan for j in range(20)])
        prof = np.array(prof)
        m, se = np.nanmean(prof, 0), np.nanstd(prof, 0) / np.sqrt(np.sum(~np.isnan(prof), 0))
        ax[0].plot(mids, m, color=col, lw=1.4, label=lab)
        ax[0].fill_between(mids, m - 1.96 * se, m + 1.96 * se, color=col, alpha=0.2, lw=0)
    ax[0].axhline(0, color='0.6', lw=0.5)
    ax[0].set_xlabel('relative position in chain')
    ax[0].set_ylabel('z of mean betweenness')
    fs.panel(ax[0], 'a', x=-0.3)
    ax[0].legend(frameon=False, loc='lower center')

    # B: residue-type means
    def per_aa(d, col):
        acc = defaultdict(list)
        for u in ids:
            seq, z, zp = d[u]
            v = z if col == 'raw' else zp
            for a, x in zip(seq, v):
                acc[a].append(x)
        return {a: float(np.mean(acc[a])) for a in AAS}
    s2raw, s2adj, s3adj = per_aa(d2, 'raw'), per_aa(d2, 'adj'), per_aa(d3, 'adj')
    order = sorted(AAS, key=lambda a: s3adj[a])
    x = np.arange(len(order))
    for a_i, a in enumerate(order):
        ax[1].plot([a_i, a_i], [s2adj[a], s3adj[a]], color='0.8', lw=0.8, zorder=1)
    ax[1].scatter(x, [s2raw[a] for a in order], s=14, marker='o', facecolor='white', edgecolor=C2, lw=0.9,
                  label='|i−j| ≥ 2, raw', zorder=3)
    ax[1].scatter(x, [s2adj[a] for a in order], s=12, color=C2, label='|i−j| ≥ 2, position removed', zorder=3)
    ax[1].scatter(x, [s3adj[a] for a in order], s=12, color=C3, label='|i−j| ≥ 3, position removed', zorder=3)
    ax[1].axhline(0, color='0.6', lw=0.5)
    ax[1].set_xticks(x)
    ax[1].set_xticklabels(order)
    fs.color_ticks(ax[1])
    ax[1].set_ylabel('mean z per residue type')
    fs.panel(ax[1], 'b', x=-0.2)
    ax[1].legend(frameon=False, loc='upper left', handletextpad=0.2)

    # C: reliability vs frames. The split-half column correlates the means of two DISJOINT random sets of n snapshots,
    # which is the reliability of an n-snapshot mean. (The other column of the summary correlates a subset with the
    # 1,000-snapshot mean that contains it, a part-whole correlation that overstates reliability; not plotted.)
    txt = open(os.path.join(EV, 'frames_reliability_summary.txt')).read()
    rows = [tuple(map(float, m)) for m in re.findall(r'^\s+(\d+)\s+([\d.]+)\s+([\d.]+)\s*$', txt, re.M)]
    fr, sh, _ = zip(*rows)
    ax[2].plot(fr, sh, 'o-', color=C3, ms=3, lw=1.2, label='two disjoint sets')
    r1000 = 2 * sh[-1] / (1 + sh[-1])                  # Spearman-Brown projection from 500 to 1,000 snapshots
    ax[2].plot([fr[-1], 1000], [sh[-1], r1000], ls='--', color=C3, lw=1.0)
    ax[2].scatter([1000], [r1000], s=14, facecolor='white', edgecolor=C3, zorder=4, label=f'projected, 1,000 ({r1000:.2f})')
    ax[2].set_xscale('log')
    ax[2].set_xticks([25, 50, 100, 200, 500, 1000])
    ax[2].set_xticklabels(['25', '50', '100', '200', '500', '1,000'])
    ax[2].set_ylim(0, 1.2)
    ax[2].set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax[2].set_xlabel('snapshots per protein ($n$)')
    ax[2].set_ylabel('median reliability')
    fs.panel(ax[2], 'c', x=-0.32)
    ax[2].legend(frameon=False, loc='upper left', borderaxespad=0.1)

    fs.save(fig, 'figS_pitfalls')
    print(f'{len(ids)} proteins; contrasts (|i-j|>=2 raw): E {s2raw["E"]:+.3f} Y {s2raw["Y"]:+.3f} M {s2raw["M"]:+.3f}; '
          f'(|i-j|>=3 adj): E {s3adj["E"]:+.3f} Y {s3adj["Y"]:+.3f} G {s3adj["G"]:+.3f}')


if __name__ == '__main__':
    main()
