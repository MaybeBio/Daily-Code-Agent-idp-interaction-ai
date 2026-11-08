"""
fig2_readouts.py — Figure 2 (2026-10-14 redesign, replaces fig2_flicker_null.py): the contact network is reproducible,
and degree more so than hub occupancy. Values from hashed results; no new statistics.
  a  the part of hub occupancy that degree does not explain, agreement between CALVADOS-2 and Mpipi-Recharged for the same
     protein against a same-length unrelated protein (hub_vs_degree Part A2; cache/a2_residual.csv)
  b  residue-type enrichment among residues in the top 10% by only one readout (Part A3; cache/a3_enrichment.csv)
  c  hub overlap of successive snapshots against random sets (hub_frames_v3_*/per_protein.csv), median and IQR
  d  split-half reliability of hub occupancy, Spearman-Brown corrected (same files)
  e  within-protein Spearman of the two readouts after removing chain position (Part A1; cache/within_per_protein.csv)
  f  agreement of each readout between CALVADOS-2 and Mpipi-Recharged for the same protein (Part A2;
     cache/crossff_per_protein.csv)
  g  agreement of each readout between independent experimental ensembles of the same protein, 16 pairs
     (hub_vs_degree/part_b.txt)
  ~/.venvs/finches/bin/python fig2_readouts.py  ->  fig2_readouts.pdf / .png
"""
import csv
import os
import re
import sys

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap

import figstyle as fs

sys.path.insert(0, fs.EV)
import hub_frames as hf  # noqa: E402

UID = 'Q9DUC6'
SNAP, BRIDGE, DENSE = 870, 52, 61          # snapshot index and 0-based residues shown in panel a
CACHE = os.path.join(fs.HERE, 'cache')
DEG, OCC = '#2ca25f', '#1b6ca8'
FFS = [('hub_frames_v3_cal', 'CALVADOS-2'), ('hub_frames_v3_mpipi', 'MPIPI-GG'), ('hub_frames_v3_recharged', 'Mpipi-Recharged')]
NAME = {'MPIPI-GG': 'Mpipi-GG'}


def degree_raster():
    cp = os.path.join(CACHE, f'degree_raster_{UID}.npy')
    if os.path.exists(cp):
        return np.load(cp)
    X = np.load(os.path.join(CACHE, f'xyz_{UID}_1000.npy'))
    n = X.shape[1]
    D = np.zeros((n, len(X)), dtype=np.uint8)
    for f, x in enumerate(X):
        adj = hf.contact_lists([x[:, 0], x[:, 1], x[:, 2]], n)
        D[:, f] = [len(a) for a in adj]
    np.save(cp, D)
    return D


def readout_pair(ax, a, b, labels, ylabel, note=None):
    fs.half_violins(ax, [a, b], [OCC, DEG], labels, max_points=150)
    ax.set_ylabel(ylabel, linespacing=1.1)
    if note:
        ax.text(0.5, 1.0, note, transform=ax.transAxes, ha='center', va='top', fontsize=5.5, color='#333333')


def main():
    D = degree_raster()
    H = np.load(os.path.join(CACHE, f'raster_{UID}_ties_random.npy'))
    n, nf = H.shape
    fig = plt.figure(figsize=(7.2, 2.75))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.12, 1], hspace=0.62, wspace=0.46, left=0.07, right=0.98, top=0.93,
                          bottom=0.12)

    # a: the part of hub occupancy that degree does not explain is shared between force fields (hub_vs_degree Part A2)
    import csv as _csv
    top = gs[0, :2].subgridspec(1, 3, width_ratios=[0.85, 1.3, 0.12], wspace=0.5)
    ax = fig.add_subplot(top[0, 0])
    res = {}
    for r in _csv.DictReader(open(os.path.join(CACHE, 'a2_residual.csv'))):
        if r['rho'] not in ('', 'nan'):
            res.setdefault(r['kind'], []).append(float(r['rho']))
    ax.axhline(0, color='0.6', lw=0.5, ls=(0, (3, 2)), zorder=0)
    fs.half_violins(ax, [res['same protein'], res['control']], [OCC, '#a6a6a6'], ['same\nprotein', 'unrelated\nprotein'],
                    max_points=150)
    for k, key in enumerate(('same protein', 'control')):
        ax.text(k + 0.14, 1.0, f'{np.median(res[key]):.2f}', ha='center', va='bottom', fontsize=5.6)
    ax.set_ylim(-0.75, 1.12)
    ax.set_ylabel('ρ between\nforce fields', labelpad=1, linespacing=1.1)
    fs.panel(ax, 'a', x=-0.42, y=1.06, title='the difference is shared')

    # b: residue types among the residues that only one readout ranks in its top 10% (hub_vs_degree Part A3)
    ax = fig.add_subplot(top[0, 1])
    enr = {}
    for r in _csv.DictReader(open(os.path.join(CACHE, 'a3_enrichment.csv'))):
        enr[(r['force_field'], r['set'], r['aa'])] = float(r['enrichment'])
    types = ['G', 'S', 'A', 'N', 'Y', 'F', 'W', 'R', 'K', 'E']
    x = np.arange(len(types))
    for j, (key, col, lab) in enumerate((('deg_only', DEG, 'degree only'), ('occ_only', OCC, 'occupancy only'))):
        xo = x + (j - 0.5) * 0.36
        ax.bar(xo, [np.log2(enr[('CALVADOS-2', key, t)]) for t in types], 0.34, color=col, label=lab)
        ax.scatter(xo, [np.log2(enr[('Mpipi-Recharged', key, t)]) for t in types], s=10, marker='D', color='white',
                   edgecolor='0.25', lw=0.6, zorder=3, label='Recharged' if j == 1 else None)
    ax.axhline(0, color='0.5', lw=0.5)
    ax.set_xticks(x, types)
    for t in ax.get_xticklabels():
        t.set_color('#C0682F' if t.get_text() in 'FWYR' else '#3E8FC6' if t.get_text() in 'GS' else '#555555')
        t.set_fontweight('bold' if t.get_text() in 'FWYRGS' else 'normal')
    ax.set_yticks([-3, -2, -1, 0, 1, 2], ['1/8', '1/4', '1/2', '1', '2', '4'])
    ax.set_ylabel('enrichment', labelpad=1)
    ax.set_ylim(-3.4, 3.6)
    ax.legend(frameon=False, fontsize=5.2, loc='upper center', handletextpad=0.3, borderaxespad=0.0, ncol=3,
              columnspacing=0.8, bbox_to_anchor=(0.5, 1.04))
    fs.panel(ax, 'b', x=-0.2, y=1.06, title='and maps onto spacers and stickers')

    # b, c: flicker and reliability, three force fields
    SHORT = {'CALVADOS-2': 'CAL', 'MPIPI-GG': 'GG', 'Mpipi-GG': 'GG', 'Mpipi-Recharged': 'Rech'}   # as in panel e
    subb = gs[0, 2].subgridspec(2, 1, hspace=1.6)
    axb1, axb2 = fig.add_subplot(subb[0]), fig.add_subplot(subb[1])
    stats = []
    for d, ff in FFS:
        rows = [r for r in csv.DictReader(open(os.path.join(fs.EV, d, 'per_protein.csv'))) if r['jaccard_successive']]
        js = np.array([float(r['jaccard_successive']) for r in rows])
        jr = np.median([float(r['jaccard_random_approx']) for r in rows])
        so = np.array([float(r['splithalf_r_occupancy']) for r in rows if r['splithalf_r_occupancy']])
        so = 2 * so / (1 + so)
        stats.append((ff, js, jr, so))
    ys = np.arange(len(stats))[::-1]
    for y_, (ff, js, jr, so) in zip(ys, stats):
        col = fs.FF_COLORS[ff]
        axb1.plot([jr, np.median(js)], [y_, y_], color='0.75', lw=0.6, zorder=1)
        axb1.plot(np.percentile(js, [25, 75]), [y_, y_], color=col, lw=1.4, solid_capstyle='butt', zorder=2)
        axb1.scatter([np.median(js)], [y_], s=14, color=col, zorder=3)
        axb1.scatter([jr], [y_], s=14, facecolor='white', edgecolor='0.45', lw=0.7, zorder=3)
        axb2.plot(np.percentile(so, [25, 75]), [y_, y_], color=col, lw=1.4, solid_capstyle='butt', zorder=2)
        axb2.scatter([np.median(so)], [y_], s=14, color=col, zorder=3)
    for ax in (axb1, axb2):
        ax.set_yticks(ys)
        ax.set_yticklabels([SHORT.get(f, f) for f, *_ in stats], fontsize=5.5)
        for t, (ff, *_) in zip(ax.get_yticklabels(), stats):
            t.set_color(fs.FF_COLORS[ff])
        ax.set_ylim(-0.55, len(stats) - 0.45)
        ax.tick_params(axis='y', length=0)
        ax.spines['left'].set_visible(False)
    axb1.set_xlim(0.03, 0.085)
    axb1.set_xticks([0.04, 0.06, 0.08])
    axb1.set_xlabel('hub overlap, successive snapshots', labelpad=1, fontsize=6.3)
    axb1.text(stats[0][2], ys[0] + 0.32, 'random', fontsize=5.5, color='0.45', ha='center', va='bottom')
    fs.panel(axb1, 'c', x=-0.5, y=1.12, title='individual hubs flicker')
    axb2.set_xlim(0.8, 1.0)
    axb2.set_xticks([0.8, 0.9, 1.0])
    axb2.set_xlabel('split-half reliability of hub occupancy', labelpad=1, fontsize=6.3)
    fs.panel(axb2, 'd', x=-0.5, y=1.12, title='hub occupancy is stable')

    # d: the two readouts within proteins
    ax = fig.add_subplot(gs[1, 0])
    w = {}
    for r in csv.DictReader(open(os.path.join(CACHE, 'within_per_protein.csv'))):
        if r['rho_occ_deg'] not in ('', 'nan'):
            w.setdefault(r['force_field'], []).append(float(r['rho_occ_deg']))
    order = ['CALVADOS-2', 'MPIPI-GG', 'Mpipi-Recharged']
    fs.half_violins(ax, [w[f] for f in order], [fs.FF_COLORS[f] for f in order], ['CAL', 'GG', 'Rech'], max_points=150)
    for k, f in enumerate(order):
        ax.text(k + 0.14, 1.005, f'{np.median(w[f]):.2f}', ha='center', va='bottom', fontsize=5.6)
    ax.set_ylim(0.5, 1.04)
    ax.set_ylabel('ρ within protein')
    fs.panel(ax, 'e', x=-0.3, y=1.1, title='but the two readouts largely agree')

    # e: between force fields
    ax = fig.add_subplot(gs[1, 1])
    cf = [r for r in csv.DictReader(open(os.path.join(CACHE, 'crossff_per_protein.csv')))]
    cf = [r for r in cf if r['rho_occupancy'] not in ('', 'nan') and r['rho_degree'] not in ('', 'nan')]
    ro = np.array([float(r['rho_occupancy']) for r in cf])
    rd = np.array([float(r['rho_degree']) for r in cf])
    readout_pair(ax, ro, rd, ['hub\noccupancy', 'degree'], 'ρ between\nforce fields',
                 note=f'{len(cf):,} proteins, degree higher in {100 * np.mean(rd > ro):.0f}%')
    ax.set_ylim(-0.6, 1.3)
    fs.panel(ax, 'f', x=-0.3, y=1.1, title='degree agrees better between force fields')

    # f: between independent experimental ensembles
    ax = fig.add_subplot(gs[1, 2])
    po, pd_ = [], []
    for line in open(os.path.join(fs.EV, 'hub_vs_degree', 'part_b.txt')):
        m = re.match(r'\s+\S+\s+PED\d+-PED\d+\s+\d+\s+\|\s+([+-][0-9.]+)\s+[+-][0-9.]+\s+\|\s+([+-][0-9.]+)', line)
        if m:
            po.append(float(m.group(1)))
            pd_.append(float(m.group(2)))
    readout_pair(ax, po, pd_, ['hub\noccupancy', 'degree'], 'ρ between\nensembles',
                 note=f'{len(po)} pairs, degree higher in {sum(d > o for o, d in zip(po, pd_))}')
    ax.set_ylim(-0.1, 1.15)
    fs.panel(ax, 'g', x=-0.3, y=1.1, title='and between experimental ensembles')
    for a_ in fig.axes:
        a_.tick_params(labelsize=6)
    fs.save(fig, 'fig2_readouts')
    print('pairs', len(po), 'median occ', np.median(po), 'deg', np.median(pd_), '; crossff', np.median(ro), np.median(rd))


if __name__ == '__main__':
    main()
