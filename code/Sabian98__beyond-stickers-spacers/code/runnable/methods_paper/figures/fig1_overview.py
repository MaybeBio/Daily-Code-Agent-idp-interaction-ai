"""
fig1_overview.py — Figure 1 (2026-10-14 redesign, replaces fig1_pipeline from fig1_split.py): from ensembles to the two
readouts of the contact network, in the ALBATROSS / IDRome overview style. Descriptive, no new statistics.
  a  the sequence of Q9DUC6 with stickers (F, W, Y, R) and spacers (G, S) coloured, and the three force fields with\n     the number of ensembles of each (as in Methods)
  b  CALVADOS-2 ensemble of Q9DUC6 (87 residues): 60 snapshots of the production half aligned onto one highlighted
     snapshot and projected onto its two principal axes; the highlighted chain shows the stickers F, W, Y, R and the
     spacers G, S
  c  the contact network of the highlighted snapshot drawn on the chain (C-alpha within 0.8 nm, |i-j| >= 3), node area
     proportional to betweenness, the 10% most central residues (hubs) ringed
  d  degree and hub occupancy of every residue over 1,000 snapshots (redo_evidence/hub_frames_v3_cal/per_residue.csv)
Snapshot choice (illustrative): among the 1,000 analysed snapshots with a radius of gyration inside the interquartile
range and a contact count in the top quarter, the one whose second principal axis is longest relative to the first.
  ~/.venvs/finches/bin/python fig1_overview.py  ->  fig1_overview.pdf / .png
"""
import csv
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection

import figstyle as fs

sys.path.insert(0, fs.EV)
import hub_frames as hf  # noqa: E402

UID, GROUP = 'Q9DUC6', 'Viruses'
CACHE = os.path.join(fs.HERE, 'cache')
STICK, SPACE, OTHER = '#d95f02', '#9ecae1', '#aeaeb2'
DEG = '#2ca25f'


def cls_color(a):
    return STICK if a in 'FWYR' else SPACE if a in 'GS' else OTHER


def frames():
    """1,000 evenly spaced production snapshots (Å), cached."""
    os.makedirs(CACHE, exist_ok=True)
    cp = os.path.join(CACHE, f'xyz_{UID}_1000.npy')
    if os.path.exists(cp):
        return np.load(cp)
    d = hf.find_traj(UID, GROUP)
    fh, off, frame, natom, has_cell, n_frames = hf.dcd_open(os.path.join(d, f'{UID}.dcd'))
    start, nf = n_frames // 2, 1000
    out = np.zeros((nf, natom, 3))
    for f in range(nf):
        xyz, box = hf.dcd_frame(fh, off, frame, natom, has_cell, start + (f * (n_frames - start)) // nf)
        hf.unwrap(xyz, box)
        out[f] = np.array(xyz).T
    fh.close()
    np.save(cp, out)
    return out


def kabsch(P, Q):
    """Rotate P (centred) onto Q (centred)."""
    H = P.T @ Q
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1, 1, d])
    return P @ (U @ D @ Vt)


def main():
    rows = [r for r in csv.DictReader(open(os.path.join(fs.EV, 'hub_frames_v3_cal', 'per_residue.csv')))
            if r['UniProt_ID'] == UID]
    rows.sort(key=lambda r: int(r['pos']))
    seq = ''.join(r['aa'] for r in rows)
    occ = np.array([float(r['occupancy']) for r in rows])
    deg = np.array([float(r['deg_mean']) for r in rows])
    n = len(seq)
    X = frames()
    assert X.shape[1] == n
    Xc = X - X.mean(1, keepdims=True)
    rg = np.sqrt((Xc ** 2).sum(2).mean(1))
    ncon = np.array([sum(len(a) for a in hf.contact_lists([x[:, 0], x[:, 1], x[:, 2]], n)) // 2 for x in X])
    lo, hi = np.percentile(rg, [25, 75])
    ok = np.flatnonzero((rg >= lo) & (rg <= hi) & (ncon >= np.percentile(ncon, 75)))
    sv = np.array([np.linalg.svd(Xc[f], compute_uv=False) for f in ok])
    f0 = int(ok[np.argmax(sv[:, 1] / sv[:, 0])])          # the roundest projection, so the drawing fills the panel
    ref = Xc[f0]
    _, _, vt = np.linalg.svd(ref, full_matrices=False)
    proj = lambda P: P @ vt[:2].T  # noqa: E731
    R = proj(ref)
    if R[0, 0] > R[-1, 0]:                                  # N terminus on the left
        vt[0] *= -1
        R = proj(ref)
    adj = hf.contact_lists([X[f0][:, 0], X[f0][:, 1], X[f0][:, 2]], n)
    bc = np.array(hf.betweenness(adj, n))
    hubs = set(hf.top_hubs(list(bc), max(1, round(0.1 * n)), hf.tie_rng(UID)))
    cols = [cls_color(a) for a in seq]

    fig = plt.figure(figsize=(7.2, 1.75))
    gs = fig.add_gridspec(1, 4, width_ratios=[0.95, 1.1, 1.15, 1.5], wspace=0.28, left=0.035, right=0.99,
                          top=0.86, bottom=0.2)

    # a: the sequence, stickers and spacers coloured, and the three force fields that simulate it
    ax = fig.add_subplot(gs[0, 0])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_axis_off()
    per = 15
    for i, aa in enumerate(seq):
        r_, c_ = divmod(i, per)
        ax.text(0.02 + c_ * 0.064, 0.9 - r_ * 0.085, aa, family='DejaVu Sans Mono', fontsize=6.2, ha='left', va='top',
                color=STICK if aa in 'FWYR' else '#3E8FC6' if aa in 'GS' else '#8a8a8a',
                fontweight='bold' if aa in 'FWYRGS' else 'normal')
    ytag = 0.9 - (len(seq) // per + 1) * 0.085 - 0.07
    for k, (name, cnt, c) in enumerate([('CALVADOS-2', '10,281', '#1b6ca8'), ('Mpipi-Recharged', '10,266', '#7570b3'),
                                         ('Mpipi-GG', '908', '#d95f02')]):
        y_ = ytag - k * 0.105
        ax.text(0.02, y_, name, fontsize=5.8, va='center', ha='left', color=c,
                bbox=dict(boxstyle='round,pad=0.25,rounding_size=0.3', fc='white', ec=c, lw=0.6))
        ax.text(0.98, y_, f'{cnt} IDRs', fontsize=5.4, va='center', ha='right', color='#444444')
    ax.text(0.02, 0.9 + 0.035, 'sticker', fontsize=5.2, color=STICK, fontweight='bold', va='bottom')
    ax.text(0.27, 0.9 + 0.035, 'spacer', fontsize=5.2, color='#3E8FC6', fontweight='bold', va='bottom')
    fs.panel(ax, 'a', x=0.0, y=1.06, title='sequence, three force fields')
    ax.annotate('', xy=(1.1, 0.55), xytext=(1.0, 0.55), xycoords='axes fraction',
                arrowprops=dict(arrowstyle='-|>', color='0.55', lw=0.8, mutation_scale=7))

    # b: ensemble
    ax = fig.add_subplot(gs[0, 1])
    rng = np.random.default_rng(3)
    for f in rng.choice(len(X), 60, replace=False):
        P = proj(kabsch(Xc[f], ref))
        ax.plot(P[:, 0], P[:, 1], color='0.55', lw=0.45, alpha=0.13, solid_capstyle='round', zorder=1)
    ax.plot(R[:, 0], R[:, 1], color='0.35', lw=0.9, zorder=2)
    ax.scatter(R[:, 0], R[:, 1], s=9, c=cols, edgecolor='0.3', lw=0.25, zorder=3)
    ax.set_aspect('equal')
    ax.set_axis_off()
    fs.panel(ax, 'b', x=0.0, y=1.04, title=f'one ensemble ({UID})')

    # c: network of the highlighted snapshot
    ax = fig.add_subplot(gs[0, 2])
    ax.plot(R[:, 0], R[:, 1], color='0.82', lw=0.8, zorder=1)
    segs = [(R[i], R[j]) for i in range(n) for j in adj[i] if j > i]
    ax.add_collection(LineCollection(segs, colors='#2ca25f', linewidths=0.5, alpha=0.7, zorder=2))
    size = 5 + 80 * bc / bc.max()
    ax.scatter(R[:, 0], R[:, 1], s=size, c=cols, edgecolor='0.3', lw=0.25, zorder=3)
    hh = sorted(hubs)
    ax.scatter(R[hh, 0], R[hh, 1], s=size[hh] + 28, facecolor='none', edgecolor='#b2182b', lw=0.8, zorder=4)
    ax.set_aspect('equal')
    ax.set_axis_off()
    fs.panel(ax, 'c', x=0.0, y=1.04, title='network of one snapshot')
    for a_ in (fig.axes[1],):
        a_.annotate('', xy=(1.06, 0.5), xytext=(0.96, 0.5), xycoords='axes fraction',
                    arrowprops=dict(arrowstyle='-|>', color='0.55', lw=0.8, mutation_scale=7))

    # d: readouts along the chain
    sub = gs[0, 3].subgridspec(3, 1, height_ratios=[0.16, 1, 1], hspace=0.12)
    a0 = fig.add_subplot(sub[0])
    a1 = fig.add_subplot(sub[1], sharex=a0)
    a2 = fig.add_subplot(sub[2], sharex=a0)
    x = np.arange(1, n + 1)
    a0.bar(x, np.ones(n), width=1.0, color=cols, lw=0)
    a0.set_ylim(0, 1)
    a0.set_axis_off()
    a1.fill_between(x, deg, color=DEG, alpha=0.35, lw=0, step='mid')
    a1.plot(x, deg, color=DEG, lw=0.8, drawstyle='steps-mid')
    a1.set_ylabel('degree', fontsize=6, labelpad=1)
    a1.set_ylim(0, deg.max() * 1.15)
    a2.bar(x, occ, width=1.0, color=cols, lw=0)
    a2.axhline(0.1, color='0.4', lw=0.5, ls=(0, (3, 2)))
    a2.set_ylabel('hub\noccupancy', fontsize=6, labelpad=1)
    a2.set_xlim(0.5, n + 0.5)
    a2.set_xticks([1, 40, n])
    a2.set_xlabel('residue', labelpad=1)
    for a_ in (a1, a2):
        a_.tick_params(labelsize=5.5)
    a1.tick_params(labelbottom=False)
    a1.text(0.02, 0.95, 'how densely contacted', transform=a1.transAxes, fontsize=5.4, va='top', color='#333333')
    a2.text(0.02, 0.95, 'how often a bridge', transform=a2.transAxes, fontsize=5.4, va='top', color='#333333')
    fs.panel(a0, 'd', x=-0.2, y=2.2, title='two readouts over 1,000 snapshots')
    fig.axes[2].annotate('', xy=(1.08, 0.5), xytext=(0.98, 0.5), xycoords='axes fraction',
                         arrowprops=dict(arrowstyle='-|>', color='0.55', lw=0.8, mutation_scale=7))
    fs.save(fig, 'fig1_overview')
    print(f'snapshot {f0}: Rg {rg[f0] / 10:.2f} nm (IQR {lo / 10:.2f}-{hi / 10:.2f}), contacts {ncon[f0]}, hubs {hh}')


if __name__ == '__main__':
    main()
