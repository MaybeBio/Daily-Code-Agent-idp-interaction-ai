"""
fig_overview.py — Fig. 1 of the manuscript (2026-10-30), the larger scheme of the paper, drawn in slide units of 0.1 in
and scaled to the text width in LaTeX, with fonts set so that they print at 7 to 8.5 pt. Journal style since 2026-10-31
(white background, bold panel letters, no cards; the bolder poster version with a number band and coloured card headers
is fig_overview_poster.py). The headline numbers (as stated in the text) sit in the panel each belongs to:
  a     10,282 proteins, ring of the 2,156 organisms, one spoke per organism (a base length plus the log number of its
        proteins), in the taxonomic colours of Fig. 4d (fig_extra.KCOL)
  b     the four perturbations (7,528 phosphorylation sites, 96,672 mutants) feeding one ensemble, drawn as 30 real
        stored frames of the example chain of Figs. 2 and 3 (68 million stored conformations in all)
  c, d  force fields (67,415 trajectories), measurement techniques
  e     total response of all 10,282 proteins in CALVADOS-2 and Mpipi-Recharged, the first proteome-wide result
  f     one real stored trajectory through to its susceptibility, in three numbered steps
The group counts, the organism ring, the ensemble, the distribution and the trajectory are data, the rest is a cartoon. Slide 00 of the
talk deck (fig_slides_deck.overview) is unchanged.
    BENDER_REANALYSIS_DIR=~/Downloads/Scripts/bender_reanalysis python fig_overview.py      (openmpipi env)
Output: BENDER-MT/manuscript/figs/fig_overview.pdf
"""
import sys

sys.dont_write_bytecode = True

import os  # noqa: E402

import matplotlib  # noqa: E402

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402
from matplotlib.patches import Circle, Rectangle  # noqa: E402

HERE = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, os.environ.get('BENDER_REANALYSIS_DIR', os.path.expanduser('~/Downloads/Scripts/bender_reanalysis')))
sys.path.insert(0, HERE)
import fig_extra as FX  # noqa: E402
import fig_overview_poster as P  # noqa: E402  (organism ring and ensemble cloud)
import fig_slides_deck as S  # noqa: E402

OUT = os.path.join(os.path.dirname(HERE), 'manuscript', 'figs', 'fig_overview.pdf')
X0, X1 = 3.2, 96.3                                # crop, in slide units of 0.1 in
XL, XR = 3.5, 96.0
INK, MUTED, ACC = S.INK, S.MUTED, '#2F6E9E'       # ACC, the one accent (not a figfont role)
GREY = '#9A9A9A'


def head(ax, x, y, letter, title, fs):
    ax.text(x, y, letter, fontsize=fs(13), fontweight='bold', color=INK, va='center', ha='left')
    ax.text(x + 2.6, y, title, fontsize=fs(11.5), fontweight='bold', color=INK, va='center', ha='left')


def rule(ax, y):
    ax.plot([XL, XR], [y, y], color='#DDDDDD', lw=0.6, zorder=0)


def draw(ax, D, bio, tr, f, floor):
    fs = lambda s: max(s * f, floor)  # noqa: E731
    yt = 74.0

    # a  proteins from the tree of life
    head(ax, XL, yt - 1.6, 'a', 'Sequences from the tree of life', fs)
    ax.text(XL + 13.5, yt - 5.4, f'{len(bio):,} highly disordered proteins', ha='center', va='center',
            fontsize=fs(10.5), color=INK)
    order = P.organism_ring(ax, XL + 13.5, yt - 18.3, bio, fs)
    n = bio.kingdom.value_counts()
    for k, g in enumerate(order):
        xx, yy = XL + 0.6 + 13.8 * (k % 2), yt - 30.8 - 2.7 * (k // 2)
        ax.add_patch(Rectangle((xx, yy - 0.7), 1.4, 1.4, fc=FX.KCOL[g], ec='none', zorder=4))
        ax.text(xx + 2.1, yy, g, fontsize=fs(10), color=INK, va='center')
        ax.text(xx + 12.4, yy, f'{n[g]:,}', fontsize=fs(10), color=MUTED, va='center', ha='right')
    S.arrow(ax, (XL + 28.4, yt - 18.3), (XL + 31.0, yt - 18.3), col=GREY, lw=1.8, ms=12)

    # b  one ensemble answers many perturbations
    xb = 35.5
    head(ax, xb, yt - 1.6, 'b', 'One ensemble answers many perturbations', fs)
    hub, R = (83.5, yt - 18.6), 9.6
    ax.add_patch(Circle(hub, R, fc='white', ec='#BFBFBF', lw=0.9, zorder=2))
    P.ensemble_cloud(ax, hub, R, tr)
    ax.text(hub[0], hub[1] + R + 1.7, 'One ensemble at 300 K', ha='center', va='center', fontsize=fs(10.5),
            fontweight='bold', color=INK)
    cb0, cb1, cby = hub[0] - 5.0, hub[0] + 5.0, hub[1] - R - 1.8
    cmap = LinearSegmentedColormap.from_list('rg', ['#1E3A5C', '#3C7DB0', '#A8CBE8'])
    ax.imshow(np.linspace(0, 1, 64)[None, :], extent=(cb0, cb1, cby - 0.4, cby + 0.4), cmap=cmap, aspect='auto',
              zorder=3)
    ax.text(cb0 - 0.6, cby, 'compact', ha='right', va='center', fontsize=fs(9.5), color=MUTED)
    ax.text(cb1 + 0.6, cby, 'expanded', ha='left', va='center', fontsize=fs(9.5), color=MUTED)
    ys = [yt - 7.4, yt - 14.4, yt - 21.4, yt - 28.4]
    labs = ['Heat', 'Salt', 'Phosphorylation', 'Substitution']
    cols = [S.HYD, S.BASIC, S.PHOS, '#2E8C6F']
    how = ['hydration and electrostatics', 'electrostatic screening', '7,528 sites', '96,672 mutants']
    xi = xb + 3.4
    for yy, lab, col, hw in zip(ys, labs, cols, how):
        ax.add_patch(Circle((xi, yy), 2.45, fc='white', ec=col, lw=1.5, zorder=4))
        ax.text(xi + 3.5, yy + 0.9, lab, fontsize=fs(11), color=INK, va='center', fontweight='bold')
        ax.text(xi + 3.5, yy - 1.3, hw, fontsize=fs(10), color=MUTED, va='center')
        ang = np.arctan2(yy - hub[1], -2 * R)
        tip = (hub[0] + (R + 0.8) * np.cos(ang), hub[1] + (R + 0.8) * np.sin(ang))
        S.arrow(ax, (xb + 26.0, yy), tip, col=col, lw=1.6, ms=11)
    S.thermometer(ax, xi, ys[0] + 0.25, 0.6)
    for dx, c, s in ((-0.95, S.BASIC, '+'), (0.95, S.ACID, '−')):
        ax.add_patch(Circle((xi + dx, ys[1]), 0.86, fc=c, ec='white', lw=0.6, zorder=6))
        ax.text(xi + dx, ys[1] - 0.05, s, ha='center', va='center', color='white', fontsize=9 * f, fontweight='bold',
                zorder=7)
    ax.add_patch(Circle((xi, ys[2]), 1.25, fc=S.PHOS, ec='white', lw=0.6, zorder=6))
    ax.text(xi, ys[2] - 0.05, 'P', ha='center', va='center', color='white', fontsize=10 * f, fontweight='bold', zorder=7)
    for dx, c in ((-1.0, S.BASIC), (1.0, S.ACID)):
        ax.add_patch(Circle((xi + dx, ys[3]), 0.78, fc=c, ec='white', lw=0.6, zorder=6))
    ax.text(xi, ys[3] + 1.5, '→', ha='center', va='center', fontsize=10 * f, color=INK, zorder=7)
    ax.text((xb + XR) / 2, yt - 35.6, 'Each response is computed from the same stored conformations, 68 million in all',
            ha='center', va='center', fontsize=fs(10), color=INK, fontstyle='italic')

    # row 2
    y2 = yt - 39.4
    rule(ax, y2 + 1.4)
    S.arrow(ax, (hub[0], yt - 36.8), (hub[0], y2 - 0.6), col=GREY, lw=1.8, ms=12)
    # c  force fields
    head(ax, XL, y2 - 1.6, 'c', 'Three force fields', fs)
    for k, (m, txt) in enumerate((('CALVADOS-2', 'all proteins, five temperatures'),
                                  ('Mpipi-GG', 'subset, five temperatures'),
                                  ('Mpipi-Recharged', 'all proteins, 300 K'))):
        yy = y2 - 6.4 - 5.5 * k
        col = FX.COL[m]
        for j in range(3):
            ax.add_patch(Circle((XL + 1.4 + 1.5 * j, yy), 0.68, fc=col, ec='white', lw=0.6, zorder=4))
        ax.plot([XL + 1.4, XL + 4.4], [yy, yy], color=S.BOND, lw=1.3, zorder=3)
        ax.text(XL + 6.4, yy, m, fontsize=fs(11), fontweight='bold', color=col, va='center')
        ax.text(XL + 0.6, yy - 2.3, txt, fontsize=fs(10), color=MUTED, va='center')
    ax.text(XL + 0.6, y2 - 23.0, '67,415 trajectories in all', fontsize=fs(10), color=INK, va='center',
            fontstyle='italic')
    # d  measurements
    xd = 30.5
    head(ax, xd, y2 - 1.6, 'd', 'Tested on 10 proteins', fs)
    S._techniques(ax, xd - 0.6, y2 + 2.2, fs)
    # e  total response of every protein, the first proteome-wide result (Fig. 4, Supplementary Fig. 7)
    xe = 61.0
    head(ax, xe, y2 - 1.6, 'e', 'Nearly every protein compacts on heating', fs)
    px0, px1, py0, py1 = xe + 4.0, XR - 1.0, y2 - 19.0, y2 - 5.0
    lo, hi = -1.5, 0.25                            # 99% of proteins in either model
    X = lambda v: px0 + (px1 - px0) * (np.asarray(v) - lo) / (hi - lo)  # noqa: E731
    edges = np.linspace(lo, hi, 49)
    hs = {m: np.histogram(D[k].total.values, edges)[0] / len(D[k]) for m, k in (('CALVADOS-2', 'C'), ('Mpipi-Recharged', 'R'))}
    top = max(h.max() for h in hs.values())
    ex = np.repeat(X(edges), 2)[1:-1]
    for m, h in hs.items():
        yy = py0 + (py1 - 3.8 - py0) * np.repeat(h, 2) / top     # legend above the peak
        if m == 'CALVADOS-2':
            ax.fill_between(ex, py0, yy, color=FX.COL[m], alpha=0.35, lw=0, zorder=2)
        ax.plot(ex, yy, color=FX.COL[m], lw=1.1, zorder=3)
    ax.plot([px0, px1], [py0, py0], color=INK, lw=0.7, zorder=4)
    for v in (-1.2, -0.8, -0.4, 0.0):
        ax.plot([X(v)] * 2, [py0, py0 - 0.5], color=INK, lw=0.7, zorder=4)
        ax.text(X(v), py0 - 1.5, FX.minus(f'{v:.1f}'), ha='center', va='center', fontsize=fs(9.5), color=INK)
    ax.plot([X(0)] * 2, [py0, py1 - 4.2], color=INK, lw=0.7, ls=(0, (2, 1.5)), zorder=4)
    ax.text((px0 + px1) / 2, py0 - 3.6, 'Total response (% per K)', ha='center', va='center', fontsize=fs(10),
            color=INK)
    for k2, (m, key) in enumerate((('CALVADOS-2', 'C'), ('Mpipi-Recharged', 'R'))):
        pct = 100 * (D[key].total.values < 0).mean()
        assert round(pct, 1) == {'C': 98.5, 'R': 98.2}[key]          # as stated in the Results
        yy = py1 + 0.2 - 2.2 * k2
        ax.plot([px0 + 0.2, px0 + 1.8], [yy, yy], color=FX.COL[m], lw=1.6)
        ax.text(px0 + 2.4, yy, f'{m}, {pct:.1f}% compact on heating', fontsize=fs(9.5), color=INK, va='center')

    # row 3, f
    y3 = y2 - 26.4
    rule(ax, y3 + 1.4)
    head(ax, XL, y3 - 1.6, 'f', 'From one trajectory to a susceptibility', fs)
    steps = [(XL + 6.0, 'Store the conformations'), (XL + 42.6, 'Score each frame'), (XL + 66.1, 'Reweight for heating')]
    for k, (xs, t) in enumerate(steps):
        yy = y3 - 5.2
        ax.add_patch(Circle((xs, yy), 1.05, fc=ACC, ec='none', zorder=5))
        ax.text(xs, yy - 0.05, str(k + 1), ha='center', va='center', color='white', fontsize=fs(9.5),
                fontweight='bold', zorder=6)
        ax.text(xs + 1.7, yy, t, ha='left', va='center', fontsize=fs(10.5), fontweight='bold', color=INK)
    S._trajectory(ax, XL, y3 - 3.9, XR - XL, fs, tr)
    return y3 - 3.9 - 17.4                         # bottom of the drawing


def main():
    FX.house()
    k = FX.WIDTH / ((X1 - X0) / 10)               # LaTeX scale factor (the figure is included at the text width)
    D, bio = FX._load_totals()
    import fig_forcing as F                        # the stored 300 K run of the example chain of Figs. 2 and 3
    c = F.compute()
    _, _, _, g = F.check(c)
    val, _ = F.response_se(c, g)
    # total thermal force of each frame (fixed solvent plus both hydration terms), as fig_forcing.response_se
    RT = F.TF.fdt.R_GAS * F.TF.T0
    G = (np.load(os.path.join(F.TF.NPZ, f'{F.UID}.npz'))['G_T'].astype(float) + g['tea_C'] * c['D_tea'] / RT
         + g['C'] * c['D_hyd'] / RT)
    assert abs(-100 * F.TF.fdt._cov(c['rg'], G) / c['rg'].mean() - val['total']) < 1e-9
    tr = dict(xyz=c['xyz'], rg=c['rg'], G=G, seq=c['seq'], total=val['total'])
    f, floor = 8.5 / (12.5 * k), 7.0 / k           # 7 to 8.5 pt
    probe = plt.figure()
    y0 = draw(probe.add_axes([0, 0, 1, 1]), D, bio, tr, f, floor)
    plt.close(probe)
    Y0, Y1 = y0, 74.3
    fig = plt.figure(figsize=((X1 - X0) / 10, (Y1 - Y0) / 10))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(X0, X1)
    ax.set_ylim(Y0, Y1)
    ax.axis('off')
    draw(ax, D, bio, tr, f, floor)
    fig.savefig(OUT)
    print('written', OUT, f'(printed at {k:.3f} of drawn size, {(Y1 - Y0) / 10 * k:.2f} in tall)')


if __name__ == '__main__':
    main()
