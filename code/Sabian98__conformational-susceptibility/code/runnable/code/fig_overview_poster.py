"""
fig_overview_poster.py — poster version of Fig. 1 (2026-10-31), for talks and posters, not used in the manuscript (Fig. 1 is
fig_overview.py, in journal style). The larger scheme of the paper, drawn in slide units of 0.1 in
and scaled to the text width in LaTeX, with fonts set so that they print at 7 to 8.5 pt. Since 2026-10-31 a bolder
layout of its own (slide 00 of the talk deck, fig_slides_deck.overview, is unchanged):
  band  headline numbers (as stated in the text) with line icons on a dark gradient
  a     ring of the 2,156 organisms, one spoke per organism (a base length plus the log number of its proteins), coloured by taxonomic
        group in the colours of Fig. 4d (fig_extra.KCOL)
  b     the four perturbations feeding one ensemble, drawn as 30 real stored frames of the example chain of Figs. 2
        and 3 (aligned, coloured by Rg)
  c-e   force fields, measurement techniques, what it reveals
  f     one real stored trajectory through to its susceptibility, in three numbered steps
Card headers are coloured by their role in the story (data slate, method blue, findings terracotta). The group counts,
the organism ring, the ensemble and the trajectory are data, the rest is a cartoon.
    BENDER_REANALYSIS_DIR=$CS_REPO/bender_reanalysis python fig_overview_poster.py      (openmpipi env)
Output: BENDER-MT/slides/fig_overview_poster.pdf
"""
import sys

sys.dont_write_bytecode = True

import os  # noqa: E402

import matplotlib  # noqa: E402

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch, Rectangle  # noqa: E402

HERE = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, os.environ.get('BENDER_REANALYSIS_DIR', os.path.expanduser((__import__('os').environ['CS_REPO'] + '/bender_reanalysis'))))
sys.path.insert(0, HERE)
import fig_extra as FX  # noqa: E402
import fig_slides_deck as S  # noqa: E402

OUT = os.path.join(os.path.dirname(HERE), 'slides', 'fig_overview_poster.pdf')
X0, X1 = 3.2, 96.3                                # crop, in slide units of 0.1 in
XL, XR = 3.5, 96.0                                # outer edges of the cards
# story colours (not figfont roles, so they print as written)
NAVY, DATA, METHOD, FIND = '#1E3A5C', '#3A5A78', '#2F6E9E', '#A8574A'
TINT = {DATA: '#F3F5F8', METHOD: '#F1F6FB', FIND: '#FBF4F2'}
INK, MUTED = S.INK, S.MUTED
HH = 4.6                                          # card header height


def panel(ax, x, y, w, h, letter, title, col, fs):
    """Card with a soft shadow, a coloured header carrying the white panel letter and title, and a tinted body."""
    r = 1.4
    for dx, dy, a in ((0.25, -0.35, 0.10), (0.5, -0.7, 0.05)):
        ax.add_patch(FancyBboxPatch((x + dx, y + dy), w, h, boxstyle=f'round,pad=0,rounding_size={r}', fc='black',
                                    ec='none', alpha=a, zorder=0))
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle=f'round,pad=0,rounding_size={r}', fc=col, ec='none', zorder=0.5))
    ax.add_patch(FancyBboxPatch((x, y), w, h - HH, boxstyle=f'round,pad=0,rounding_size={r}', fc=TINT[col], ec='none',
                                zorder=0.6))
    ax.add_patch(Rectangle((x, y + h - HH - r), w, r, fc=TINT[col], ec='none', zorder=0.6))
    yc = y + h - HH / 2
    ax.add_patch(Circle((x + 2.6, yc), 1.3, fc='white', ec='none', zorder=5))
    ax.text(x + 2.6, yc - 0.05, letter, ha='center', va='center', color=col, fontsize=fs(12), fontweight='bold', zorder=6)
    ax.text(x + 4.9, yc, title, fontsize=fs(12.5), fontweight='bold', color='white', va='center', zorder=6)
    return y + h - HH                             # top of the body


def _icon(ax, kind, x, y, s=1.0):
    """Small white line icons of the band (lines and circles only)."""
    wc = dict(color='white', lw=1.1, zorder=3, solid_capstyle='round')
    if kind == 'protein':                         # a chain of beads
        t = np.linspace(-1.8, 1.8, 40)
        ax.plot(x + s * t, y + s * 0.5 * np.sin(2.4 * t), **wc)
        for u in np.linspace(-1.8, 1.8, 6):
            ax.add_patch(Circle((x + s * u, y + s * 0.5 * np.sin(2.4 * u)), 0.32 * s, fc='white', ec='none', zorder=4))
    elif kind == 'traj':                          # three stacked time series
        rng = np.random.default_rng(3)
        t = np.linspace(-1.9, 1.9, 22)
        for k in range(3):
            ax.plot(x + s * t, y + s * (0.75 * (k - 1) + 0.13 * rng.standard_normal(len(t))), color='white',
                    lw=0.7, alpha=1.0 - 0.3 * k, zorder=3)
    elif kind == 'frames':                        # a stack of frames
        for k in range(3):
            ax.add_patch(Rectangle((x - 1.5 * s + 0.5 * s * k, y - 1.1 * s + 0.35 * s * k), 2.0 * s, 1.5 * s,
                                   fc='#2C5A86', ec='white', lw=0.9, zorder=3 + k))
    elif kind == 'mutant':                        # one bead turning into another
        ax.add_patch(Circle((x - 1.3 * s, y), 0.65 * s, fc='none', ec='white', lw=1.1, zorder=3))
        ax.add_patch(Circle((x + 1.3 * s, y), 0.65 * s, fc='white', ec='none', zorder=3))
        ax.add_patch(FancyArrowPatch((x - 0.55 * s, y), (x + 0.55 * s, y), arrowstyle='-|>', mutation_scale=7,
                                     lw=1.0, color='white', zorder=3, shrinkA=0, shrinkB=0))
    elif kind == 'phospho':                       # a phosphate on a chain
        t = np.linspace(-1.8, 1.8, 30)
        ax.plot(x + s * t, y - 0.6 * s + s * 0.25 * np.sin(2.4 * t), **wc)
        ax.plot([x, x], [y - 0.5 * s, y + 0.1 * s], **wc)
        ax.add_patch(Circle((x, y + 0.75 * s), 0.7 * s, fc='white', ec='none', zorder=4))
        ax.text(x, y + 0.73 * s, 'P', ha='center', va='center', color=NAVY, fontsize=5.5, fontweight='bold', zorder=5)


def band(ax, y, h, fs, tr):
    """Headline numbers on a dark gradient, each as stated in the manuscript text (Methods and Results), with a line
    icon above each and faint real conformations of the example chain drifting behind."""
    box = FancyBboxPatch((XL, y), XR - XL, h, boxstyle='round,pad=0,rounding_size=1.6', fc='none', ec='none')
    ax.add_patch(box)
    cmap = LinearSegmentedColormap.from_list('band', [NAVY, '#2C5A86', '#3C7DB0'])
    im = ax.imshow(np.linspace(0, 1, 256)[None, :], extent=(XL, XR, y, y + h), cmap=cmap, aspect='auto', zorder=0.5,
                   interpolation='bicubic')
    im.set_clip_path(box)
    xyz = np.asarray(tr['xyz'], float)            # faint conformations behind the numbers
    rng = np.random.default_rng(7)
    for k, xc in enumerate(np.linspace(XL - 2, XR + 2, 14)):
        p = xyz[rng.integers(len(xyz))]
        p = p - p.mean(0)
        _, _, vt = np.linalg.svd(p, full_matrices=False)
        q = p @ vt[:2].T
        q = q * (0.42 * h / np.abs(q).max())
        ln, = ax.plot(xc + q[:, 0], y + h / 2 + rng.uniform(-0.25, 0.25) * h + q[:, 1], color='white', lw=0.5,
                      alpha=0.10, zorder=1)
        ln.set_clip_path(box)
    items = [('protein', '10,282', 'highly disordered proteins'), ('traj', '67,415', 'simulated trajectories'),
             ('frames', '68 million', 'stored conformations'), ('mutant', '96,672', 'mutants evaluated\nwithout simulation'),
             ('phospho', '7,528', 'phosphorylation sites\ntested')]
    for k, (ic, num, lab) in enumerate(items):
        xc = XL + (XR - XL) * (k + 0.5) / len(items)
        _icon(ax, ic, xc, y + h - 2.6)
        ax.text(xc, y + h - 6.6, num, ha='center', va='center', fontsize=fs(21), fontweight='bold', color='white',
                zorder=3)
        ax.text(xc, y + 2.6, lab, ha='center', va='center', fontsize=fs(10.5), color='#D6E4F0', linespacing=1.15,
                zorder=3)
        if k:
            xs = XL + (XR - XL) * k / len(items)
            ax.plot([xs, xs], [y + 1.6, y + h - 1.6], color='white', alpha=0.25, lw=0.8, zorder=1)


def organism_ring(ax, cx, cy, bio, fs, r0=5.6, L=4.6):
    """One spoke per organism, grouped by taxonomic group, length rising with the log number of its proteins."""
    order = ['Bacteria', 'Protists', 'Plants', 'Fungi', 'Mammals', 'Viruses']
    per = bio.groupby(['kingdom', 'taxonomy_id']).size()
    gap = np.radians(4.0)
    span = 2 * np.pi - gap * len(order)
    n_org = per.groupby(level=0).size()
    th = np.pi / 2 - gap / 2
    nmax = np.log1p(per.max())
    for g in order:
        v = per.loc[g].sort_values(ascending=False).values
        ang = th - span * (np.arange(len(v)) + 0.5) / n_org.sum()
        rr = r0 + L * (0.3 + 0.7 * np.log1p(v) / nmax)    # a base length, so single-protein organisms show
        segs = [[(cx + r0 * np.cos(a), cy + r0 * np.sin(a)), (cx + q * np.cos(a), cy + q * np.sin(a))]
                for a, q in zip(ang, rr)]
        from matplotlib.collections import LineCollection
        ax.add_collection(LineCollection(segs, colors=FX.KCOL[g], linewidths=0.4, zorder=3))
        a0, a1 = th, th - span * len(v) / n_org.sum()
        t = np.linspace(a0, a1, 60)
        ax.plot(cx + (r0 - 0.8) * np.cos(t), cy + (r0 - 0.8) * np.sin(t), color=FX.KCOL[g], lw=2.6, zorder=3,
                solid_capstyle='butt')
        th = a1 - gap
    ax.text(cx, cy + 0.9, f'{bio.taxonomy_id.nunique():,}', ha='center', va='center', fontsize=fs(15),
            fontweight='bold', color=INK)
    ax.text(cx, cy - 1.6, 'organisms', ha='center', va='center', fontsize=fs(10), color=MUTED)
    return order


def ensemble_cloud(ax, hub, R, tr, n=30):
    """n real stored frames spanning the Rg range, aligned onto the median frame, coloured from compact (dark) to
    expanded (pale), with the median frame on top and its charged residues marked."""
    xyz, rg, seq = np.asarray(tr['xyz'], float), np.asarray(tr['rg'], float), str(tr['seq'])
    order = np.argsort(rg)
    pick = order[np.linspace(len(rg) * 0.03, len(rg) * 0.97, n).astype(int)]
    ref = xyz[order[len(rg) // 2]] - xyz[order[len(rg) // 2]].mean(0)
    _, _, vt = np.linalg.svd(ref, full_matrices=False)

    def align(p):
        p = p - p.mean(0)
        u, _, wt = np.linalg.svd(p.T @ ref)
        d = np.sign(np.linalg.det(u @ wt))
        return p @ (u @ np.diag([1, 1, d]) @ wt) @ vt[:2].T

    proj = [align(xyz[k]) for k in pick]
    rad = np.hypot(*np.concatenate(proj).T)
    sc = 0.9 * R / np.percentile(rad, 99)           # scaled to the typical frames, the far tails clipped by the circle
    clip = Circle(hub, R - 0.15, transform=ax.transData)
    cmap = LinearSegmentedColormap.from_list('rg', ['#1E3A5C', '#3C7DB0', '#A8CBE8'])
    for j, (k, q) in enumerate(zip(pick, proj)):
        ln, = ax.plot(hub[0] + sc * q[:, 0], hub[1] + sc * q[:, 1], color=cmap(j / (n - 1)), lw=0.7, alpha=0.75,
                      zorder=3, solid_capstyle='round')
        ln.set_clip_path(clip)
    q = align(ref)
    ax.plot(hub[0] + sc * q[:, 0], hub[1] + sc * q[:, 1], color=INK, lw=1.3, zorder=4, solid_capstyle='round')
    chg = np.array([{'K': 1, 'R': 1, 'D': -1, 'E': -1}.get(r, 0) for r in seq])
    for v, col in ((1, S.BASIC), (-1, S.ACID)):
        m = chg == v
        ax.scatter(hub[0] + sc * q[m, 0], hub[1] + sc * q[m, 1], s=7, color=col, ec='white', lw=0.3, zorder=5)
    return cmap


def draw(ax, bio, tr, f, floor):
    fs = lambda s: max(s * f, floor)  # noqa: E731
    yt = 74.0
    # band
    yb = yt - 13.2
    band(ax, yb, 13.2, fs, tr)
    # row 1, a and b
    y1, h1 = yb - 2.0 - 38.5, 38.5
    top = panel(ax, XL, y1, 29.0, h1, 'a', 'Proteins from the tree of life', DATA, fs)
    ring = (XL + 14.5, top - 11.6)
    order = organism_ring(ax, *ring, bio, fs)
    n = bio.kingdom.value_counts()
    for k, g in enumerate(order):
        xx, yy = XL + 2.0 + 13.8 * (k % 2), y1 + 8.3 - 2.9 * (k // 2)
        ax.add_patch(Rectangle((xx, yy - 0.75), 1.5, 1.5, fc=FX.KCOL[g], ec='none', zorder=4))
        ax.text(xx + 2.2, yy, g, fontsize=fs(10), color=INK, va='center')
        ax.text(xx + 12.6, yy, f'{n[g]:,}', fontsize=fs(10), color=MUTED, va='center', ha='right')
    ax.text(XL + 14.5, y1 + 10.9, 'Spoke length, proteins per organism', fontsize=fs(9.5), color=MUTED, ha='center',
            va='center')
    S.arrow(ax, (XL + 29.6, y1 + h1 / 2 - 2), (37.0, y1 + h1 / 2 - 2), col=NAVY, lw=2.6, ms=16)

    xb, wb = 37.5, XR - 37.5
    top = panel(ax, xb, y1, wb, h1, 'b', 'One ensemble answers many perturbations', METHOD, fs)
    hub, R = (xb + 45.0, top - 14.6), 10.0
    for rr, a in ((R + 2.2, 0.10), (R + 1.2, 0.18), (R + 0.4, 0.30)):           # soft glow
        ax.add_patch(Circle(hub, rr, fc='#A8CBE8', ec='none', alpha=a, zorder=1))
    ax.add_patch(Circle(hub, R, fc='white', ec='#3C7DB0', lw=2.0, zorder=2))
    ensemble_cloud(ax, hub, R, tr)
    ax.text(hub[0], hub[1] + R + 3.3, 'One ensemble at 300 K', ha='center', va='center', fontsize=fs(11.5),
            fontweight='bold', color=METHOD)
    cb0, cb1, cby = hub[0] - 6.0, hub[0] + 6.0, hub[1] - R - 2.4                # compact to expanded key
    cmap = LinearSegmentedColormap.from_list('rg', ['#1E3A5C', '#3C7DB0', '#A8CBE8'])
    ax.imshow(np.linspace(0, 1, 64)[None, :], extent=(cb0, cb1, cby - 0.45, cby + 0.45), cmap=cmap, aspect='auto',
              zorder=3)
    ax.text(cb0 - 0.6, cby, 'compact', ha='right', va='center', fontsize=fs(9.5), color=MUTED)
    ax.text(cb1 + 0.6, cby, 'expanded', ha='left', va='center', fontsize=fs(9.5), color=MUTED)
    ys = [top - 4.6, top - 11.4, top - 18.2, top - 25.0]
    labs = ['Heat', 'Salt', 'Phosphorylation', 'Substitution']
    cols = [S.HYD, S.BASIC, S.PHOS, '#2E8C6F']
    how = ['hydration and electrostatics', 'electrostatic screening', 'charge of the site', 'residue parameters']
    xi = xb + 5.0
    for yy, lab, col, hw in zip(ys, labs, cols, how):
        ax.add_patch(Circle((xi, yy), 2.55, fc='white', ec=col, lw=2.2, zorder=4))
        ax.text(xi + 3.7, yy + 0.9, lab, fontsize=fs(12), color=INK, va='center', fontweight='bold')
        ax.text(xi + 3.7, yy - 1.3, hw, fontsize=fs(10), color=MUTED, va='center')
        ang = np.arctan2(yy - hub[1], -2 * R)
        tip = (hub[0] + (R + 0.9) * np.cos(ang), hub[1] + (R + 0.9) * np.sin(ang))
        S.arrow(ax, (xb + 26.5, yy), tip, col=col, lw=2.6, ms=15)
    S.thermometer(ax, xi, ys[0] + 0.25, 0.62)
    for dx, c, s in ((-0.95, S.BASIC, '+'), (0.95, S.ACID, '−')):
        ax.add_patch(Circle((xi + dx, ys[1]), 0.88, fc=c, ec='white', lw=0.6, zorder=6))
        ax.text(xi + dx, ys[1] - 0.05, s, ha='center', va='center', color='white', fontsize=9 * f, fontweight='bold',
                zorder=7)
    ax.add_patch(Circle((xi, ys[2]), 1.3, fc=S.PHOS, ec='white', lw=0.6, zorder=6))
    ax.text(xi, ys[2] - 0.05, 'P', ha='center', va='center', color='white', fontsize=10 * f, fontweight='bold', zorder=7)
    for dx, c in ((-1.0, S.BASIC), (1.0, S.ACID)):
        ax.add_patch(Circle((xi + dx, ys[3]), 0.8, fc=c, ec='white', lw=0.6, zorder=6))
    ax.text(xi, ys[3] + 1.5, '→', ha='center', va='center', fontsize=10 * f, color=INK, zorder=7)
    ax.text(xb + wb / 2, y1 + 2.8, 'Each perturbation is computed from the same stored conformations', ha='center',
            va='center', fontsize=fs(11), color=INK, fontstyle='italic')

    # row 2, c d e
    h2 = 24.5
    y2 = y1 - 2.4 - h2
    S.arrow(ax, (hub[0], y1 - 0.2), (hub[0], y2 + h2 + 0.2), col=NAVY, lw=2.6, ms=14)
    top = panel(ax, XL, y2, 26.5, h2, 'c', 'Three force fields', DATA, fs)
    S._force_fields(ax, XL, top + HH + 0.6, fs)
    top = panel(ax, 33.0, y2, 28.5, h2, 'd', 'Tested on 10 proteins', DATA, fs)
    for cx in (33.0 + 1.6, 33.0 + 15.1):            # a white tile behind each technique
        for cy in (top - 4.7, top - 14.2):
            ax.add_patch(FancyBboxPatch((cx, cy - 5.0), 11.8, 8.0, boxstyle='round,pad=0,rounding_size=1.0',
                                        fc='white', ec='#D5DEE8', lw=0.8, zorder=1))
    S._techniques(ax, 33.0, top + HH + 0.9, fs)
    xe = 64.5
    top = panel(ax, xe, y2, XR - xe, h2, 'e', 'What it reveals', FIND, fs)
    rows = [('Rules', 'Nearly every protein compacts on\nheating, set by compaction and charge.', S.C2),
            ('Signals', 'Phosphorylation and salt change heat\nsensitivity in step with compaction.', S.PHOS),
            ('Evolution', 'Composition adapts to growth\ntemperature, heat sensitivity does not.', '#2E8C6F')]
    for k, (hd, t, col) in enumerate(rows):
        yy = top - 2.2 - 6.3 * k
        ax.add_patch(Rectangle((xe + 1.4, yy - 4.6), 0.6, 5.6, fc=col, ec='none', zorder=4))
        ax.text(xe + 3.0, yy, hd, fontsize=fs(12), fontweight='bold', color=INK, va='center')
        ax.text(xe + 3.0, yy - 2.8, t, fontsize=fs(10.5), color=INK, va='center', linespacing=1.2)

    # row 3, f
    h3 = 23.6
    y3 = y2 - 2.0 - h3
    top = panel(ax, XL, y3, XR - XL, h3, 'f', 'From one trajectory to a susceptibility', METHOD, fs)
    steps = [(XL + 6.0, 'Store the conformations'), (XL + 42.6, 'Score each frame'),
             (XL + 66.1, 'Reweight for heating')]
    for k, (xs, t) in enumerate(steps):
        yy = top - 2.4
        ax.add_patch(Circle((xs, yy), 1.15, fc=METHOD, ec='none', zorder=5))
        ax.text(xs, yy - 0.05, str(k + 1), ha='center', va='center', color='white', fontsize=fs(10),
                fontweight='bold', zorder=6)
        ax.text(xs + 1.8, yy, t, ha='left', va='center', fontsize=fs(10.5), fontweight='bold', color=METHOD)
    S._trajectory(ax, XL, top - 1.0, XR - XL, fs, tr)

    return y3 - 1.0                               # bottom of the drawing


def main():
    FX.house()
    k = FX.WIDTH / ((X1 - X0) / 10)               # LaTeX scale factor (the figure is included at the text width)
    _, bio = FX._load_totals()
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
    y0 = draw(probe.add_axes([0, 0, 1, 1]), bio, tr, f, floor)
    plt.close(probe)
    Y0, Y1 = y0, 74.3
    fig = plt.figure(figsize=((X1 - X0) / 10, (Y1 - Y0) / 10))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(X0, X1)
    ax.set_ylim(Y0, Y1)
    ax.axis('off')
    draw(ax, bio, tr, f, floor)
    fig.savefig(OUT)
    print('written', OUT, f'(printed at {k:.3f} of drawn size, {(Y1 - Y0) / 10 * k:.2f} in tall)')


if __name__ == '__main__':
    main()
