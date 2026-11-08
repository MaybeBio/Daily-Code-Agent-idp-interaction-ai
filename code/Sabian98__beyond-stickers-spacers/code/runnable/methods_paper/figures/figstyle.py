"""
figstyle.py — shared look for all paper figures: Arial everywhere (text, math, PDF embedding as TrueType),
consistent sizes, and one colour per force field / rule.  `import figstyle` before creating any figure.
Palette  the scripts keep their own hex colours; at save time (save() or recolor()) each colour of the original palette
         is replaced by its light counterpart after ALBATROSS (Lotthammer 2024) and IDRome (Tesei 2024): soft orange,
         teal and sky blue fills with thin outlines a shade darker, mid tones for lines, darker shades for coloured text,
         and viridis for density hexbins, as in IDRome. Set FIG_PALETTE=original to keep the old colours (2026-10-14).
"""
import os
import re

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib import colors as mcolors  # noqa: E402
from matplotlib.collections import Collection, LineCollection, PolyCollection  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
from matplotlib.text import Text  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
EV = os.path.join(HERE, '..', '..', 'redo_evidence')

plt.rcParams.update({
    'font.family': 'Arial', 'font.sans-serif': ['Arial'],
    'mathtext.fontset': 'custom', 'mathtext.rm': 'Arial', 'mathtext.it': 'Arial:italic', 'mathtext.bf': 'Arial:bold',
    'font.size': 7.5, 'axes.titlesize': 8, 'axes.labelsize': 7.5, 'legend.fontsize': 6.5,
    'xtick.labelsize': 6.5, 'ytick.labelsize': 6.5, 'axes.linewidth': 0.6,
    'xtick.major.width': 0.6, 'ytick.major.width': 0.6, 'xtick.major.size': 2.5, 'ytick.major.size': 2.5,
    'axes.spines.top': False, 'axes.spines.right': False,
    'pdf.fonttype': 42, 'ps.fonttype': 42, 'savefig.dpi': 300,
})

FF_COLORS = {'CALVADOS-2': '#1b6ca8', 'MPIPI-GG': '#d95f02', 'Mpipi-Recharged': '#7570b3'}
RULE_COLORS = {'sep2': '#d95f02', 'sep3': '#1b6ca8', 'fixed': '#1b6ca8', 'size': '#e7298a'}
CHARGED, AROMATIC = '#b2182b', '#2166ac'


def color_ticks(ax, axis='x'):
    """Charged residue letters in red, aromatic in blue."""
    for t in (ax.get_xticklabels() if axis == 'x' else ax.get_yticklabels()):
        if t.get_text() in ('D', 'E', 'K', 'R'):
            t.set_color(CHARGED)
        elif t.get_text() in ('F', 'W', 'Y'):
            t.set_color(AROMATIC)


# original colour -> role
ROLE = {'#a6a6a6': 'grey', '#1b6ca8': 'blue', '#0072b2': 'blue', '#2166ac': 'blue', '#2171b5': 'blue', '#4292c6': 'blue_pale',
        '#d95f02': 'orange', '#d55e00': 'orange', '#ff7f00': 'orange', '#e69f00': 'gold',
        '#7570b3': 'purple', '#6a3d9a': 'purple', '#7b3294': 'purple', '#b39ddb': 'purple_pale',
        '#e7298a': 'rose', '#d6006f': 'rose', '#c51b7d': 'rose', '#cc79a7': 'rose',
        '#00a651': 'teal', '#2ca25f': 'teal', '#2bb673': 'teal_pale', '#00843d': 'teal', '#b2182b': 'red',
        '#9ecae1': 'blue_pale', '#8fb4d6': 'blue_pale', '#c6dbef': 'blue_pale', '#6baed6': 'blue_mid',
        '#0096c7': 'blue_mid', '#08306b': 'blue_deep'}
# role -> (fill, line, text); fills sampled from the ALBATROSS and IDRome figures
LIGHT = {'grey': ('#D9D9D9', '#8C8C8C', '#666666'), 'blue': ('#7DB9E2', '#4F9BD0', '#2F7DB5'), 'orange': ('#E7A479', '#DB8A55', '#C0682F'),
         'gold': ('#F2C879', '#E3AA45', '#B8841F'), 'purple': ('#ABA3D9', '#8B80C6', '#6A5DAF'),
         'purple_pale': ('#D3CDEB', '#B4A9DA', '#8B80C6'), 'rose': ('#EE9FC0', '#DE6FA0', '#C24B85'),
         'teal': ('#6BB6AC', '#4E9E93', '#2F7F74'), 'red': ('#E58A7E', '#D0645A', '#B5473A'),
         'blue_pale': ('#BDDDF2', '#90CBEC', '#4F9BD0'), 'blue_mid': ('#90CBEC', '#6AB0DD', '#3E8FC6'),
         'blue_deep': ('#3A78AD', '#2C5F8C', '#24557F'),
         'teal_pale': ('#B5DED7', '#6BB6AC', '#2F7F74')}
PALETTE = os.environ.get('FIG_PALETTE', 'light')
EDGE_LW = 0.4


def _key(c):
    return tuple(round(v, 3) for v in mcolors.to_rgb(c))


_TABLE = {_key(k): LIGHT[r] for k, r in ROLE.items()}


def _new(c, which):
    """Light colour for one original colour (keeping its alpha); None when it is not in the original palette."""
    if c is None or (isinstance(c, str) and c in ('none', 'None', 'auto', 'face')):
        return None
    try:
        rgba = mcolors.to_rgba(c)
    except (ValueError, TypeError):
        return None
    if rgba[3] == 0:
        return None
    t = _TABLE.get(_key(rgba))
    return None if t is None else mcolors.to_rgba(t[which], rgba[3])


def _arr(cols, which):
    out, hit = [], False
    for c in cols:
        n = _new(tuple(c), which)
        hit |= n is not None
        out.append(n if n is not None else tuple(c))
    return out if hit else None


def recolor(fig):
    """Map the original palette onto the light one throughout a figure (0 = fill, 1 = line, 2 = text)."""
    if PALETTE != 'light':
        return
    for a in fig.findobj():
        if isinstance(a, Line2D):
            face = a.get_markerfacecolor()
            n = _new(a.get_color(), 1)
            if n is not None:
                a.set_color(n)
            for get, put, w in (('get_markeredgecolor', 'set_markeredgecolor', 1), ('get_markerfacecolor', 'set_markerfacecolor', 0)):
                v = face if put == 'set_markerfacecolor' else getattr(a, get)()
                m = _new(v, w)
                if m is not None:
                    getattr(a, put)(m)
        elif isinstance(a, Patch):
            f = _new(a.get_facecolor(), 0)
            e = _new(a.get_edgecolor(), 1)
            if f is not None:
                a.set_facecolor(f)
                if e is None and not a.get_hatch() and (a.get_linewidth() == 0 or a.get_edgecolor()[3] == 0):
                    a.set_edgecolor(_shade(f))
                    a.set_linewidth(EDGE_LW)
            if e is not None:
                a.set_edgecolor(e)
        elif isinstance(a, Collection):
            if a.get_array() is not None:
                if isinstance(a, PolyCollection) and a.get_cmap().name == 'Greys':
                    a.set_cmap('viridis_r')
                continue
            if isinstance(a, LineCollection):
                c = _arr(a.get_colors(), 1)
                if c is not None:
                    a.set_colors(c)
                continue
            fc_old = [tuple(c) for c in a.get_facecolors()]
            ec_old = [tuple(c) for c in a.get_edgecolors()]
            fc = _arr(fc_old, 0)
            if fc is not None:
                a.set_facecolors(fc)
            ec = _arr(ec_old, 1)
            thin = max(a.get_linewidths(), default=0) == 0
            sparse = len(a.get_offsets()) < 500           # outline markers of sparse scatters only, not dense clouds
            if ec is not None:
                a.set_edgecolors(ec)
            elif fc is not None and (len(ec_old) == 0 or ec_old == fc_old or thin):
                a.set_edgecolors([_shade(c) for c in fc])
            if fc is not None and thin and sparse and (ec is not None or len(ec_old) == 0 or ec_old == fc_old or thin):
                a.set_linewidths(EDGE_LW)
        elif isinstance(a, Text):
            n = _new(a.get_color(), 2)
            if n is not None:
                a.set_color(n)


def _shade(rgba):
    """Outline a shade darker than a light fill (the line tone of its role when known)."""
    for f, l, _ in LIGHT.values():
        if _key(f) == _key(rgba):
            return mcolors.to_rgba(l, 1)
    r, g, b = mcolors.to_rgb(rgba)
    return (0.8 * r, 0.8 * g, 0.8 * b, 1)


def save(fig, name):
    recolor(fig)
    for ext in ('pdf', 'png'):
        fig.savefig(os.path.join(HERE, f'{name}.{ext}'))


def half_violins(ax, groups, colors, labels, seed=0, max_points=None):
    """Distribution of each group as the left half of a violin, the values as points to its right and the median as a
    dark line across both (ALBATROSS / IDRome style; chosen by the user 2026-10-14 over box plots and median bars).
    max_points draws a random subset of the points for large groups (the violin and median use all values)."""
    rng = np.random.default_rng(seed)
    for k, (v, c) in enumerate(zip(groups, colors)):
        v = np.asarray(v, float)
        v = v[np.isfinite(v)]
        vp = ax.violinplot(v, positions=[k - 0.08], widths=0.7, showextrema=False, bw_method=0.5)
        for body in vp['bodies']:
            xy = body.get_paths()[0].vertices
            xy[:, 0] = np.minimum(xy[:, 0], k - 0.08)
            body.set_facecolor(c)
            body.set_edgecolor(c)
            body.set_linewidth(0.6)
            body.set_alpha(0.75)
        pv = v if max_points is None or len(v) <= max_points else rng.choice(v, max_points, replace=False)
        sz = 10 if len(pv) <= 60 else 3
        ax.scatter(k + 0.14 + rng.uniform(-0.05, 0.05, len(pv)), pv, s=sz, color=c, edgecolor=c, lw=0.5 if sz > 5 else 0.2,
                   zorder=3)
        ax.plot([k - 0.3, k + 0.24], [np.median(v)] * 2, color='#333333', lw=1.0, zorder=4)
    ax.set_xticks(range(len(groups)), labels)
    ax.set_xlim(-0.55, len(groups) - 0.5)


def cap(title):
    """Capital first letter for a panel header whose first word is plain lowercase (keeps hnRNPA1, p53, pSic1, alpha-Syn)."""
    w = re.match(r'[^\s,;:()]+', title)
    if w and w.group().isascii() and w.group().isalpha() and w.group().islower():
        return title[0].upper() + title[1:]
    return title


def panel(ax, letter, x=-0.16, y=1.04, title=None):
    """Nature-style panel letter (bold lowercase) at the top left of an axis, with an optional short header 11 points to
    its right (ALBATROSS and IDRome put the explanation in the caption, so headers stay a few words long)."""
    ax.annotate(letter, (x, y), xycoords='axes fraction', fontsize=9, fontweight='bold', ha='left', va='baseline')
    if title:
        ax.annotate(cap(title), (x, y), xycoords='axes fraction', xytext=(10, 0), textcoords='offset points', fontsize=6.5,
                    ha='left', va='baseline', color='#333333')
