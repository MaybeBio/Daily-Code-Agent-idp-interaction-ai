"""
figfont.py — shared style for every manuscript and note figure.
  Font     Arial for text and math labels (mathtext), TrueType embedding in PDF.
  Palette  the figure scripts keep their own hex colours; at save time every colour of the original (muted) palette is
           replaced by its counterpart in the active palette, and a few colormaps by stronger ones. The default is
           'petroff6' (Petroff 2021, matplotlib's accessible colour cycle: blue CALVADOS-2, red Mpipi-GG, purple
           Mpipi-Recharged, orange for the charged-group form and hydrophobic forcing; chosen 2026-10-24 because every pair of
           the four main colours stays distinct under simulated deuteranopia, protanopia and tritanopia; previously
           'f4p-green' after github.com/ChenLiu-1996/figures4papers, 'soft' and 'okabe-ito'); another
           can be set with set_palette(name) or FIG_PALETTE ('current' keeps the original muted colours).
  Labels   at save time every text that starts with a lowercase word is capitalised (titles, legends, annotations,
           axis and colour-bar labels, categorical tick labels); single letters (panel letters), names and symbols
           whose first word carries a capital or a digit (sfAFP, smFRET, p53, pH), d$x$/d$T$ and log are kept
           (added 2026-10-20; FIG_CAPITALIZE=0 switches it off).
Call apply() after a script's own rcParams.update; it also installs the save-time recolouring (once).
"""
import os
import re

import matplotlib
import matplotlib.image
from matplotlib import colors as mcolors
from matplotlib.collections import Collection, LineCollection, PathCollection
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.text import Text
from matplotlib.ticker import FixedFormatter, FuncFormatter

ARIAL = {'font.family': 'Arial', 'font.weight': 'normal',
         'mathtext.fontset': 'custom', 'mathtext.rm': 'Arial', 'mathtext.it': 'Arial:italic',
         'mathtext.bf': 'Arial:bold', 'mathtext.sf': 'Arial', 'mathtext.fallback': 'stixsans',
         # Greek nu: Arial's nu is indistinguishable from a Latin v, so $\mathtt{\nu}$ draws it from STIX italic
         'mathtext.tt': 'STIXGeneral:italic',
         'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'path'}

# original colour -> role, so that every palette is defined by role
ROLE = {'#2E6FA3': 'cal', '#3E6FA3': 'cal', '#A8CBE8': 'cal_pale', '#7C8BA0': 'slate', '#6F8199': 'charge',
        '#A9502A': 'gg', '#E8B9A0': 'gg_pale', '#B5473A': 'red',
        '#6B4FA0': 'rch', '#D7CDEB': 'rch_pale',
        '#8F5E14': 'gold', '#B07D2B': 'gold2', '#C9A86A': 'gold_light', '#E3C98F': 'gold_pale', '#D9CDB8': 'beige',
        '#C8792A': 'orange', '#2E8C6F': 'teal', '#3E8E7E': 'teal2',
        '#E2DCCF': 'grid', '#F7F4EE': 'bg_warm', '#EFEBE3': 'bg_warm2', '#C8C2B4': 'bond', '#B9B2A4': 'beige_grey',
        '#9AA3AE': 'grey1', '#8A8F98': 'grey2', '#B9C3CF': 'grey3'}
PALETTES = {
    'current': {},
    'vivid': {'cal': '#0A6CF1', 'cal_pale': '#8EC5FF', 'slate': '#4C6EF5', 'charge': '#339AF0',
              'gg': '#F2541B', 'gg_pale': '#FFB59A', 'red': '#E03131',
              'rch': '#8F3BE8', 'rch_pale': '#D6B8FF',
              'gold': '#F08C00', 'gold2': '#F59F00', 'gold_light': '#FFC94D', 'gold_pale': '#FFE08A', 'beige': '#DEE2E6',
              'orange': '#FF8A00', 'teal': '#00B37E', 'teal2': '#12B886',
              'grid': '#E9ECEF', 'bg_warm': '#F8F9FA', 'bg_warm2': '#F1F3F5', 'bond': '#CED4DA', 'beige_grey': '#ADB5BD',
              'grey1': '#868E96', 'grey2': '#868E96', 'grey3': '#ADB5BD'},
    'okabe-ito': {'cal': '#0072B2', 'cal_pale': '#9ED0F0', 'slate': '#3B6FB6', 'charge': '#56B4E9',
                  'gg': '#D55E00', 'gg_pale': '#F6B98C', 'red': '#C0392B',
                  'rch': '#6A3D9A', 'rch_pale': '#CAB2D6',   # violet (was magenta #C5157F / #F0B3D6 until 2026-10-19)
                  'gold': '#E69F00', 'gold2': '#E69F00', 'gold_light': '#F3C766', 'gold_pale': '#F8E0A6', 'beige': '#DEE2E6',
                  'orange': '#E69F00', 'teal': '#009E73', 'teal2': '#009E73',
                  'grid': '#E9ECEF', 'bg_warm': '#F8F9FA', 'bg_warm2': '#F1F3F5', 'bond': '#CED4DA', 'beige_grey': '#ADB5BD',
                  'grey1': '#868E96', 'grey2': '#868E96', 'grey3': '#ADB5BD'}}
PALETTES['cream-pastel'] = {'bg': '#FBF6EC', 'cal': '#6E9FD6', 'cal_pale': '#C9DCF1', 'slate': '#8EA6C8', 'charge': '#8FB8E0',
                            'gg': '#EE9D63', 'gg_pale': '#F8D5BC', 'red': '#E58A7E',
                            'rch': '#C08BCB', 'rch_pale': '#E6D2EA',
                            'gold': '#E9B949', 'gold2': '#EDC25E', 'gold_light': '#F4D78E', 'gold_pale': '#F8E7BC',
                            'beige': '#EADFCB', 'orange': '#F0A964', 'teal': '#7FC4A7', 'teal2': '#8CCBB0',
                            'grid': '#EDE3D1', 'bg_warm': '#F6EFE2', 'bg_warm2': '#F2E9D8', 'bond': '#DCCFB8',
                            'beige_grey': '#C9BCA4', 'grey1': '#A39A8B', 'grey2': '#A39A8B', 'grey3': '#C2B8A6'}
PALETTES['cream-jewel'] = {'bg': '#FBF6EC', 'cal': '#1F5FA8', 'cal_pale': '#A9C4E6', 'slate': '#3F6496', 'charge': '#4F86C6',
                           'gg': '#D2581F', 'gg_pale': '#F2BE9F', 'red': '#B8322A',
                           'rch': '#B0307A', 'rch_pale': '#E9B5D2',
                           'gold': '#D9921A', 'gold2': '#E0A22E', 'gold_light': '#F0C870', 'gold_pale': '#F6DFA8',
                           'beige': '#E6D9C2', 'orange': '#E08A1E', 'teal': '#1F8A70', 'teal2': '#2A9478',
                           'grid': '#EBDFCB', 'bg_warm': '#F6EFE2', 'bg_warm2': '#F2E9D8', 'bond': '#D8CAB2',
                           'beige_grey': '#BFB19A', 'grey1': '#968C7C', 'grey2': '#968C7C', 'grey3': '#BBB09E'}
NEUTRAL = {k: PALETTES['okabe-ito'][k] for k in ('beige', 'grid', 'bg_warm', 'bg_warm2', 'bond', 'beige_grey', 'grey1', 'grey2', 'grey3')}
# muted, mid-tone palettes in the manner of IDRome / ALBATROSS / FINCHES (2026-10-19; 'soft' is the default)
PALETTES['muted'] = {**NEUTRAL, 'cal': '#4C72B0', 'cal_pale': '#B5C7E3', 'slate': '#5A78A8', 'charge': '#7FB2D5',
                     'gg': '#DD8452', 'gg_pale': '#F2C9AD', 'red': '#C44E52',
                     'rch': '#8172B3', 'rch_pale': '#CEC6E3',
                     'gold': '#CCA43F', 'gold2': '#CCA43F', 'gold_light': '#E3CB86', 'gold_pale': '#F0E2B8',
                     'orange': '#D9A441', 'teal': '#55A868', 'teal2': '#55A868'}
PALETTES['soft'] = {**NEUTRAL, 'cal': '#3E8FA3', 'cal_pale': '#B4D8E0', 'slate': '#4F7FA0', 'charge': '#86BCD6',
                    'gg': '#E3926A', 'gg_pale': '#F5CDB8', 'red': '#C65A52',
                    'rch': '#6E6FA8', 'rch_pale': '#C9C9E2',
                    'gold': '#D2A03E', 'gold2': '#D2A03E', 'gold_light': '#E6C885', 'gold_pale': '#F2E1B6',
                    'orange': '#D9A441', 'teal': '#7FAE5E', 'teal2': '#7FAE5E'}
# candidates after figures4papers (github.com/ChenLiu-1996/figures4papers): navy anchor, brick red, teal/green (2026-10-20)
PALETTES['f4p-teal'] = {**NEUTRAL, 'cal': '#0F4D92', 'cal_pale': '#AFC6E2', 'slate': '#3775BA', 'charge': '#5B8FCB',
                        'gg': '#B64342', 'gg_pale': '#E9A6A1', 'red': '#B64342',
                        'rch': '#42949E', 'rch_pale': '#B3D7DB',
                        'gold': '#C9A13B', 'gold2': '#C9A13B', 'gold_light': '#E2C77E', 'gold_pale': '#F1E3B8',
                        'orange': '#C9A13B', 'teal': '#5DAA5D', 'teal2': '#5DAA5D'}
PALETTES['f4p-green'] = {**NEUTRAL, 'cal': '#0F4D92', 'cal_pale': '#AFC6E2', 'slate': '#3775BA', 'charge': '#5B8FCB',
                         'gg': '#B64342', 'gg_pale': '#E9A6A1', 'red': '#B64342',
                         'rch': '#5DAA5D', 'rch_pale': '#C2E3C1',
                         'gold': '#C9A13B', 'gold2': '#C9A13B', 'gold_light': '#E2C77E', 'gold_pale': '#F1E3B8',
                         'orange': '#C9A13B', 'teal': '#42949E', 'teal2': '#42949E'}
# Petroff (2021) accessible colour cycle, matplotlib style 'petroff6' (blue, red, purple, orange for the four main roles;
# pale tints mixed with white; teal kept from f4p), added and made the default 2026-10-24
PALETTES['petroff6'] = {**NEUTRAL, 'cal': '#5790FC', 'cal_pale': '#BCD3FE', 'slate': '#4673CA', 'charge': '#89B1FD', 'gg': '#E42536', 'gg_pale': '#F4A8AF', 'red': '#E42536', 'rch': '#964A8B', 'rch_pale': '#D5B7D1', 'gold': '#F89C20', 'gold2': '#F89C20', 'gold_light': '#FBC479', 'gold_pale': '#FDE1BC', 'orange': '#F89C20', 'teal': '#42949E', 'teal2': '#42949E'}
# after Adduri et al., Cell 189, 5914 (2026) ('State'): steel blue hero, their icon olive, amber and tan, their greys and
# terracotta; adapted for three force fields by lightness (blue 51, olive 30, amber 54, tan 85), colour-blind gaps 32/25/36
# (deutan/protan/tritan, CIELAB, Machado 2009); chosen by the user 2026-10-24. The tan role (gold, orange) uses #E5C27A,
# a deeper version of their #F1D18A, for legibility of text and small markers on white (gaps then 27/25/31)
PALETTES['state'] = {'cal': '#387FB8', 'cal_pale': '#AFCCE3', 'slate': '#2A5F8A', 'charge': '#74A5CD', 'gg': '#3A4C12', 'gg_pale': '#B0B7A0', 'red': '#CA534E', 'rch': '#B47327', 'rch_pale': '#E1C7A9', 'gold': '#E5C27A', 'gold2': '#E5CD90', 'gold_light': '#F6E5C3', 'gold_pale': '#FAF2E1', 'orange': '#E5C27A', 'teal': '#AC7D5C', 'teal2': '#AC7D5C', 'beige': '#E0E0E0', 'grid': '#E5E5E5', 'bg_warm': '#F5F5F5', 'bg_warm2': '#EFEFEF', 'bond': '#C7C7C7', 'beige_grey': '#B5B5B5', 'grey1': '#919191', 'grey2': '#636363', 'grey3': '#C7C7C7'}
# lighter version on the State paper's own principle (blue-yellow axis plus a grey ladder): light steel blue CALVADOS-2,
# dark grey Mpipi-GG, light amber Mpipi-Recharged, their terracotta for the fourth role; colour-blind gap 39 for the three
# force fields (all deficiencies), 27/27/27/19 with terracotta (2026-10-24). Later that night, on the user's request, Mpipi-GG
# softened to mid grey #858585 (gaps 30 to the force fields, 23 to terracotta) and the neutral greys lightened (plain #BDBDBD)
# so that the plain models stay distinct from Mpipi-GG
PALETTES['state-light'] = {'cal': '#6FA3D2', 'cal_pale': '#BED6EB', 'slate': '#537A9E', 'charge': '#9ABFE0', 'gg': '#C5E0C4', 'gg_line': '#3A6E38', 'gg_pale': '#E2F0E1', 'red': '#CA534E', 'rch': '#E3A35A', 'rch_pale': '#F2D6B5', 'gold': '#AC7D5C', 'gold2': '#B89074', 'gold_light': '#D1B8A5', 'gold_pale': '#E6D8CE', 'orange': '#AC7D5C', 'teal': '#7F9E8F', 'teal2': '#7F9E8F', 'beige': '#E0E0E0', 'grid': '#E5E5E5', 'bg_warm': '#F5F5F5', 'bg_warm2': '#EFEFEF', 'bond': '#D0D0D0', 'beige_grey': '#C2C2C2', 'grey1': '#BDBDBD', 'grey2': '#ABABAB', 'grey3': '#D0D0D0'}
CMAPS = {'current': {}, 'vivid': {'Greys': 'Blues', 'BrBG_r': 'RdBu_r', 'BrBG': 'RdBu'},
         'cream-pastel': {'Greys': 'YlOrBr'}, 'cream-jewel': {'Greys': 'YlOrBr', 'BrBG_r': 'PuOr_r', 'BrBG': 'PuOr'},
         'okabe-ito': {'Greys': 'Blues', 'BrBG_r': 'PuOr_r', 'BrBG': 'PuOr'},
         'muted': {'Greys': 'Blues', 'BrBG_r': 'PuOr_r', 'BrBG': 'PuOr'}, 'soft': {'Greys': 'Blues', 'BrBG_r': 'PuOr_r', 'BrBG': 'PuOr'},
         'f4p-teal': {'Greys': 'Blues', 'BrBG_r': 'PuOr_r', 'BrBG': 'PuOr'}, 'f4p-green': {'Greys': 'Blues', 'BrBG_r': 'PuOr_r', 'BrBG': 'PuOr'},
         'petroff6': {'Greys': 'Blues', 'BrBG_r': 'PuOr_r', 'BrBG': 'PuOr'},
         'state': {'Greys': 'Blues', 'BrBG_r': 'RdBu_r', 'BrBG': 'RdBu'},
         'state-light': {'Greys': 'Blues', 'BrBG_r': 'RdBu_r', 'BrBG': 'RdBu'}}

_active = os.environ.get('FIG_PALETTE', 'state-light')  # lighter State-paper scheme (chosen 2026-10-24 night; 'state' and
#                                                        'petroff6' earlier that day, 'f4p-green' 2026-10-20, 'soft' 2026-10-19,
#                                                        'okabe-ito' 2026-10-11)
_installed = False


def set_palette(name):
    global _active
    if name not in PALETTES:
        raise ValueError(f'unknown palette {name}')
    _active = name


def _key(c):
    return tuple(round(v, 3) for v in mcolors.to_rgb(c))


def _table(stroke=False):
    """Old colour -> active colour. stroke=True (lines, markers, outlines, text) takes a role's '<role>_line' shade where
    the palette has one (state-light Mpipi-GG: the user's pale green for fills, a deeper green for thin marks)."""
    pal = PALETTES[_active]
    return {_key(old): pal[f'{role}_line'] if stroke and f'{role}_line' in pal else pal[role]
            for old, role in ROLE.items() if role in pal}


def _swap(c, table):
    """New RGBA for one colour, keeping its alpha; None when unchanged."""
    if c is None or (isinstance(c, str) and c in ('none', 'None', 'auto')):
        return None
    try:
        rgba = mcolors.to_rgba(c)
    except (ValueError, TypeError):
        return None
    new = table.get(_key(rgba))
    return None if new is None else mcolors.to_rgba(new, rgba[3])


def _swap_array(arr, table):
    out, changed = [], False
    for c in arr:
        n = _swap(tuple(c), table)
        changed |= n is not None
        out.append(n if n is not None else tuple(c))
    return out if changed else None


def recolor(fig):
    """Map every colour of the original palette in a figure onto the active palette (and the listed colormaps)."""
    if _active == 'current':
        return
    table, stroke, cmaps = _table(), _table(stroke=True), CMAPS[_active]
    bg = PALETTES[_active].get('bg')
    if bg:
        white = _key('#FFFFFF')
        if _key(fig.get_facecolor()) == white:
            fig.set_facecolor(bg)
        for ax in fig.axes:
            if _key(ax.get_facecolor()) == white:
                ax.set_facecolor(bg)
    for a in fig.findobj():
        if isinstance(a, Line2D):
            for get, put in (('get_color', 'set_color'), ('get_markerfacecolor', 'set_markerfacecolor'),
                             ('get_markeredgecolor', 'set_markeredgecolor')):
                n = _swap(getattr(a, get)(), stroke)
                if n is not None:
                    getattr(a, put)(n)
        elif isinstance(a, Patch):
            fc = a.get_facecolor()                                           # an edge in the face colour of a filled patch stays a
            same = a.get_fill() and fc[3] > 0 and _key(fc) == _key(a.get_edgecolor())   # fill; outlines (step histograms) are strokes
            for get, put, tb in (('get_facecolor', 'set_facecolor', table),
                                 ('get_edgecolor', 'set_edgecolor', table if same else stroke)):
                n = _swap(getattr(a, get)(), tb)
                if n is not None:
                    getattr(a, put)(n)
        elif isinstance(a, Collection):
            cm = a.get_cmap() if a.get_array() is not None else None
            if cm is not None and cm.name in cmaps:
                a.set_cmap(cmaps[cm.name])
            if a.get_array() is None:                      # scatter markers are thin marks; areas (violins, fills) are fills
                fc = _swap_array(a.get_facecolors(), stroke if isinstance(a, PathCollection) else table)
                if fc is not None:
                    a.set_facecolors(fc)
            ec = _swap_array(a.get_edgecolors(), stroke)
            if ec is not None:
                (a.set_colors if isinstance(a, LineCollection) else a.set_edgecolors)(ec)
        elif isinstance(a, Text):
            n = _swap(a.get_color(), stroke)
            if n is not None:
                a.set_color(n)
        elif isinstance(a, matplotlib.image.AxesImage):
            cm = a.get_cmap()
            if cm.name in cmaps:
                a.set_cmap(cmaps[cm.name])


_CAPITALIZE = os.environ.get('FIG_CAPITALIZE', '1') != '0'


def _cap(t):
    """Capitalise a label that starts with a lowercase word; keep panel letters, mixed-case names and symbols whose
    first word carries a capital or a digit (sfAFP, smFRET, p53, pH), derivatives d$x$/d$T$, and log."""
    if not isinstance(t, str) or len(t) < 2 or not ('a' <= t[0] <= 'z') or not ('a' <= t[1] <= 'z'):
        return t
    first = re.split(r'[\s,(/+]', t, maxsplit=1)[0]
    if t.startswith('log') or any(c.isupper() or c.isdigit() for c in first):
        return t
    return t[0].upper() + t[1:]


def capitalize(fig):
    """Capitalise every label of a figure in place (texts, legends, axis labels and fixed tick labels)."""
    for a in fig.findobj(Text):
        n = _cap(a.get_text())
        if n != a.get_text():
            a.set_text(n)
    for ax in fig.axes:
        for axis in (ax.xaxis, ax.yaxis):
            f = axis.get_major_formatter()
            if isinstance(f, FixedFormatter):
                f.seq = [_cap(x) for x in f.seq]
            elif isinstance(f, FuncFormatter) and not getattr(f, '_capitalized', False):   # string tick labels (mpl >= 3.10)
                f.func = (lambda g: (lambda x, pos=None: _cap(g(x, pos))))(f.func)
                f._capitalized = True


def apply():
    global _installed
    matplotlib.rcParams.update(ARIAL)
    if not _installed:
        orig = Figure.savefig

        def savefig(self, *args, **kwargs):
            recolor(self)
            if _CAPITALIZE:
                capitalize(self)
            bg = PALETTES[_active].get('bg')
            if bg and ('facecolor' not in kwargs or _key(kwargs['facecolor']) == _key('#FFFFFF')):
                kwargs['facecolor'] = bg
            return orig(self, *args, **kwargs)
        Figure.savefig = savefig
        _installed = True
