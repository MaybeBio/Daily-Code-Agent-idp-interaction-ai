"""
slides_pptx.py — turn a matplotlib slide figure into an editable PowerPoint slide (python-pptx, 2026-10-30).
The figure has one full-slide background axes in slide units of 0.1 in (fig_slides_deck.new_slide) and any number of
plot axes. On export
  * rounded boxes (FancyBboxPatch) on the background become native rounded rectangles,
  * step badges (Circle of radius 1.25 filled with ink) become native circles,
  * every text on the background and every axes title becomes a native text box, placed on the rendered extent of
    the matplotlib text, with its font size, weight, style, colour and alignment (text with mathtext stays drawn),
  * each plot axes (with its labels, ticks and legend) becomes its own transparent picture, and what is left on the
    background (cartoon chains, arrows, figure legends) one full-slide transparent picture.
Stacking: boxes, background picture, plot pictures, badges, text. Colours are read after the figure has been saved
once, so figfont's palette and capitalisation are already applied.
"""
import io

import matplotlib.colors as mcolors
import numpy as np
from matplotlib.patches import Circle, FancyBboxPatch
from matplotlib.transforms import Bbox
from pptx import Presentation
from pptx.dml.color import RGBColor
from pptx.enum.shapes import MSO_SHAPE
from pptx.enum.text import MSO_ANCHOR, PP_ALIGN
from pptx.oxml.ns import qn
from pptx.util import Emu, Inches, Pt

DPI = 250
EMU_IN = 914400
ALIGN = {'left': PP_ALIGN.LEFT, 'center': PP_ALIGN.CENTER, 'right': PP_ALIGN.RIGHT}


def _rgb(c):
    return RGBColor.from_string(mcolors.to_hex(c, keep_alpha=False)[1:].upper())


def _plain(sh):
    """Drop the theme style reference of an autoshape (no shadow or theme outline in any viewer)."""
    st = sh._element.find(qn('p:style'))
    if st is not None:
        sh._element.remove(st)
    sh.line.fill.background()


def _emu(inches):
    return Emu(int(round(inches * EMU_IN)))


class Deck:
    def __init__(self, w, h):
        self.prs = Presentation()
        self.prs.slide_width, self.prs.slide_height = Inches(w), Inches(h)
        self.w, self.h = w, h
        self.blank = self.prs.slide_layouts[6]

    # ------------------------------------------------------------------ collect the native items
    def _collect(self, fig, bg, r):
        ux = self.w / bg.get_xlim()[1]                     # inches per background unit
        uy = self.h / bg.get_ylim()[1]
        boxes, badges, texts = [], [], []
        for p in list(bg.patches):
            if isinstance(p, FancyBboxPatch):
                bs = p.get_boxstyle()
                rad = getattr(bs, 'rounding_size', None) or 0
                boxes.append(dict(x=p.get_x() * ux, y=self.h - (p.get_y() + p.get_height()) * uy, w=p.get_width() * ux,
                                  h=p.get_height() * uy, r=rad * ux, fc=p.get_facecolor()))
                p.remove()
            elif isinstance(p, Circle) and abs(p.get_radius() - 1.25) < 1e-9:
                cx, cy = p.center
                badges.append(dict(x=(cx - 1.25) * ux, y=self.h - (cy + 1.25) * uy, d=2.5 * ux, fc=p.get_facecolor()))
                p.remove()
        tx = [t for t in bg.texts if t.get_visible() and t.get_text() and '$' not in t.get_text()]
        for a in fig.axes[1:]:
            for t in (a.title, a._left_title, a._right_title):
                if t.get_text() and '$' not in t.get_text():
                    tx.append(t)
        for t in tx:
            bb = t.get_window_extent(r)
            texts.append(dict(x=bb.x0 / fig.dpi, y=self.h - bb.y1 / fig.dpi, w=bb.width / fig.dpi, h=bb.height / fig.dpi,
                              text=t.get_text(), size=t.get_fontsize(), bold=t.get_fontweight() in ('bold', 700, 'heavy'),
                              italic=t.get_fontstyle() == 'italic', color=t.get_color(), ha=t.get_horizontalalignment(),
                              ls=t._linespacing if isinstance(t._linespacing, (int, float)) else 1.2))
            t.set_text('')
        return boxes, badges, texts

    # ------------------------------------------------------------------ pictures
    def _png(self, fig, bbox=None):
        buf = io.BytesIO()
        kw = dict(dpi=DPI, transparent=True)
        if bbox is not None:
            kw.update(bbox_inches=bbox, pad_inches=0)
        fig.savefig(buf, **kw)
        buf.seek(0)
        return buf

    def add(self, fig, bg, name=''):
        fig.set_dpi(100)
        fig.canvas.draw()
        r = fig.canvas.get_renderer()
        boxes, badges, texts = self._collect(fig, bg, r)
        s = self.prs.slides.add_slide(self.blank)
        for b in boxes:
            sh = s.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE, _emu(b['x']), _emu(b['y']), _emu(b['w']), _emu(b['h']))
            sh.adjustments[0] = min(0.5, b['r'] / max(1e-6, min(b['w'], b['h'])))
            sh.fill.solid()
            sh.fill.fore_color.rgb = _rgb(b['fc'])
            _plain(sh)
        others = fig.axes[1:]
        legends = list(fig.legends)
        # background picture: cartoons, arrows and figure legends
        for a in others:
            a.set_visible(False)
        s.shapes.add_picture(self._png(fig), 0, 0, self.prs.slide_width, self.prs.slide_height).name = 'cartoons'
        for lg in legends:
            lg.set_visible(False)
        bg.set_visible(False)
        # one picture per plot axes, cropped to its tight extent
        fig.canvas.draw()
        for k, a in enumerate(others):
            for o in others:
                o.set_visible(o is a)
            fig.canvas.draw()
            bb = a.get_tightbbox(fig.canvas.get_renderer())
            ib = Bbox.from_extents(bb.x0 / fig.dpi - 0.03, bb.y0 / fig.dpi - 0.03, bb.x1 / fig.dpi + 0.03,
                                   bb.y1 / fig.dpi + 0.03)
            pic = s.shapes.add_picture(self._png(fig, ib), _emu(ib.x0), _emu(self.h - ib.y1), _emu(ib.width),
                                       _emu(ib.height))
            pic.name = f'plot {k + 1}'
        for a in others + [bg]:
            a.set_visible(True)
        for lg in legends:
            lg.set_visible(True)
        for b in badges:
            sh = s.shapes.add_shape(MSO_SHAPE.OVAL, _emu(b['x']), _emu(b['y']), _emu(b['d']), _emu(b['d']))
            sh.fill.solid()
            sh.fill.fore_color.rgb = _rgb(b['fc'])
            _plain(sh)
        for t in texts:
            pad = 0.04 * t['size'] / 12
            tb = s.shapes.add_textbox(_emu(t['x'] - pad), _emu(t['y'] - pad), _emu(t['w'] + 2 * pad),
                                      _emu(t['h'] + 2 * pad))
            tf = tb.text_frame
            tf.word_wrap = False
            tf.auto_size = None
            tf.margin_left = tf.margin_right = tf.margin_top = tf.margin_bottom = 0
            tf.vertical_anchor = MSO_ANCHOR.MIDDLE
            for i, line in enumerate(t['text'].split('\n')):
                p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
                p.alignment = ALIGN.get(t['ha'], PP_ALIGN.LEFT)
                p.line_spacing = float(np.clip(t['ls'] / 1.2, 0.8, 1.6))
                run = p.add_run()
                run.text = line
                f = run.font
                f.name = 'Arial'
                f.size = Pt(t['size'])
                f.bold = t['bold']
                f.italic = t['italic']
                f.color.rgb = _rgb(t['color'])
        return s

    def save(self, path):
        self.prs.save(path)
        return path
