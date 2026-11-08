"""
Scene 24 -- "steady-state plots showing sharp transition" (storyboard note).

Panels A and B of the paper's Fig. 2, one over the other as published:
histograms, one column per background concentration c, of the probe-enhancer
distance (A) and of the local concentration (B), with their means (dots).
Around c = 35 uM the distance collapses to ~3 sites and the local
concentration jumps. Above A, a legend for the dots; right of each panel, a
cartoon of what it plots:
    A  scene 20's chain, its promoter and enhancer joined by a red arrow,
       "distance"
    B  scene 21's locus, a droplet of TFs condensed at its promoter, ringed by
       a dashed circle

The two panels are cropped together from the published figure (slides/figures/
scene24/fig2.pdf, a copy of the paper's arXiv/fig2.pdf) and placed by pdflatex,
as in scene 23, so they stay vector and keep their own lettering. The panel
letters are left out. The legend and cartoons are drawn with matplotlib on a
transparent slide that pdflatex lays over them.

Output: slides/figures/scene24/scene24.pdf, a full 128 x 96 mm slide.

Usage (from slides/):
    python3 scripts/scene24/scene24.py
Requires the `pdflatex` command-line tool.
"""

import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # slides/scripts, for `common`

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Circle

from common.drawing import (FIGURES_DIR, LABEL_PT, SLIDE_MM, SLIDE_UNITS, TEXT_COLOR, UNIT_MM,
                            draw_arrow, new_slide, zoomed_axes)
from common.molecules import DNA_LINE_PT, draw_dna, draw_legend
from scene20.scene20 import DIM_DNA, SEGMENTS, STRETCHED, locus
from scene21.scene21 import DROPLET_CENTER, draw_locus
from scene22.scene22 import MASK, NODE, PT_MM, typeset

FIG2 = FIGURES_DIR / "scene24" / "fig2.pdf"

# Panels A and B, a viewport of fig2.pdf (x0, y0, x1, y1 in points, y up): from
# the y-labels to the colour bars' labels, and from the top of A to the
# x-label under B. A's letter, above its y-label, is outside it.
PANELS_BOX = (1.6, 2.4, 403.5, 355.4)
# B's letter stands between the panels, above B's y-label: it is painted white
# (box in points, as PANELS_BOX; nothing else is drawn there). The box clears
# the letter by 2 pt and overhangs the viewport's left edge, or rasterizers
# leave a speck of its anti-aliased edges.
B_LETTER_MASK = (0.6, 181.4, 13.0, 202.4)
# The plot areas (points, as PANELS_BOX): A and B share their left and right
# spines; bottom and top spines of each.
SPINES_X = (49.0, 356.0)
SPINES_Y = dict(a=(206.6, 350.4), b=(36.4, 180.1))
# The means' marker as Fig. 2 draws it: a black dot 4.4 pt across with a 1-pt
# white edge. Fig. 2's legends are set in 12 pt of Computer Modern, which 10 pt
# of FreeSans matches.
MEAN_MARKER = dict(size=4.4, edge=1.0)
LEGEND_PT = 10

# The panels are drawn at SCALE; right of them, after GAP, is a column COLUMN
# wide for the cartoons. The row is centred, the panels and LEGEND_ROOM above
# them centred vertically; the legend sits there, at the right end of A's plot
# area (slide units).
SCALE = 0.65
GAP = 0.3
COLUMN = 4.6
LEGEND_ROOM = 0.6

# Each cartoon is centred in the column beside its panel's plot area: drawn
# `zoom` times the size of its drawing (in scene-4 coordinates), the drawing's
# point `center` lands at the middle of its slot.
DISTANCE = dict(zoom=0.4, center=(6.9, 13.3))      # scene 20's stretched chain
CONDENSATE = dict(zoom=0.42, center=(5.4, 9.9),    # scene 21's locus,
                  lift=0.75)                       # raised to leave room for the legend
# The arrow from promoter to enhancer, stopping `gap` (scene-4 units) short of
# each, in the red of scenes 2-3; its label beside its middle.
ARROW = dict(color="#D93A45", lw=1.5, head=9, gap=0.4)
# Around the droplet (scene-4 units, as scene 21's DROPLET).
RING = dict(radius=1.1, lw=1.0, dashes=(3, 2))
LEGEND_AT = (14.75, 2.2)  # in the corner under the cartoons (slide units)


def draw_distance(fig, center):
    """Scene 20's chain, its promoter and enhancer joined by the red arrow,
    centred at `center` (slide units)."""
    zoom = DISTANCE["zoom"]
    ax = zoomed_axes(fig, np.asarray(center) - zoom * np.asarray(DISTANCE["center"]), zoom)
    dna = locus(**STRETCHED).dna
    draw_dna(ax, dna, SEGMENTS, z=1, line_pt=tuple(zoom * w for w in DNA_LINE_PT), **DIM_DNA)
    ends = [dna.at(sum(dna.span(kind)) / 2)[0] for kind in ("promoter", "enhancer")]
    middle = (ends[0] + ends[1]) / 2
    for end in ends:
        tip = end + ARROW["gap"] * (middle - end) / np.linalg.norm(middle - end)
        draw_arrow(ax, middle, tip, ARROW["color"], z=10, lw=ARROW["lw"], head=ARROW["head"])
    ax.text(*(middle - (0.35, 0)), "distance", ha="right", va="center", fontsize=LABEL_PT,
            color=ARROW["color"], zorder=11)


def draw_condensate(fig, center):
    """Scene 21's locus, its droplet ringed, centred at `center` (slide units)."""
    zoom = CONDENSATE["zoom"]
    ax = draw_locus(fig, np.asarray(center) - zoom * np.asarray(CONDENSATE["center"]), zoom)
    ax.add_patch(Circle(DROPLET_CENTER, RING["radius"], fill=False, ec=TEXT_COLOR,
                        lw=RING["lw"], ls=(0, RING["dashes"]), zorder=7))


if __name__ == "__main__":
    x0, y0, x1, y1 = PANELS_BOX
    pt = SCALE * PT_MM / UNIT_MM  # slide units per point of fig2.pdf
    left = (SLIDE_UNITS[0] - (x1 - x0) * pt - GAP - COLUMN) / 2
    top = (SLIDE_UNITS[1] + (y1 - y0) * pt + LEGEND_ROOM) / 2 - LEGEND_ROOM

    def to_slide(x, y):
        """Slide position (slide units) of a point of fig2.pdf (points)."""
        return left + (x - x0) * pt, top - (y1 - y) * pt

    fig, ax = new_slide()
    column_x = left + (x1 - x0) * pt + GAP + COLUMN / 2
    draw_distance(fig, (column_x, to_slide(0, sum(SPINES_Y["a"]) / 2)[1]))
    draw_condensate(fig, (column_x, to_slide(0, sum(SPINES_Y["b"]) / 2)[1] + CONDENSATE["lift"]))
    draw_legend(ax, ["DNA", "Promoter", "Gene", "TF"], LEGEND_AT)
    mean = Line2D([], [], ls="none", marker="o", ms=MEAN_MARKER["size"] * SCALE,
                  mfc="black", mec="white", mew=MEAN_MARKER["edge"] * SCALE)
    ax.legend([mean], ["mean"], loc="center right", frameon=False,
              bbox_to_anchor=(to_slide(SPINES_X[1], 0)[0], top + LEGEND_ROOM / 2),
              bbox_transform=ax.transData, fontsize=LEGEND_PT * SCALE, borderpad=0,
              borderaxespad=0, handlelength=1.0, handletextpad=0.4)

    with tempfile.TemporaryDirectory() as tmp:
        overlay = Path(tmp) / "overlay.pdf"
        fig.savefig(overlay, transparent=True)
        plt.close(fig)
        mm = lambda x, y: (x * UNIT_MM, y * UNIT_MM)
        typeset([NODE % (*mm(left, top), SCALE, *PANELS_BOX, FIG2),
                 MASK % (*mm(*to_slide(*B_LETTER_MASK[:2])), *mm(*to_slide(*B_LETTER_MASK[2:]))),
                 NODE % (0, SLIDE_MM[1], 1, 0, 0, SLIDE_MM[0] / PT_MM, SLIDE_MM[1] / PT_MM,
                         overlay)],
                "scene24", "scene24")
