"""
Scene 23 -- "multi-scale temporal information in transcription activity"
(storyboard note).

Scene 22's two 3D plots stay at the top, where scene 22 has them. Below them,
left: panel F of the paper's Fig. 1, its two transcription time series put on
top of each other as published: the polymer's probe (top), switching with the
polymer between closed (red) and open (green) states, and the free receptor
(bottom). The second slide reveals panel E right of them: one stretch of the
polymer's trajectory, its smoothed transcription activity (black) against the
probe-enhancer distance (orange).

F and E are cropped from the published figure (slides/figures/scene23/
fig1.pdf) and placed by pdflatex, so they stay vector and keep the figure's own
lettering. The panel letters are left out.

Output: slides/figures/scene23/scene23_{1,2}.pdf, full 128 x 96 mm slides.

Usage (from slides/):
    python3 scripts/scene23/scene23.py
Requires the `gs` and `pdflatex` command-line tools.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # slides/scripts, for `common`

from common.drawing import FIGURES_DIR, SLIDE_UNITS, UNIT_MM
from scene22.scene22 import MASK, NODE, PT_MM, plots, typeset

FIG1 = FIGURES_DIR / "scene23" / "fig1.pdf"

# Panel F's plots share one geometry, the free receptor's 243 pt right of the
# polymer's, so one viewport of fig1.pdf (x0, y0, x1, y1 in points, y up),
# shifted, crops each: from the y-label to the right spine, and from the line
# of labels above the axes to the T's under them. The panel letter, above the
# polymer's y-label, is outside it.
POLYMER_BOX = (4.0, 4.6, 230.8, 73.2)
RECEPTOR_SHIFT = 243.0

# Panel E, a viewport of fig1.pdf: from one y-label to the other, and from the
# top tick labels to the x-label.
E_BOX = (295.3, 79.9, 476.6, 192.0)
# E's letter stands left of the y-label, its lowest stroke level with the top
# of the axes, and the end of that stroke falls inside the viewport: it is
# painted white (box in points, as E_BOX; nothing else is drawn there). The box
# clears the stroke's end by 0.8 pt and its underside by 2 pt, and overhangs
# the viewport's top and left edges, or rasterizers leave a speck of its
# anti-aliased edges.
E_LETTER_MASK = (294.3, 187.3, 297.0, 193.0)

# Fig. 1's lettering is 9 pt; drawn a little smaller (7.6 pt), F and E fit
# side by side across the slide. They are centred as a row, and on each other
# vertically, in the space under the 3D plots.
FIG_SCALE = 0.84
GAP = 0.3       # between F's two plots (slide units)
ROW_GAP = 0.4   # between F and E


def size(box):
    """Width and height (slide units) of a viewport of fig1.pdf at FIG_SCALE."""
    return ((box[2] - box[0]) * FIG_SCALE * PT_MM / UNIT_MM,
            (box[3] - box[1]) * FIG_SCALE * PT_MM / UNIT_MM)


def place(box, left, top):
    """TikZ node: viewport `box` of fig1.pdf, its top-left corner at (left,
    top) (slide units)."""
    return NODE % (left * UNIT_MM, top * UNIT_MM, FIG_SCALE, *box, FIG1)


def mask(box, left, top, letter):
    """TikZ: white over the box `letter` (points) of viewport `box`, placed at
    (left, top)."""
    def to_mm(x, y):
        """Slide position (mm) of a point of fig1.pdf (points)."""
        return (left * UNIT_MM + (x - box[0]) * FIG_SCALE * PT_MM,
                top * UNIT_MM - (box[3] - y) * FIG_SCALE * PT_MM)

    return MASK % (*to_mm(*letter[:2]), *to_mm(*letter[2:]))


if __name__ == "__main__":
    plot_nodes, plots_box = plots()
    f_width, f_height = size(POLYMER_BOX)
    e_width, e_height = size(E_BOX)
    stack = 2 * f_height + GAP
    left = (SLIDE_UNITS[0] - f_width - ROW_GAP - e_width) / 2
    top = (plots_box[1] * PT_MM / UNIT_MM + stack) / 2

    # F: the polymer's plot, then the receptor's.
    x0, y0, x1, y1 = POLYMER_BOX
    f_nodes = [place((x0 + shift, y0, x1 + shift, y1), left, top - i * (f_height + GAP))
               for i, shift in enumerate((0, RECEPTOR_SHIFT))]
    e_left, e_top = left + f_width + ROW_GAP, top - (stack - e_height) / 2

    typeset(plot_nodes + f_nodes, "scene23", "scene23_1")
    typeset(plot_nodes + f_nodes + [place(E_BOX, e_left, e_top),
                                    mask(E_BOX, e_left, e_top, E_LETTER_MASK)],
            "scene23", "scene23_2")
