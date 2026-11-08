"""
Scene 27 -- "However something more interesting goes on when we think about the
distal enhancer. When two attractive surfaces come together, they form a
capillary bridge. [...] So we solved the Cahn-Hilliard equation with two
spheres, and from there we could compute exactly the potential and the force
that arise from this capillary bridge, and crucially, their dependence on
concentration."

Left: panel A of the paper's Fig. 5, the bulk-mediated attraction between two
spheres through a capillary bridge, and the attractive potential it makes,
deeper as c increases. Right, one above the other: panels A and B of Fig. 6,
the ligands' volume fraction at the probe and the attractive force against the
distance between the spheres, at c/c_sat = 0.1 and 0.9.

The panels are cropped from the published figures (slides/figures/scene27/
fig5.pdf, fig6.pdf) and placed by pdflatex as in scene 22, so they stay vector
and keep their own lettering. The crops leave the panel letters out.

Output: slides/figures/scene27/scene27.pdf, a full 128 x 96 mm slide.

Usage (from slides/):
    python3 scripts/scene27/scene27.py
Requires the `pdflatex` command-line tool.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # slides/scripts, for `common`

from common.drawing import FIGURES_DIR, SLIDE_UNITS, UNIT_MM
from scene22.scene22 import NODE, PT_MM, typeset

FIG5 = FIGURES_DIR / "scene27" / "fig5.pdf"
FIG6 = FIGURES_DIR / "scene27" / "fig6.pdf"

# Viewports (x0, y0, x1, y1 in points, y up). Fig. 5's A: from the faint left
# edge of the spheres' glow to the tip of the distance axis, and from the titles
# to its label; the panel letter, left of the glow, is outside. Fig. 6's A and
# B: from the y-label to the last tick label (which overhangs the axes), and
# from the top tick label to the x-label; the panel letters, above the axes,
# are outside.
FIG5_A = (16.5, 85.7, 197.8, 200.0)
FIG6_A = (4.5, 182.7, 124.6, 280.2)
FIG6_B = (128.7, 182.7, 245.1, 280.2)

# The figures' lettering is 9 pt; at this scale the panels fill the slide's
# width, Fig. 5's A on the left, Fig. 6's A over B on the right.
FIG_SCALE = 1.12
MARGIN = 0.45   # left of Fig. 5's A, right of Fig. 6's panels (slide units)
GAP = 0.35      # between Fig. 6's A and B


def size(box):
    """Width and height (slide units) of a viewport at FIG_SCALE."""
    return ((box[2] - box[0]) * FIG_SCALE * PT_MM / UNIT_MM,
            (box[3] - box[1]) * FIG_SCALE * PT_MM / UNIT_MM)


def place(pdf, box, left, top):
    """TikZ node: viewport `box` of `pdf`, its top-left corner at (left, top)
    (slide units)."""
    return NODE % (left * UNIT_MM, top * UNIT_MM, FIG_SCALE, *box, pdf)


if __name__ == "__main__":
    # Fig. 5's A, centred vertically on the left.
    nodes = [place(FIG5, FIG5_A, MARGIN, (SLIDE_UNITS[1] + size(FIG5_A)[1]) / 2)]
    # Fig. 6's A over B, the pair centred vertically, each flush right.
    top = (SLIDE_UNITS[1] + size(FIG6_A)[1] + GAP + size(FIG6_B)[1]) / 2
    for box in (FIG6_A, FIG6_B):
        width, height = size(box)
        nodes.append(place(FIG6, box, SLIDE_UNITS[0] - MARGIN - width, top))
        top -= height + GAP
    typeset(nodes, "scene27", "scene27")
