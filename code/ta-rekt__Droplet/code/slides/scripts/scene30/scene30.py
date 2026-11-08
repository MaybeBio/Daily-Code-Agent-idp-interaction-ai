"""
Scene 30 -- "One thing I'll tell you so that you appreciate this result. [...]
Still, comparing the Fisher information rates in the same units, we found that
the polymer can still give you more sensitivity, per unit time, than a single
receptor. [...]"

Panel E of the paper's Fig. 6 as published: the Fisher information rate against
the concentration c / c_sat, for the proximal promoter alone, the polymer alone,
both together, and a free receptor.

The panel is cropped from the published figure (slides/figures/scene30/
fig6.pdf) and placed by pdflatex as in scenes 22 and 27, so it stays vector and
keeps its own lettering. The crop leaves out the panel letter.

Output: slides/figures/scene30/scene30.pdf, a full 128 x 96 mm slide.

Usage (from slides/):
    python3 scripts/scene30/scene30.py
Requires the `pdflatex` command-line tool.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # slides/scripts, for `common`

from common.drawing import FIGURES_DIR, SLIDE_UNITS, UNIT_MM
from scene22.scene22 import NODE, PT_MM, typeset

FIG6 = FIGURES_DIR / "scene30" / "fig6.pdf"

# Panel E, a viewport of fig6.pdf (x0, y0, x1, y1 in points, y up): from the
# y-label to the right spine, and from the top tick label to the x-label. The
# panel letter, left of the y-label, and panel F's, right of the spine, are
# outside.
PANEL_E = (14.8, 1.0, 168.0, 162.5)

# The figure's lettering is 9 pt; at this scale the panel, centred, fills the
# slide's height with a margin.
FIG_SCALE = 1.4


if __name__ == "__main__":
    x0, y0, x1, y1 = PANEL_E
    width = (x1 - x0) * FIG_SCALE * PT_MM / UNIT_MM
    height = (y1 - y0) * FIG_SCALE * PT_MM / UNIT_MM
    node = NODE % ((SLIDE_UNITS[0] - width) / 2 * UNIT_MM,
                   (SLIDE_UNITS[1] + height) / 2 * UNIT_MM, FIG_SCALE, *PANEL_E, FIG6)
    typeset([node], "scene30", "scene30")
