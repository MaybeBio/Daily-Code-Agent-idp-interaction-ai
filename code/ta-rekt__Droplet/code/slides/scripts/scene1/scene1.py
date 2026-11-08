"""
Scene 1 -- "Different cells in your body do different things."

Two cell cartoons traced from the storyboard sketch (storyboard/scene1.svg):
a neuron and a rounded cell. Their geometry and palettes live in
scripts/common/cells.py, shared with the following scenes.

Outputs, written to slides/figures/scene1/ as .pdf (for beamer) and .svg (for
Inkscape), all with a transparent background:
    neuron.*      the neuron alone
    round_cell.*  the rounded cell alone
    scene1.*      both cells laid out as in the storyboard, on a full 128 x 96 mm
                  (beamer 4:3) slide

Every asset uses the same physical scale (UNIT_MM per drawing unit), so
\\includegraphics without a width keeps the cells at their relative sketch sizes.

Usage (from slides/):
    python3 scripts/scene1/scene1.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # slides/scripts, for `common`

from common.cells import draw_neuron, draw_round_cell
from common.drawing import export_single, new_slide, save

if __name__ == "__main__":
    export_single(draw_neuron, "scene1", "neuron")
    export_single(draw_round_cell, "scene1", "round_cell")

    fig, ax = new_slide()
    draw_neuron(ax)
    draw_round_cell(ax)
    save(fig, "scene1", "scene1")
