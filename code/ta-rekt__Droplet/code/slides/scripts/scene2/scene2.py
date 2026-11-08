"""
Scene 2 -- "...even when they contain the exact same DNA. How is that?"

The two cells of scene 1, unchanged and in the same place, plus the chromosomes
they both carry, with a line from each nucleus to the chromosomes.

Output: slides/figures/scene2/scene2.{pdf,svg}, a full 128 x 96 mm slide.

Usage (from slides/):
    python3 scripts/scene2/scene2.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # slides/scripts, for `common`

import numpy as np
from matplotlib.patches import Circle

from common.cells import draw_neuron, draw_round_cell
from common.chromosomes import CHROMOSOMES, draw_chromosome
from common.drawing import new_slide, save

LINK_COLOR = "#4A4A55"
LINK_PT = 1.2
LINK_DOT_RADIUS = 0.08

# Where each nucleus' line ends, next to the chromosomes (traced from the sketch).
NEURON_LINK_END = (5.67, 8.59)
ROUND_CELL_LINK_END = (7.4, 7.0)


def draw_link(ax, nucleus, end, z):
    """Line from a cell's nucleus (marked with a dot) to the chromosomes."""
    ax.plot(*np.transpose([nucleus, end]), color=LINK_COLOR, lw=LINK_PT,
            solid_capstyle="round", zorder=z)
    ax.add_patch(Circle(nucleus, LINK_DOT_RADIUS, fc=LINK_COLOR, ec="none", zorder=z))


if __name__ == "__main__":
    fig, ax = new_slide()
    neuron = draw_neuron(ax)
    round_cell = draw_round_cell(ax)
    for chromosome in CHROMOSOMES:
        draw_chromosome(ax, chromosome, z=5)
    draw_link(ax, neuron["nucleus"], NEURON_LINK_END, z=10)
    draw_link(ax, round_cell["nucleus"], ROUND_CELL_LINK_END, z=10)
    save(fig, "scene2", "scene2")
