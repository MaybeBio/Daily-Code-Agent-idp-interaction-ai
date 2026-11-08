"""
Scene 3 -- "genes are regulated: they are read out and used selectively,
depending on the cell's context"

The scene-2 layout with the DNA coloured by use: red stretches are genes used by
the neuron, blue by the rounded cell, purple by both, and black is non-coding
DNA. Each cell takes the colour of its genes, and an arrow runs from the
chromosomes to it.

Output: slides/figures/scene3/scene3.{pdf,svg}, a full 128 x 96 mm slide.

Usage (from slides/):
    python3 scripts/scene3/scene3.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # slides/scripts, for `common`

from common.cells import draw_neuron, draw_round_cell
from common.chromosomes import CHROMOSOMES, GENE_COLORS, draw_chromosome
from common.drawing import draw_arrow, new_slide, save

NEURON_RED = dict(fill="#F8CFCF", edge=GENE_COLORS["neuron"], shine="#FFFFFF")
ROUND_CELL_BLUE = dict(fill="#CEDBF5", edge=GENE_COLORS["round"], shine="#FFFFFF",
                       organelle="#93AFE6")

# Arrows from the chromosomes to each cell (start, end), traced from the sketch.
NEURON_ARROW = ((5.94, 8.45), (5.94, 9.95))
ROUND_CELL_ARROW = ((7.6, 6.85), (10.8, 6.85))


if __name__ == "__main__":
    fig, ax = new_slide()
    draw_neuron(ax, colors=NEURON_RED)
    draw_round_cell(ax, colors=ROUND_CELL_BLUE)
    for chromosome in CHROMOSOMES:
        draw_chromosome(ax, chromosome, banded=True, z=5)
    draw_arrow(ax, *NEURON_ARROW, GENE_COLORS["neuron"], z=10)
    draw_arrow(ax, *ROUND_CELL_ARROW, GENE_COLORS["round"], z=10)
    save(fig, "scene3", "scene3")
