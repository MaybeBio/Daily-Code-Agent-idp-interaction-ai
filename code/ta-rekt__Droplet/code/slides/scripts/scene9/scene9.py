"""
Scene 9 -- "The story I'll be telling you today is about how this mechanism of
transcriptional control, despite being constantly bombarded on all sides by
randomly moving molecules [...] could still exert fine control on gene
expression [...]"

The locus of scene 6, transcribed by one polymerase, with transcription factors
gathered at the start of the gene and a few more scattered around (in 3D, as in
scene 6).

Output: slides/figures/scene9/scene9.{pdf,svg}, a full 128 x 96 mm slide.

Usage (from slides/):
    python3 scripts/scene9/scene9.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # slides/scripts, for `common`

import numpy as np

from common.drawing import SLIDE_UNITS, new_slide, save
from common.locus import Locus
from common.molecules import (TF_RADIUS, circle_points, draw_legend, draw_scattered_tfs,
                              legend_area, scatter_proteins)

# Where scene 6 has it (the storyboard's copies differ by < 0.3 units).
LOCUS = Locus(offset=(4.81, -1.34))
POLYMERASES = [0.75]      # fraction of the gene transcribed
TRANSCRIPT_SEED = 9
# TFs bound around the start of the gene: (nodes from the gene start, side).
BOUND_SITES = [(n, side) for n in (-11, -7.5, -4, -0.5) for side in (-1, 1)]
N_FREE_TFS = 10
TF_SEED = 90
MARGIN = 0.5
REGION = (MARGIN, SLIDE_UNITS[0] - MARGIN, MARGIN, SLIDE_UNITS[1] - MARGIN)
LEGEND = ["DNA", "Promoter", "Gene", "TF", "RNAP", "mRNA"]


if __name__ == "__main__":
    tfs = LOCUS.bound_tfs(BOUND_SITES)
    states = LOCUS.transcribe(POLYMERASES, TRANSCRIPT_SEED,
                              [circle_points(c, TF_RADIUS) for c in tfs])
    free = scatter_proteins(N_FREE_TFS, REGION, np.random.default_rng(TF_SEED),
                            LOCUS.obstacles(states, tfs) + [legend_area(LEGEND)])

    fig, ax = new_slide()
    LOCUS.draw(ax, states, tfs)
    draw_scattered_tfs(ax, free, z_back=0.5, z_front=6)
    draw_legend(ax, LEGEND)
    save(fig, "scene9", "scene9")
