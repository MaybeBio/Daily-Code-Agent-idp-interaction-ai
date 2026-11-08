"""
Scene 10 -- "To give you a feeling for how surprising this is, let's take a
closer look. The segment of DNA highlighted in orange is called a promoter [...]
Transcription typically happens in bursts of many mRNA molecules elongating in
rapid succession"

The locus of scene 9 with its promoter in orange (labelled), transcription
factors bound on it, and two polymerases on the gene, each with its elongating
transcript.

Output: slides/figures/scene10/scene10.{pdf,svg}, a full 128 x 96 mm slide.

Usage (from slides/):
    python3 scripts/scene10/scene10.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # slides/scripts, for `common`

import numpy as np

from common.drawing import SLIDE_UNITS, draw_label, new_slide, save
from common.locus import Locus
from common.molecules import (TF_RADIUS, circle_points, draw_legend, draw_scattered_tfs,
                              legend_area, scatter_proteins)

LOCUS = Locus(offset=(4.81, -1.34))  # where scenes 6 and 9 have it
POLYMERASES = [0.95, 0.35]  # the polymerase of scene 9, further on, and a new one
TRANSCRIPT_SEED = 9
# TFs bound on the promoter: (nodes from the gene start, side).
BOUND_SITES = [(n, side) for n in (-10, -6.5, -3) for side in (-1, 1)]
N_FREE_TFS = 14
TF_SEED = 100
MARGIN = 0.5
REGION = (MARGIN, SLIDE_UNITS[0] - MARGIN, MARGIN, SLIDE_UNITS[1] - MARGIN)

# "Promoter" sits above the orange segment, over the arc of DNA that leads
# into it (the nearest clear space): offset from the promoter's upstream end.
PROMOTER_LABEL_OFFSET = (0.4, 1.0)
PROMOTER_LABEL_HALF_WIDTH = 1.0
LEGEND = ["DNA", "Promoter", "Gene", "TF", "RNAP", "mRNA"]


if __name__ == "__main__":
    tfs = LOCUS.bound_tfs(BOUND_SITES)
    states = LOCUS.transcribe(POLYMERASES, TRANSCRIPT_SEED,
                              [circle_points(c, TF_RADIUS) for c in tfs])
    promoter_top = LOCUS.dna.at(LOCUS.dna.span("promoter")[0])[0]
    label = promoter_top + np.array(PROMOTER_LABEL_OFFSET)
    label_box = np.linspace(label - (PROMOTER_LABEL_HALF_WIDTH, 0),
                            label + (PROMOTER_LABEL_HALF_WIDTH, 0), 9)
    free = scatter_proteins(N_FREE_TFS, REGION, np.random.default_rng(TF_SEED),
                            LOCUS.obstacles(states, tfs)
                            + [(label_box, 0.25), legend_area(LEGEND)])

    fig, ax = new_slide()
    LOCUS.draw(ax, states, tfs)
    draw_label(ax, label, "Promoter")
    draw_scattered_tfs(ax, free, z_back=0.5, z_front=6)
    draw_legend(ax, LEGEND)
    save(fig, "scene10", "scene10")
