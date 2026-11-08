"""
Scene 4 -- "reading out a gene's information to do something in the cell begins
with transcription. this is when an enzyme called RNA polymerase runs along the
gene and produces a messenger RNA molecule. this molecule later gets translated
into proteins which then carry out the cell's needs"

Three slides, each a full 128 x 96 mm slide in slides/figures/scene4/:
    scene4_1  transcription panel: RNA polymerase at the start of the gene, no mRNA yet
    scene4_2  transcription panel: polymerase at the end of the gene with its mRNA
              (as in the storyboard)
    scene4_3  scene4_2 plus the translation panel (mRNA -> protein, drawn as the
              standard transcription factor in the blue of its gene)

DNA and RNA are chunked polymers (see common/molecules.py), ready to be animated.

Usage (from slides/):
    python3 scripts/scene4/scene4.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # slides/scripts, for `common`

import numpy as np

from common.drawing import draw_arrow, draw_label, new_slide, save
from common.locus import Locus
from common.molecules import (BLUE, RNA_BOND, Polymer, draw_dna, draw_legend, draw_rna, draw_rnap,
                              draw_tf)

# ---------------------------------------------------------------------------
# Geometry traced from the storyboard (slide units); the locus (DNA, polymerase
# positions, finished mRNA) is in common/locus.py, shared with later scenes.
# ---------------------------------------------------------------------------

# From the end by the arrow; the mRNA is cut to the length of the gene.
MRNA_TRANSLATED = [(14.95, 9.3), (14.76, 9.23), (14.49, 9.2), (14.2, 9.14), (14.05, 9.06),
                   (13.88, 8.81), (13.65, 8.64), (13.5, 8.6), (13.2, 8.6), (13.06, 8.57),
                   (12.93, 8.46), (12.82, 8.24)]

TRANSLATION_ARROW = ((15.22, 8.99), (16.68, 8.99))
PROTEIN_CENTER = (17.67, 8.91)
PROTEIN_SCALE = 1.75

LABELS = {
    "DNA": (2.27, 6.49), "mRNA": (6.2, 7.3), "Protein": (17.7, 7.94),
    "Transcription": (5.27, 5.33), "Translation": (15.58, 6.13),
}
# "RNA Pol" is placed relative to the polymerase: to its right when it is at the
# end of the gene (as in the sketch); below it at the start, where the mRNA
# will later be, clear of the DNA.
RNAP_LABEL_AT_END = dict(offset=(0.53, 0.03), ha="left")
RNAP_LABEL_AT_START = dict(offset=(0.0, -0.6), ha="center")

LOCUS = Locus()


# ---------------------------------------------------------------------------
# Panels
# ---------------------------------------------------------------------------

def draw_transcription(ax, finished):
    """DNA with RNA polymerase at the start of the gene, or (finished) at its end
    with the mRNA it made."""
    draw_dna(ax, LOCUS.dna, z=1)
    if finished:
        u = LOCUS.end
        label = RNAP_LABEL_AT_END
        draw_rna(ax, LOCUS.full_transcript(), z=3)
        draw_label(ax, LABELS["mRNA"], "mRNA")
    else:
        u = LOCUS.start
        label = RNAP_LABEL_AT_START
    polymerase = draw_rnap(ax, LOCUS.dna, u, z=4)
    draw_label(ax, polymerase + np.array(label["offset"]), "RNA Pol", ha=label["ha"])
    draw_label(ax, LABELS["DNA"], "DNA")
    draw_label(ax, LABELS["Transcription"], "Transcription", bold=True)


def draw_translation(ax):
    draw_rna(ax, Polymer.from_trace([("rna", MRNA_TRANSLATED)], RNA_BOND).cut(LOCUS.gene_length),
             z=3)
    draw_arrow(ax, *TRANSLATION_ARROW, color="#3A3A3A", z=3)
    draw_tf(ax, PROTEIN_CENTER, scale=PROTEIN_SCALE, colors=BLUE, z=3)
    draw_label(ax, LABELS["Protein"], "Protein")
    draw_label(ax, LABELS["Translation"], "Translation", bold=True)


if __name__ == "__main__":
    for name, finished, translation in [("scene4_1", False, False),
                                        ("scene4_2", True, False),
                                        ("scene4_3", True, True)]:
        fig, ax = new_slide()
        draw_transcription(ax, finished)
        if translation:
            draw_translation(ax)
        draw_legend(ax, ["DNA", "Promoter", "Gene", "RNAP"] + ["mRNA"] * finished)
        save(fig, "scene4", name)
