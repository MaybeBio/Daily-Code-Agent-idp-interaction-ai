"""
Scene 5 -- "whether a given gene gets expressed can depend indirectly on many
things [...] Most often, it is directly controlled by transcription factors.
these are proteins that bind near the gene and influence the transcription
process, either by enabling or suppressing it."

The gene of a transcription factor (green) is transcribed and translated; the
factor binds next to the blue gene of scene 4, near a free RNA polymerase.
Transcription factors are copies of slides/figures/tf.svg.

Output: slides/figures/scene5/scene5.{pdf,svg}, a full 128 x 96 mm slide.

Usage (from slides/):
    python3 scripts/scene5/scene5.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # slides/scripts, for `common`

from common.drawing import draw_arrow, new_slide, save
from common.locus import Locus
from common.molecules import (DNA_BOND, GREEN, RNA_BOND, SEGMENT_COLORS, Polymer, draw_dna,
                              draw_legend, draw_rna, draw_rnap_at, draw_tf)

# ---------------------------------------------------------------------------
# Geometry traced from the storyboard (slide units)
# ---------------------------------------------------------------------------

# The transcription factor's own gene, from its lower end to its upper end.
TF_DNA_TRACE = [
    ("coil", [(4.5, 8.94), (4.72, 9.18), (4.64, 9.46), (4.44, 9.73), (4.28, 10.01),
              (4.37, 10.28), (4.59, 10.45), (4.86, 10.56), (5.14, 10.65), (5.36, 10.65),
              (5.53, 10.59)]),
    ("promoter", [(5.53, 10.59), (5.75, 10.56), (5.91, 10.43), (6.08, 10.26)]),
    ("gene", [(6.08, 10.26), (6.35, 10.17), (6.63, 10.23), (6.91, 10.28), (7.09, 10.43),
              (7.15, 10.62), (7.04, 11.0)]),
    ("linker", [(7.04, 11.0), (6.96, 11.17), (6.85, 11.39), (6.69, 11.5), (6.46, 11.53),
                (6.24, 11.5), (6.02, 11.37), (5.91, 11.41), (5.85, 11.55), (5.8, 11.72)]),
    ("coil", [(5.8, 11.72), (5.64, 11.89), (5.41, 11.97), (5.19, 12.11), (5.03, 12.27),
              (5.05, 12.49), (5.19, 12.66), (5.41, 12.75), (5.64, 12.77)]),
]
# Between the two arrows; it is cut to the length of the gene.
TF_MRNA_TRACE = [
    (9.5, 9.83), (9.64, 9.85), (9.78, 9.89), (9.76, 10.04), (9.75, 10.17), (9.88, 10.22),
    (10.02, 10.25), (10.11, 10.35), (10.02, 10.47), (9.89, 10.55), (9.78, 10.62),
    (9.64, 10.68), (9.5, 10.76), (9.47, 10.9), (9.58, 11.0)]
TF_PROTEINS = [((12.42, 10.14), 15), ((12.85, 9.92), 300)]  # (centre, angle)

# gene -> mRNA -> protein -> binds near the blue gene
ARROWS = [((7.79, 10.54), (8.59, 10.54)), ((10.8, 10.2), (11.62, 10.2)),
          ((12.8, 9.35), (12.8, 8.16))]
ARROW_COLOR = "#3A3A3A"

LOCUS = Locus(offset=(6.96, -1.60))  # the scene-4 locus, where the sketch copies it
FREE_POLYMERASE = (13.31, 5.7)


if __name__ == "__main__":
    fig, ax = new_slide()

    tf_dna = Polymer.from_trace(TF_DNA_TRACE, DNA_BOND)
    gene_first, gene_last = tf_dna.span("gene")
    draw_dna(ax, tf_dna, segment_colors={**SEGMENT_COLORS, "gene": GREEN}, z=1)
    draw_rna(ax, Polymer.from_trace([("rna", TF_MRNA_TRACE)], RNA_BOND)
             .cut((gene_last - gene_first) * DNA_BOND), color=GREEN["fill"], z=3)
    for center, angle in TF_PROTEINS:
        draw_tf(ax, center, angle=angle, z=3)
    for start, end in ARROWS:
        draw_arrow(ax, start, end, ARROW_COLOR, z=3)

    draw_dna(ax, LOCUS.dna, z=1)
    draw_tf(ax, LOCUS.dna.at(LOCUS.promoter)[0], angle=20, z=5)
    draw_rnap_at(ax, FREE_POLYMERASE, angle=10, z=4)
    draw_legend(ax, ["DNA", "Promoter", "Gene", "TF", "RNAP", "mRNA"], mrna=GREEN["fill"])

    save(fig, "scene5", "scene5")
