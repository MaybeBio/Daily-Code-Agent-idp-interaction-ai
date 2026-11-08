"""
Scene 20 -- "Let's now step back and take a better look at how biology actually
carries out transcription. [...] In eukaryotes [...] distant parts of DNA are
known to play a central role in transcription initiation. Like promoters,
these are non-coding sequences, known as enhancers. However, they can be
hundreds of thousands of base pairs away from a gene. [...]"

The scene-4 locus, zoomed out to show more of the chain beyond it: promoter
(orange), gene (blue) and, far along the chain, an enhancer (cyan). The
chromatin is dimmed so that these stand out, as in the storyboard.

    scene20_1  the chain stretched out, labelled (left)
    scene20_2  + the same chain folded so that the enhancer touches the
               promoter, with two polymerases transcribing the gene (right)

Each panel is drawn in the scene-4 frame (common/locus.py) on axes that cover
the slide, zoomed out by ZOOM: a scene-4 point p lands at shift + ZOOM * p.
The chain beyond the scene-4 drawing is traced from the storyboard in slide
units.

Outputs: slides/figures/scene20/scene20_<k>.{pdf,svg}, full 128 x 96 mm slides.

Usage (from slides/):
    python3 scripts/scene20/scene20.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # slides/scripts, for `common`

from common.drawing import draw_label, new_slide, save, zoomed_axes
from common.locus import DNA_TRACE, Locus
from common.molecules import DNA_LINE_PT, SEGMENT_COLORS, draw_dna, draw_legend, draw_rna, draw_rnap

ZOOM = 0.61                    # as the storyboard copies the scene-4 drawing
# Chromatin: dimmed, and its lines thinned with the zoom so the helix stays open.
DIM_DNA = dict(line_colors=dict(edge="#9A9A9A", fill="white"),
               nucleosome_colors=dict(fill="#B4BAC4", shade="#E9ECF1"))
ENHANCER = dict(fill="#3CC3E0", edge="#1E7F96")
SEGMENTS = {**SEGMENT_COLORS, "enhancer": ENHANCER}

# ---------------------------------------------------------------------------
# The chain beyond the scene-4 drawing, traced from the storyboard (slide
# units). It continues from the free upper end of the scene-4 DNA.
# ---------------------------------------------------------------------------

STRETCHED = dict(
    shift=(1.35, 1.55),
    continuation=[
        ("coil", [(5.25, 9.72), (5.12, 10.0), (5.15, 10.35), (5.4, 10.62), (5.8, 10.75),
                  (6.2, 10.7), (6.45, 10.45), (6.55, 10.1), (6.75, 9.85), (7.0, 9.9),
                  (7.12, 10.25), (7.1, 10.65), (6.98, 11.0)]),
        ("linker", [(6.98, 11.0), (6.85, 11.25), (6.62, 11.42), (6.38, 11.52)]),
        ("enhancer", [(6.38, 11.52), (6.17, 11.66), (6.05, 11.88), (6.05, 12.12),
                      (6.18, 12.32)]),
        ("linker", [(6.18, 12.32), (6.38, 12.45), (6.62, 12.5), (6.85, 12.6)]),
        ("coil", [(6.85, 12.6), (7.1, 12.85), (7.4, 13.05), (7.7, 13.0), (7.9, 12.75),
                  (8.15, 12.6), (8.42, 12.72), (8.6, 13.0), (8.72, 13.35)]),
    ],
)
# (text, left end or centre, alignment): beside the enhancer and the gene, under
# the promoter.
LABELS = [("enhancer", (5.82, 11.82), "right"), ("promoter", (4.85, 6.2), "center"),
          ("gene", (6.18, 6.85), "left")]

# Folded: the chain loops round, up and back down into the bend of the locus,
# where the enhancer lies along the inside of the promoter; the rest of the
# chain leaves to the upper left.
FOLDED = dict(
    shift=(10.72, 2.84),
    continuation=[
        ("coil", [(14.45, 10.98), (14.1, 11.2), (13.8, 11.55), (13.45, 11.75), (13.1, 11.6),
                  (12.9, 11.25), (12.92, 10.8), (13.1, 10.35), (13.4, 9.9), (13.75, 9.6),
                  (14.1, 9.42), (14.45, 9.3)]),
        ("linker", [(14.45, 9.3), (14.57, 9.05), (14.58, 8.75), (14.55, 8.45), (14.49, 8.2)]),
        ("enhancer", [(14.49, 8.2), (14.39, 8.11), (14.28, 8.16), (14.22, 8.32), (14.21, 8.52),
                      (14.22, 8.72)]),
        ("linker", [(14.24, 8.75), (14.15, 8.95), (13.92, 9.08)]),
        ("coil", [(13.92, 9.08), (13.55, 9.15), (13.15, 9.27), (12.75, 9.45), (12.35, 9.7),
                  (11.95, 9.95), (11.6, 10.2)]),
    ],
    polymerases=[1.0, 0.3],    # fraction of the gene transcribed by each
    seed=20,
)


def locus(shift, continuation):
    """The scene-4 locus (with its promoter), its DNA continued by `continuation`."""
    pieces = [(kind, [((x - shift[0]) / ZOOM, (y - shift[1]) / ZOOM) for x, y in points])
              for kind, points in continuation]
    first_kind, first_points = pieces[0]
    pieces[0] = (first_kind, [DNA_TRACE[-1][1][-1], *first_points])
    return Locus(trace=DNA_TRACE + pieces)


def draw_chain(ax, dna):
    draw_dna(ax, dna, SEGMENTS, z=1, line_pt=tuple(ZOOM * w for w in DNA_LINE_PT), **DIM_DNA)


def draw_stretched(fig, ax):
    draw_chain(zoomed_axes(fig, STRETCHED["shift"], ZOOM), locus(**STRETCHED).dna)
    for text, xy, ha in LABELS:
        draw_label(ax, xy, text, ha=ha, va="baseline")


def draw_folded(fig):
    folded = locus(FOLDED["shift"], FOLDED["continuation"])
    pax = zoomed_axes(fig, FOLDED["shift"], ZOOM)
    draw_chain(pax, folded.dna)
    for u, rna in folded.transcribe(FOLDED["polymerases"], FOLDED["seed"]):
        if rna is not None:
            draw_rna(pax, rna, z=3)
        draw_rnap(pax, folded.dna, u, z=4)


if __name__ == "__main__":
    for name, folded in [("scene20_1", False), ("scene20_2", True)]:
        fig, ax = new_slide()
        draw_stretched(fig, ax)
        if folded:
            draw_folded(fig)
        draw_legend(ax, ["DNA", "Promoter", "Gene"] + ["RNAP", "mRNA"] * folded)
        save(fig, "scene20", name)
