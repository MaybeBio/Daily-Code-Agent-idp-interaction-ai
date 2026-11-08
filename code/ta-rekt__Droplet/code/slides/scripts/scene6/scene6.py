"""
Scene 6 -- "This type of control is not all-or-nothing: the degree of expression
of a gene depends in a continuous way on the concentration of proteins that
regulate it. This control mechanism is however very noisy [...]"

Eight slides around the scene-4 locus, with transcription factors scattered in
3D (smaller and fainter = further back):

  Expression rises continuously with TF concentration. One population of TFs
  grows: each slide keeps the TFs of the previous one and adds new ones.
    scene6_1  few TFs, no polymerase
    scene6_2  more TFs, one polymerase on the gene with its nascent mRNA
    scene6_3  more TFs still, two polymerases transcribing the gene at once
  ...but noisily, at one fixed concentration (30 TFs), each slide a new sample:
    scene6_4..6  three independent random configurations, nothing bound
    scene6_7  a TF binds the promoter (highlighted, "Binding") and a polymerase
              initiates (no mRNA yet)
    scene6_8  the polymerase terminates, with its full mRNA; the bound TF is the
              only one that has not moved since 6_7

Transcripts trail behind their polymerase, opposite to its direction of travel.

Outputs: slides/figures/scene6/scene6_<k>.{pdf,svg}, full 128 x 96 mm slides.
Change a slide's seed below to draw a different random configuration.

Usage (from slides/):
    python3 scripts/scene6/scene6.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # slides/scripts, for `common`

import numpy as np
from matplotlib.patches import Circle

from common.drawing import SLIDE_UNITS, draw_label, new_slide, save
from common.locus import POLYMERASE_RADIUS, Locus
from common.molecules import (RNA_BEAD_RADIUS, TF_RADIUS, circle_points, dna_obstacles, draw_dna,
                              draw_legend, draw_rna, draw_rnap, draw_scattered_tfs, draw_tf,
                              legend_area, scatter_proteins)

LOCUS = Locus(offset=(4.81, -1.34))  # the scene-4 locus, where the sketch copies it
MARGIN = 0.5
REGION = (MARGIN, SLIDE_UNITS[0] - MARGIN, MARGIN, SLIDE_UNITS[1] - MARGIN)
BOUND_TF = LOCUS.dna.at(LOCUS.promoter)[0]
BOUND_TF_ANGLE = 20

# Highlight of the binding event in scene6_7, and its label.
HIGHLIGHT = dict(radius=0.5, color="#FFD54A", alpha=0.45)
BINDING_LABEL_OFFSET = (0.2, 1.45)  # from the bound TF; clear of the DNA around it
BINDING_LABEL_HALF_WIDTH = 1.0

# Each slide: total number of TFs (including a bound one), random seed,
# polymerases given as the fraction of the gene each has transcribed
# (0 = initiating, 1 = terminating), whether a TF is bound at the promoter,
# whether this slide extends the TFs of the previous "growing" slide, and
# whether the binding is highlighted.
SLIDES = {
    "scene6_1": dict(n_tfs=16, seed=11, polymerases=[], bound=False, growing=True),
    "scene6_2": dict(n_tfs=30, seed=12, polymerases=[0.55], bound=True, growing=True),
    "scene6_3": dict(n_tfs=52, seed=13, polymerases=[0.95, 0.3], bound=True, growing=True),
    "scene6_4": dict(n_tfs=30, seed=14, polymerases=[], bound=False),
    "scene6_5": dict(n_tfs=30, seed=15, polymerases=[], bound=False),
    "scene6_6": dict(n_tfs=30, seed=16, polymerases=[], bound=False),
    "scene6_7": dict(n_tfs=30, seed=17, polymerases=[0.0], bound=True, highlight=True),
    "scene6_8": dict(n_tfs=30, seed=18, polymerases=[1.0], bound=True),
}
LEGEND = ["DNA", "Promoter", "Gene", "TF", "RNAP", "mRNA"]


def transcription(polymerases, bound, seed):
    """(node index, transcript or None) for each polymerase on the gene; their
    transcripts keep clear of the TF bound at the promoter, if any."""
    blocked = [circle_points(BOUND_TF, TF_RADIUS)] if bound else []
    return LOCUS.transcribe(polymerases, seed, blocked)


def obstacles(spec):
    """What the scattered TFs of a slide must keep clear of, as (points, radius) pairs."""
    avoid = dna_obstacles(LOCUS.dna) + [legend_area(LEGEND)]
    for u, rna in transcription(spec["polymerases"], spec["bound"], spec["seed"]):
        avoid.append((LOCUS.dna.at(u)[0], POLYMERASE_RADIUS))
        if rna is not None:
            avoid.append((rna.nodes, RNA_BEAD_RADIUS))
    if spec["bound"]:
        avoid.append((BOUND_TF, TF_RADIUS))
    if spec.get("highlight"):
        avoid.append((BOUND_TF, HIGHLIGHT["radius"]))
        label = BOUND_TF + np.array(BINDING_LABEL_OFFSET)
        dx = BINDING_LABEL_HALF_WIDTH
        avoid.append((np.linspace(label - (dx, 0), label + (dx, 0), 9), 0.25))
    return avoid


def free_tfs():
    """The scattered (unbound) TFs of every slide.

    The growing slides share one population, kept clear of everything any of
    them draws, so that TFs placed early never collide with later transcripts.
    """
    growing = [s for s in SLIDES.values() if s.get("growing")]
    growing_avoid = [a for s in growing for a in obstacles(s)]
    tfs, previous = {}, []
    for name, spec in SLIDES.items():
        n_free = spec["n_tfs"] - spec["bound"]
        rng = np.random.default_rng(spec["seed"] + 1000)
        if spec.get("growing"):
            previous = tfs[name] = scatter_proteins(n_free, REGION, rng, growing_avoid,
                                                    existing=previous)
        else:
            tfs[name] = scatter_proteins(n_free, REGION, rng, obstacles(spec))
    return tfs


def draw_slide(name, spec, tfs):
    fig, ax = new_slide()
    draw_dna(ax, LOCUS.dna, z=1)
    states = transcription(spec["polymerases"], spec["bound"], spec["seed"])
    for u, rna in states:
        if rna is not None:
            draw_rna(ax, rna, z=3)
        draw_rnap(ax, LOCUS.dna, u, z=4)
    if spec["bound"]:
        draw_tf(ax, BOUND_TF, angle=BOUND_TF_ANGLE, z=5)
    if spec.get("highlight"):
        ax.add_patch(Circle(BOUND_TF, HIGHLIGHT["radius"], fc=HIGHLIGHT["color"], ec="none",
                            alpha=HIGHLIGHT["alpha"], zorder=0.8))
        draw_label(ax, BOUND_TF + np.array(BINDING_LABEL_OFFSET), "Binding")
    draw_scattered_tfs(ax, tfs, z_back=0.5, z_front=6)
    # Polymerases, then transcripts, are the last rows: the others stay put.
    has_mrna = any(rna is not None for _, rna in states)
    draw_legend(ax, LEGEND[:4 + bool(states) + has_mrna])
    save(fig, "scene6", name)


if __name__ == "__main__":
    tfs = free_tfs()
    for name, spec in SLIDES.items():
        draw_slide(name, spec, tfs[name])
