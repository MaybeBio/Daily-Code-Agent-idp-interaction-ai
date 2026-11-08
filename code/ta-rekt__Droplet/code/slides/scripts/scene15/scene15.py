"""
Scene 15 -- "Up till now, we have used the simplest possible model of
transcription initiation [...] We can capture this in the following model
with 4 states instead of 2. You can notice the delay in activating
transcription, governed by this rate kon."

The four-state model (TF unbound / bound x transcription OFF / ON) in the middle
of the slide, going round clockwise: binding, activation (kon), unbinding,
deactivation. Beside each state, a simplified stretch of double helix shows it
on two nearby sites: the binding site (orange), holding a TF or not, and the
gene (blue), with or without a polymerase and its elongating mRNA.

    scene15_1  nothing bound, no transcription; top-left state highlighted
    scene15_2  a TF binds (highlighted, "Binding"), no transcription yet; top right
    scene15_3  bound and transcribing; bottom right
    scene15_4  transcribing, the TF gone; bottom left

fisher_off, below, is the Fisher information of an OFF time (TF binding at rate
x = c kb, then activation at rate a = kon), with the counter-clockwise rates all
zero (SI of the SPRT paper, eq. fisher_off); scene 18 uses it:

    F(c) = 1 - alpha^2 + (2a/x) alpha^4 zeta(3, beta),
    alpha = x / (x - a),  beta = alpha if x > a, 1 - alpha if x < a.

Outputs: slides/figures/scene15/scene15_<k>.{pdf,svg}, full 128 x 96 mm slides.

Usage (from slides/):
    python3 scripts/scene15/scene15.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # slides/scripts, for `common`

import numpy as np
from matplotlib.patches import Circle
from scipy.special import zeta

from common.drawing import draw_label, new_slide, save, zoomed_axes
from common.molecules import (DNA_BOND, DNA_LINE_PT, Polymer, draw_dna, draw_legend, draw_rna,
                              draw_rnap, draw_tf, grow_rna)
from common.telegraph import draw_state_graph

# ---------------------------------------------------------------------------
# Fisher information of an OFF time
# ---------------------------------------------------------------------------

SINGULAR_WINDOW = 1e-3  # relative distance from x = a below which we interpolate


def fisher_off(x, a):
    """Fisher information on ln c of one OFF time: binding at rate x, then
    activation at rate a."""
    def closed_form(x):
        alpha = x / (x - a)
        beta = alpha if x > a else 1 - alpha
        return 1 - alpha ** 2 + 2 * a / x * alpha ** 4 * zeta(3, beta)

    if abs(x - a) < SINGULAR_WINDOW * a:
        # x = a is a removable singularity, near which the closed form loses
        # precision: interpolate across it.
        low, high = a * (1 - SINGULAR_WINDOW), a * (1 + SINGULAR_WINDOW)
        return np.interp(x, [low, high], [closed_form(low), closed_form(high)])
    return closed_form(x)


# ---------------------------------------------------------------------------
# Slides
# ---------------------------------------------------------------------------

# Highlighted state; "Binding" marks the TF that has just arrived.
SLIDES = {
    "scene15_1": dict(state="TL"),
    "scene15_2": dict(state="TR", binding=True),
    "scene15_3": dict(state="BR"),
    "scene15_4": dict(state="BL"),
}

# ---------------------------------------------------------------------------
# Layout (slide units)
# ---------------------------------------------------------------------------

# Four-state graph, in the middle of the slide. Both ways between neighbours
# (clockwise outside). Rate labels: binding and unbinding (rows), activation and
# deactivation of the bound promoter (right column).
STATES = {"TL": (8.35, 9.6), "TR": (11.65, 9.6), "BL": (8.35, 6.1), "BR": (11.65, 6.1)}
CYCLE = ["TL", "TR", "BR", "BL"]
ARROWS = [arrow for a, b in zip(CYCLE, CYCLE[1:] + CYCLE[:1]) for arrow in ((a, b), (b, a))]
RATE_LABELS = {("TL", "TR"): r"$c\,k_b$", ("TR", "TL"): r"$k_u$",
               ("BL", "BR"): r"$c\,k_b$", ("BR", "BL"): r"$k_u$",
               ("TR", "BR"): r"$k_{\rm on}$", ("BR", "TR"): r"$k_{\rm off}$"}
ROW_LABELS = [("Transcription OFF", 11.05), ("Transcription ON", 4.65)]

# Beside each state, its cartoon: where its DNA is centred, and what is on it.
CARTOONS = {"TL": dict(center=(4.3, 9.6), bound=False, transcribing=False),
            "TR": dict(center=(15.7, 9.6), bound=True, transcribing=False),
            "BR": dict(center=(15.7, 6.1), bound=True, transcribing=True),
            "BL": dict(center=(4.3, 6.1), bound=False, transcribing=True)}
ZOOM = 1.45          # the cartoons are drawn this many times their usual size
# Stretches of the DNA, from left to right along it (drawing units from its middle).
STRETCHES = [("linker", -1.6, -1.05), ("promoter", -1.05, -0.4), ("linker", -0.4, 0.05),
             ("gene", 0.05, 1.2), ("linker", 1.2, 1.6)]
TF_AT = (-0.72, 0.27)  # from the middle of the DNA: sitting on the binding site
POLYMERASE_AT = 0.75   # along the DNA: on the gene
# The elongating mRNA trails behind the polymerase, away from the DNA, and is
# not longer than the gene.
MRNA = dict(beads=8, direction=(-1.0, -0.7), seed=15, wiggle=0.3, pull=0.3, clearance=0.2)
BINDING = dict(radius=0.34, color="#FFD54A", alpha=0.45)  # as in scene 6
BINDING_LABEL_OFFSET = (0.0, 0.6)  # from the bound TF


def draw_four_state_graph(ax, highlight=None):
    draw_state_graph(ax, STATES, ARROWS, highlight, labels=RATE_LABELS)
    for text, y in ROW_LABELS:
        draw_label(ax, (10.0, y), text)


def draw_cartoon(ax, center, bound, transcribing, binding=False):
    """A stretch of DNA centred at `center` (slide units), on axes zoomed by ZOOM."""
    x, y = np.asarray(center) / ZOOM
    dna = Polymer.from_trace([(kind, [(x + a, y), (x + b, y)]) for kind, a, b in STRETCHES],
                             DNA_BOND)
    draw_dna(ax, dna, z=1, line_pt=tuple(ZOOM * w for w in DNA_LINE_PT))
    if bound:
        tf = np.array((x + TF_AT[0], y + TF_AT[1]))
        draw_tf(ax, tf, z=5)
        if binding:
            ax.add_patch(Circle(tf, BINDING["radius"], fc=BINDING["color"], ec="none",
                                alpha=BINDING["alpha"], zorder=0.8))
            draw_label(ax, tf + BINDING_LABEL_OFFSET, "Binding")
    if transcribing:
        polymerase = draw_rnap(ax, dna, dna.nearest((x + POLYMERASE_AT, y)), z=4)
        draw_rna(ax, grow_rna(polymerase, MRNA["direction"], MRNA["beads"],
                              np.random.default_rng(MRNA["seed"]), avoid=dna.nodes,
                              clearance=MRNA["clearance"], wiggle=MRNA["wiggle"],
                              pull=MRNA["pull"]), z=3)


def draw_slide(name, state, binding=False):
    fig, ax = new_slide()
    cartoons = zoomed_axes(fig, (0, 0), ZOOM)
    for key, spec in CARTOONS.items():
        draw_cartoon(cartoons, **spec, binding=binding and key == state)
    draw_four_state_graph(ax, state)
    draw_legend(ax, ["DNA", "Promoter", "Gene", "TF", "RNAP", "mRNA"])
    save(fig, "scene15", name)


if __name__ == "__main__":
    for name, spec in SLIDES.items():
        draw_slide(name, **spec)
