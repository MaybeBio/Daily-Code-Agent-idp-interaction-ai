"""
Scene 16 -- "This plot shows the maximum possible Fisher rate for each given
concentration. You can clearly see the two regimes: diffusion-limited, and
activation-time limited."

Figure 2a of the SPRT paper, re-drawn from its source (ConstrainedPromoters,
figures.ipynb, "Figure 2 / A / Plot") in the slides' font: the maximum Fisher
rate over four-state circuits ("full" decoding), the diffusion limit 4aD c
and the activation limit 0.16 kmax. Left out: the "acc" curve, its diffusion
limit (8aD c/3) and the vertical line at the crossover.

Under the curve, on either side, the circuits optimal at low and high c, as the
paper's Fig. 3a shows its own (slides/figures/scene16/fig3a.png), in its colour
code. Each arrow's colour is read off the data: black if binding is at the
diffusion limit (kb = 4aD), red if the rate is at kmax, light grey if it
saturates neither constraint. Rates at the floor kmin are not drawn (white on the
paper's colour bar). States are darker the more they are occupied.

At high c, binding is not at the diffusion limit: the optimum is c kb = 0.5 kmax,
below both constraints, so that arrow is grey there (the paper's inset, for two
sites, has it black).

Data: fig2a_full.json, a copy of ConstrainedPromoters/data/fig2/a1/
ir_full_vs_ep_pareto_pid1.json (commit dea9418); one entry per c0: the eight
log10 rates of the optimal circuit, its drift V at dc/c = 0.1, and c0.

Output: slides/figures/scene16/scene16.{pdf,svg}, a full 128 x 96 mm slide.

Usage (from slides/):
    python3 scripts/scene16/scene16.py
"""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # slides/scripts, for `common`

import numpy as np
from matplotlib.patches import Circle

from common.drawing import (LABEL_PT, SLIDE_UNITS, TEXT_COLOR, draw_arrow, draw_label, new_slide,
                            save, zoomed_axes)

# ---------------------------------------------------------------------------
# The paper's parameters and data
# ---------------------------------------------------------------------------

FOUR_AD = 0.0888       # 4 a D: a = 3e-3 um, D = 7.4 um^2/s
K_MAX = 5.0            # maximum activation and unbinding rates (1/s)
K_MIN = 1e-5           # lower bound on every rate (1/s)
ACTIVATION_LIMIT = 1 / (1.0 + 2.5 ** 2 / 25)  # 0.8 = 0.16 kmax, as plotted in the paper
DRIFT_TO_FISHER = 2 / 0.1 ** 2  # V = 1/2 (dc/c)^2 Fdot, at dc/c = 0.1

DATA = Path(__file__).with_name("fig2a_full.json")
# The eight rates of an entry, as edges between the paper's states, the binding
# ones (rate proportional to c) marked; and where those states go on the square.
PAPER_EDGES = [(1, 3), (1, 2), (2, 4), (2, 1), (3, 1), (3, 4), (4, 2), (4, 3)]
PAPER_BINDING = {(1, 3), (2, 4)}
PAPER_STATES = {1: (-1, 1), 3: (1, 1), 4: (1, -1), 2: (-1, -1)}  # clockwise from the top left
C0_ENTRY = 9


def paper_rows():
    rows = json.loads(DATA.read_text())
    return [rows[str(i)] for i in range(1, len(rows) + 1)]


def paper_curve():
    """(c0, maximum Fisher rate) as plotted in the paper."""
    rows = paper_rows()
    return (np.array([r[C0_ENTRY] for r in rows]),
            DRIFT_TO_FISHER * np.array([r[8] for r in rows]))


# ---------------------------------------------------------------------------
# Optimal circuits
# ---------------------------------------------------------------------------

def occupancy(rates):
    """Steady-state probability of each state, given {edge: rate}."""
    generator = np.zeros((4, 4))
    for (u, v), rate in rates.items():
        generator[u - 1, v - 1] = rate
    generator -= np.diag(generator.sum(axis=1))
    values, vectors = np.linalg.eig(generator.T)
    probability = np.real(vectors[:, np.argmin(np.abs(values))])
    return probability / probability.sum()


def constraint(edge, log_rate, c):
    """What limits an edge of the optimal circuit at c: "diffusion", "activation",
    "none", or None for a rate at the floor (no arrow)."""
    k = 10 ** log_rate
    binding = edge in PAPER_BINDING
    if k < 1.01 * K_MIN:
        return None
    if (c * k if binding else k) > 0.99 * K_MAX:
        return "activation"
    if binding and k > 0.99 * FOUR_AD:
        return "diffusion"
    return "none"


# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------

AXES = (3.2, 2.9, 15.6, 9.6)  # left, bottom, width, height (slide units)
X_LIMITS, Y_LIMITS = (8e-2, 1.2e3), (6e-3, 16)  # as in the paper
TICK_PT = 8.4
Y_LABEL = r"$\max_{\rm rates} \dot{ℱ}(c)$"  # the best circuit at each c
LIMIT_STYLE = dict(color="black", alpha=0.5, lw=1.2)
DIFFUSION_LABEL_C = 0.25   # where the diffusion limit is labelled...
LABEL_LIFT = 2.0           # ...this factor above the line
ACTIVATION_LABEL_C = 0.1

# The circuits: entries of the data shown (from 1), and the middle of each square
# under the curve (slide units), as in Fig. 3a.
LOW_ROW, HIGH_ROW = 35, 70  # c0 = 1.15 and 152
CIRCUITS = {LOW_ROW: (8.3, 4.35), HIGH_ROW: (15.6, 6.9)}
CIRCUIT_HALF_SIDE = 0.8
NODE_RADIUS = 0.13
NODE_LIGHTEST = 0.72        # grey level of an empty state; a full one is black
ARROW_CLEARANCE = 0.13      # between arrow ends and states
ARROW_PT, ARROW_HEAD = 2.2, 11
ARROW_COLORS = {"diffusion": "black", "activation": "#FF0000", "none": "#B4B4B4"}
LEGEND_LABELS = {"diffusion": "Diffusion limited", "activation": "Activation-limited",
                 "none": "No constraint saturated"}
LEGEND_AT, LEGEND_ROW, LEGEND_ARROW = (11.5, 4.6), 0.62, 0.75


def paper_axes(fig, box, x_limits, y_limits, y_label, tick_pt=TICK_PT):
    """Log-log axes at `box` (left, bottom, width, height in slide units), in
    the style of the paper's figures but the slides' font, with the y label
    horizontal on top of the axis."""
    left, bottom, width, height = box
    ax = fig.add_axes([left / SLIDE_UNITS[0], bottom / SLIDE_UNITS[1],
                       width / SLIDE_UNITS[0], height / SLIDE_UNITS[1]])
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(*x_limits)
    ax.set_ylim(*y_limits)
    # Ticks at the powers of ten in range.
    for limits, set_ticks in ((x_limits, ax.set_xticks), (y_limits, ax.set_yticks)):
        set_ticks(10.0 ** np.arange(np.ceil(np.log10(limits[0])),
                                    np.floor(np.log10(limits[1])) + 1))
    ax.tick_params(labelsize=tick_pt, colors=TEXT_COLOR)
    ax.set_xlabel("Concentration", fontsize=LABEL_PT, color=TEXT_COLOR)
    ax.text(0, 1.03, y_label, transform=ax.transAxes, ha="left", va="bottom",
            fontsize=LABEL_PT, color=TEXT_COLOR)
    return ax


def label_along(ax, c, y, text, slope, size=LABEL_PT, below=False):
    """Text starting at (c, y), above (or `below`) a line of the given log-log
    slope, rotated along it."""
    (x0, y0), (x1, y1) = ax.transData.transform([(c, y), (10 * c, y * 10 ** slope)])
    angle = np.degrees(np.arctan2(y1 - y0, x1 - x0))
    ax.text(c, y, text, rotation=angle, rotation_mode="anchor", ha="left",
            va="top" if below else "bottom", fontsize=size, color=TEXT_COLOR)


def draw_circuit(ax, row, center):
    """The optimal circuit of an entry of the data, as a square around `center`."""
    c = row[C0_ENTRY]
    rates = {e: 10 ** w * (c if e in PAPER_BINDING else 1) for e, w in zip(PAPER_EDGES, row)}
    where = {state: np.asarray(center) + CIRCUIT_HALF_SIDE * np.array(corner)
             for state, corner in PAPER_STATES.items()}
    for state, probability in zip(sorted(where), occupancy(rates)):
        ax.add_patch(Circle(where[state], NODE_RADIUS, fc=str(NODE_LIGHTEST * (1 - probability)),
                            ec="none", zorder=5))
    for (u, v), log_rate in zip(PAPER_EDGES, row):
        limit = constraint((u, v), log_rate, c)
        if limit:
            direction = (where[v] - where[u]) / np.linalg.norm(where[v] - where[u])
            reach = NODE_RADIUS + ARROW_CLEARANCE
            draw_arrow(ax, where[u] + reach * direction, where[v] - reach * direction,
                       ARROW_COLORS[limit], z=4, lw=ARROW_PT, head=ARROW_HEAD)


def draw_circuit_legend(ax):
    x, y = LEGEND_AT
    for k, (limit, text) in enumerate(LEGEND_LABELS.items()):
        row_y = y - k * LEGEND_ROW
        draw_arrow(ax, (x, row_y), (x + LEGEND_ARROW, row_y), ARROW_COLORS[limit], z=4,
                   lw=ARROW_PT, head=ARROW_HEAD)
        draw_label(ax, (x + LEGEND_ARROW + 0.25, row_y), text, ha="left")


if __name__ == "__main__":
    c, maximum = paper_curve()

    fig, _ = new_slide()
    ax = paper_axes(fig, AXES, X_LIMITS, Y_LIMITS, Y_LABEL)
    ax.plot(c, maximum, color="black", lw=1.8, zorder=4)
    ax.plot(c, FOUR_AD * c, **LIMIT_STYLE)
    ax.axhline(ACTIVATION_LIMIT, linestyle="dotted", **LIMIT_STYLE)
    label_along(ax, DIFFUSION_LABEL_C, LABEL_LIFT * FOUR_AD * DIFFUSION_LABEL_C,
                "Diffusion limit", slope=1)
    label_along(ax, ACTIVATION_LABEL_C, 1.15 * ACTIVATION_LIMIT, "Activation limit", slope=0)

    overlay = zoomed_axes(fig, (0, 0), 1)  # on top of the plot, in slide units
    rows = paper_rows()
    for row, center in CIRCUITS.items():
        draw_circuit(overlay, rows[row - 1], center)
    draw_circuit_legend(overlay)
    save(fig, "scene16", "scene16")
