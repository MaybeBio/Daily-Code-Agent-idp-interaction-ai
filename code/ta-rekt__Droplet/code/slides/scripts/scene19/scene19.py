"""
Scene 19 -- "However, we've learned something more interesting: the best
promoters don't necessarily saturate the limits of their performance allowable
by general principles. They are always constrained by physical laws, but those
might come up as system-dependent chemical rates. [...]"

Takeaway. Left: the landscape of a function to optimize over two parameters,
in oblique projection. The constrained parameter space is a rectangle on the
floor; the landscape's peaks lie outside it, so the best a system can do is on
its boundary: the dots are the highest points of the landscape above the left
(blue) and back (red) edges. Right: each optimum is reached by a different
mechanism, a four-state cycle (red) and a six-state one (blue).

Output: slides/figures/scene19/scene19.{pdf,svg}, a full 128 x 96 mm slide.

Usage (from slides/):
    python3 scripts/scene19/scene19.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # slides/scripts, for `common`

import numpy as np
from matplotlib.collections import PolyCollection
from matplotlib.colors import to_rgb
from matplotlib.patches import Circle, Polygon

from common.drawing import (AXIS_COLOR, LABEL_PT, TEXT_COLOR, draw_arrow, draw_label,
                            new_slide, save)
from common.molecules import BLUE

# As in the sketch, the title is 1.95 times the size of the labels.
TITLE_PT = LABEL_PT * 0.686485 / 0.352777
TITLE = dict(text="Takeaway", xy=(10.0, 13.03))
RED = dict(fill="#D93A45", edge="#8E1F2A")  # the red of scenes 2-3

# ---------------------------------------------------------------------------
# Landscape (left): parameters and function values are lengths in slide units
# ---------------------------------------------------------------------------

# Oblique projection: param 1 runs right, param 2 into the page (up and right,
# at DEPTH_ANGLE), the function up.
ORIGIN = np.array([1.75, 6.55])
DEPTH_ANGLE = np.radians(20)
AXIS_LENGTHS = dict(param1=5.3, param2=3.6, function=4.2)

# The function: a base level plus two Gaussian hills ((param 1, param 2) of the
# peak, height, width). The tall one is behind the constrained region, the
# other to its left.
BASE = 1.4
HILLS = [((2.55, 3.55), 2.2, 0.75), ((0.45, 1.9), 1.2, 0.7)]

SURFACE_RANGE = dict(param1=(0.4, 4.3), param2=(0.4, 4.0))
MESH = (11, 10)            # cells along param 1, param 2
MESH_SUBSTEPS = 6          # points per cell edge, so cells follow the surface
SURFACE_COLORS = dict(lit="#FBFBFB", shade="#A7A9B6", line="#5C5E68")
LIGHT = np.array([-0.45, -0.55, 1.0])  # from the upper left, in front

# Constrained parameter space: a slab on the floor, with its left edge (low
# param 1) blue and its back edge (high param 2) red.
REGION = dict(param1=(1.3, 3.6), param2=(1.2, 3.0))
SLAB_THICKNESS = 0.12
SLAB_COLORS = dict(top="#D5D6DE", side="#9A9CA8", edge=AXIS_COLOR)
EDGE_PT = 3.2
DOT_RADIUS = 0.13

LABELS = {
    "function": [("function", (1.75, 11.73)), ("to optimize", (1.75, 11.0))],
    "param 1": [("param 1", (6.3, 6.05))],
    "param 2": [("param 2", (3.03, 7.72))],
    "region": [("constrained", (10.3, 7.8)), ("parameter space", (10.3, 7.07))],
}

# ---------------------------------------------------------------------------
# Mechanisms (right)
# ---------------------------------------------------------------------------

HEADING = [("same optimized function", (16.4, 11.75)), ("different mechanisms", (16.4, 11.02))]
OPTIMUM_X = 14.4           # the optima again, left of their mechanism
MAPS_TO = (14.8, 15.8)     # x of the arrow from optimum to mechanism
STATE_RADIUS = 0.1
STATE_CLEARANCE = 0.06     # between transition arrows and states
PAIR_GAP = 0.11           # half-distance between the arrows of a reversible step

# States (slide units) and transitions, going round clockwise; the six-state
# cycle's bottom steps are reversible.
FOUR_STATE = dict(
    y=9.85,
    states={"TL": (16.55, 10.45), "TR": (17.95, 10.45), "BL": (16.55, 9.25),
            "BR": (17.95, 9.25)},
    steps=[("TL", "TR"), ("TR", "BR"), ("BR", "BL"), ("BL", "TL")],
)
SIX_STATE = dict(
    y=7.85,
    states={"T0": (16.55, 8.45), "T1": (17.4, 8.45), "T2": (18.25, 8.45),
            "B0": (16.55, 7.25), "B1": (17.4, 7.25), "B2": (18.25, 7.25)},
    steps=[("T0", "T1"), ("T1", "T2"), ("T2", "B2"), ("B2", "B1"), ("B1", "B2"),
           ("B1", "B0"), ("B0", "B1"), ("B0", "T0")],
)


# ---------------------------------------------------------------------------
# Landscape
# ---------------------------------------------------------------------------

def project(p1, p2, f):
    """Slide position of the point (param 1, param 2, function value)."""
    return (ORIGIN[0] + p1 + p2 * np.cos(DEPTH_ANGLE),
            ORIGIN[1] + p2 * np.sin(DEPTH_ANGLE) + f)


def landscape(p1, p2):
    return BASE + sum(h * np.exp(-((p1 - c1) ** 2 + (p2 - c2) ** 2) / (2 * w ** 2))
                      for (c1, c2), h, w in HILLS)


def nearness(p1, p2, f):
    """How near a point is to the viewer: points along the viewing direction
    (-cos, 1, -sin) all project to one spot; the nearest is drawn on top."""
    return p1 * np.cos(DEPTH_ANGLE) - p2 + f * np.sin(DEPTH_ANGLE)


def shade(p1, p2, eps=1e-4):
    """Fill colour from the surface normal, lit from LIGHT."""
    d1 = (landscape(p1 + eps, p2) - landscape(p1 - eps, p2)) / (2 * eps)
    d2 = (landscape(p1, p2 + eps) - landscape(p1, p2 - eps)) / (2 * eps)
    normal = np.array([-d1, -d2, 1.0])
    light = max(normal @ LIGHT / np.linalg.norm(normal) / np.linalg.norm(LIGHT), 0)
    lit, dark = np.array(to_rgb(SURFACE_COLORS["lit"])), np.array(to_rgb(SURFACE_COLORS["shade"]))
    return dark + (lit - dark) * light


def draw_surface(ax, z=2):
    """The landscape as mesh cells, filled and outlined, far ones first (so
    nearer cells hide what is behind them)."""
    p1s = np.linspace(*SURFACE_RANGE["param1"], MESH[0] + 1)
    p2s = np.linspace(*SURFACE_RANGE["param2"], MESH[1] + 1)
    t = np.linspace(0, 1, MESH_SUBSTEPS, endpoint=False)
    cells = []
    for a1, b1 in zip(p1s[:-1], p1s[1:]):
        for a2, b2 in zip(p2s[:-1], p2s[1:]):
            # Round the cell: front edge, right edge, back edge, left edge.
            q1 = np.concatenate([a1 + (b1 - a1) * t, np.full_like(t, b1),
                                 b1 - (b1 - a1) * t, np.full_like(t, a1)])
            q2 = np.concatenate([np.full_like(t, a2), a2 + (b2 - a2) * t,
                                 np.full_like(t, b2), b2 - (b2 - a2) * t])
            c1, c2 = (a1 + b1) / 2, (a2 + b2) / 2
            cells.append((nearness(c1, c2, landscape(c1, c2)),
                          np.column_stack(project(q1, q2, landscape(q1, q2))),
                          shade(c1, c2)))
    cells.sort(key=lambda cell: cell[0])
    ax.add_collection(PolyCollection([c[1] for c in cells], facecolors=[c[2] for c in cells],
                                     edgecolors=SURFACE_COLORS["line"], linewidths=0.5,
                                     joinstyle="round", zorder=z))


def draw_region(ax, z=1):
    """The constrained parameter space: a thin slab whose top is on the floor."""
    (a1, b1), (a2, b2) = REGION["param1"], REGION["param2"]
    top = np.column_stack(project(np.array([a1, b1, b1, a1]), np.array([a2, a2, b2, b2]), 0))
    down = np.array([0, -SLAB_THICKNESS])
    front_left, front_right, back_right, back_left = top
    for face in ([front_left, front_right, front_right + down, front_left + down],
                 [front_right, back_right, back_right + down, front_right + down]):
        ax.add_patch(Polygon(face, closed=True, fc=SLAB_COLORS["side"], ec=SLAB_COLORS["edge"],
                             lw=0.8, joinstyle="round", zorder=z))
    ax.add_patch(Polygon(top, closed=True, fc=SLAB_COLORS["top"], ec=SLAB_COLORS["edge"],
                         lw=0.8, joinstyle="round", zorder=z))
    for (start, end), colors in (((front_left, back_left), BLUE),
                                 ((back_left, back_right), RED)):
        ax.plot(*np.array([start, end]).T, color=colors["fill"], lw=EDGE_PT,
                solid_capstyle="round", zorder=z + 0.5)


def optima():
    """Highest point of the landscape along the blue and the red edge,
    as (param 1, param 2)."""
    (a1, b1), (a2, b2) = REGION["param1"], REGION["param2"]
    along_blue, along_red = np.linspace(a2, b2, 2001), np.linspace(a1, b1, 2001)
    return {"blue": (a1, along_blue[np.argmax(landscape(a1, along_blue))]),
            "red": (along_red[np.argmax(landscape(along_red, b2))], b2)}


def draw_dot(ax, xy, colors, z=6):
    ax.add_patch(Circle(xy, DOT_RADIUS, fc=colors["fill"], ec=colors["edge"], lw=0.8,
                        zorder=z))


def draw_landscape(ax):
    x0, y0 = ORIGIN
    draw_arrow(ax, (x0, y0), project(AXIS_LENGTHS["param1"], 0, 0), AXIS_COLOR, z=3,
               lw=1.3, head=9)
    draw_arrow(ax, (x0, y0), project(0, AXIS_LENGTHS["param2"], 0), AXIS_COLOR, z=3,
               lw=1.3, head=9)
    draw_arrow(ax, (x0, y0), project(0, 0, AXIS_LENGTHS["function"]), AXIS_COLOR, z=3,
               lw=1.3, head=9)
    draw_region(ax)
    draw_surface(ax)
    for name, (p1, p2) in optima().items():
        top = project(p1, p2, landscape(p1, p2))
        ax.plot(*np.array([project(p1, p2, 0), top]).T, color=AXIS_COLOR, lw=1.0,
                dashes=(3, 2), zorder=5)
        draw_dot(ax, top, BLUE if name == "blue" else RED)
    for lines in LABELS.values():
        for text, xy in lines:
            draw_label(ax, xy, text, va="baseline")


# ---------------------------------------------------------------------------
# Mechanisms
# ---------------------------------------------------------------------------

def draw_mechanism(ax, mechanism, colors):
    y = mechanism["y"]
    draw_dot(ax, (OPTIMUM_X, y), colors)
    draw_arrow(ax, (MAPS_TO[0], y), (MAPS_TO[1], y), AXIS_COLOR, z=3, lw=2.2, head=10)
    states = {name: np.array(xy) for name, xy in mechanism["states"].items()}
    for center in states.values():
        ax.add_patch(Circle(center, STATE_RADIUS, fc=colors["fill"], ec=colors["edge"], lw=0.8,
                            zorder=4))
    steps = set(mechanism["steps"])
    for source, target in steps:
        a, b = states[source], states[target]
        along = (b - a) / np.linalg.norm(b - a)
        reach = STATE_RADIUS + STATE_CLEARANCE
        # The two arrows of a reversible step sit side by side, each on its
        # own right.
        shift = PAIR_GAP * np.array([along[1], -along[0]]) if (target, source) in steps else 0
        draw_arrow(ax, a + along * reach + shift, b - along * reach + shift, AXIS_COLOR, z=3,
                   lw=1.0, head=6)


if __name__ == "__main__":
    fig, ax = new_slide()
    ax.text(*TITLE["xy"], TITLE["text"], ha="center", va="baseline", fontsize=TITLE_PT,
            color=TEXT_COLOR)
    draw_landscape(ax)
    for text, xy in HEADING:
        draw_label(ax, xy, text, va="baseline")
    draw_mechanism(ax, FOUR_STATE, RED)
    draw_mechanism(ax, SIX_STATE, BLUE)
    save(fig, "scene19", "scene19")
