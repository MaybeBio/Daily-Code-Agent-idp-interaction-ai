"""
Cell cartoons traced from the storyboard (scene 1): a neuron and a rounded cell.

Each draw_* function takes a placement (origin, scale, angle) and a palette, and
returns {"outlines": [polygons], "nucleus": centre of the nucleus}, both in axes
coordinates, so scenes can crop to the cell or attach things to its nucleus.
"""

import numpy as np
from matplotlib.patches import Polygon

from .drawing import (EDGE_PT, LIGHT_DIR, draw_shine, draw_silhouette, draw_tube_shine,
                      placer, ribbon, smooth_closed, smooth_open, tapered_tube)

# ---------------------------------------------------------------------------
# Palettes and placement
# ---------------------------------------------------------------------------

NEURON_COLORS = dict(fill="#F9D9A8", edge="#C47A2C", shine="#FFFFFF")
ROUND_CELL_COLORS = dict(fill="#CFE8E4", edge="#2E8B84", shine="#FFFFFF",
                         organelle="#8CC7BF")
# Both nuclei share a colour on purpose: same DNA in every cell (scene 2).
NUCLEUS_COLORS = dict(fill="#CDCDCD", edge="#9D9D9D", nucleolus="#696969")

# Where the cells sit on the slide in scenes 1-3. Kept identical across those
# scenes (the sketches differ by < 1 mm) so the slides change without jumps.
NEURON_HOME = (6.4, 11.0)
ROUND_CELL_HOME = (14.95, 7.5)


def draw_nucleus(ax, outline, nucleolus_center, nucleolus_radius, z):
    draw_silhouette(ax, [outline], NUCLEUS_COLORS["fill"], NUCLEUS_COLORS["edge"], z,
                    edge_pt=0.8 * EDGE_PT)
    center = outline.mean(axis=0)
    draw_shine(ax, outline, center + 0.25 * (outline.max(0) - outline.min(0)) * LIGHT_DIR,
               "#FFFFFF", z + 0.2, layers=8, alpha=0.06)
    angle = np.linspace(0, 2 * np.pi, 40)
    blob = nucleolus_center + nucleolus_radius * np.column_stack(
        [np.cos(angle), 0.8 * np.sin(angle)])
    ax.add_patch(Polygon(blob, closed=True, fc=NUCLEUS_COLORS["nucleolus"], ec="none",
                         zorder=z + 0.3))


# ---------------------------------------------------------------------------
# Neuron
# ---------------------------------------------------------------------------
# Local frame: origin at the soma, 1 unit = 100 sketch px, y up.

NEURON_SOMA = [(-1.25, -0.42), (-0.9, 0.02), (-0.3, 0.25), (0.4, 0.42), (1.1, 0.45),
               (1.35, 0.2), (1.1, -0.2), (0.5, -0.33), (-0.2, -0.3), (-0.9, -0.48)]

NEURON_NUCLEUS = [(-0.7, -0.19), (-0.4, 0.0), (-0.2, 0.15), (0.12, 0.25), (0.45, 0.13),
                  (0.4, -0.05), (0.2, -0.08), (-0.2, -0.06), (-0.4, -0.07)]

# Each process: centerline from its base (inside the soma or a parent process,
# so that its rounded base end stays hidden) out to its tip, and half-widths at
# base and tip. `shine_from` is where its highlight stripe starts (fraction of
# its length), so stripes of side branches don't cross the parent process.
NEURON_PROCESSES = {
    "left_main":    dict(path=[(0.0, -0.05), (-1.4, -0.53), (-2.4, -0.85),
                               (-3.4, -1.08), (-4.35, -1.33)], base=0.15, tip=0.045,
                         shine_from=0.1),
    "left_upper":   dict(path=[(-0.75, -0.08), (-1.3, -0.15), (-1.9, -0.15),
                               (-2.62, -0.13)], base=0.08, tip=0.035, shine_from=0.25),
    "left_lower":   dict(path=[(-1.4, -0.55), (-2.2, -0.55), (-3.15, -0.65)],
                         base=0.08, tip=0.035, shine_from=0.12),
    "right_upper":  dict(path=[(0.6, 0.25), (1.0, 0.43), (1.6, 0.54), (2.2, 0.62),
                               (2.67, 0.67)], base=0.11, tip=0.04, shine_from=0.15),
    "right_middle": dict(path=[(1.0, 0.2), (1.6, 0.25), (2.3, 0.25), (2.65, 0.2)],
                         base=0.09, tip=0.045, shine_from=0.15),
    "right_long":   dict(path=[(0.8, -0.05), (1.6, -0.01), (2.4, -0.07), (3.2, -0.15),
                               (4.15, -0.21)], base=0.09, tip=0.04, shine_from=0.1),
    "right_lower":  dict(path=[(0.7, -0.1), (1.4, -0.33), (1.8, -0.48), (2.1, -0.65)],
                         base=0.1, tip=0.04, shine_from=0.2),
}


def draw_neuron(ax, origin=NEURON_HOME, scale=1.0, angle=0.0, colors=NEURON_COLORS, z=0):
    place = placer(origin, scale, angle)
    soma = place(smooth_closed(np.array(NEURON_SOMA)))
    processes = list(NEURON_PROCESSES.values())
    tubes = [tapered_tube(place(p["path"]), scale * p["base"], scale * p["tip"])
             for p in processes]
    shapes = [soma] + [poly for poly, *_ in tubes]

    draw_silhouette(ax, shapes, colors["fill"], colors["edge"], z)

    # A thin shine along the lit side of every process, and a soft glow on the soma.
    for p, (_, curve, s, core) in zip(processes, tubes):
        fade_in = np.clip((s - p["shine_from"]) / 0.15, 0, 1)
        draw_tube_shine(ax, curve, core, fade_in, z + 0.2, color=colors["shine"])
    draw_shine(ax, soma, place([0.2, 0.15]), colors["shine"], z + 0.2)

    nucleus = place(smooth_closed(np.array(NEURON_NUCLEUS)))
    draw_nucleus(ax, nucleus, place([0.05, 0.07]), 0.07 * scale, z + 1)
    return {"outlines": shapes, "nucleus": nucleus.mean(axis=0)}


# ---------------------------------------------------------------------------
# Rounded cell
# ---------------------------------------------------------------------------
# Local frame: origin at the cell centre, 1 unit = 100 sketch px, y up.

ROUND_CELL_OUTLINE = [
    (0.25, 1.8), (0.7, 1.88), (1.05, 1.81), (1.55, 1.65), (2.05, 1.45), (2.35, 1.35),
    (2.55, 1.2), (2.67, 0.95), (2.6, 0.7), (2.45, 0.35), (2.35, -0.05), (2.2, -0.35),
    (2.12, -0.65), (2.12, -0.95), (2.05, -1.2), (1.85, -1.35), (1.55, -1.5),
    (1.05, -1.67), (0.55, -1.77), (0.05, -1.8), (-0.35, -1.85), (-0.7, -1.82),
    (-1.2, -1.8), (-1.7, -1.7), (-2.05, -1.55), (-2.35, -1.3), (-2.6, -0.95),
    (-2.68, -0.55), (-2.65, -0.2), (-2.45, 0.15), (-2.05, 0.45), (-1.7, 0.65),
    (-1.2, 0.88), (-0.85, 1.15), (-0.65, 1.3), (-0.35, 1.38), (-0.15, 1.5), (0.05, 1.65),
]

ROUND_CELL_NUCLEUS = [
    (-0.75, 0.63), (-0.3, 0.58), (-0.05, 0.45), (0.04, 0.18), (-0.05, 0.03),
    (-0.3, -0.07), (-0.5, -0.2), (-0.6, -0.45), (-0.75, -0.6), (-0.85, -0.85),
    (-1.0, -0.97), (-1.25, -1.07), (-1.5, -1.12), (-1.7, -0.95), (-1.9, -0.77),
    (-1.95, -0.45), (-1.9, -0.2), (-1.65, 0.0), (-1.35, 0.25), (-1.1, 0.43),
    (-0.95, 0.55),
]

# The short strokes scattered in the sketch's cytoplasm, drawn as small capsules.
ROUND_CELL_ORGANELLES = [
    [(0.2, 1.2), (0.37, 1.39)], [(1.15, 1.33), (1.32, 1.15)], [(-0.68, 0.98), (-0.5, 1.03)],
    [(-1.74, 0.31), (-1.76, 0.15)], [(-2.34, -0.14), (-2.39, -0.2)],
    [(1.09, 0.44), (1.01, 0.26)], [(1.95, 0.66), (1.72, 0.23)], [(0.44, 0.09), (0.47, -0.17)],
    [(1.18, -0.35), (1.26, -0.36)], [(0.12, -0.8), (0.06, -0.86)],
    [(1.17, -0.86), (1.32, -0.92), (1.45, -0.86)], [(1.64, -1.05), (1.77, -1.06), (1.85, -0.91)],
    [(-2.0, -1.1), (-1.89, -1.15)], [(-1.51, -1.41), (-1.29, -1.47)],
    [(-0.6, -1.07), (-0.38, -1.01)], [(0.24, -1.19), (0.39, -1.22)],
    [(-0.45, -1.55), (-0.34, -1.55)],
]
ORGANELLE_HALF_WIDTH = 0.055


def draw_round_cell(ax, origin=ROUND_CELL_HOME, scale=1.0, angle=0.0,
                    colors=ROUND_CELL_COLORS, z=0):
    place = placer(origin, scale, angle)
    outline = place(smooth_closed(np.array(ROUND_CELL_OUTLINE)))
    draw_silhouette(ax, [outline], colors["fill"], colors["edge"], z)
    draw_shine(ax, outline, place([-0.3, 0.6]), colors["shine"], z + 0.2)

    hw = scale * ORGANELLE_HALF_WIDTH
    for stroke in ROUND_CELL_ORGANELLES:
        curve, _ = smooth_open(place(stroke), n=30)
        capsule = ribbon(curve, np.full(len(curve), hw))
        draw_silhouette(ax, [capsule], colors["organelle"], colors["edge"], z + 0.5,
                        edge_pt=0.4 * EDGE_PT)

    nucleus = place(smooth_closed(np.array(ROUND_CELL_NUCLEUS)))
    draw_nucleus(ax, nucleus, place([-1.2, -0.35]), 0.15 * scale, z + 1)
    return {"outlines": [outline], "nucleus": nucleus.mean(axis=0)}
