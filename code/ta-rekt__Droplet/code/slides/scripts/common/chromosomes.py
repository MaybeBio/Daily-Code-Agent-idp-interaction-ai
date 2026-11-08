"""
Chromosome cartoons traced from the storyboard (scene 2): three replicated
chromosomes, each made of two sister chromatids joined at the centromere.

Scene 2 draws them plain; scene 3 colours their DNA by which cell type uses it.
"""

import numpy as np
from matplotlib.patches import Polygon

from .drawing import EDGE_PT, draw_silhouette, placer, ribbon, smooth_open, unit_normals

# ---------------------------------------------------------------------------
# Palettes
# ---------------------------------------------------------------------------

GENE_COLORS = {
    "neuron":    "#D93A45",  # genes used by the neuron (red)
    "round":     "#2F62C8",  # genes used by the rounded cell (blue)
    "shared":    "#8E4FC0",  # genes used by both cells (purple = red + blue)
    "noncoding": "#3A3A3A",  # non-coding DNA (black)
}
PLAIN_COLORS = dict(fill="#DCDCE4", edge="#4A4A55")  # uncoloured DNA (scene 2)

CHROMATID_HALF_WIDTH = 0.12
CHROMOSOME_EDGE_PT = 0.7 * EDGE_PT
CENTROMERE_PINCH = 0.35  # chromatids narrow by this fraction at the centromere

# ---------------------------------------------------------------------------
# Geometry
# ---------------------------------------------------------------------------
# Each chromosome: centromere position on the slide; its two sister chromatids,
# each a path from its top end, through the centromere at (0, 0), to its bottom
# end (local units, y up); and its bands for scene 3, as (kind, end) pairs along
# a chromatid from top (0) to bottom (1), each band running from the previous
# end. Sister chromatids are copies of each other, so both get the same bands;
# centromeres and chromosome ends are non-coding.

CHROMOSOMES = [
    dict(centromere=(4.62, 6.73),
         chromatids=[
             [(-0.18, 1.31), (-0.22, 0.96), (-0.24, 0.63), (-0.22, 0.33), (-0.14, 0.15),
              (0.0, 0.0), (-0.07, -0.14), (-0.18, -0.4), (-0.25, -0.62), (-0.31, -0.88),
              (-0.35, -1.14), (-0.35, -1.39), (-0.31, -1.58), (-0.27, -1.83)],
             [(0.29, 1.24), (0.19, 0.93), (0.17, 0.66), (0.15, 0.38), (0.08, 0.15),
              (0.0, 0.0), (0.14, -0.14), (0.2, -0.33), (0.19, -0.62), (0.24, -0.88),
              (0.27, -1.17), (0.29, -1.5), (0.32, -1.83)],
         ],
         bands=[("noncoding", 0.04), ("neuron", 0.17), ("noncoding", 0.2), ("round", 0.33),
                ("shared", 0.38), ("noncoding", 0.47), ("neuron", 0.6), ("shared", 0.65),
                ("round", 0.8), ("noncoding", 0.84), ("neuron", 0.96), ("noncoding", 1.0)]),
    dict(centromere=(5.89, 6.35),
         chromatids=[
             [(-0.32, 0.6), (-0.29, 0.34), (-0.23, 0.12), (-0.12, 0.02), (0.0, 0.0),
              (-0.09, -0.24), (-0.22, -0.35), (-0.33, -0.49), (-0.4, -0.69), (-0.44, -1.06),
              (-0.39, -1.45)],
             [(0.0, 0.56), (0.02, 0.34), (0.01, 0.12), (0.0, 0.0), (0.1, -0.07),
              (0.15, -0.4), (0.19, -0.76), (0.22, -1.06), (0.21, -1.45)],
         ],
         bands=[("noncoding", 0.05), ("round", 0.2), ("noncoding", 0.36), ("neuron", 0.5),
                ("noncoding", 0.54), ("shared", 0.62), ("round", 0.78), ("neuron", 0.94),
                ("noncoding", 1.0)]),
    dict(centromere=(6.69, 6.89),
         chromatids=[
             [(-0.17, 0.8), (-0.16, 0.43), (-0.12, 0.17), (-0.06, 0.04), (0.0, 0.0),
              (-0.11, -0.12), (-0.22, -0.31), (-0.27, -0.56), (-0.25, -0.94), (-0.22, -1.3),
              (-0.25, -1.61)],
             [(0.2, 0.85), (0.15, 0.5), (0.08, 0.17), (0.0, 0.0), (0.09, -0.12),
              (0.17, -0.39), (0.19, -0.72), (0.14, -1.05), (0.13, -1.23), (0.19, -1.38)],
         ],
         bands=[("noncoding", 0.05), ("neuron", 0.2), ("shared", 0.26), ("noncoding", 0.4),
                ("round", 0.56), ("noncoding", 0.6), ("neuron", 0.74), ("shared", 0.8),
                ("round", 0.95), ("noncoding", 1.0)]),
]


def band_strip(curve, s, start, end, half_width):
    """Strip slightly wider than the tube over its stretch [start, end].

    It is clipped to the tube afterwards; keeping it only a little wider than the
    tube stops it from spilling onto the other side of sharp bends. Bands touching
    an end of the tube overshoot it, so the round cap is covered, and every band
    overlaps the next one slightly so no hairline seam shows between them.
    """
    ss = np.linspace(start, min(end + 0.005, 1.0), 25)
    piece = np.column_stack([np.interp(ss, s, curve[:, 0]), np.interp(ss, s, curve[:, 1])])
    if start == 0:
        t = (piece[1] - piece[0]) / np.linalg.norm(piece[1] - piece[0])
        piece = np.vstack([piece[0] - t * 1.5 * half_width, piece])
    if end >= 1:
        t = (piece[-1] - piece[-2]) / np.linalg.norm(piece[-1] - piece[-2])
        piece = np.vstack([piece, piece[-1] + t * 1.5 * half_width])
    n = unit_normals(piece) * 1.2 * half_width
    return np.vstack([piece + n, (piece - n)[::-1]])


def draw_chromosome(ax, chromosome, banded=False, z=0):
    """Draw one chromosome, plain (scene 2) or with coloured gene bands (scene 3)."""
    center = np.asarray(chromosome["centromere"])
    place = placer(center, 1.0, 0.0)
    tubes = []
    for path in chromosome["chromatids"]:
        curve, s = smooth_open(place(path))
        s_c = s[np.argmin(np.linalg.norm(curve - center, axis=1))]
        hw = CHROMATID_HALF_WIDTH * (1 - CENTROMERE_PINCH * np.exp(-((s - s_c) / 0.05) ** 2))
        tubes.append((ribbon(curve, hw), curve, s))
    outlines = [poly for poly, *_ in tubes]

    draw_silhouette(ax, outlines, PLAIN_COLORS["fill"], PLAIN_COLORS["edge"], z,
                    edge_pt=CHROMOSOME_EDGE_PT)

    if banded:
        for poly, curve, s, *_ in tubes:
            clip = Polygon(poly, closed=True, transform=ax.transData)
            start = 0.0
            for kind, end in chromosome["bands"]:
                band = Polygon(band_strip(curve, s, start, end, CHROMATID_HALF_WIDTH),
                               closed=True, fc=GENE_COLORS[kind], ec="none", zorder=z + 0.2)
                ax.add_patch(band)
                band.set_clip_path(clip)
                start = end
    return {"outlines": outlines}
