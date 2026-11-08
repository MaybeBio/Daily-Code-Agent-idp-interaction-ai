"""
Shared style, geometry helpers, drawing primitives and export for the slide
illustrations (scripts/scene*/scene*.py).

Coordinates are in "drawing units", y up. One unit = 100 px of a 2000-px-wide
render of the scene-1 sketch; a 4:3 slide is 20 x 15 units, with (10, 7.5) at
its centre.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import ArrowStyle, FancyArrowPatch, PathPatch, Polygon
from matplotlib.path import Path as DrawnPath
from scipy.interpolate import splev, splprep

FIGURES_DIR = Path(__file__).resolve().parents[2] / "figures"


# ---------------------------------------------------------------------------
# Style
# ---------------------------------------------------------------------------

UNIT_MM = 6.4                       # 20 units across a 128 mm (beamer 4:3) slide
SLIDE_MM = (128, 96)
SLIDE_UNITS = (SLIDE_MM[0] / UNIT_MM, SLIDE_MM[1] / UNIT_MM)

EDGE_PT = 1.0                       # visible outline width of shapes, in points
LIGHT_DIR = np.array([-0.5, 0.85])  # highlights are lit from the upper left
ARROW_PT = 3.0

# Labels are in FreeSans, like the slide text: matplotlib sets them itself,
# equations included (mathtext, so no LaTeX commands beyond math), and they stay
# real text in the SVGs, editable in Inkscape. FreeSans looks larger than
# Computer Modern Sans at the same size, so 9.2 pt reads like the old 11 pt.
LABEL_PT = 9.2
TEXT_COLOR = "#222222"
AXIS_COLOR = "#3A3A3A"
BOLD = 600                          # FreeSans's bold declares itself semi-bold
plt.rcParams.update({
    "font.family": "FreeSans",
    "mathtext.fontset": "custom",
    "mathtext.rm": "FreeSans",
    "mathtext.it": "FreeSans:italic",
    "mathtext.bf": f"FreeSans:weight={BOLD}",
    "mathtext.sf": "FreeSans",
    "mathtext.cal": "FreeSans:italic",
    "mathtext.fallback": "stixsans",   # for what FreeSans lacks: tall brackets
    "pdf.fonttype": 42,                # TrueType in the PDFs...
    "svg.fonttype": "none",            # ...and <text> in the SVGs
})


# ---------------------------------------------------------------------------
# Geometry helpers
# ---------------------------------------------------------------------------

def smooth_closed(points, n=300, smoothing=0.002):
    """Smooth periodic spline through hand-traced outline points."""
    pts = np.vstack([points, points[:1]])  # splprep(per=True) wants the loop closed
    tck, _ = splprep(pts.T, s=smoothing * len(pts), per=True)
    return np.column_stack(splev(np.linspace(0, 1, n), tck))


def smooth_open(points, n=150):
    """Smooth curve through points; returns (curve, arclength fraction 0..1)."""
    pts = np.asarray(points, float)
    tck, _ = splprep(pts.T, s=0, k=min(3, len(pts) - 1))
    curve = np.column_stack(splev(np.linspace(0, 1, n), tck))
    step = np.linalg.norm(np.diff(curve, axis=0), axis=1)
    s = np.concatenate([[0], np.cumsum(step)])
    return curve, s / s[-1]


def unit_normals(curve):
    d = np.gradient(curve, axis=0)
    t = d / np.linalg.norm(d, axis=1, keepdims=True)
    return np.column_stack([-t[:, 1], t[:, 0]])


def ribbon(curve, half_width):
    """Closed polygon of a band of varying half-width along a curve, with round ends."""
    n = unit_normals(curve)
    left = curve + n * half_width[:, None]
    right = curve - n * half_width[:, None]

    def half_circle(center, start_normal, r):
        # from +start_normal round to -start_normal, bulging 90 degrees clockwise of it
        a = np.arctan2(start_normal[1], start_normal[0]) - np.linspace(0, np.pi, 16)
        return center + r * np.column_stack([np.cos(a), np.sin(a)])

    tip_cap = half_circle(curve[-1], n[-1], half_width[-1])
    base_cap = half_circle(curve[0], -n[0], half_width[0])
    return np.vstack([left, tip_cap, right[::-1], base_cap])


def tapered_tube(path, base, tip, flare=0.3):
    """Tube along `path` tapering from half-width `base` to `tip`.

    `flare` widens the first ~10% so the tube blends into what it grows from.
    Returns (polygon, centerline, arclength fraction, half-width without flare).
    """
    curve, s = smooth_open(path)
    core = tip + (base - tip) * (1 - s)
    return ribbon(curve, core + flare * base * (1 - s) ** 12), curve, s, core


def placer(origin, scale, angle):
    """Map local coordinates to axes coordinates (rotate, scale, translate)."""
    c, s = np.cos(np.radians(angle)), np.sin(np.radians(angle))
    rot = np.array([[c, -s], [s, c]])
    return lambda pts: np.asarray(pts) @ rot.T * scale + np.asarray(origin)


# ---------------------------------------------------------------------------
# Drawing primitives
# ---------------------------------------------------------------------------

def draw_silhouette(ax, shapes, fill, edge, z, edge_pt=EDGE_PT):
    """Fill overlapping shapes as one silhouette with a single outer outline.

    Every shape is first stroked at twice the outline width; filling all shapes
    on top then hides the inner half of each stroke and every stroke that runs
    inside another shape, so only the outer edge of the union remains.
    """
    for pts in shapes:
        ax.add_patch(Polygon(pts, closed=True, fc=edge, ec=edge, lw=2 * edge_pt,
                             joinstyle="round", zorder=z))
    for pts in shapes:
        ax.add_patch(Polygon(pts, closed=True, fc=fill, ec="none", zorder=z + 0.1))


def draw_shine(ax, shape, focus, color, z, layers=14, alpha=0.05, smallest=0.1):
    """Soft highlight: stacked translucent copies of `shape` shrinking toward `focus`.

    Gives a radial-gradient look while staying pure vector (no embedded bitmap).
    """
    for f in np.linspace(0.97, smallest, layers):
        ax.add_patch(Polygon(focus + (shape - focus) * f, closed=True, fc=color,
                             ec="none", alpha=alpha, zorder=z))


def draw_tube_shine(ax, curve, half_width, weight, z, color="#FFFFFF", alpha=0.55):
    """Thin highlight stripe along the lit side of a tube, for a rounded 3D look.

    `weight` (0..1 per centerline point) scales the stripe along the tube, e.g. to
    fade it in from a base; wherever it is 0 the stripe is interrupted.
    """
    n = unit_normals(curve)
    side = np.sign((n @ LIGHT_DIR).mean())
    center = curve + side * n * 0.35 * half_width[:, None]
    width = 0.3 * half_width * weight
    on = np.r_[0, (width > 0).astype(int), 0]
    starts, stops = np.flatnonzero(np.diff(on) == 1), np.flatnonzero(np.diff(on) == -1)
    for a, b in zip(starts, stops):
        if b - a > 1:
            ax.add_patch(Polygon(ribbon(center[a:b], width[a:b]), closed=True, fc=color,
                                 ec="none", alpha=alpha, zorder=z))


def draw_arrow(ax, start, end, color, z, lw=ARROW_PT, head=12):
    ax.add_patch(FancyArrowPatch(
        start, end, arrowstyle=ArrowStyle("-|>", head_length=0.5, head_width=0.3),
        mutation_scale=head, lw=lw, color=color, capstyle="round",
        joinstyle="miter", shrinkA=0, shrinkB=0, zorder=z))


def draw_label(ax, xy, text, bold=False, ha="center", va="center", size=LABEL_PT, z=20):
    """Text label; `xy` is its centre (ha="center") or its left end (ha="left")."""
    ax.text(*xy, text, ha=ha, va=va, fontsize=size, fontweight=BOLD if bold else "normal",
            color=TEXT_COLOR, zorder=z)


def draw_underbrace(ax, x_range, y, label, color=TEXT_COLOR, depth=0.26, size=LABEL_PT - 1, z=20):
    """A curly brace under `x_range`, its arms up at `y`, with `label` under its tip."""
    x0, x1 = x_range
    middle, bend, low = (x0 + x1) / 2, depth / 2, y - depth / 2
    points = [(x0, y), (x0, low), (x0 + bend, low), (middle - bend, low), (middle, low),
              (middle, y - depth), (middle, low), (middle + bend, low), (x1 - bend, low),
              (x1, low), (x1, y)]
    codes = [DrawnPath.MOVETO, DrawnPath.CURVE3, DrawnPath.CURVE3, DrawnPath.LINETO,
             DrawnPath.CURVE3, DrawnPath.CURVE3, DrawnPath.CURVE3, DrawnPath.CURVE3,
             DrawnPath.LINETO, DrawnPath.CURVE3, DrawnPath.CURVE3]
    ax.add_patch(PathPatch(DrawnPath(points, codes), fc="none", ec=color, lw=0.8, zorder=z))
    ax.text(middle, y - depth - 0.08, label, ha="center", va="top", fontsize=size, color=color,
            zorder=z)


def draw_plot_axes(ax, origin, tips, y_label=None, x_label=None):
    """Axes of a sketched plot: two arrows from `origin` to `tips` (x of the
    horizontal arrow's tip, y of the vertical one's), with `y_label` right of
    the vertical tip and `x_label` right of the horizontal one."""
    (x0, y0), (x_tip, y_tip) = origin, tips
    draw_arrow(ax, (x0, y0), (x_tip, y0), AXIS_COLOR, z=3, lw=1.3, head=9)
    draw_arrow(ax, (x0, y0), (x0, y_tip), AXIS_COLOR, z=3, lw=1.3, head=9)
    if y_label:
        draw_label(ax, (x0 + 0.27, y_tip), y_label, ha="left")
    if x_label:
        draw_label(ax, (x_tip + 0.25, y0), x_label, ha="left")


# ---------------------------------------------------------------------------
# Export
# ---------------------------------------------------------------------------

def new_canvas(width_units, height_units):
    """Figure whose axes span it exactly, at UNIT_MM per drawing unit."""
    fig = plt.figure(figsize=(width_units * UNIT_MM / 25.4, height_units * UNIT_MM / 25.4))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_axis_off()
    ax.set_xlim(0, width_units)
    ax.set_ylim(0, height_units)
    return fig, ax


def new_slide():
    """Canvas the size of a full 4:3 slide (128 x 96 mm)."""
    return new_canvas(*SLIDE_UNITS)


def zoomed_axes(fig, shift, zoom):
    """Axes over the whole slide on which a point p of a drawing lands at
    shift + zoom * p (slide units): the drawing, details included (helix,
    nucleosomes, proteins), is drawn `zoom` times its size. Line widths, in points,
    stay."""
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_axis_off()
    ax.set_xlim(-shift[0] / zoom, (SLIDE_UNITS[0] - shift[0]) / zoom)
    ax.set_ylim(-shift[1] / zoom, (SLIDE_UNITS[1] - shift[1]) / zoom)
    return ax


def save(fig, scene, name):
    """Write figures/<scene>/<name>.pdf (for beamer) and .svg (for Inkscape)."""
    out_dir = FIGURES_DIR / scene
    out_dir.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "svg"):
        # No bbox_inches="tight": the page size is set exactly by new_canvas.
        fig.savefig(out_dir / f"{name}.{ext}", transparent=True)
    plt.close(fig)
    print(f"wrote {out_dir / name}.pdf/.svg")


def export_single(draw, scene, name, pad=0.15):
    """One object on its own page, cropped to its outlines plus `pad` units.

    `draw(ax)` must return a dict with the drawn outline polygons under "outlines".
    """
    fig, ax = new_canvas(1, 1)
    shapes = draw(ax)["outlines"]
    lo = np.min([s.min(axis=0) for s in shapes], axis=0) - pad
    hi = np.max([s.max(axis=0) for s in shapes], axis=0) + pad
    ax.set_xlim(lo[0], hi[0])
    ax.set_ylim(lo[1], hi[1])
    fig.set_size_inches(*((hi - lo) * UNIT_MM / 25.4))
    save(fig, scene, name)
