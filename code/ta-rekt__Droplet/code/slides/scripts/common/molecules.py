"""
Molecules for the gene-expression scenes: DNA (chromatin), RNA, RNA polymerase,
proteins.

DNA and RNA are `Polymer`s: chains of backbone nodes. Chunk i is the short
stretch of chain between nodes i and i+1, and every chunk has a kind ("coil"
for compacted chromatin, "linker" for open DNA, "gene", "rna", ...).

Everything drawn along a chain -- the double helix, nucleosomes, RNA beads, a
polymerase riding on the DNA -- is computed from the node positions and placed
by node index (the chain's material coordinate), never by absolute position.
To animate thermal motion later, move the nodes (e.g. a small random-walk or
Rouse update per frame) and redraw: every detail follows its chunk.

DNA is drawn after slides/figures/dna.png, and transcription factors are exact
copies of slides/figures/tf.svg.
"""

import re
import xml.etree.ElementTree as ET

import numpy as np
from matplotlib.colors import to_rgba
from matplotlib.patches import Circle, PathPatch, Polygon
from matplotlib.path import Path
from scipy.interpolate import splev, splprep
from scipy.ndimage import gaussian_filter1d

from .drawing import (EDGE_PT, FIGURES_DIR, UNIT_MM, draw_label, draw_silhouette, placer,
                      smooth_closed, smooth_open, unit_normals)

# ---------------------------------------------------------------------------
# Palettes and sizes (drawing units unless stated)
# ---------------------------------------------------------------------------

BLUE = dict(fill="#2F62C8", edge="#1D3F85")  # same blue as the rounded cell's genes
GREEN = dict(fill="#07B503", edge="#2B6E31")  # transcription factor, its gene and mRNA
ORANGE = dict(fill="#F39A2B", edge="#B5651D")  # promoter
RNAP_COLORS = dict(fill="#D957DB", edge="#8E1F90")
MOLECULE_EDGE_PT = 0.7 * EDGE_PT

DNA_BOND = 0.11            # node spacing along DNA
# The double helix is a ribbon waving across the chain, an outline around a fill.
DNA_LINE = dict(edge="#2B2B2B", fill="white")
DNA_LINE_PT = (1.3, 0.5)   # width (points) of the ribbon, and of the fill inside its outline
HELIX_TURN = 0.2           # length of chain per turn...
HELIX_RADIUS = 0.1         # ...and how far the ribbon swings to either side of it
HELIX_SHARPNESS = 0.93     # 0: a sine wave; towards 1: straight legs, tight turns
SEGMENT_COLORS = {"gene": BLUE, "promoter": ORANGE}  # chunk kinds drawn in colour

# Nucleosomes: neutral, so that colour along the chain only ever means a gene or
# a promoter. "fill" is the histone core, "shade" what lies behind the DNA.
NUCLEOSOME_COLORS = dict(fill="#56657E", shade="#CFD6E2")
NUCLEOSOME_SPACING = 1.25  # along "coil" stretches, between nucleosomes
WRAP_RADIUS = 0.3          # seen face-on: the DNA runs round the core at this radius...
CORE_RADIUS = 0.2          # ...which has this radius
CYLINDER = dict(half_width=0.31, half_height=0.23, cap=0.11)  # seen from the side (see side_on)
STRAIGHTEST_WITHIN = 0.4   # a side view looks this far along the chain for a straight stretch
MIN_STRAIGHTNESS = 0.9     # chord / chain length it needs there; else, a face-on view
NUCLEOSOME_EXTENT = WRAP_RADIUS + HELIX_RADIUS  # how far a nucleosome reaches from its centre

RNA_BOND = 0.15            # node spacing along RNA = distance between beads
RNA_BEAD_RADIUS = 0.052

TF_HEIGHT = 0.38           # transcription factor, drawn at scale 1
TF_RADIUS = 0.2            # outer radius of the same


# ---------------------------------------------------------------------------
# Polymer chains
# ---------------------------------------------------------------------------

def resample(points, spacing):
    """Points evenly spaced (as close to `spacing` as fits) along a smooth curve
    through `points`, including both ends."""
    curve, s = smooth_open(points, n=400)
    length = np.linalg.norm(np.diff(curve, axis=0), axis=1).sum()
    ss = np.linspace(0, 1, max(1, round(length / spacing)) + 1)
    return np.column_stack([np.interp(ss, s, curve[:, 0]), np.interp(ss, s, curve[:, 1])])


def runs(mask):
    """(start, stop) index pairs of the contiguous True stretches of a boolean array."""
    edges = np.diff(np.r_[0, np.asarray(mask, int), 0])
    return list(zip(np.flatnonzero(edges == 1), np.flatnonzero(edges == -1)))


def smoothstep(x):
    return x * x * (3 - 2 * x)


def circle_points(center, radius, n=24):
    """Points around a circle, e.g. to keep growing chains clear of a protein."""
    a = np.linspace(0, 2 * np.pi, n, endpoint=False)
    return np.asarray(center) + radius * np.column_stack([np.cos(a), np.sin(a)])


class Polymer:
    """A chain of backbone nodes; chunk i joins nodes i and i+1 and has kind kinds[i]."""

    def __init__(self, nodes, kinds):
        self.nodes = np.asarray(nodes, float)        # (N, 2): the state to move when animating
        self.kinds = np.asarray(kinds, dtype=object)  # (N - 1,); editable, e.g. to mark a promoter

    @classmethod
    def from_trace(cls, pieces, bond):
        """Chain from traced pieces [(kind, points), ...] laid end to end.

        Each piece is resampled to nodes `bond` apart; consecutive pieces share
        their joining node.
        """
        nodes, kinds = [], []
        for kind, points in pieces:
            pts = resample(points, bond)
            nodes.extend(pts if not nodes else pts[1:])
            kinds.extend([kind] * (len(pts) - 1))
        return cls(nodes, kinds)

    def smoothed(self, nodes=1.5):
        """The same chain with kinks of its hand-traced backbone ironed out, over
        a few `nodes`, so that the helix doesn't knot around them."""
        return Polymer(gaussian_filter1d(self.nodes, nodes, axis=0, mode="nearest"), self.kinds)

    def spline(self):
        u = np.arange(len(self.nodes), dtype=float)
        return splprep(self.nodes.T, u=u, s=0, k=min(3, len(self.nodes) - 1))[0]

    def curve(self, per_chunk=16):
        """Smooth curve through the nodes: (points, u), u = fractional node index."""
        u = np.linspace(0, len(self.nodes) - 1, (len(self.nodes) - 1) * per_chunk + 1)
        return np.column_stack(splev(u, self.spline())), u

    def at(self, u):
        """Position and unit tangent at fractional node index u."""
        tck = self.spline()
        d = np.array(splev(u, tck, der=1))
        return np.array(splev(u, tck)), d / np.linalg.norm(d)

    def kind_of(self, u):
        """Kind of the chunk containing each fractional node index in u."""
        return self.kinds[np.minimum(np.asarray(u).astype(int), len(self.kinds) - 1)]

    def span(self, kind):
        """(first, last) node index of the stretch of chunks of this kind."""
        idx = np.flatnonzero(self.kinds == kind)
        return float(idx[0]), float(idx[-1] + 1)

    def nearest(self, point):
        """Fractional node index of the point of the chain closest to `point`."""
        pts, u = self.curve()
        return u[np.argmin(np.linalg.norm(pts - np.asarray(point), axis=1))]

    def cut(self, length):
        """The first nodes of the chain, as many as stay within `length` along it."""
        along = np.r_[0, np.cumsum(np.linalg.norm(np.diff(self.nodes, axis=0), axis=1))]
        n = np.searchsorted(along, length, side="right")
        return Polymer(self.nodes[:n], self.kinds[:n - 1])


# ---------------------------------------------------------------------------
# DNA
# ---------------------------------------------------------------------------

def bridge(p, p_tangent, q, q_tangent, n=24):
    """Smooth curve from p to q, leaving and arriving along the given unit tangents."""
    s = np.linspace(0, 1, n)[:, None]
    reach = np.linalg.norm(q - p)
    return ((2 * s ** 3 - 3 * s ** 2 + 1) * p + (s ** 3 - 2 * s ** 2 + s) * reach * p_tangent
            + (3 * s ** 2 - 2 * s ** 3) * q + (s ** 3 - s ** 2) * reach * q_tangent)


def face_on(dna, at):
    """A nucleosome seen along its axis at node index `at`: the chain runs once
    round the core, on the side of the chain with more room, and carries on."""
    p, t = dna.at(at)
    left = np.array([-t[1], t[0]])
    chain = dna.curve()[0]
    others = chain[np.linalg.norm(chain - p, axis=1) > 2 * NUCLEOSOME_EXTENT]

    def room(side):
        return np.linalg.norm(others - (p + side * WRAP_RADIUS * left), axis=1).min()

    side = max((1, -1), key=room)
    center = p + side * WRAP_RADIUS * left
    start = np.arctan2(*(p - center)[::-1])
    angle = start + side * np.linspace(0, 2 * np.pi, 160)
    wrap = center + WRAP_RADIUS * np.column_stack([np.cos(angle), np.sin(angle)])
    return dict(center=center, enter=at, leave=at, wraps=[wrap],
                shapes=[(circle_points(center, NUCLEOSOME_EXTENT, 48), "shade"),
                        (circle_points(center, CORE_RADIUS, 48), "fill")])


def side_on_span(dna):
    """Nodes the chain gives up on each side of a nucleosome seen from the side."""
    return (CYLINDER["half_width"] + 0.2) / DNA_BOND


def straightness(dna, at):
    """Chord / chain length over the stretch a side view would replace."""
    span = side_on_span(dna)
    if at - span < 0 or at + span > len(dna.nodes) - 1:
        return 0.0
    chord = np.linalg.norm(dna.at(at + span)[0] - dna.at(at - span)[0])
    return chord / (2 * span * DNA_BOND)


def side_on(dna, at, tilt, cap_side):
    """A nucleosome seen from the side at node index `at`: a cylinder across the
    chain, turned by `tilt` degrees, showing its end cap on `cap_side` (+1 or -1).

    The chain leaves its course, runs across the front of the cylinder, passes
    behind it (out of sight) and crosses the front again, a little further up the
    cylinder, before rejoining its course.
    """
    w, h, e = CYLINDER["half_width"], CYLINDER["half_height"], CYLINDER["cap"]
    span = side_on_span(dna)
    (p, p_tangent), (q, q_tangent) = dna.at(at - span), dna.at(at + span)
    center = (p + q) / 2
    place = placer(center, 1.0, np.degrees(np.arctan2(*(q - p)[::-1])) + tilt)

    def local(x, y):
        return place(np.column_stack([x, cap_side * y]))

    x = np.linspace(-w, w, 60)
    arc = 1 - (x / w) ** 2                 # the front of a circle round the cylinder
    lower, upper = local(x, -0.19 - 0.7 * e * arc), local(x, 0.05 - 0.7 * e * arc)
    along = place([1.0, 0.0]) - center
    wraps = [np.vstack([bridge(p, p_tangent, lower[0], along), lower]),
             np.vstack([upper, bridge(upper[-1], along, q, q_tangent)])]

    around = np.linspace(0, 2 * np.pi, 48)
    cap = local(w * np.cos(around), h + e * np.sin(around))
    body = np.vstack([cap[:24], local(x, -h - e * arc)])  # the far half of the cap, then round
    return dict(center=center, enter=at - span, leave=at + span, wraps=wraps,
                shapes=[(body, "shade"), (cap, "fill")])


# The nucleosomes along a chain take these views in turn; side views show their
# cylinders at various angles, as in dna.png.
NUCLEOSOME_VIEWS = [None, dict(tilt=0, cap_side=1), None, dict(tilt=-28, cap_side=-1),
                    None, dict(tilt=24, cap_side=1)]


def nucleosomes(dna):
    """Nucleosomes spread evenly over the "coil" stretches of a chain, in node
    order, each a dict: its `center`, the node indices where the chain `enter`s
    and `leave`s it, the `wraps` of DNA around it and its `shapes` (polygons
    with their tone in NUCLEOSOME_COLORS)."""
    _, u = dna.curve()
    sites = []
    for a, b in runs(dna.kind_of(u) == "coil"):
        nodes = u[b - 1] - u[a]
        count = int(nodes * DNA_BOND / NUCLEOSOME_SPACING)
        for k in range(count):
            at = u[a] + (k + 0.5) * nodes / count
            view = NUCLEOSOME_VIEWS[len(sites) % len(NUCLEOSOME_VIEWS)]
            if view:
                reach = STRAIGHTEST_WITHIN / DNA_BOND
                best = max(np.linspace(at - reach, at + reach, 9),
                           key=lambda v: straightness(dna, v))
                if straightness(dna, best) >= MIN_STRAIGHTNESS:
                    sites.append(side_on(dna, best, **view))
                    continue
            sites.append(face_on(dna, at))
    return sites


def dna_pieces(dna, sites):
    """The axis of the helix as drawn: the chain, detouring by the wraps of the
    nucleosomes `sites`. A list of unbroken pieces (points, kind of each point);
    the chain breaks where it passes behind a nucleosome."""
    center, u = dna.curve()
    kinds = dna.kind_of(u)
    pieces, points, labels, last = [], [], [], -1.0
    for site in sites:
        keep = (u > last) & (u < site["enter"])
        first, *others = site["wraps"]
        points += [center[keep], first]
        labels += [kinds[keep], np.full(len(first), "coil", dtype=object)]
        for wrap in others:
            pieces.append((np.vstack(points), np.concatenate(labels)))
            points, labels = [wrap], [np.full(len(wrap), "coil", dtype=object)]
        last = site["leave"]
    keep = u > last
    pieces.append((np.vstack(points + [center[keep]]), np.concatenate(labels + [kinds[keep]])))
    return pieces


def helix(axis):
    """The ribbon of the double helix around a curve `axis`: it swings from one
    side of the axis to the other, one full swing per HELIX_TURN."""
    s = np.r_[0, np.cumsum(np.linalg.norm(np.diff(axis, axis=0), axis=1))]
    swing = np.arcsin(HELIX_SHARPNESS * np.sin(2 * np.pi * s / HELIX_TURN)) / np.arcsin(HELIX_SHARPNESS)
    return axis + HELIX_RADIUS * swing[:, None] * unit_normals(axis)


def draw_ribbon(ax, points, colors, line_pt, z):
    """A line as wide as line_pt[0], its outline in colors["edge"] and the
    line_pt[1] wide middle in colors["fill"]."""
    for color, width in zip((colors["edge"], colors["fill"]), line_pt):
        ax.plot(*np.transpose(points), color=color, lw=width, solid_capstyle="round",
                solid_joinstyle="round", zorder=z)


def draw_dna(ax, dna, segment_colors=SEGMENT_COLORS, z=0, line_colors=DNA_LINE, line_pt=DNA_LINE_PT,
             nucleosome_colors=NUCLEOSOME_COLORS):
    """A double helix along the chain, in `line_colors` except over chunks whose
    kind is in `segment_colors` (e.g. the gene), with nucleosomes on its "coil"
    stretches. Lighter `line_colors` and `nucleosome_colors` dim the chromatin;
    narrower `line_pt` keep the helix open on a zoomed-out drawing.
    """
    dna = dna.smoothed()
    sites = nucleosomes(dna)
    for site in sites:
        for shape, tone in site["shapes"]:
            ax.add_patch(Polygon(shape, closed=True, fc=nucleosome_colors[tone], ec="none",
                                 zorder=z - 0.1))
    for axis, kinds in dna_pieces(dna, sites):
        strand = helix(axis)
        for kind in set(kinds):
            for a, b in runs(kinds == kind):
                draw_ribbon(ax, strand[a:b + 1], segment_colors.get(kind, line_colors), line_pt, z)


def dna_obstacles(dna, clearance=0.2):
    """What scattered molecules must keep clear of, as (points, radius) pairs."""
    return [(dna.curve()[0], clearance)] + [(site["center"], NUCLEOSOME_EXTENT)
                                            for site in nucleosomes(dna.smoothed())]


# ---------------------------------------------------------------------------
# RNA
# ---------------------------------------------------------------------------

def draw_rna(ax, rna, color=BLUE["fill"], z=0):
    """RNA as a string of beads, one bead per node."""
    for p in rna.nodes:
        ax.add_patch(Circle(p, RNA_BEAD_RADIUS, fc=color, ec="none", zorder=z))


def grow_rna(start, direction, n_beads, rng, avoid=(), clearance=0.25, wiggle=0.45,
             pull=0.12, tries=40):
    """A transcript as a random, self-avoiding worm-like chain of `n_beads` beads,
    growing from `start` (the polymerase) roughly along `direction`.

    Each step turns by a random angle (`wiggle`) with a weak `pull` back to the
    initial heading; steps that land within `clearance` of the `avoid` points
    (e.g. the DNA) or on the chain itself are redrawn, with wider and wider
    turns, so a blocked chain bends around the obstacle.
    """
    avoid = np.asarray(avoid, float).reshape(-1, 2)
    goal = heading = np.arctan2(direction[1], direction[0])
    nodes = [np.asarray(start, float)]
    for i in range(n_beads - 1):
        for attempt in range(tries):
            spread = wiggle * (1 + attempt / 4)
            turn = heading + spread * rng.standard_normal() + pull * (goal - heading)
            p = nodes[-1] + RNA_BOND * np.array([np.cos(turn), np.sin(turn)])
            clear = (i < 3 or len(avoid) == 0
                     or np.linalg.norm(avoid - p, axis=1).min() > clearance)
            self_clear = (len(nodes) < 3
                          or np.linalg.norm(np.array(nodes[:-2]) - p, axis=1).min() > 1.5 * RNA_BOND)
            if clear and self_clear:
                break
        heading = turn
        nodes.append(p)
    return Polymer(nodes, ["rna"] * (len(nodes) - 1))


# ---------------------------------------------------------------------------
# Proteins
# ---------------------------------------------------------------------------

# RNA polymerase: a lumpy globule around its centre, in a frame where the DNA
# under it runs along +x; it turns with the DNA as it moves along it.
RNAP_OUTLINE = [(0.37, 0.0), (0.26, 0.19), (0.09, 0.29), (-0.11, 0.33), (-0.27, 0.19),
                (-0.32, 0.0), (-0.29, -0.21), (-0.1, -0.31), (0.09, -0.28), (0.28, -0.2)]


def draw_rnap_at(ax, center, angle=0.0, colors=RNAP_COLORS, scale=1.0, alpha=1.0, z=0):
    """RNA polymerase at `center`, turned by `angle` degrees (e.g. floating free)."""
    blob = placer(center, scale, angle)(smooth_closed(np.array(RNAP_OUTLINE)))
    if alpha < 1:  # one translucent patch (the silhouette trick would show its layers)
        ax.add_patch(Polygon(blob, closed=True, fc=colors["fill"], ec=colors["edge"],
                             lw=MOLECULE_EDGE_PT, alpha=alpha, zorder=z))
    else:
        draw_silhouette(ax, [blob], colors["fill"], colors["edge"], z, edge_pt=MOLECULE_EDGE_PT)


def draw_rnap(ax, dna, u, colors=RNAP_COLORS, z=0):
    """RNA polymerase sitting on the DNA at fractional node index u."""
    p, t = dna.at(u)
    draw_rnap_at(ax, p, np.degrees(np.arctan2(t[1], t[0])), colors=colors, z=z)
    return p


# ---------------------------------------------------------------------------
# Transcription factors: exact copies of figures/tf.svg
# ---------------------------------------------------------------------------

def svg_path(d):
    """Vertices and codes of an SVG path of moveto, lineto and cubic Bezier
    commands, kept as curves. Lower-case commands are relative."""
    tokens = re.findall(r"[a-zA-Z]|-?\d*\.?\d+(?:e-?\d+)?", d)
    vertices, codes, at, i, command = [], [], np.zeros(2), 0, None
    while i < len(tokens):
        if tokens[i].isalpha():
            command = tokens[i]
            i += 1
        n = 6 if command in "cC" else 2
        points = np.array(tokens[i:i + n], float).reshape(-1, 2)
        i += n
        points = points + at if command.islower() else points
        vertices.extend(points)
        codes.extend({"m": [Path.MOVETO], "l": [Path.LINETO], "c": [Path.CURVE4] * 3}[
            command.lower()])
        at = points[-1]
        if command in "mM":  # further pairs after a moveto are linetos
            command = "l" if command == "m" else "L"
    return np.array(vertices), codes


def load_tf():
    """The outline of tf.svg, y up, centred, 1 tall, and its colours."""
    shape = ET.parse(FIGURES_DIR / "tf.svg").getroot().find(".//{http://www.w3.org/2000/svg}path")
    style = dict(item.split(":") for item in shape.get("style").split(";"))
    vertices, codes = svg_path(shape.get("d"))
    vertices = vertices * (1, -1)
    low, high = vertices.min(axis=0), vertices.max(axis=0)
    height = high[1] - low[1]
    return dict(vertices=(vertices - (low + high) / 2) / height, codes=codes,
                fill=to_rgba(style["fill"], float(style["fill-opacity"])),
                edge=to_rgba(style["stroke"], float(style["stroke-opacity"])),
                edge_width=float(style["stroke-width"]) / height)


TF = load_tf()
POINTS_PER_UNIT = UNIT_MM / 25.4 * 72


def draw_tf(ax, center, scale=1.0, angle=0.0, alpha=1.0, z=0, colors=None):
    """A transcription factor at `center`, turned by `angle` degrees; `colors`
    (fill and edge) repaint it, e.g. as a protein of another gene."""
    height = scale * TF_HEIGHT
    fill, edge = (to_rgba(colors[key]) if colors else TF[key] for key in ("fill", "edge"))
    fill, edge = (color[:3] + (color[3] * alpha,) for color in (fill, edge))
    ax.add_patch(PathPatch(Path(placer(center, height, angle)(TF["vertices"]), TF["codes"]),
                           fc=fill, ec=edge, lw=TF["edge_width"] * height * POINTS_PER_UNIT,
                           joinstyle="round", zorder=z))


# ---------------------------------------------------------------------------
# Proteins scattered in 3D
# ---------------------------------------------------------------------------
# Each protein gets a depth d in [0, 1] (0 = far back, 1 = front): far proteins
# are drawn smaller and fainter, and behind the DNA.

def depth_scale(d):
    return 0.6 + 0.5 * d


def depth_alpha(d):
    return 0.3 + 0.7 * d


def tf_alpha(d):
    """As depth_alpha, for TFs: they are small, so far ones fade less."""
    return 0.65 + 0.35 * d


def scatter_proteins(n, region, rng, avoid=(), existing=(), gap=0.1, tries=20000,
                     radius=TF_RADIUS):
    """`n` proteins at random positions, depths and orientations.

    `region` is (x0, x1, y0, y1); `avoid` is a list of (points, radius) pairs
    that proteins must stay clear of (DNA, polymerases, transcripts...).
    `existing` proteins (an earlier scatter) are kept, and new ones added to them.
    `radius` is the size of one protein at depth scale 1.
    Returns a list of dicts (center, depth, angle), ordered from back to front.
    """
    x0, x1, y0, y1 = region
    placed = list(existing)
    for _ in range(tries):
        if len(placed) == n:
            break
        d = rng.uniform()
        r = radius * depth_scale(d)
        c = np.array([rng.uniform(x0 + r, x1 - r), rng.uniform(y0 + r, y1 - r)])
        if all(np.linalg.norm(np.asarray(pts).reshape(-1, 2) - c, axis=1).min() > rad + r + gap
               for pts, rad in avoid) and \
           all(np.linalg.norm(q["center"] - c) > r + radius * depth_scale(q["depth"]) + gap
               for q in placed):
            placed.append(dict(center=c, depth=d, angle=rng.uniform(0, 360)))
    return sorted(placed, key=lambda q: q["depth"])


def draw_scattered_tfs(ax, tfs, z_back=0, z_front=10):
    """Draw scatter_proteins() output as TFs: the far half behind z_back, the near
    half at z_front."""
    for q in tfs:
        draw_tf(ax, q["center"], scale=depth_scale(q["depth"]), angle=q["angle"],
                alpha=tf_alpha(q["depth"]), z=z_back if q["depth"] < 0.5 else z_front)


# ---------------------------------------------------------------------------
# Legend
# ---------------------------------------------------------------------------

# The legend of scene 11, whose rows are these in this order.
LEGEND_NAMES = {"DNA": "DNA", "Promoter": "Promoter", "Gene": "Gene",
                "TF": "Transcription factor", "RNAP": "RNA Polymerase", "mRNA": "mRNA"}
LEGEND_HOME = (1.28, 4.44)  # the middle of the first row's icon, in most scenes
LEGEND_ROW = 0.6            # distance between rows
LEGEND_TEXT_X = 0.61        # from the middle of the icons to the labels
LEGEND_ICON = 0.34          # half-width of the icons of DNA, promoter, gene and mRNA
LEGEND_HELIX = 0.85         # the helix of those three is drawn this many times its size
LEGEND_BEADS = dict(n=9, radius=0.028)
LEGEND_PT = 8.5
LEGEND_WIDTH = 5.0          # room taken, from the left of the icons to the labels' end


def draw_legend(ax, items, at=LEGEND_HOME, mrna=BLUE["fill"], z=20):
    """A small legend of `items` (keys of LEGEND_NAMES, in the order wanted), going
    down from the row at `at`. Rows don't move when later ones come and go.
    `mrna` is the colour of its transcript."""
    strands = {"DNA": DNA_LINE, "Promoter": ORANGE, "Gene": BLUE}
    for row, item in enumerate(items):
        x, y = at[0], at[1] - row * LEGEND_ROW
        ends = np.array([(x - LEGEND_ICON, y), (x + LEGEND_ICON, y)])
        if item in strands:
            center = np.array((x, y))
            axis = center + (np.linspace(*ends, 60) - center) / LEGEND_HELIX
            draw_ribbon(ax, center + LEGEND_HELIX * (helix(axis) - center), strands[item],
                        tuple(LEGEND_HELIX * width for width in DNA_LINE_PT), z)
        elif item == "TF":
            draw_tf(ax, (x, y), scale=1.1, z=z)
        elif item == "RNAP":
            draw_rnap_at(ax, (x, y), scale=0.72, z=z)
        elif item == "mRNA":
            for bead in np.linspace(*ends, LEGEND_BEADS["n"]):
                ax.add_patch(Circle(bead, LEGEND_BEADS["radius"], fc=mrna, ec="none", zorder=z))
        draw_label(ax, (x + LEGEND_TEXT_X, y), LEGEND_NAMES[item], ha="left", size=LEGEND_PT, z=z)


def legend_area(items, at=LEGEND_HOME):
    """What scattered molecules must keep clear of around a legend, as a
    (points, radius) pair."""
    x, y = at
    xs = x - LEGEND_ICON + np.linspace(0, LEGEND_WIDTH, 14)
    ys = y - LEGEND_ROW * np.arange(len(items))
    return np.array([(a, b) for a in xs for b in ys]), 0.32
