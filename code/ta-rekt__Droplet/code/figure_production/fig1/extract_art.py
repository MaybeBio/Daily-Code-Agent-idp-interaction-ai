#!/usr/bin/env python3
"""Extract the hand-drawn artwork of Fig. 1 (panels A and D) into standalone SVGs.

One-off provenance step, as for Fig. 5. The source is the paper's canvas
(Tarek_Droplets/fig_sources/fig1.svg, an Illustrator export later touched
in Inkscape), which is not in this repository; the art_*.svg it writes are what
fig1.ipynb uses.

    python figure_production/fig1/extract_art.py /path/to/fig_sources/fig1.svg

Only the drawing is kept -- beads, arrows, interaction icons, the two-state
cartoon, the polymers and their clouds, and the two 3D lattice plots B and C,
which are pasted as published (matplotlib's own rendering of chain, panes,
grid and depth-shaded ligands; re-plotting them depends on the matplotlib
version). All lettering (glyph outlines in the source) is dropped and re-set by
the notebook at the caption size; for B and C the position and text of every
tick label is written to art_bc_labels.json. E and F are regenerated.

Each piece is also converted to PDF here, once, with Inkscape (art_*.pdf), and
the notebook stamps those PDFs: SVG renderers disagree on this file's gradients
(cairosvg on another machine drew D's glows differently), so no conversion
happens at render time.

Elements are selected by id. The outputs keep the source canvas (viewBox
0 0 350.94 340.33, y down), so the notebook places them with one canvas -> page
map per piece.

One repair: several circles still carry sodipodi:type="arc" attributes from an
earlier Inkscape session, with centres that no longer match their path data
(Illustrator rewrote `d`). Inkscape regenerates such shapes from the sodipodi
attributes, which drew D's three clouds and its probe dots in the wrong place
(off canvas, and one stray dot beside A's "Off" state). All sodipodi:*
attributes are stripped, so every renderer uses `d`.
"""
import copy, csv, json, os, shutil, subprocess, sys, tempfile
import xml.etree.ElementTree as ET

SVG = "http://www.w3.org/2000/svg"
SOD = "{http://sodipodi.sourceforge.net/DTD/sodipodi-0.dtd}"
for p, u in (("", SVG), ("xlink", "http://www.w3.org/1999/xlink"),
             ("inkscape", "http://www.inkscape.org/namespaces/inkscape"),
             ("sodipodi", SOD[1:-1])):
    ET.register_namespace(p, u)

HERE = os.path.dirname(os.path.abspath(__file__))

# ---- panel A: the drawing, by id (direct children of svg1) --------------------
A_DRAWING = [
    # particle-interaction icons and the brace
    "g1684", "g1685", "g58", "path3446-6-8-7-6-4-8", "path3446-6-8-7-6-4-8-8", "g50",
    "surface1-1-1-4-7",
    # two-state probe: Off / On circles, their dots, the two arcs and heads, leader
    "path3060-9", "path3075-1", "path3506-2-6", "g31",
    "path3110-2", "path3111-7", "path3115-0", "path3116-9", "path3406-6-9",
    # the chain: beads, probe, ligand, ellipsis dots
    "path3446-41-9", "path3446-1-0", "path3446-7-5", "path3446-4-8",
    "path3460-2", "path3461-4", "path3465-8", "path3470-4", "path3471-8", "path3475-0",
    "path3466-1", "path3446-6", "path3450-4", "path3451-8", "path3495-7", "path3496-7",
    "path3500-2", "path3455-7", "path3476-6", "path3480-3", "path3481-6", "path3486-6",
    "path3490-5", "path3491-4", "path3485-7", "path3505-5", "path3506-2",
    "path3446-4-3-3", "path3446-7-3-0", "path3446-1-8-8", "path3446-41-6-0",
    "path3446-5-0-9", "path3446-18-4-9", "path1262-92-6-9",
    # the three length arrows under the chain
    "path3510-6", "path3511-9", "path3515-6", "path3516-7", "path3520-6", "path3521-4",
    "path3406-0", "path3410-9", "path3411-0",
]
# ---- panel D: group g33 without its three labels ------------------------------
D_GROUP = "g33"
D_TEXT = ["g1796-8", "surface1-1-1-1-1-1-1-1-1-1-1-1-1-1-1-1-1-1-1-1-1-1-1-1-0", "g1798-7"]
D_STRAY = []     # (path1786-4 looks stray only when drawn from its stale sodipodi
                 # attributes; from its path data it is the closed polymer's probe dot)
# The three small glows at the contact sites (radial gradients with a mirrored,
# partly rotated gradientTransform) land ~1 unit away from where the published
# Illustrator PDF draws them, in every SVG renderer tried (Inkscape, cairosvg).
# Size, shape and tilt agree; only the centres differ. Each glow element is
# shifted -- shape and gradient together -- by the offset measured between the
# PDF and a render of the extracted art (blue-excess centroids at 500 dpi).
GLOW_SHIFT = {"path1782-25": (-0.21, -1.36),       # closed polymer, at the contact
              "path1782-2-4": (0.62, -0.94),       # open polymer, at the probe
              "path1782-2-2-1": (-0.61, 0.51)}     # open polymer, at the enhancer

# D's clouds and glows are radial gradients that fade to transparent. Apple's PDF
# renderer (Preview, Quick Look) fills such soft masks opaquely when they come
# from Inkscape or cairo -- the clouds print as solid discs -- although poppler and
# Illustrator's own PDF are fine. So the transparency is removed: each gradient
# is pre-composited over what lies under it (white for a cloud, the flat cloud
# tint for a glow) into opaque stops, which every renderer draws alike.
CLOUD_TINT_ALPHA, CLOUD_COLOR = 0.23, "#9d9df3"
CLOUDS = ["radial-gradient", "radial-gradient2", "radial-gradient3"]
# glow gradient -> (does the cloud lie on top of it in the source?)
GLOWS = {"radial-gradient1": True, "radial-gradient4": False, "radial-gradient5": False}
GLOW_UNDER_CLOUD = ("path1782-25", "path1782-6-9")   # closed: glow drawn before its cloud
# The closed polymer's probe dot lies under that cloud and glow too, and shows
# through them; once they are opaque it is drawn on top, in its printed colour
# (sampled from the published PDF at its centre).
DOT_UNDER_CLOUD, DOT_PRINTED = "path1786-4", "#5555a4"


# ---- panels B, C: the matplotlib groups, without their tick-label glyphs ----------
# Each label is one unclassed, black-filled path. Their text, by position (canvas
# units, measured): z ticks bottom-up, y ticks front to back, x ticks left to
# right, then the axis letter.
BC = {
    "b": ("g1-0", {"path54-9": "15", "path56-6": "20", "path58-02": "25", "path60-47": "30",
                   "path62-16": "35", "path64-5": "40", "path65-7": "z",
                   "path40-0": "15", "path42-2": "20", "path44-93": "25", "path46-2": "30",
                   "path48-9": "35", "path50-2": "40", "path51-08": "y",
                   "path26-76": "15", "path28-10": "20", "path30-11": "25", "path32-64": "30",
                   "path34-4": "35", "path36-0": "40", "path37-7": "x"}),
    "c": ("g1-9", {"path54-3": "20", "path56-8": "25", "path58-6": "30", "path60-3": "35",
                   "path62-3": "40", "path64-6": "45", "path65-0": "z",
                   "path40-3": "20", "path42-59": "25", "path44-2": "30", "path46-85": "35",
                   "path48-4": "40", "path50-59": "45", "path51-7": "y",
                   "path26-6": "20", "path28-5": "25", "path30-8": "30", "path32-8": "35",
                   "path34-3": "40", "path36-44": "45", "path37-4": "x"}),
}
INKSCAPE = "/Applications/Inkscape.app/Contents/MacOS/inkscape"


def _rgb(h):
    h = h.lstrip("#")
    return [int(h[i:i + 2], 16) / 255 for i in (0, 2, 4)]


def _hex(c):
    return "#%02x%02x%02x" % tuple(int(round(255 * min(max(v, 0), 1))) for v in c)


def _stops(grad):
    """(offset, rgb, alpha) of a gradient's stops."""
    out = []
    for st in grad:
        out.append((float(st.get("offset")), _rgb(st.get("stop-color")),
                    float(st.get("stop-opacity", 1))))
    return out


def _sample(stops, t):
    """Colour and alpha at t, interpolated separately (as SVG does), padded."""
    if t <= stops[0][0]:
        return stops[0][1], stops[0][2]
    for (t0, c0, a0), (t1, c1, a1) in zip(stops, stops[1:]):
        if t <= t1:
            f = (t - t0) / (t1 - t0)
            return [x + f * (y - x) for x, y in zip(c0, c1)], a0 + f * (a1 - a0)
    return stops[-1][1], stops[-1][2]


def opaque_gradients(doc):
    """Pre-composite D's transparent gradients (see CLOUDS / GLOWS above)."""
    white, cloud = [1, 1, 1], _rgb(CLOUD_COLOR)
    tint = [CLOUD_TINT_ALPHA * c + (1 - CLOUD_TINT_ALPHA) * w for c, w in zip(cloud, white)]
    grads = {g.get("id"): g for g in doc.iter() if g.tag.endswith("radialGradient")}

    def rewrite(gid, under, over_alpha=0.0):
        g = grads[gid]
        st = _stops(g)
        for e in list(g):
            g.remove(e)
        for t in [i / 20 for i in range(21)]:
            c, a = _sample(st, t)
            col = [a * x + (1 - a) * u for x, u in zip(c, under)]
            if over_alpha:                       # the cloud's flat tint on top
                col = [over_alpha * x + (1 - over_alpha) * y for x, y in zip(cloud, col)]
            ET.SubElement(g, "{%s}stop" % SVG, {"offset": "%.3f" % t, "stop-color": _hex(col)})

    for gid in CLOUDS:
        rewrite(gid, white)
    for gid, under_cloud in GLOWS.items():
        # a glow under the cloud: composite over white, then the tint on top; it is
        # then drawn after the cloud, opaque, so the result is the same
        rewrite(gid, white if under_cloud else tint, CLOUD_TINT_ALPHA if under_cloud else 0.0)
    glow, cloud_el = [next(e for e in doc.iter() if e.get("id") == i) for i in GLOW_UNDER_CLOUD]
    parent = {c: p for p in doc.iter() for c in p}
    par = parent[glow]
    assert parent[cloud_el] is par
    par.remove(glow)
    par.insert(list(par).index(cloud_el) + 1, glow)
    dot = next(e for e in doc.iter() if e.get("id") == DOT_UNDER_CLOUD)
    parent[dot].remove(dot)
    dot.set("style", "fill:%s" % DOT_PRINTED)
    par.insert(list(par).index(glow) + 1, dot)


def strip_sodipodi(root):
    for e in root.iter():
        for k in [k for k in e.attrib if k.startswith(SOD)]:
            del e.attrib[k]


def document(root, elements):
    """A copy of the source canvas holding only `elements` (plus all defs/styles)."""
    out = ET.Element(root.tag, {k: v for k, v in root.attrib.items() if k in ("viewBox", "version")})
    for c in root:
        if c.tag.split("}")[1] == "defs":
            out.append(copy.deepcopy(c))
    g = ET.SubElement(out, "{%s}g" % SVG, {"id": "art"})
    for e in elements:
        g.append(copy.deepcopy(e))
    return out


def main(src):
    root = ET.parse(src).getroot()
    strip_sodipodi(root)
    svg1 = [c for c in root if c.get("id") == "svg1"][0]
    kids = {c.get("id"): c for c in svg1}
    missing = [i for i in A_DRAWING + [D_GROUP] if i not in kids]
    assert not missing, "ids not found: %s" % missing

    a = document(root, [c for c in svg1 if c.get("id") in set(A_DRAWING)])   # document order
    d_el = copy.deepcopy(kids[D_GROUP])
    parent = {c: p for p in d_el.iter() for c in p}
    for e in list(d_el.iter()):
        if e.get("id") in set(D_TEXT + D_STRAY):
            parent[e].remove(e)
        elif e.get("id") in GLOW_SHIFT:
            assert not e.get("transform"), e.get("id")
            e.set("transform", "translate(%g %g)" % GLOW_SHIFT[e.get("id")])
    d = document(root, [d_el])
    opaque_gradients(d)

    pieces = {"a": a, "d": d}
    glyph_ids = {}
    for name, (gid, labels) in BC.items():
        g = copy.deepcopy(kids[gid])
        parent = {c: p for p in g.iter() for c in p}
        glyphs = [e for e in g.iter() if e.tag.endswith("path") and not e.get("class")]
        assert sorted(e.get("id") for e in glyphs) == sorted(labels), name
        for e in glyphs:
            parent[e].remove(e)
        pieces[name] = document(root, [g])
        glyph_ids[name] = labels

    for name, doc in pieces.items():
        path = os.path.join(HERE, "art_%s.svg" % name)
        ET.ElementTree(doc).write(path, xml_declaration=True, encoding="utf-8")
        subprocess.run([INKSCAPE, path, "--export-type=pdf", "--export-area-page",
                        "--export-filename=" + path[:-4] + ".pdf"], check=True, capture_output=True)
        print("wrote", os.path.relpath(path), "and .pdf (%d elements)" % sum(1 for _ in doc.iter()))

    # label positions: glyph bounding-box centres, from Inkscape on the stripped source
    tmp = tempfile.mkdtemp()
    src_clean = os.path.join(tmp, "src.svg")
    ET.ElementTree(root).write(src_clean, xml_declaration=True, encoding="utf-8")
    q = subprocess.run([INKSCAPE, src_clean, "--query-all"], check=True, capture_output=True, text=True).stdout
    bbox = {row[0]: [float(v) for v in row[1:5]] for row in csv.reader(q.splitlines()) if len(row) == 5}
    shutil.rmtree(tmp)
    out = {}
    for name, labels in glyph_ids.items():
        out[name] = [{"text": t, "u": round(bbox[i][0] + bbox[i][2] / 2, 2),
                      "v": round(bbox[i][1] + bbox[i][3] / 2, 2)} for i, t in labels.items()]
    json.dump(out, open(os.path.join(HERE, "art_bc_labels.json"), "w"), indent=1)
    print("wrote art_bc_labels.json")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else
         os.path.join(HERE, "../../../Tarek_Droplets/fig_sources/fig1.svg"))
