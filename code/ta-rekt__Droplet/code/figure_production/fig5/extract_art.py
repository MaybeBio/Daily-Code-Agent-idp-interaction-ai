#!/usr/bin/env python3
"""Extract the hand-drawn artwork of the old Fig. 4 (now Fig. 5) into standalone SVGs.

One-off provenance step. The source is Tarek's Inkscape canvas, which is not in
this repository; the art_*.svg it writes are what fig5.ipynb uses.

    python figure_production/fig5/extract_art.py /path/to/figs/fig4.svg

What is kept is only the drawing: droplets, chain, arrows, sketch curves. Every
label and every plot is dropped, because the notebook re-sets all text in
Computer Modern at the caption size and regenerates all plots from code.

Elements are selected by id, not by position: the canvas also holds scraps of
other figures, and an SVG parser that splits <text>/<tspan> numbers children
differently from the XML.

Two changes to the drawing itself. The d arrow of panel (a) runs edge to edge
in the source; d is the distance between the centres (the paper's definition),
so it is redrawn centre to centre (see D_ARROW). And in this SVG the capillary bridge in panel (a)
is three flat purple shapes, whereas the published PDF shows a soft glow like
the halos around the other droplets. Those three shapes are rendered, Gaussian
blurred and embedded as a raster under the vector circles -- a glow is soft by
nature, so nothing visible is lost, and it renders identically under any
converter (cairosvg ignores SVG blur filters, Inkscape does not).
"""
import base64, copy, io, os, re, sys
import xml.etree.ElementTree as ET

import numpy as np
import cairosvg
from PIL import Image
from scipy.ndimage import gaussian_filter

SVG = "http://www.w3.org/2000/svg"
NS = "{%s}" % SVG
XLINK = "{http://www.w3.org/1999/xlink}"
for p, u in (("", SVG), ("svg", SVG), ("xlink", XLINK[1:-1]),
             ("inkscape", "http://www.inkscape.org/namespaces/inkscape"),
             ("sodipodi", "http://sodipodi.sourceforge.net/DTD/sodipodi-0.dtd")):
    ET.register_namespace(p, u)

HERE = os.path.dirname(os.path.abspath(__file__))

# ---- the artwork, by element id (direct children of layer1) -----------------
INSET_STROKE_SCALE = 0.5 / 0.85     # see main(): inset curves to 0.5 pt as printed

PIECES = {
    # (a) two droplets with halos, probe dot, and the double arrow between them
    "a_pair":    ["g11-9", "g40", "path4-61", "path4-61-7", "path4-61-77",
                  "path4-61-7-6", "path89", "path3510-8"],
    # (a) the bridged pair, its probe dot, and the leader to the "Capillary bridge" label
    "a_bridge":  ["g2467", "path89-3", "path3510-8-7-1"],
    # (b) the Rouse chain: halos, springs, beads, enhancer, anchor + probe,
    #     the ellipsis, the s and a arrows and the ticks on s
    "b_chain":   ["path4-61-77-0", "path4-61-77-0-2"]
                 + ["path3812-0", "path3814-9", "path3812-5-6", "path3814-6-0",
                    "path3812-7", "path3814-92", "path3812-5-2", "path3814-6-3",
                    "path4", "path3812-0-5", "path3814-9-9", "path3812-0-4",
                    "path3814-9-5", "path3812-0-4-8", "path3814-9-5-6",
                    "path3812-5-6-3", "path3814-6-0-6", "path3812-7-0",
                    "path3814-92-6", "path3812-5-2-2", "path3814-6-3-0",
                    "path3812-0-5-1", "path3814-9-9-5"]
                 + ["surface1", "path3510", "path3510-0-5-6", "path4-4",
                    "path89-3-4", "path17", "path17-1", "path17-36", "path17-3"],
    # (e) the polymer around the anchor, with the d^o double arrow
    "e_polymer": ["g79", "path3510-0-5-6-4"],
}
# (c) the three sketch curves of the "U_tot = Rouse + surface attraction" inset.
# They sit inside group g6 next to the glyphs for "=", "+", "Rouse", "Surface
# attraction", which are dropped: only these children of g6 are kept.
INSET_PARENT, INSET_KEEP = "g6", ["path47-3-7", "path47-3-7-9", "path47-3-7-0"]

# the bridge shapes that get the glow, inside g2467
GLOW_PARENT, GLOW_IDS, GLOW_SIGMA_MM = "g2467", ["path4-61-6-9", "path4-61-6-9-4", "path112"], 0.55

# panel (a)'s d arrow, redrawn between the centres of the two droplets
D_ARROW, D_ENDS = "path3510-8", ("path4-61-7-6", "path4-61-7")

PX_PER_MM = 24          # ink-bbox resolution
MARGIN_MM = 0.3


def url_refs(el):
    out = set()
    for e in el.iter():
        for k, v in e.attrib.items():
            out |= set(re.findall(r"url\(#([^)]+)\)", v))
            if k == XLINK + "href" and v.startswith("#"):
                out.add(v[1:])
    return out


def pruned_defs(defs, seeds):
    """The defs the selected elements need, followed transitively."""
    byid = {d.get("id"): d for d in defs}
    need, todo = set(), list(seeds)
    while todo:
        i = todo.pop()
        if i in need or i not in byid:
            continue
        need.add(i)
        todo.extend(url_refs(byid[i]))
    new = ET.Element(NS + "defs")
    for d in defs:
        if d.get("id") in need:
            new.append(copy.deepcopy(d))
    return new


def document(root, defs, elements):
    out = ET.Element(root.tag, {k: v for k, v in root.attrib.items()
                                if k in ("width", "height", "viewBox", "version")})
    seeds = set()
    for e in elements:
        seeds |= url_refs(e)
    out.append(pruned_defs(defs, seeds))
    g = ET.SubElement(out, NS + "g", {"id": "art"})
    for e in elements:
        g.append(copy.deepcopy(e))
    return out


def ink_bbox_mm(doc):
    """Ink extent in page millimetres, from a render (alpha > 2/255)."""
    w_mm = float(doc.get("viewBox").split()[2])
    data = ET.tostring(doc, encoding="unicode").encode()
    png = cairosvg.svg2png(bytestring=data, output_width=int(w_mm * PX_PER_MM))
    a = np.asarray(Image.open(io.BytesIO(png)).convert("RGBA"))[..., 3]
    ys, xs = np.where(a > 2)
    vb = [float(v) for v in doc.get("viewBox").split()]
    k = vb[2] / a.shape[1]
    return (vb[0] + xs.min() * k, vb[1] + ys.min() * k,
            vb[0] + (xs.max() + 1) * k, vb[1] + (ys.max() + 1) * k)


def crop(doc, bb, margin=MARGIN_MM):
    x0, y0, x1, y1 = bb[0] - margin, bb[1] - margin, bb[2] + margin, bb[3] + margin
    doc.set("viewBox", "%.4f %.4f %.4f %.4f" % (x0, y0, x1 - x0, y1 - y0))
    doc.set("width", "%.4fmm" % (x1 - x0))
    doc.set("height", "%.4fmm" % (y1 - y0))
    return doc


def add_glow(doc, root, defs, bb):
    """Replace the flat bridge shapes by a blurred raster of themselves."""
    parent = next(e for e in doc.iter() if e.get("id") == GLOW_PARENT)
    shapes = [c for c in list(parent) if c.get("id") in GLOW_IDS]
    assert len(shapes) == len(GLOW_IDS), "bridge shapes not found"
    # render just the shapes, in page coordinates, over a region padded by 4 sigma
    pad = 4 * GLOW_SIGMA_MM
    x0, y0, x1, y1 = bb[0] - pad, bb[1] - pad, bb[2] + pad, bb[3] + pad
    tmp = document(root, defs, [])
    holder = ET.SubElement(tmp.find(NS + "g"), NS + "g",
                           {k: v for k, v in parent.attrib.items() if k == "transform"})
    for s in shapes:
        holder.append(copy.deepcopy(s))
    crop(tmp, (x0, y0, x1, y1), margin=0)
    px = 40                                            # ~1000 dpi
    png = cairosvg.svg2png(bytestring=ET.tostring(tmp, encoding="unicode").encode(),
                           output_width=int((x1 - x0) * px))
    rgba = np.asarray(Image.open(io.BytesIO(png)).convert("RGBA")).astype(float)
    a = rgba[..., 3] / 255.0
    pre = rgba[..., :3] * a[..., None]                 # premultiplied, so edges stay purple
    s = GLOW_SIGMA_MM * px
    a_b = gaussian_filter(a, s)
    pre_b = np.stack([gaussian_filter(pre[..., c], s) for c in range(3)], -1)
    rgb = np.where(a_b[..., None] > 1e-6, pre_b / np.maximum(a_b[..., None], 1e-6), 0)
    out = np.dstack([np.clip(rgb, 0, 255), np.clip(a_b * 255, 0, 255)]).astype(np.uint8)
    buf = io.BytesIO(); Image.fromarray(out, "RGBA").save(buf, "PNG", optimize=True)
    img = ET.Element(NS + "image", {
        "x": "%.4f" % x0, "y": "%.4f" % y0,
        "width": "%.4f" % (x1 - x0), "height": "%.4f" % (y1 - y0),
        "preserveAspectRatio": "none",
        XLINK + "href": "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode()})
    for s_ in shapes:
        parent.remove(s_)
    # the raster sits at page coordinates, under everything else in the piece
    doc.find(NS + "g").insert(0, img)
    return (x0, y0, x1, y1)


def centre_to_centre(doc):
    """Point the d arrow from the centre of one droplet to the other's."""
    byid = {e.get("id"): e for e in doc.iter()}
    cx = "{http://sodipodi.sourceforge.net/DTD/sodipodi-0.dtd}cx"
    cy = "{http://sodipodi.sourceforge.net/DTD/sodipodi-0.dtd}cy"
    (x0, y0), (x1, y1) = [(float(byid[i].get(cx)), float(byid[i].get(cy))) for i in D_ENDS]
    byid[D_ARROW].set("d", "M %.6f,%.6f %.6f,%.6f" % (x0, y0, x1, y1))


def main(src):
    root = ET.parse(src).getroot()
    defs = root.find(NS + "defs")
    layer = root.find(".//%sg[@id='layer1']" % NS)
    assert not layer.get("transform"), "layer1 has a transform -- handle it"
    kids = {c.get("id"): c for c in layer}

    report = []
    for name, ids in PIECES.items():
        missing = [i for i in ids if i not in kids]
        assert not missing, "%s: ids not found: %s" % (name, missing)
        order = [c for c in layer if c.get("id") in set(ids)]     # document order
        doc = document(root, defs, order)
        if name == "a_pair":
            centre_to_centre(doc)
        bb = ink_bbox_mm(doc)
        if name == "a_bridge":
            gb = add_glow(doc, root, defs, bb)
            bb = (min(bb[0], gb[0] + 2 * GLOW_SIGMA_MM), min(bb[1], gb[1] + 2 * GLOW_SIGMA_MM),
                  max(bb[2], gb[2] - 2 * GLOW_SIGMA_MM), max(bb[3], gb[3] - 2 * GLOW_SIGMA_MM))
            bb = ink_bbox_mm(doc)
        report.append((name, crop(doc, bb), bb))

    # (c) inset: g6 with only the three curves kept
    g6 = copy.deepcopy(kids[INSET_PARENT])
    for c in list(g6):
        if c.get("id") not in INSET_KEEP:
            g6.remove(c)
    assert len(list(g6)) == 3, "inset curves not found"
    doc = document(root, defs, [g6])
    bb = ink_bbox_mm(doc)          # crop on the original ink, so the anchor stays put
    # The notebook draws this group 0.92 x 1.45 = 1.33x larger than drawn, which
    # would print the curves at 0.85 pt; thin them so they print at the house 0.5 pt.
    for el in doc.iter():
        st = el.get("style")
        if st and "stroke-width" in st:
            el.set("style", re.sub(r"stroke-width:([0-9.eE+-]+)",
                                   lambda m: "stroke-width:%.6g" % (float(m.group(1)) * INSET_STROKE_SCALE), st))
    report.append(("c_inset", crop(doc, bb), bb))

    for name, doc, bb in report:
        path = os.path.join(HERE, "art_%s.svg" % name)
        ET.ElementTree(doc).write(path, xml_declaration=True, encoding="utf-8")
        print("art_%-10s %6.2f x %5.2f mm  at page (%6.2f, %6.2f) mm   %6.1f KB"
              % (name + ".svg", bb[2] - bb[0], bb[3] - bb[1], bb[0], bb[1],
                 os.path.getsize(path) / 1024))


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else
         os.path.expanduser("~/NextCloudCNRS/Code/droplet/Tarek_Droplets/figs/fig4.svg"))
