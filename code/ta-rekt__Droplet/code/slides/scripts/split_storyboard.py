"""
Split the master storyboard into one SVG per scene.

The storyboard (storyboard/defense-storyboard.svg) is a single Inkscape layer laid
out on a grid: columns 1..11, rows a..k, read left-to-right then top-to-bottom
("slide order" arrow). Each non-empty grid cell is a scene. This script writes
storyboard/scene<N>.svg for every non-empty cell, numbered in reading order.

Each scene file keeps the original elements verbatim (same coordinates, styles,
ids); only the page is changed so that it frames the cell, and only the <defs>
entries that the scene actually references (gradients, path effects) are copied.

Usage (from slides/):
    python3 scripts/split_storyboard.py
Requires the `inkscape` CLI (used to get the bounding box of every element).
"""

import csv
import io
import re
import subprocess
import xml.etree.ElementTree as ET
from pathlib import Path

SLIDES_DIR = Path(__file__).resolve().parent.parent
MASTER = SLIDES_DIR / "storyboard" / "defense-storyboard.svg"
OUT_DIR = SLIDES_DIR / "storyboard"

PX_PER_UNIT = 96 / 25.4  # inkscape --query-all reports px; the document unit is mm

NS = {
    "": "http://www.w3.org/2000/svg",
    "svg": "http://www.w3.org/2000/svg",
    "inkscape": "http://www.inkscape.org/namespaces/inkscape",
    "sodipodi": "http://sodipodi.sourceforge.net/DTD/sodipodi-0.dtd",
    "xlink": "http://www.w3.org/1999/xlink",
}
for prefix, uri in NS.items():
    if prefix != "svg":
        ET.register_namespace(prefix, uri)
SVG = "{%s}" % NS[""]
INK = "{%s}" % NS["inkscape"]


# ---------------------------------------------------------------------------
# Geometry: bounding boxes and the grid
# ---------------------------------------------------------------------------

def query_bboxes(svg_path):
    """Return {id: (x, y, w, h)} in document units (mm) for every element."""
    out = subprocess.run(
        ["inkscape", "--query-all", str(svg_path)],
        capture_output=True, text=True, check=True,
    ).stdout
    boxes = {}
    for row in csv.reader(io.StringIO(out)):
        if len(row) == 5:
            boxes[row[0]] = tuple(float(v) / PX_PER_UNIT for v in row[1:])
    return boxes


def find_grid(layer, boxes):
    """The grid is drawn with long straight lines spanning the whole board.

    Returns (grid line ids, sorted x positions of verticals, sorted y of horizontals).
    """
    xs, ys, ids = [], [], set()
    full_w = max(boxes[c.get("id")][2] for c in layer if c.get("id") in boxes)
    full_h = max(boxes[c.get("id")][3] for c in layer if c.get("id") in boxes)
    for child in layer:
        x, y, w, h = boxes.get(child.get("id"), (0, 0, 0, 0))
        if child.tag == SVG + "path" and h > 0.9 * full_h and w < 1:
            xs.append(x + w / 2)
            ids.add(child.get("id"))
        elif child.tag == SVG + "path" and w > 0.9 * full_w and h < 1:
            ys.append(y + h / 2)
            ids.add(child.get("id"))
    return ids, sorted(xs), sorted(ys)


def anchor_point(elem, boxes):
    """Point used to decide which cell an element belongs to.

    Normally the centre of its bounding box; empty <text> stubs have no box, so
    fall back to their x/y attribute (on the element or its first tspan).
    """
    if elem.get("id") in boxes:
        x, y, w, h = boxes[elem.get("id")]
        return x + w / 2, y + h / 2
    for e in [elem, *elem]:
        if e.get("x") is not None and e.get("y") is not None:
            return float(e.get("x")), float(e.get("y"))
    return None


def cell_index(value, edges):
    for i in range(len(edges) - 1):
        if edges[i] <= value < edges[i + 1]:
            return i
    return None


# ---------------------------------------------------------------------------
# References into <defs>
# ---------------------------------------------------------------------------

REF = re.compile(r"#([A-Za-z_][\w.\-]*)")


def referenced_ids(elem):
    """All '#id' targets mentioned in attributes of elem and its descendants."""
    found = set()
    for e in elem.iter():
        for key, value in e.attrib.items():
            if "url(#" in value or key.endswith("href") or key == INK + "path-effect":
                found.update(REF.findall(value))
    return found


def defs_closure(scene_elems, defs):
    """Ids of every <defs> child needed by the scene, following chains (href)."""
    by_id = {d.get("id"): d for d in defs}
    needed, todo = set(), set()
    for e in scene_elems:
        todo |= referenced_ids(e)
    while todo:
        ref = todo.pop()
        if ref in by_id and ref not in needed:
            needed.add(ref)
            todo |= referenced_ids(by_id[ref])
    return needed


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    tree = ET.parse(MASTER)
    root = tree.getroot()
    defs = root.find(SVG + "defs")
    namedview = root.find("{%s}namedview" % NS["sodipodi"])
    layer = root.find(SVG + "g")
    boxes = query_bboxes(MASTER)
    grid_ids, xs, ys = find_grid(layer, boxes)

    # Assign every non-grid element to a (row, col) cell, in document order.
    cells, skipped = {}, []
    for child in layer:
        if child.get("id") in grid_ids:
            continue
        p = anchor_point(child, boxes)
        col = cell_index(p[0], xs) if p else None
        row = cell_index(p[1], ys) if p else None
        if row is None or col is None:
            skipped.append(child.get("id"))  # row/column labels, "slide order" arrow
            continue
        cells.setdefault((row, col), []).append(child)

    # Write one file per non-empty cell, in reading order (row-major).
    for n, (row, col) in enumerate(sorted(cells), start=1):
        x0, x1, y0, y1 = xs[col], xs[col + 1], ys[row], ys[row + 1]
        elems = cells[(row, col)]
        needed = defs_closure(elems, defs)

        new_root = ET.Element(root.tag, {
            k: v for k, v in root.attrib.items()
            if k not in ("width", "height", "viewBox", "{%s}docname" % NS["sodipodi"])
        })
        new_root.set("width", f"{x1 - x0:.4f}mm")
        new_root.set("height", f"{y1 - y0:.4f}mm")
        new_root.set("viewBox", f"{x0:.4f} {y0:.4f} {x1 - x0:.4f} {y1 - y0:.4f}")
        new_root.set("{%s}docname" % NS["sodipodi"], f"scene{n}.svg")

        # Drop the saved view so Inkscape opens the scene fitted to its page.
        view_keys = {INK + k for k in ("zoom", "cx", "cy", "window-width", "window-height",
                                       "window-x", "window-y", "window-maximized")}
        nv = ET.SubElement(new_root, namedview.tag, {
            k: v for k, v in namedview.attrib.items() if k not in view_keys
        })
        new_defs = ET.SubElement(new_root, defs.tag, defs.attrib)
        new_defs.extend(d for d in defs if d.get("id") in needed)
        new_layer = ET.SubElement(new_root, layer.tag, layer.attrib)
        new_layer.extend(elems)
        # Copied elements keep their own whitespace; ET.indent() is avoided on
        # purpose because it would inject spaces into xml:space="preserve" text.
        new_root.text = nv.tail = new_defs.tail = "\n  "
        new_defs.text = new_layer.text = "\n    "
        new_layer.tail = "\n"

        out = OUT_DIR / f"scene{n}.svg"
        ET.ElementTree(new_root).write(out, encoding="UTF-8", xml_declaration=True)
        print(f"scene{n:>2}  cell {'abcdefghijk'[row]}{col + 1:<2}  "
              f"{len(elems):>4} elements  {len(needed):>4} defs  -> {out.name}")

    print(f"not in any cell (labels etc.): {skipped}")


if __name__ == "__main__":
    main()
