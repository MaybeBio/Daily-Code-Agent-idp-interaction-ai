"""
Slides that reuse storyboard artwork as drawn (scenes 7, 8, ...).

The drawing comes straight from storyboard/scene<N>.svg; only its text is
redone, in the slides' font like every other label:

  1. the scene's SVG, minus its text, is placed on the 128 x 96 mm slide with
     the same mapping as the traced scenes (sketch frame centre -> slide
     centre, scaled by SKETCH_SCALE) and exported to PDF by Inkscape;
  2. its labels are drawn by matplotlib at the same places, on a transparent
     slide-sized page;
  3. qpdf overlays the labels on the artwork.

Optionally (`raster_dpi`), the line art is flattened into one JPEG and only the
gradient fills stay vector, laid on top in multiply mode: over white they look
as usual and black lines stay black, so the lines show "through" them as if
drawn on top. Hundreds of hand-drawn paths otherwise make PDF viewers lag.

Requires the `inkscape` and `qpdf` command-line tools.
"""

import copy
import hashlib
import re
import subprocess
import tempfile
import xml.etree.ElementTree as ET
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from PIL import Image

from .drawing import LABEL_PT, SLIDE_MM, SLIDE_UNITS, UNIT_MM, draw_label, new_slide

STORYBOARD_DIR = Path(__file__).resolve().parents[2] / "storyboard"

# Slide mm per sketch mm: the scene-1 sketch frame (16.2138 mm) spans 20 units.
SKETCH_SCALE = 20 * UNIT_MM / 16.2138
# Storyboard captions are 0.3528 mm tall; they become LABEL_PT labels, and
# other text sizes scale in proportion.
SKETCH_LABEL_SIZE = 0.352777

NS = {
    "": "http://www.w3.org/2000/svg",
    "inkscape": "http://www.inkscape.org/namespaces/inkscape",
    "sodipodi": "http://sodipodi.sourceforge.net/DTD/sodipodi-0.dtd",
    "xlink": "http://www.w3.org/1999/xlink",
}
for prefix, uri in NS.items():
    ET.register_namespace(prefix, uri)
SVG = "{%s}" % NS[""]


def style_value(style, key):
    m = re.search(rf"(?:^|;){key}:([^;]+)", style or "")
    return m.group(1) if m else None


class StoryboardScene:
    """A storyboard scene's artwork and labels, placed on a slide.

    `narration` is the id of the scene's narration text, which is left out;
    `shift` moves artwork and labels together (slide units), e.g. to make room
    for labels that come out wider than in the sketch.
    """

    def __init__(self, number, narration, shift=(0.0, 0.0)):
        self.tree = ET.parse(STORYBOARD_DIR / f"scene{number}.svg")
        root = self.tree.getroot()
        x, y, w, h = map(float, root.get("viewBox").split())
        self.center = (x + w / 2, y + h / 2)
        self.shift = shift
        self.layer = root.find(SVG + "g")
        texts = [t for t in self.layer if t.tag == SVG + "text"]
        self.labels = [label for t in texts if t.get("id") != narration
                       for label in self.text_lines(t)]
        for t in texts:
            self.layer.remove(t)

    def to_slide(self, x, y):
        """Sketch coordinates (mm, y down) -> slide units (y up)."""
        k = SKETCH_SCALE / UNIT_MM
        return (SLIDE_UNITS[0] / 2 + self.shift[0] + (x - self.center[0]) * k,
                SLIDE_UNITS[1] / 2 + self.shift[1] - (y - self.center[1]) * k)

    def text_lines(self, text):
        """One label per line (tspan) of a storyboard text element."""
        style = text.get("style")
        size = float(style_value(style, "font-size").rstrip("px"))
        # Text drawn as outlines (no fill, black stroke) reads as bold.
        outlined = style_value(style, "fill") == "none"
        anchor = {"start": "left", "middle": "center", "end": "right"}[
            style_value(style, "text-anchor") or "start"]
        size_pt = LABEL_PT * size / SKETCH_LABEL_SIZE
        # Lines keep the sketch's first baseline, but are spaced for the label
        # size (the sketch's text was smaller): 1.2 x font size, in slide units.
        line_step = 1.2 * size_pt * 25.4 / 72 / UNIT_MM
        lines = []
        for tspan in text.iter(SVG + "tspan"):
            content = "".join(tspan.itertext()).strip()
            if content:
                bold = outlined or style_value(tspan.get("style"), "font-weight") == "bold"
                x, y = self.to_slide(float(tspan.get("x", text.get("x"))),
                                     float(tspan.get("y", text.get("y"))))
                if lines:
                    y = lines[0]["xy"][1] - len(lines) * line_step
                lines.append(dict(text=content, bold=bold, ha=anchor, xy=(x, y), size=size_pt))
        return lines

    def placed_svg(self, edit=None):
        """The artwork on a slide-sized page; `edit(root)` may modify a copy first."""
        root = copy.deepcopy(self.tree.getroot())
        if edit:
            edit(root)
        layer = root.find(SVG + "g")
        tx = SLIDE_MM[0] / 2 + self.shift[0] * UNIT_MM
        ty = SLIDE_MM[1] / 2 - self.shift[1] * UNIT_MM
        layer.set("transform", f"translate({tx},{ty}) scale({SKETCH_SCALE}) "
                               f"translate({-self.center[0]},{-self.center[1]})")
        root.set("width", f"{SLIDE_MM[0]}mm")
        root.set("height", f"{SLIDE_MM[1]}mm")
        root.set("viewBox", f"0 0 {SLIDE_MM[0]} {SLIDE_MM[1]}")
        return ET.ElementTree(root)

    def export(self, out_pdf, edits=(None,), labels=True, raster_dpi=None):
        """Write the slide (artwork + labels) to `out_pdf`, one page per edit.

        Several edits make a multi-page PDF, e.g. the frames of an animation.
        With `raster_dpi`, the line art is flattened to a JPEG at that resolution
        (see the module docstring).
        """
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            svgs = []
            for i, edit in enumerate(edits):
                tree = self.placed_svg(edit)
                if raster_dpi:
                    tree = flattened(tree, tmp, raster_dpi)
                svgs.append(tmp / f"art{i:03d}.svg")
                tree.write(svgs[-1], encoding="UTF-8", xml_declaration=True)
            subprocess.run(["inkscape", "--export-type=pdf", "--export-area-page",
                            *map(str, svgs)], check=True, capture_output=True)
            art = tmp / "art.pdf"
            subprocess.run(["qpdf", "--empty", "--pages", *[str(s.with_suffix(".pdf"))
                            for s in svgs], "--", str(art)], check=True)

            fig, ax = new_slide()
            for label in self.labels if labels else []:
                draw_label(ax, label["xy"], label["text"], bold=label["bold"],
                           ha=label["ha"], va="baseline", size=label["size"])
            label_pdf = tmp / "labels.pdf"
            fig.savefig(label_pdf, transparent=True)
            plt.close(fig)

            Path(out_pdf).parent.mkdir(parents=True, exist_ok=True)
            subprocess.run(["qpdf", str(art), "--overlay", str(label_pdf), "--from=",
                            "--repeat=1", "--", str(out_pdf)], check=True)
        print(f"wrote {out_pdf} ({len(edits)} page{'s' * (len(edits) > 1)})")


# ---------------------------------------------------------------------------
# Flattening heavy line art to a JPEG
# ---------------------------------------------------------------------------

JPEG_QUALITY = 90
JPEG_PAD_PX = 12


def is_gradient_filled(e):
    return "fill:url(#" in (e.get("style") or "")


def line_art_jpeg(root, tmp, dpi):
    """Render the slide-sized SVG `root` on white, crop it to its drawing, and
    save it as a greyscale JPEG. Returns (path, (x, y, width, height) in mm).

    Identical line art (e.g. the same drawing under different overlays) is
    rendered once per export.
    """
    text = ET.tostring(root)
    name = hashlib.sha1(text).hexdigest()[:12]
    jpg = tmp / f"lines_{name}.jpg"
    box_file = jpg.with_suffix(".box")
    if not jpg.exists():
        svg, png = tmp / f"lines_{name}.svg", tmp / f"lines_{name}.png"
        svg.write_bytes(text)
        subprocess.run(["inkscape", str(svg), "--export-type=png", f"--export-dpi={dpi}",
                        "--export-area-page", "--export-background=white",
                        "--export-background-opacity=1", f"--export-filename={png}"],
                       check=True, capture_output=True)
        image = Image.open(png).convert("L")
        ys, xs = np.nonzero(np.asarray(image) < 250)
        left, top = max(xs.min() - JPEG_PAD_PX, 0), max(ys.min() - JPEG_PAD_PX, 0)
        right = min(xs.max() + JPEG_PAD_PX, image.width)
        bottom = min(ys.max() + JPEG_PAD_PX, image.height)
        image.crop((left, top, right, bottom)).save(jpg, quality=JPEG_QUALITY)
        mm = 25.4 / dpi
        box_file.write_text(" ".join(str(v * mm) for v in
                                     (left, top, right - left, bottom - top)))
    return jpg, tuple(map(float, box_file.read_text().split()))


def flattened(tree, tmp, dpi):
    """The placed artwork `tree` with its line art as one JPEG and its gradient
    fills (top-level shapes) kept as vector overlays in multiply mode."""
    root = tree.getroot()
    lines = copy.deepcopy(root)
    for e in lines.iter():
        if is_gradient_filled(e):
            e.set("style", re.sub(r"fill:url\(#[^)]+\)", "fill:none", e.get("style")))
    jpg, (x, y, w, h) = line_art_jpeg(lines, tmp, dpi)

    layer = root.find(SVG + "g")
    for e in list(layer):
        if not is_gradient_filled(e):
            layer.remove(e)
        else:
            # No stroke; and no object opacity (the storyboard has 0.994 on
            # everything): with it, Inkscape exports the shape as an isolated
            # group, inside which the multiply blend would not see the lines.
            style = re.sub(r"(^|;)(stroke|opacity):[^;]*", "", e.get("style"))
            e.set("style", f"{style};stroke:none;mix-blend-mode:multiply")
    # The image must share the overlays' group (a blend only sees its own group's
    # backdrop), so it goes first in the layer, undoing the layer's placement.
    image = ET.Element(SVG + "image", {
        "x": str(x), "y": str(y), "width": str(w), "height": str(h),
        "preserveAspectRatio": "none", "{%s}href" % NS["xlink"]: str(jpg),
        "transform": inverse_transform(layer.get("transform"))})
    layer.insert(0, image)
    return tree


def inverse_transform(transform):
    """Inverse of a "translate(a,b) scale(s) translate(c,d)" placement."""
    (a, b), (s,), (c, d) = [tuple(map(float, re.split(r"[ ,]+", args.strip())))
                            for args in re.findall(r"\(([^)]*)\)", transform)]
    return f"translate({-c},{-d}) scale({1 / s}) translate({-a},{-b})"
