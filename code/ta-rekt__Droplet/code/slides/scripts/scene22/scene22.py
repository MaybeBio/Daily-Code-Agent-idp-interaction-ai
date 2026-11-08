"""
Scene 22 -- "We wondered what a gene regulation mechanism that accounts for
this biological reality would look like. I started by building intuition with
lattice simulations of polymers in 3D space interacting with dilute molecules,
representing regulatory proteins. [...]"

Panels B and C of the paper's Fig. 1, side by side: the lattice polymer (green)
with its promoter (orange) and enhancer (red), among ligands (blue) that are
dilute (B, c = 0.003), or condensed on the regulatory regions, holding promoter
and enhancer together (C, c = 0.017). Below them, panel A: the model (the
chain's regions, the particle interactions, the transcription model).

B and C are the published plots as extracted by figure_production/fig1/
extract_art.py (art_b.pdf, art_c.pdf), which carry no lettering: the axis and
tick labels that Fig. 1 sets on top of them are left out. Each is cropped to its
ink (Ghostscript's bounding box). A is cropped from the published figure
(slides/figures/scene22/fig1.pdf) with its own lettering, drawn at the slides'
label size; the panel letter is left out. All are placed by pdflatex, so they
stay vector.

Output: slides/figures/scene22/scene22.pdf, a full 128 x 96 mm slide.

Usage (from slides/):
    python3 scripts/scene22/scene22.py
Requires the `gs` and `pdflatex` command-line tools.
"""

import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # slides/scripts, for `common`

from common.drawing import FIGURES_DIR, LABEL_PT, SLIDE_MM, SLIDE_UNITS, UNIT_MM

FIG1_DIR = Path(__file__).resolve().parents[3] / "figure_production" / "fig1"
FIG1 = FIGURES_DIR / "scene22" / "fig1.pdf"
PT_MM = 25.4 / 72

# The 3D plots, B then C, side by side and centred at the top: top edge, B's
# width (C is drawn at the same scale) and the gap between them (slide units).
PANELS = ["art_b.pdf", "art_c.pdf"]
PLOTS = dict(top=14.65, width=6.5, gap=0.9)

# Panel A, a viewport of fig1.pdf (x0, y0, x1, y1 in points, y up): from the
# leftmost bead to "Ligand", and from the labels over the chain to the last
# interaction. Centred under the plots, A_GAP below them.
A_BOX = (1.9, 216.5, 188.6, 330.8)
A_GAP = 0.45
# The panel letter stands above the leftmost bead, and its lowest strokes fall
# inside the viewport: they are painted white (box in points, as A_BOX;
# nothing else is drawn there). The box overhangs the viewport's top and left
# edges, or rasterizers leave a speck of the clipped letter along them.
A_LETTER_MASK = (1.0, 329.0, 7.8, 331.8)
# Fig. 1's lettering is 9 pt; this scale draws it at the slides' label size.
FIG_SCALE = LABEL_PT / 9

TEX = r"""\documentclass[border=0pt]{standalone}
\usepackage{graphicx}
\usepackage{tikz}
\begin{document}
\begin{tikzpicture}[x=1mm, y=1mm]
\useasboundingbox (0, 0) rectangle (%(width)s, %(height)s);
%(nodes)s
\end{tikzpicture}
\end{document}
"""
NODE = (r"\node[anchor=north west, inner sep=0] at (%.3f, %.3f) "
        r"{\includegraphics[scale=%.5f, viewport=%.3f %.3f %.3f %.3f, clip]{%s}};")
MASK = r"\fill[white] (%.3f, %.3f) rectangle (%.3f, %.3f);"


def ink_box(pdf):
    """(x0, y0, x1, y1), in points, of everything drawn on the page."""
    out = subprocess.run(["gs", "-q", "-dBATCH", "-dNOPAUSE", "-sDEVICE=bbox", str(pdf)],
                         capture_output=True, text=True, check=True).stderr
    return tuple(map(float, re.search(r"%%HiResBoundingBox: (.*)", out).group(1).split()))


def typeset(nodes, scene, name):
    """Compile TikZ `nodes` on a full slide into figures/<scene>/<name>.pdf."""
    out = FIGURES_DIR / scene / f"{name}.pdf"
    out.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp:
        tex = Path(tmp) / f"{name}.tex"
        tex.write_text(TEX % dict(width=SLIDE_MM[0], height=SLIDE_MM[1], nodes="\n".join(nodes)))
        subprocess.run(["pdflatex", "-interaction=nonstopmode", "-halt-on-error", tex.name],
                       cwd=tmp, check=True, capture_output=True)
        shutil.copy(tex.with_suffix(".pdf"), out)
    print(f"wrote {out}")


def plots():
    """TikZ nodes of the 3D plots, and the box they fill on the slide (x0, y0,
    x1, y1 in points, y up, as ink_box()): scene 23 has them in the same place."""
    boxes = [ink_box(FIG1_DIR / name) for name in PANELS]
    scale = PLOTS["width"] * UNIT_MM / ((boxes[0][2] - boxes[0][0]) * PT_MM)
    widths = [(b[2] - b[0]) * scale * PT_MM / UNIT_MM for b in boxes]
    height = max(b[3] - b[1] for b in boxes) * scale * PT_MM / UNIT_MM
    left = (SLIDE_UNITS[0] - sum(widths) - PLOTS["gap"]) / 2
    nodes, x = [], left
    for name, box, width in zip(PANELS, boxes, widths):
        nodes.append(NODE % (x * UNIT_MM, PLOTS["top"] * UNIT_MM, scale, *box, FIG1_DIR / name))
        x += width + PLOTS["gap"]
    to_pt = UNIT_MM / PT_MM
    return nodes, (left * to_pt, (PLOTS["top"] - height) * to_pt,
                   (x - PLOTS["gap"]) * to_pt, PLOTS["top"] * to_pt)


def panel_a(top):
    """TikZ nodes of panel A, centred, its top edge at `top` (slide units)."""
    x0, y0, x1, y1 = A_BOX
    left = (SLIDE_MM[0] - (x1 - x0) * FIG_SCALE * PT_MM) / 2

    def to_mm(x, y):
        """Slide position (mm) of a point of fig1.pdf (points)."""
        return left + (x - x0) * FIG_SCALE * PT_MM, top * UNIT_MM - (y1 - y) * FIG_SCALE * PT_MM

    return [NODE % (left, top * UNIT_MM, FIG_SCALE, *A_BOX, FIG1),
            MASK % (*to_mm(*A_LETTER_MASK[:2]), *to_mm(*A_LETTER_MASK[2:]))]


if __name__ == "__main__":
    nodes, box = plots()
    typeset(nodes + panel_a(box[1] * PT_MM / UNIT_MM - A_GAP), "scene22", "scene22")
