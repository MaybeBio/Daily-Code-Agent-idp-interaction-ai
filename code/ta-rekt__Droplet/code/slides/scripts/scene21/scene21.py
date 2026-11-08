"""
Scene 21 -- "Another ingredient has been emerging in this recipe, building on
recent experimental evidence: surface condensation. [...] condensation of a
liquid-like droplet on the surface of DNA. this is a physical phenomenon that
has been studied since the 1970s, called pre-wetting. [...] below the
concentration at which it saturates."

Left: the scene-4 locus, with transcription factors dilute in the bulk
and a dense droplet of them condensed on the DNA at the promoter. Right: clips
from papers (screenshots in slides/figures/scene21/), each with its citation
below:
    Morin et al., 2022 (Klf4 condensates on DNA; Klf4 droplets vs [Klf4], [DNA])
    Rouches & Machta, 2025 (polymer with strongly interacting stretches)

Output: slides/figures/scene21/scene21.{pdf,svg}, a full 128 x 96 mm slide.

Usage (from slides/):
    python3 scripts/scene21/scene21.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # slides/scripts, for `common`

import numpy as np
from PIL import Image
from matplotlib.patches import Rectangle

from common.drawing import FIGURES_DIR, TEXT_COLOR, draw_shine, new_slide, save, zoomed_axes
from common.locus import Locus
from common.molecules import (DNA_LINE_PT, TF, circle_points, draw_dna, draw_legend,
                              draw_scattered_tfs, scatter_proteins)

CLIPS_DIR = FIGURES_DIR / "scene21"

# ---------------------------------------------------------------------------
# Left: the locus, a droplet at its promoter, dilute TFs
# ---------------------------------------------------------------------------

LOCUS = Locus(offset=(1.16, -0.34))  # the scene-4 locus, where the sketch copies it
DROPLET_CENTER = LOCUS.dna.at(LOCUS.promoter)[0]
# TFs of the droplet are packed within DROPLET_RADIUS of its centre (closer than
# in the bulk: they overlap in projection); a soft glow marks the dense phase.
DROPLET = dict(n_tfs=32, radius=0.55, spacing=0.11, seed=21)
GLOW = dict(radius=1.0, alpha=0.027, layers=14)
FREE_TFS = dict(n_tfs=20, seed=1021)
REGION = (0.5, 10.4, 5.0, 14.4)  # where the free TFs are scattered
# The cartoon above is drawn `zoom` times its size, a point p of it landing at
# shift + zoom * p, to leave room for the clips; its DNA's lines are thinned
# with it, as in scene 20.
CARTOON = dict(zoom=0.8, shift=(0.0, 1.3))

# ---------------------------------------------------------------------------
# Right: clips from papers, in two columns, each with its citation below
# ---------------------------------------------------------------------------

# Crop boxes (left, top, right, bottom, in pixels) trim the screenshots' white
# margins; the Rouches clip's last row holds a stray piece of the next panel.
CLIPS = {
    "morin_dna": dict(file="morin1.png", crop=(20, 22, 810, 232)),
    "morin_droplets": dict(file="morin2.png", crop=(6, 8, 868, 746)),
    "rouches": dict(file="rouches.png", crop=(37, 31, 1115, 1401)),
}
# Columns: left edge and width (slide units), clips from top to bottom. The
# Rouches column is as wide as its author line.
COLUMNS = [
    dict(x=8.61, width=5.25, clips=["morin_dna", "morin_droplets"],
         title="Sequence-dependent surface\ncondensation of a pioneer\n"
               "transcription factor on DNA",
         authors="Morin et al., 2022"),
    dict(x=14.31, width=5.44, clips=["rouches"],
         title="Protein-DNA Co-condensation\nis Prewetting to a Collapsed\nPolymer",
         authors="Rouches & Machta, 2025"),
]
CLIPS_TOP = 14.2
CLIP_GAP = 0.19            # between clips of a column
CITATION_GAP = 0.27        # between a column's last clip and its citation
CITATION_PT = 7.4
CITATION_LINE = 0.46       # from one line of a citation to the next
FRAME = dict(color="#C8C8C8", lw=0.5)


def droplet_tfs():
    """TFs packed in a disc around the droplet's centre, at random depths and
    orientations, ordered from back to front (as scatter_proteins())."""
    rng = np.random.default_rng(DROPLET["seed"])
    placed = []
    while len(placed) < DROPLET["n_tfs"]:
        r, a = DROPLET["radius"] * np.sqrt(rng.uniform()), rng.uniform(0, 2 * np.pi)
        c = DROPLET_CENTER + r * np.array([np.cos(a), np.sin(a)])
        if all(np.linalg.norm(q["center"] - c) > DROPLET["spacing"] for q in placed):
            placed.append(dict(center=c, depth=rng.uniform(), angle=rng.uniform(0, 360)))
    return sorted(placed, key=lambda q: q["depth"])


def draw_locus(fig, shift=CARTOON["shift"], zoom=CARTOON["zoom"]):
    """The cartoon, placed as CARTOON says (or by `shift` and `zoom`: scene 24
    copies it); returns the axes it is drawn on, in scene-4 coordinates."""
    ax = zoomed_axes(fig, shift, zoom)
    draw_shine(ax, circle_points(DROPLET_CENTER, GLOW["radius"], n=64), DROPLET_CENTER,
               TF["fill"][:3], z=0.3, layers=GLOW["layers"], alpha=GLOW["alpha"])
    draw_dna(ax, LOCUS.dna, z=1, line_pt=tuple(zoom * w for w in DNA_LINE_PT))
    draw_scattered_tfs(ax, droplet_tfs(), z_back=0.5, z_front=6)
    avoid = LOCUS.obstacles() + [(DROPLET_CENTER, GLOW["radius"])]
    free = scatter_proteins(FREE_TFS["n_tfs"], REGION, np.random.default_rng(FREE_TFS["seed"]),
                            avoid)
    draw_scattered_tfs(ax, free, z_back=0.5, z_front=6)
    return ax


def draw_clips(ax):
    for column in COLUMNS:
        x0, width, top = column["x"], column["width"], CLIPS_TOP
        for name in column["clips"]:
            clip = CLIPS[name]
            image = Image.open(CLIPS_DIR / clip["file"]).convert("RGB").crop(clip["crop"])
            height = width * image.height / image.width
            # interpolation="none": the PDF keeps the screenshot's own pixels.
            ax.imshow(np.asarray(image), extent=(x0, x0 + width, top - height, top),
                      interpolation="none", aspect="auto", zorder=2)
            ax.add_patch(Rectangle((x0, top - height), width, height, fc="none",
                                   ec=FRAME["color"], lw=FRAME["lw"], zorder=3))
            top -= height + CLIP_GAP
        top += CLIP_GAP
        title_top = top - CITATION_GAP
        ax.text(x0 + width / 2, title_top, column["title"], style="italic", ha="center", va="top",
                fontsize=CITATION_PT, color=TEXT_COLOR)
        authors_top = title_top - (column["title"].count("\n") + 1) * CITATION_LINE
        ax.text(x0 + width / 2, authors_top, column["authors"], ha="center", va="top",
                fontsize=CITATION_PT, color=TEXT_COLOR)


if __name__ == "__main__":
    fig, ax = new_slide()
    draw_locus(fig)
    draw_clips(ax)
    draw_legend(ax, ["DNA", "Promoter", "Gene", "TF"])
    save(fig, "scene21", "scene21")
