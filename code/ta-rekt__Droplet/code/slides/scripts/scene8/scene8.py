"""
Scene 8 -- "Whether to express a gene or not is a decision every nucleus in this
embryo is making. And it's making it quickly: in this example, in just 4
minutes [...]"

The storyboard drawing as is (embryo through division cycles 1-14), with its
labels typeset in the slides' font, and the blue expression gradient of cycle
14 fading in: the cycle-14 embryo starts identical to cycle 13 (same nuclei, no
colour) and turns blue.

Outputs, in slides/figures/scene8/ (full 128 x 96 mm pages):
    scene8_1..4.pdf  four stills with the blue at opacity 0, 1/3, 2/3 and 1, to
                     step through by hand (what the deck uses for now); their
                     line art is a JPEG, the blue a vector overlay
    scene8_fade.pdf  the animation, one frame per page, for the `animate`
                     package: only the outline of the cycle-14 embryo, with its
                     blue fill at increasing opacity
    scene8_base.pdf  everything else (drawing and labels, transparent
                     background), laid over the animation in the beamer frame
                     so the nuclei show on top of the blue, as in the storyboard

Usage (from slides/):
    python3 scripts/scene8/scene8.py
"""

import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # slides/scripts, for `common`

import numpy as np

from common.artwork import SVG, StoryboardScene
from common.drawing import FIGURES_DIR

# Shifted left so the labels, larger than in the sketch, fit to the right.
SCENE = StoryboardScene(8, narration="text498", shift=(-0.5, 0.0))
OUT_DIR = FIGURES_DIR / "scene8"

# In the storyboard the cycle-14 embryo's outline is filled with the blue
# gradient, and so is each of its nuclei (copies of the cycle-13 ones) -- the
# same gradient over the same spot, so the nuclei's own fills are redundant
# and are dropped (checked: renders pixel-identical).
CYCLE14_OUTLINE = "path197-7-1-0"

STILL_OPACITIES = [0, 1 / 3, 2 / 3, 1]
# In the stills, the embryo drawings (hundreds of paths) are flattened to a JPEG
# at this resolution, with the blue gradient kept as a vector overlay.
RASTER_DPI = 500

FPS = 12          # must match the frame rate given to \animategraphics
HOLD_FRAMES = 6   # cycle 14 shown uncoloured first (0.5 s)...
FADE_FRAMES = 24  # ...then the blue fades in (2 s)


def blue_at(opacity, part="all"):
    """Edit setting the opacity of the blue fill of cycle 14, keeping `part` of
    the drawing: "all", only the cycle-14 "outline" (the animated part), or the
    "rest". The outline is painted stroke first, then fill (as in the
    storyboard), so it must stay whole in the frames rather than in the base."""
    def edit(root):
        layer = root.find(SVG + "g")
        for e in list(layer):
            style = e.get("style") or ""
            is_outline = e.get("id") == CYCLE14_OUTLINE
            if (part == "outline" and not is_outline) or (part == "rest" and is_outline):
                layer.remove(e)
            elif is_outline:
                style = re.sub(r"(^|;)fill-opacity:[^;]*", "", style)
                e.set("style", f"{style};fill-opacity:{opacity:.4f}")
            elif "fill:url(#" in style:  # a cycle-14 nucleus
                e.set("style", re.sub(r"fill:url\(#[^)]+\)", "fill:none", style))
    return edit


if __name__ == "__main__":
    for k, opacity in enumerate(STILL_OPACITIES, start=1):
        SCENE.export(OUT_DIR / f"scene8_{k}.pdf", [blue_at(opacity)], raster_dpi=RASTER_DPI)

    t = np.linspace(0, 1, FADE_FRAMES)
    opacities = [0.0] * HOLD_FRAMES + list(t * t * (3 - 2 * t))  # ease in and out
    SCENE.export(OUT_DIR / "scene8_fade.pdf", [blue_at(a, "outline") for a in opacities],
                 labels=False)
    SCENE.export(OUT_DIR / "scene8_base.pdf", [blue_at(0.0, "rest")])
