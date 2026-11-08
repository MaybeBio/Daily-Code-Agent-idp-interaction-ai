"""
Scene 7 -- "Surprisingly, cell nuclei still manage to control the expression of
genes very precisely using the concentration of transcription factors [...]"

The storyboard drawing as is (fly embryos: TF gradient above, expression of its
target gene below), with its labels typeset in the slides' font.

Output: slides/figures/scene7/scene7.pdf, a full 128 x 96 mm slide.

Usage (from slides/):
    python3 scripts/scene7/scene7.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # slides/scripts, for `common`

from common.artwork import StoryboardScene
from common.drawing import FIGURES_DIR

# Shifted left so the labels, larger than in the sketch, fit to the right.
SCENE = StoryboardScene(7, narration="text200", shift=(-0.6, 0.0))
# The embryo drawings (hundreds of paths) are flattened to a JPEG at this
# resolution, with their colour gradients kept as vector overlays.
RASTER_DPI = 500

if __name__ == "__main__":
    SCENE.export(FIGURES_DIR / "scene7" / "scene7.pdf", raster_dpi=RASTER_DPI)
