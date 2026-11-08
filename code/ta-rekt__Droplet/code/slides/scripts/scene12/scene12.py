"""
Scene 12 -- "Let's think about what you were actually doing when you attempted
this problem. You looked at the transcription trace, and asked yourself: how
likely is it that this trace came from the low concentration versus the high
concentration? Which is the same as asking whether this ratio is larger than 1"

The storyboard's equations, in the slides' font, with the
sketched little traces as a pictogram (\trace, see common/telegraph.py):

    scene12_1  P(trace | c_low) >? P(trace | c_high), and the ratio of the two > 1
    scene12_2  the same, naming the ratio R(t) =, and the drift-diffusion
               equation log R(t) = V t + sqrt(2D) xi

Outputs: slides/figures/scene12/scene12_<k>.{pdf,svg}, full 128 x 96 mm slides.

Usage (from slides/):
    python3 scripts/scene12/scene12.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # slides/scripts, for `common`

from common.drawing import TEXT_COLOR, new_slide, placer, save
from common.molecules import POINTS_PER_UNIT
from common.telegraph import (draw_likelihood, draw_likelihood_ratio, likelihood_width,
                              text_width)

EQUATION_PT = 15
CENTER_X = 10.0
GAP = 0.35 * EQUATION_PT  # between the parts of an equation (points)
# Baselines of the three lines (slide units).
QUESTION_Y = 11.6
RATIO_Y = 8.0
DRIFT_Y = 4.4

COMPARE = r"$\overset{?}{>}$"
NAME, TAIL = "$R(t) =$", "$> 1$"
DRIFT = r"$\log R(t) = V t + \sqrt{2D}\,\xi$"


def equation(ax, xy, text, ha="center"):
    return ax.text(*xy, text, ha=ha, va="baseline", fontsize=EQUATION_PT, color=TEXT_COLOR)


def draw_question(ax):
    """P(trace | c_low) ?> P(trace | c_high), centred."""
    place = placer((CENTER_X, QUESTION_Y), 1 / POINTS_PER_UNIT, 0)
    compare = text_width(COMPARE, EQUATION_PT)
    width = (likelihood_width("low", EQUATION_PT) + likelihood_width("high", EQUATION_PT)
             + compare + 2 * GAP)
    x = draw_likelihood(ax, place, (-width / 2, 0.0), "low", EQUATION_PT)
    equation(ax, place((x + GAP, 0.0)), COMPARE, ha="left")
    draw_likelihood(ax, place, (x + 2 * GAP + compare, 0.0), "high", EQUATION_PT)


def draw_ratio(ax, named):
    """The ratio > 1, with "R(t) =" in front if `named`. The fraction stays put
    when "R(t) =" appears; the line with it is centred."""
    fraction = max(likelihood_width(which, EQUATION_PT) for which in ("low", "high"))
    name, tail = (text_width(text, EQUATION_PT) for text in (NAME, TAIL))
    middle = (name - tail) / 2  # of the fraction, in points from the centre of the line
    place = placer((CENTER_X, RATIO_Y), 1 / POINTS_PER_UNIT, 0)
    draw_likelihood_ratio(ax, place((middle, 0.0)), EQUATION_PT, log=False)
    equation(ax, place((middle + fraction / 2 + GAP, -0.33 * EQUATION_PT)), TAIL, ha="left")
    if named:
        equation(ax, place((middle - fraction / 2 - GAP, -0.33 * EQUATION_PT)), NAME, ha="right")


def draw_slide(name, drift):
    fig, ax = new_slide()
    draw_question(ax)
    draw_ratio(ax, named=drift)
    if drift:
        equation(ax, (CENTER_X, DRIFT_Y), DRIFT)
    save(fig, "scene12", name)


if __name__ == "__main__":
    draw_slide("scene12_1", drift=False)
    draw_slide("scene12_2", drift=True)
