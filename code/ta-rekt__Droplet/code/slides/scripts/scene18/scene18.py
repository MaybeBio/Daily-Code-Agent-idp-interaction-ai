"""
Scene 18 -- "How well are we doing in our objective of making a fast decision
in under 4 minutes? Here is a comparison of the two promoter models we have
seen so far in terms of decision time."

Mean decision time <T_dec> = (K/V) tanh(K/2), V = 1/2 (dc/c)^2 Fdot(c),
against concentration, for the three models of transcription seen so far,
all with the same, generic rates (not optimized): binding at the diffusion
limit, kb = 4aD per site (the SPRT paper's value, in its concentration
units), and ku = kon = koff = RATE. Counter-clockwise rates are zero, as in
the Fisher rate of scene 15:

    2 states   OFF -(c kb)-> ON -(ku)-> OFF:  Fdot = c kb ku / (c kb + ku)
    4 states   OFF time: binding (c kb), then activation (kon); ON time:
               unbinding (ku), then deactivation (koff)
    6 states   two sites, transcription only with both bound: OFF time:
               binding (2 c kb, then c kb), then activation (kon); ON time:
               unbinding (2 ku, then ku), then deactivation (koff)

Only OFF times depend on c; Fdot = F_off(c) / mean cycle time, with F_off the
Fisher information of an OFF time on ln c: the closed form of scene 15 for 4
states, numerical integration for 6 (checked against the closed form when
the script runs). Decisions between c and 1.1 c (dc/c = 0.1, as in scene 11)
at the paper's error level, 0.32.

Output: slides/figures/scene18/scene18.{pdf,svg}, a full 128 x 96 mm slide.

Usage (from slides/):
    python3 scripts/scene18/scene18.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # slides/scripts, for `common`

import numpy as np
from matplotlib.ticker import FuncFormatter
from scipy.linalg import expm

from common.drawing import AXIS_COLOR, LABEL_PT, TEXT_COLOR, new_slide, save
from scene15.scene15 import fisher_off
from scene16.scene16 import FOUR_AD, paper_axes

# ---------------------------------------------------------------------------
# Models and decision times
# ---------------------------------------------------------------------------

RATE = 1.0          # ku = kon = koff (1/s)
DC_OVER_C = 0.1
ERROR = 0.32
K = np.log((1 - ERROR) / ERROR)


def chain_fisher(rates_of_c, c, h=1e-4, points=800):
    """Fisher information on ln c of the time to go through a chain of
    exponential stages, whose rates at concentration c are rates_of_c(c):
    the score d ln f / d ln c by finite differences, the expectation by
    quadrature on a log-spaced grid of times."""
    def density(cc, t):
        rates = np.asarray(rates_of_c(cc), float)
        generator = np.diag(-rates) + np.diag(rates[:-1], 1)
        return expm(generator[None] * t[:, None, None])[:, 0, -1] * rates[-1]

    rates = np.asarray(rates_of_c(c), float)
    t = np.geomspace(1e-4 / rates.max(), 60 / rates.min(), points)
    f = density(c, t)
    score = (np.log(density(c * np.exp(h), t)) - np.log(density(c * np.exp(-h), t))) / (2 * h)
    return np.trapezoid(f * score ** 2 * t, np.log(t))


def two_state(c):
    x = FOUR_AD * c
    return x * RATE / (x + RATE)


def four_state(c):
    x = FOUR_AD * c
    return fisher_off(x, RATE) / (1 / x + 1 / RATE + 1 / RATE + 1 / RATE)


def six_state(c):
    x = FOUR_AD * c
    off = chain_fisher(lambda cc: [2 * FOUR_AD * cc, FOUR_AD * cc, RATE], c)
    cycle = 1 / (2 * x) + 1 / x + 1 / RATE + 1 / (2 * RATE) + 1 / RATE + 1 / RATE
    return off / cycle


def decision_minutes(fisher_rate):
    return K * np.tanh(K / 2) / (0.5 * DC_OVER_C ** 2 * fisher_rate) / 60


def check_numerics(cs, tolerance=1e-4):
    """The numerical Fisher information against the closed form, for 4 states."""
    worst = max(abs(chain_fisher(lambda cc: [FOUR_AD * cc, RATE], c)
                    / fisher_off(FOUR_AD * c, RATE) - 1) for c in cs)
    print(f"numerical vs closed-form Fisher information (4 states): {worst:.1e}")
    assert worst < tolerance


# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------

AXES = (2.6, 2.5, 11.4, 9.4)  # left, bottom, width, height (slide units)
X_LIMITS, Y_LIMITS = (1e-1, 1e3), (0.5, 1e3)
Y_LABEL = "Decision time (A.U.)"
TARGET_MINUTES = 4
# Model, colour, label (lines) and its placement at the right end of the curve.
MODELS = [
    (two_state, "#222222", ["2 state model"], "center"),
    (four_state, "#2F62C8", ["4 state model"], "bottom"),
    (six_state, "#D9532B", ["6 state model", "with cooperativity"], "top"),
]
LINE_STEP = 1.25  # between label lines, in units of the font size


def label_lines(ax, xy, lines, color, va):
    """Lines of text left-aligned at xy (display offset), stacked from va."""
    n = len(lines)
    first = {"bottom": n - 1, "center": (n - 1) / 2, "top": 0}[va]
    for i, line in enumerate(lines):
        ax.annotate(line, xy, xytext=(8, (first - i) * LINE_STEP * LABEL_PT),
                    textcoords="offset points", ha="left", va=va if n == 1 else "center",
                    fontsize=LABEL_PT, color=color, annotation_clip=False)


if __name__ == "__main__":
    cs = np.geomspace(*X_LIMITS, 160)
    check_numerics(cs[::8])

    fig, _ = new_slide()
    ax = paper_axes(fig, AXES, X_LIMITS, Y_LIMITS, Y_LABEL)
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"${v:g}$"))
    for model, color, lines, va in MODELS:
        minutes = decision_minutes(np.array([model(c) for c in cs]))
        ax.plot(cs, minutes, color=color, lw=1.8, zorder=4)
        label_lines(ax, (cs[-1], minutes[-1]), lines, color, va)
        print(f"{lines[0]}: fastest {minutes.min():.2f} min at c = {cs[minutes.argmin()]:.3g}")
    # ax.axhline(TARGET_MINUTES, color=AXIS_COLOR, lw=1.2, dashes=(5, 4), zorder=3)
    # ax.annotate(f"{TARGET_MINUTES} min", (X_LIMITS[0], TARGET_MINUTES), xytext=(6, 3),
    #             textcoords="offset points", ha="left", va="bottom", fontsize=LABEL_PT,
    #             color=TEXT_COLOR)
    save(fig, "scene18", "scene18")
