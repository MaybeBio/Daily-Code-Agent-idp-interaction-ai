"""
Scene 13 -- The sequential probability ratio test (SPRT) of scene 12, watched
in real time on simulated transcription.

A transcription time series at the true concentration (top); below it, on the
same time axis, its log-likelihood ratio, c_low vs c_high, on a linear scale
from -2K to 2K. The test stops as soon as the ratio leaves (-K, K), the dashed
horizontal lines: +K decides c_low, -K decides c_high. The vertical dashed line
marks that decision time, through both plots, and the box of the decision taken
lights up. The ratio is drawn on past it, until it leaves the plot at +-2K.

    scene13_1  true c = c_low: a typical trace, decided correctly
    scene13_2  true c = c_high: a typical trace, decided correctly
    scene13_3  true c = c_low, decided incorrectly; K = log (1 - eps) / eps is
               added on the side, with what eps is
    scene13_4  the trace of 13_1 again, with c_high further from c_low: the
               ratio climbs faster and the decision comes sooner

Each trace is the first random seed meeting the criteria of typical_run() or
misled_run(), so they are reproducible and easy to re-draw.

Outputs: slides/figures/scene13/scene13_<k>.{pdf,svg}, full 128 x 96 mm slides.

Usage (from slides/):
    python3 scripts/scene13/scene13.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # slides/scripts, for `common`

import numpy as np
from matplotlib.colors import to_rgba

from common.drawing import AXIS_COLOR, LABEL_PT, TEXT_COLOR, draw_label, new_slide, save
from common.telegraph import (HIGHLIGHT, TimeAxes, draw_likelihood_ratio, draw_time_series,
                              first_seed, log_likelihood_ratio, telegraph, typical)

# ---------------------------------------------------------------------------
# The test
# ---------------------------------------------------------------------------

C_LOW = 1.0
EPSILON = 0.05                         # error level: simulated errors come out at 3-5%
K = np.log((1 - EPSILON) / EPSILON)    # the test stops at log ratio +-K
DURATION = 36.0                        # time shown, on every slide

# dc_over_c sets c_high = c_low (1 + dc/c), `truth` the concentration simulated.
BASE_DC, LARGE_DC = 1.0, 2.0
SLIDES = {
    "scene13_1": dict(truth="low", dc_over_c=BASE_DC, explain_k=False),
    "scene13_2": dict(truth="high", dc_over_c=BASE_DC, explain_k=False),
    "scene13_3": dict(truth="low", dc_over_c=BASE_DC, explain_k=True),
    "scene13_4": dict(truth="low", dc_over_c=LARGE_DC, explain_k=True),
}


def concentrations(dc_over_c):
    return {"low": C_LOW, "high": C_LOW * (1 + dc_over_c)}


def ratio(seed, truth, dc_over_c, duration=DURATION):
    """A trace at the true concentration, and its log ratio as a polyline (t, log R)."""
    c = concentrations(dc_over_c)
    trace = telegraph(c[truth], duration, seed)
    return trace, log_likelihood_ratio(trace, c["low"], c["high"])


def stopped(t, log_r, level):
    """The polyline (t, log R), cut where it first reaches +-level (all of it
    if it never does). log R rises steadily or jumps down, so that point is on
    a straight piece, or at a jump."""
    beyond = np.flatnonzero(np.abs(log_r) >= level)
    if not len(beyond):
        return t, log_r
    i = beyond[0]
    boundary = np.sign(log_r[i]) * level
    fraction = (boundary - log_r[i - 1]) / (log_r[i] - log_r[i - 1])
    return (np.append(t[:i], t[i - 1] + fraction * (t[i] - t[i - 1])),
            np.append(log_r[:i], boundary))


def decision(t, log_r):
    """The outcome of the test: (decision time, "low" or "high"), or None if
    the log ratio never leaves (-K, K)."""
    t_stop, log_r_stop = stopped(t, log_r, K)
    if abs(log_r_stop[-1]) < K:
        return None
    return t_stop[-1], "low" if log_r_stop[-1] > 0 else "high"


# ---------------------------------------------------------------------------
# Choosing the traces
# ---------------------------------------------------------------------------

TYPICAL_TIME_TOLERANCE = 0.15   # decision time within this fraction of the median
GRAZING = 0.9                   # before deciding, the ratio stays within this fraction of K...
CARRIES_ON = 1.4                # ...and it ends beyond this many K, towards the truth
MEDIAN_RUNS = 2000
LONG_ENOUGH = 400.0             # every simulated run reaches a decision by then
MISLED_TIMES = (10.0, 20.0)     # when the wrong decision is taken, on 13_3
MISLED_ON_TOLERANCE = 0.05      # its trace looks like c_low overall


def median_decision_time(truth, dc_over_c):
    return np.median([decision(*ratio(seed, truth, dc_over_c, LONG_ENOUGH)[1])[0]
                      for seed in range(MEDIAN_RUNS)])


def stays_on_plot(log_r):
    return np.abs(log_r).max() < 2 * K


def carries_on(log_r):
    return stays_on_plot(log_r) and abs(log_r[-1]) > CARRIES_ON * K


def no_near_miss(t, log_r):
    """Does the ratio stay clear of +-K until it crosses, so that the decision
    time cannot be misread on the plot?"""
    return np.abs(stopped(t, log_r, K)[1][:-1]).max() < GRAZING * K


def typical_run(seed, truth, dc_over_c, median):
    """Is this a typical trace, decided correctly in about the median time?"""
    trace, (t, log_r) = ratio(seed, truth, dc_over_c)
    outcome = decision(t, log_r)
    return (outcome is not None and outcome[1] == truth
            and abs(outcome[0] / median - 1) < TYPICAL_TIME_TOLERANCE
            and no_near_miss(t, log_r)
            and typical(concentrations(dc_over_c)[truth])(trace))


def misled_run(seed):
    """At c_low: is the wrong decision taken midway, the ratio then recovering
    on the side of the truth?"""
    trace, (t, log_r) = ratio(seed, "low", BASE_DC)
    outcome = decision(t, log_r)
    return (outcome is not None and outcome[1] == "high"
            and MISLED_TIMES[0] < outcome[0] < MISLED_TIMES[1]
            and no_near_miss(t, log_r) and stays_on_plot(log_r) and log_r[-1] > 0
            and typical(C_LOW, MISLED_ON_TOLERANCE)(trace))


def choose_seeds():
    """Seeds of the slides. 13_4 shows the trace of 13_1, so that seed has to
    suit both concentration differences."""
    median_low = median_decision_time("low", BASE_DC)
    median_high = median_decision_time("high", BASE_DC)
    median_large = median_decision_time("low", LARGE_DC)
    print(f"median decision times: {median_low:.1f} at c_low, {median_high:.1f} at c_high, "
          f"{median_large:.1f} at c_low with the larger difference")

    def low_and_large(seed):
        return (typical_run(seed, "low", BASE_DC, median_low)
                and carries_on(ratio(seed, "low", BASE_DC)[1][1])
                and typical_run(seed, "low", LARGE_DC, median_large))

    def high(seed):
        return (typical_run(seed, "high", BASE_DC, median_high)
                and carries_on(ratio(seed, "high", BASE_DC)[1][1]))

    first = first_seed(low_and_large, 1)
    return {"scene13_1": first, "scene13_2": first_seed(high, 1),
            "scene13_3": first_seed(misled_run, 1), "scene13_4": first}


# ---------------------------------------------------------------------------
# Layout (slide units)
# ---------------------------------------------------------------------------

ACCENT = "#D9532B"  # the decision time
TRACE_PT = 0.6
RATIO_PT = 1.0
DASHED = dict(color=AXIS_COLOR, lw=1.1, dashes=(4, 3), zorder=3)

# Left: the time series, and the log ratio below it, on the same time axis.
PLOT_LEFT, PLOT_RIGHT = 3.0, 15.5
SPAN = (PLOT_LEFT + 0.07, PLOT_RIGHT - 0.4)
TIME_SERIES_AXES = TimeAxes(origin=(PLOT_LEFT, 11.3), tips=(PLOT_RIGHT, 13.67), span=SPAN)
RATIO_ZERO = 5.55     # height of log ratio 0
K_UNIT = 1.95         # slide units per K of log ratio
RATIO_AXES = TimeAxes(origin=(PLOT_LEFT, RATIO_ZERO - 2 * K_UNIT),
                      tips=(PLOT_RIGHT, RATIO_ZERO + 2 * K_UNIT + 0.45), span=SPAN)
RATIO_TICKS = [(2, "$+2K$"), (1, "$+K$"), (0, "$0$"), (-1, "$-K$"), (-2, "$-2K$")]
RATIO_LABEL_X = 0.95

# Right: what is true, the decisions at the ends of +-K, and what K is.
SIDE_X = 15.95
TRUTH_Y = [12.95, 12.25]
K_LABEL_Y = RATIO_ZERO + 0.45
EPSILON_LABEL_Y = RATIO_ZERO - 0.6


# ---------------------------------------------------------------------------
# Panels
# ---------------------------------------------------------------------------

def ratio_y(log_r):
    return RATIO_ZERO + np.asarray(log_r) / K * K_UNIT


def draw_ratio(ax, t, log_r):
    """The log ratio from -2K to 2K, with +-K dashed."""
    RATIO_AXES.draw(ax)
    x0 = RATIO_AXES.origin[0]
    for multiple, text in RATIO_TICKS:
        y = ratio_y(multiple * K)
        ax.plot([x0 - 0.1, x0], [y, y], color=AXIS_COLOR, lw=1.3, zorder=3)
        draw_label(ax, (x0 - 0.2, y), text, ha="right")
    for level in (K, -K):
        ax.plot([x0, RATIO_AXES.tips[0]], [ratio_y(level)] * 2, **DASHED)
    draw_likelihood_ratio(ax, (RATIO_LABEL_X, RATIO_ZERO), rotation=90)

    t_shown, log_r_shown = stopped(t, log_r, 2 * K)
    ax.plot(RATIO_AXES.x(t_shown, DURATION), ratio_y(log_r_shown), color="#222222", lw=RATIO_PT,
            solid_joinstyle="miter", zorder=4)


def draw_decisions(ax, taken):
    """The two decisions, at the ends of +-K; the one `taken` lit up."""
    for level, which in ((K, "low"), (-K, "high")):
        ax.text(SIDE_X, ratio_y(level), rf"Decision: $c_{{\mathrm{{{which}}}}}$", ha="left",
                va="center", fontsize=LABEL_PT - 1, color=TEXT_COLOR, zorder=20,
                bbox=dict(boxstyle="square,pad=0.35", ec=AXIS_COLOR, lw=0.8,
                          fc=to_rgba(HIGHLIGHT["color"], HIGHLIGHT["alpha"]) if which == taken
                          else "white"))


def draw_decision_time(ax, t):
    """A vertical dashed line at t, through both plots."""
    x = RATIO_AXES.x(t, DURATION)
    bottom = RATIO_AXES.origin[1]
    ax.plot([x, x], [bottom, TIME_SERIES_AXES.tips[1] - 0.35], color=ACCENT, lw=1.4,
            dashes=(3, 2), zorder=5)
    ax.text(x, bottom - 0.5, r"$T_{\rm dec}$", ha="center", va="center", fontsize=LABEL_PT,
            color=ACCENT)


def draw_truth(ax, truth, dc_over_c):
    """The concentration simulated, and the difference tested for."""
    lines = (rf"True $c = c_{{\mathrm{{{truth}}}}}$",
             rf"$c_{{\mathrm{{high}}}} = {1 + dc_over_c:.0f}\,c_{{\mathrm{{low}}}}$")
    for text, y in zip(lines, TRUTH_Y):
        draw_label(ax, (SIDE_X, y), text, ha="left")


def draw_k_explained(ax):
    draw_label(ax, (SIDE_X, K_LABEL_Y), r"$K = \log\dfrac{1-\epsilon}{\epsilon}$", ha="left")
    draw_label(ax, (SIDE_X, EPSILON_LABEL_Y), r"$\epsilon$: error level", ha="left")


def draw_slide(name, truth, dc_over_c, explain_k, seed):
    fig, ax = new_slide()
    trace, (t, log_r) = ratio(seed, truth, dc_over_c)
    t_decision, taken = decision(t, log_r)
    print(f"{name}: seed {seed}, true c_{truth}, decided c_{taken} at t = {t_decision:.1f}")
    draw_time_series(ax, trace, DURATION, TIME_SERIES_AXES, lw=TRACE_PT, time_label=False)
    draw_ratio(ax, t, log_r)
    draw_decisions(ax, taken)
    draw_decision_time(ax, t_decision)
    draw_truth(ax, truth, dc_over_c)
    if explain_k:
        draw_k_explained(ax)
    save(fig, "scene13", name)


if __name__ == "__main__":
    seeds = choose_seeds()
    for name, spec in SLIDES.items():
        draw_slide(name, seed=seeds[name], **spec)
