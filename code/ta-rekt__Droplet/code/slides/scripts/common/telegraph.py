"""
The two-state ("telegraph") model of transcription, shared by scenes 11-13.

The gene switches OFF -> ON at rate c kb (a TF binds; c is its concentration)
and ON -> OFF at rate ku (it unbinds). This module simulates the model,
computes the likelihood ratio of a trace under two concentrations, and draws
the panels that show them: the two-state graph, the concentration slider, plots
against time, and a little trace pictogram for equations.

Panels are drawn where scene 11 has them; other scenes move them with `offset`
(slide units) or give their own axes.
"""

from dataclasses import dataclass

import numpy as np
from matplotlib.patches import Ellipse
from matplotlib.textpath import TextPath

from .drawing import AXIS_COLOR, LABEL_PT, TEXT_COLOR, draw_arrow, draw_label, placer
from .molecules import GREEN, POINTS_PER_UNIT

# ---------------------------------------------------------------------------
# Model and simulated time series
# ---------------------------------------------------------------------------

K_UNBIND = 1.0  # ON -> OFF rate, ku
K_BIND = 1.0    # OFF -> ON rate per unit concentration, kb (the rate is c kb)


def telegraph(c, duration, seed):
    """Two-state process over [0, duration], started in its steady state.

    Returns (times, states): the process is in states[i] (True = ON) from
    times[i] to times[i + 1]; the last time is `duration`.
    """
    rng = np.random.default_rng(seed)
    state = rng.uniform() < c * K_BIND / (c * K_BIND + K_UNBIND)
    times, states, t = [0.0], [], 0.0
    while t < duration:
        states.append(state)
        t += rng.exponential(1 / (K_UNBIND if state else c * K_BIND))
        times.append(min(t, duration))
        state = not state
    return np.array(times), np.array(states)


def truncated(trace, duration):
    times, states = trace
    keep = times[:-1] < duration
    return np.append(np.minimum(times[:-1][keep], duration), duration), states[keep]


def stats(trace):
    """Number of ON bursts, complete ON bursts, mean complete OFF time, ON fraction."""
    times, states = trace
    lengths = np.diff(times)
    inside = np.zeros(len(states), bool)
    inside[1:-1] = True  # not cut by either end of the window
    off = lengths[~states & inside]
    return dict(bursts=int(states.sum()), complete=int((states & inside).sum()),
                mean_off=off.mean() if len(off) else np.inf,
                on_fraction=lengths[states].sum() / times[-1])


def first_seed(ok, seed):
    """The first seed from `seed` on for which ok(seed) holds."""
    while not ok(seed):
        seed += 1
    return seed


def first_trace(c, duration, ok, seed):
    """The first trace (by seed, from `seed` on) for which ok(trace) holds."""
    return telegraph(c, duration, first_seed(lambda s: ok(telegraph(c, duration, s)), seed))


def typical(c, tolerance=0.02):
    """Criterion: ON fraction within `tolerance` of its expected value c kb / (c kb + ku)."""
    expected = c * K_BIND / (c * K_BIND + K_UNBIND)
    return lambda tr: abs(stats(tr)["on_fraction"] - expected) < tolerance


# ---------------------------------------------------------------------------
# Telling two concentrations apart
# ---------------------------------------------------------------------------

def binding_rate(c):
    """Binding events per unit time in steady state, c kb ku / (c kb + ku):
    one per OFF -> ON -> OFF cycle."""
    return c * K_BIND * K_UNBIND / (c * K_BIND + K_UNBIND)


def log_likelihood_ratio(trace, c_low, c_high):
    """log R(t) = log P(trace until t | c_low) - log P(trace until t | c_high),
    given the initial state; returned as a polyline (t, log R).

    ON stretches are equally likely at both concentrations (same ku), so log R
    stays flat. OFF stretches raise it by (c_high - c_low) kb per unit time
    (waiting long for a TF favours the low concentration), and each switch ON
    lowers it by log(c_high / c_low) (a binding favours the high one).
    """
    times, states = trace
    t, log_r = [0.0], [0.0]
    for start, end, on in zip(times[:-1], times[1:], states):
        if not on:
            t.append(end)
            log_r.append(log_r[-1] + (c_high - c_low) * K_BIND * (end - start))
            if end < times[-1]:  # not cut by the end of the window: switched ON
                t.append(end)
                log_r.append(log_r[-1] - np.log(c_high / c_low))
        else:
            t.append(end)
            log_r.append(log_r[-1])
    return np.array(t), np.array(log_r)


# ---------------------------------------------------------------------------
# Panels
# ---------------------------------------------------------------------------

HIGHLIGHT = dict(color="#FFD54A", alpha=0.55)

# Two-state graph, as placed in scene 11.
STATES = {"OFF": (14.3, 13.22), "ON": (14.3, 10.10)}
STATE_SIZE = (1.0, 0.7)
STATE_LABELS = {"OFF": [("Transcription", 13.25), ("OFF", 12.57)],
                "ON": [("Transcription", 10.40), ("ON", 9.72)]}
STATE_LABEL_X = 16.85
UNBIND = dict(x=14.0, label=r"$k_u$", color="#9B3DB8")
BIND = dict(x=14.6, label=r"$c\,k_b$", color=GREEN["edge"])
RATE_ARROW_SPAN = (10.5, 12.82)
RATE_LABEL_Y = 11.66

# Concentration slider (log scale: 1/2x, 1x, 2x evenly spaced), as placed in
# scene 11.
SLIDER_X = 1.57
SLIDER_1X = 11.2
SLIDER_DOUBLING = 0.73
SLIDER_ENDS = (9.7, 12.65)
SLIDER_TICKS = [(2.0, r"$2\times$"), (1.0, r"$1\times$"), (0.5, r"$\frac{1}{2}\times$")]
SLIDER_COLOR = GREEN["fill"]
SLIDER_TITLE = [("Conc.", 13.45), ("difference", 12.85)]


def draw_two_state_graph(ax, highlight=None, offset=(0.0, 0.0)):
    """OFF and ON states with the rates between them; `highlight` is "OFF" or "ON"."""
    dx, dy = offset
    for name, (x, y) in STATES.items():
        if name == highlight:
            ax.add_patch(Ellipse((x + dx, y + dy), STATE_SIZE[0] + 0.3, STATE_SIZE[1] + 0.3,
                                 fc=HIGHLIGHT["color"], ec="none", alpha=HIGHLIGHT["alpha"],
                                 zorder=2))
        ax.add_patch(Ellipse((x + dx, y + dy), *STATE_SIZE, fc="none", ec=AXIS_COLOR, lw=1.3,
                             zorder=3))
        for text, label_y in STATE_LABELS[name]:
            draw_label(ax, (STATE_LABEL_X + dx, label_y + dy), text)
    low, high = RATE_ARROW_SPAN[0] + dy, RATE_ARROW_SPAN[1] + dy
    draw_arrow(ax, (UNBIND["x"] + dx, low), (UNBIND["x"] + dx, high), UNBIND["color"], z=3,
               lw=1.8)
    draw_arrow(ax, (BIND["x"] + dx, high), (BIND["x"] + dx, low), BIND["color"], z=3, lw=1.8)
    ax.text(UNBIND["x"] + dx - 0.22, RATE_LABEL_Y + dy, UNBIND["label"], ha="right",
            va="center", fontsize=LABEL_PT + 2, color=UNBIND["color"])
    ax.text(BIND["x"] + dx + 0.22, RATE_LABEL_Y + dy, BIND["label"], ha="left", va="center",
            fontsize=LABEL_PT + 2, color=BIND["color"])


PAIR_GAP = 0.15         # half-distance between the two arrows of a pair
ARROW_CLEARANCE = 0.12  # between arrow ends and states
FAINT = "#C4C4C4"       # greyed-out states and arrows
ARROW_HIGHLIGHT_PT = 5


def draw_state_graph(ax, states, arrows, highlight=None, faded=(), labels=None, styles=None):
    """A graph of states (name -> centre), drawn as ellipses -- `highlight`
    (a name) behind a yellow halo, `faded` ones greyed out -- with arrows
    (source, target) between them.

    Each arrow is shifted to its left, so the two arrows between a pair of
    states run side by side in opposite directions (round a square, the
    clockwise ones outside). `labels` {(source, target): text} go beside their
    arrow, on its outer side; `styles` {(source, target): "faint" or
    "highlight"} greys an arrow out or highlights it. Arrows to or from a
    faded state are faint.
    """
    labels, styles = labels or {}, styles or {}
    for name, center in states.items():
        if name == highlight:
            ax.add_patch(Ellipse(center, STATE_SIZE[0] + 0.3, STATE_SIZE[1] + 0.3,
                                 fc=HIGHLIGHT["color"], ec="none", alpha=HIGHLIGHT["alpha"],
                                 zorder=2))
        ax.add_patch(Ellipse(center, *STATE_SIZE, fc="none", lw=1.3, zorder=3,
                             ec=FAINT if name in faded else AXIS_COLOR))
    for source, target in arrows:
        p, q = np.array(states[source]), np.array(states[target])
        direction = (q - p) / np.linalg.norm(q - p)
        left = np.array([-direction[1], direction[0]])
        # Ends clear of the ellipses (semi-axis along the arrow).
        reach = np.abs(direction) @ np.array(STATE_SIZE) / 2 + ARROW_CLEARANCE
        start, end = p + direction * reach + PAIR_GAP * left, q - direction * reach + PAIR_GAP * left
        style = "faint" if {source, target} & set(faded) else styles.get((source, target))
        color = FAINT if style == "faint" else AXIS_COLOR
        if style == "highlight":
            ax.plot(*np.transpose([start, end]), color=HIGHLIGHT["color"], alpha=HIGHLIGHT["alpha"],
                    lw=ARROW_HIGHLIGHT_PT, solid_capstyle="round", zorder=2)
        draw_arrow(ax, start, end, color, z=3, lw=1.4, head=8)
        if (source, target) in labels:
            xy = (p + q) / 2 + (PAIR_GAP + 0.2) * left
            ax.text(*xy, labels[source, target], fontsize=LABEL_PT + 2, color=color,
                    ha={1: "left", -1: "right", 0: "center"}[int(round(left[0]))],
                    va={1: "bottom", -1: "top", 0: "center"}[int(round(left[1]))])


def slider_y(c):
    return SLIDER_1X + np.log2(c) * SLIDER_DOUBLING


def draw_slider(ax, c, offset=(0.0, 0.0)):
    """Concentration slider set at `c` (relative to 1x)."""
    dx, dy = offset
    x = SLIDER_X + dx
    for end in SLIDER_ENDS:
        draw_arrow(ax, (x, SLIDER_1X + dy), (x, end + dy), SLIDER_COLOR, z=3, lw=1.5, head=8)
    for value, text in SLIDER_TICKS:
        y = slider_y(value) + dy
        ax.plot([x - 0.13, x + 0.13], [y, y], color=SLIDER_COLOR, lw=1.5, zorder=3)
        draw_label(ax, (x + 0.3, y), text, ha="left")
    ax.plot([x - 0.45, x + 0.15], [slider_y(c) + dy] * 2, color="black", lw=5,
            solid_capstyle="round", zorder=4)
    for text, y in SLIDER_TITLE:
        draw_label(ax, (x, y + dy), text)


@dataclass(frozen=True)
class TimeAxes:
    """Axes of a plot against time: two arrows from `origin` to the `tips`
    (x of the time axis' tip, y of the vertical axis' tip); the time window
    is drawn over the x range `span`. All in slide units."""
    origin: tuple
    tips: tuple
    span: tuple

    def moved(self, dy):
        (x0, y0), (x_tip, y_tip) = self.origin, self.tips
        return TimeAxes((x0, y0 + dy), (x_tip, y_tip + dy), self.span)

    def x(self, t, duration):
        return self.span[0] + np.asarray(t) / duration * (self.span[1] - self.span[0])

    def draw(self, ax, time_label=True):
        (x0, y0), (x_tip, y_tip) = self.origin, self.tips
        draw_arrow(ax, (x0, y0), (x_tip, y0), AXIS_COLOR, z=3, lw=1.3, head=9)
        draw_arrow(ax, (x0, y0), (x0, y_tip), AXIS_COLOR, z=3, lw=1.3, head=9)
        if time_label:
            draw_label(ax, (x_tip - 0.15, y0 - 0.5), "time")


# Time series, as placed in scene 11; its ON and OFF levels, above the time axis.
TIME_SERIES_AXES = TimeAxes(origin=(12.75, 5.96), tips=(18.05, 8.33), span=(12.82, 17.65))
LEVELS = {True: 1.47, False: 0.40}


def draw_time_series(ax, trace, duration, axes=TIME_SERIES_AXES, lw=1.3, highlight=None,
                     title=True, time_label=True, label=None):
    """Plot of a trace (ON/OFF against time) over `duration`, with its ON or
    OFF stretches highlighted if `highlight` is True or False; `label` is
    written on its right."""
    axes.draw(ax, time_label)
    (x0, y0), (x_tip, y_tip) = axes.origin, axes.tips
    levels = {s: y0 + h for s, h in LEVELS.items()}
    if title:
        draw_label(ax, (x0 + 0.27, y_tip), "transcription", ha="left")
    for state, text in ((True, "ON"), (False, "OFF")):
        draw_label(ax, (x0 - 0.2, levels[state]), text, ha="right")
    if label:
        draw_label(ax, (x_tip + 0.3, (levels[True] + levels[False]) / 2), label, ha="left",
                   size=LABEL_PT + 2)

    times, states = trace
    x = axes.x(times, duration)
    y = np.array([levels[s] for s in states])
    ax.plot(np.repeat(x, 2)[1:-1], np.repeat(y, 2), color="#222222", lw=lw,
            solid_joinstyle="miter", zorder=4)
    if highlight is not None:
        for i in np.flatnonzero(states == highlight):
            ax.plot(x[i:i + 2], [levels[highlight]] * 2, color=HIGHLIGHT["color"], lw=8,
                    alpha=HIGHLIGHT["alpha"], solid_capstyle="butt", zorder=2)


# ---------------------------------------------------------------------------
# The trace pictogram, for P(trace | c) in equations
# ---------------------------------------------------------------------------

# A little ON/OFF trace, as sketched in the storyboard: lengths of its OFF and
# ON stretches, alternating from OFF, in em (it scales with the text).
GLYPH_STRETCHES = [0.18, 0.26, 0.28, 0.72, 0.15, 0.14, 0.15, 0.14, 0.55, 0.3, 0.14]
GLYPH_HEIGHT = 0.62  # em; about the height of capitals
GLYPH_LINE = 0.05    # em
GLYPH_GAP = 0.14     # em, between the pictogram and the text around it

LIKELIHOOD_OPEN = "$P($"


def likelihood_close(which):
    return rf"$\mid c_{{\mathrm{{{which}}}}})$"


def text_width(text, size):
    """Width of `text` as set at `size` (both in points)."""
    return TextPath((0, 0), text, size=size).get_extents().width


def likelihood_width(which, size):
    """Width of P(trace | c_<which>) at `size`, in points."""
    return (text_width(LIKELIHOOD_OPEN, size) + text_width(likelihood_close(which), size)
            + (sum(GLYPH_STRETCHES) + 2 * GLYPH_GAP) * size)


def trace_glyph():
    """The pictogram as a polyline (em), sitting on the baseline at the origin."""
    x, level, points = 0.0, 0.0, [(0.0, 0.0)]
    for length in GLYPH_STRETCHES:
        x += length
        points.append((x, level))
        level = GLYPH_HEIGHT - level  # the switch
        points.append((x, level))
    return np.array(points[:-1])


def draw_likelihood(ax, place, start, which, size, rotation=0, color=TEXT_COLOR):
    """P(trace | c_<which>), the likelihood of the trace at concentration c_low
    or c_high (`which` = "low" or "high"). `place` maps points from a local frame
    to slide units (see placer), `start` is the left end of the baseline there,
    and `rotation` is the frame's angle. Returns the right end of the line."""
    x, y = start

    def write(text, x):
        ax.text(*place((x, y)), text, rotation=rotation, rotation_mode="anchor", ha="left",
                va="baseline", fontsize=size, color=color)
        return x + text_width(text, size)

    x = write(LIKELIHOOD_OPEN, x) + GLYPH_GAP * size
    glyph = trace_glyph() * size + (x, y)
    ax.plot(*np.transpose(place(glyph)), color=color, lw=GLYPH_LINE * size,
            solid_joinstyle="miter", solid_capstyle="butt")
    return write(likelihood_close(which), glyph[-1, 0] + GLYPH_GAP * size)


def draw_likelihood_ratio(ax, center, size=LABEL_PT, rotation=0, log=True):
    """The fraction P(trace | c_low) / P(trace | c_high), with "log" in front if
    `log`, centred at `center` and turned by `rotation` degrees."""
    place = placer(center, 1 / POINTS_PER_UNIT, rotation)
    log_width = text_width("log", size) + 0.3 * size if log else 0.0
    widths = {which: likelihood_width(which, size) for which in ("low", "high")}
    bar_width = max(widths.values()) + 0.2 * size
    left = -(log_width + bar_width) / 2
    bar_start = left + log_width
    if log:
        ax.text(*place((left, -0.33 * size)), "log", rotation=rotation, rotation_mode="anchor",
                ha="left", va="baseline", fontsize=size, color=TEXT_COLOR)
    ax.plot(*np.transpose(place([(bar_start, 0.0), (bar_start + bar_width, 0.0)])),
            color=TEXT_COLOR, lw=0.04 * size)
    for which, y in (("low", 0.33 * size), ("high", -1.0 * size)):
        draw_likelihood(ax, place, (bar_start + (bar_width - widths[which]) / 2, y), which, size,
                        rotation)
