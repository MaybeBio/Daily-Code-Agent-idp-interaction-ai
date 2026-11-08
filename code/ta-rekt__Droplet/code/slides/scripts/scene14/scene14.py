"""
Scene 14 -- "This slope is just the derivative of the average of a log
likelihood ratio. We can approximate this with a quantity commonly used in
statistical inference, called the Fisher information rate. [...] In this
model, faster binding = faster decisions, up to the limit allowed by
diffusion. [...]"

Equations on the left, appearing one by one; on the right, the Fisher rate of
the two-state model, Fdot(c) = 1/tau_burst = c kb ku / (c kb + ku), saturating
at ku (dashed), and the mean decision time <T_dec> = (K/V) tanh(K/2), down to
its own floor (dashed).

    scene14_1  V = 1/2 (dc/c)^2 Fdot(c), Fdot labelled "Fisher information rate"
    scene14_2  + the definition of Fdot(c)
    scene14_3  + "= 1/tau_burst"; the plot of Fdot(c) and the note on c kb
    scene14_4  + <T_dec> = (K/V) tanh(K/2); the plot of <T_dec>

Outputs: slides/figures/scene14/scene14_<k>.{pdf,svg}, full 128 x 96 mm slides.

Usage (from slides/):
    python3 scripts/scene14/scene14.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # slides/scripts, for `common`

import numpy as np

from common.drawing import (AXIS_COLOR, TEXT_COLOR, draw_label, draw_plot_axes, draw_underbrace,
                            new_slide, save)
from common.telegraph import K_BIND, K_UNBIND, binding_rate

# ---------------------------------------------------------------------------
# Equations (left half), in display style
# ---------------------------------------------------------------------------

EQUATION_PT = 10
EQUATIONS_X = 0.7  # left end of the equations
DRIFT = r"V = \dfrac{1}{2}\left(\dfrac{\delta c}{c}\right)^2"
FISHER_RATE = r"\dot{ℱ}(c)"  # after DRIFT, with a brace under it
FISHER_LHS = FISHER_RATE + " ="
FISHER = (FISHER_LHS + r"-\lim_{t\rightarrow\infty} \dfrac{1}{t}"
          r"\left\langle\dfrac{\partial^2 \ln P(\mathcal{O}(t)\,|\,c)}{(\partial \ln c)^2}"
          r"\right\rangle_c")
BURST = r"= \dfrac{1}{\tau_{\rm burst}}"  # on the next line, under the "="
DECISION = r"\langle T_{\rm dec}\rangle = \dfrac{K}{V}\tanh\left(\dfrac{K}{2}\right)"
# Baselines (slide units).
DRIFT_Y = 12.6
BRACE_DROP = 0.3   # from the baseline to the top of the brace
FISHER_Y = 9.3
BURST_Y = 7.3
DECISION_Y = 4.4

NOTES = [r"$c\,k_b$: diffusion-limited arrival of molecules", r"$k_b = 4aD$"]
NOTES_Y = [1.55, 0.8]

# ---------------------------------------------------------------------------
# Plots (right half): c from 0 to C_MAX, drawn up to C_RIGHT
# ---------------------------------------------------------------------------

C_MAX = 4.0
FISHER_ORIGIN, FISHER_TIPS = (12.3, 8.9), (19.2, 13.5)
DECISION_ORIGIN, DECISION_TIPS = (12.3, 3.3), (19.2, 7.9)
C_RIGHT = 18.7
FISHER_UNIT = 3.5    # slide units per unit of Fdot (ku = 1)
# <T_dec> is proportional to 1/Fdot(c); in units of its floor (c -> infinity),
# it is 1 + ku / (c kb). It is drawn from where it equals DECISION_TOP.
DECISION_TOP = 6.0
DECISION_UNIT = 0.62  # slide units per unit of the floor

SLIDES = {
    "scene14_1": dict(fisher=False, burst=False, decision=False),
    "scene14_2": dict(fisher=True, burst=False, decision=False),
    "scene14_3": dict(fisher=True, burst=True, decision=False),
    "scene14_4": dict(fisher=True, burst=True, decision=True),
}


def equation(ax, xy, tex):
    return ax.text(*xy, rf"${tex}$", ha="left", va="baseline", fontsize=EQUATION_PT,
                   color=TEXT_COLOR)


def ink_x(ax, tex):
    """Left and right ends of the ink of an equation set at EQUATIONS_X."""
    text = equation(ax, (EQUATIONS_X, 0), tex)
    box = text.get_window_extent().transformed(ax.transData.inverted())
    text.remove()
    return box.x0, box.x1


def draw_drift(ax):
    """V = ..., then the Fisher information rate, with its name under a brace."""
    equation(ax, (EQUATIONS_X, DRIFT_Y), DRIFT)
    rate = equation(ax, (ink_x(ax, DRIFT)[1] + 0.1, DRIFT_Y), FISHER_RATE)
    box = rate.get_window_extent().transformed(ax.transData.inverted())
    draw_underbrace(ax, (box.x0, box.x1), DRIFT_Y - BRACE_DROP, "Fisher information rate")


def draw_burst(ax):
    """ "= 1/tau_burst" under the definition of Fdot, its "=" under the other.
    (Text is placed by its ink, so this is measured.)"""
    left, right = ink_x(ax, FISHER_LHS)
    equals_left, equals_right = ink_x(ax, "=")
    equation(ax, (right - (equals_right - equals_left), BURST_Y), BURST)


def plot_curve(ax, origin, c, y, unit, asymptote):
    """Curve y(c) on plot axes at `origin`, with its high-c asymptote dashed."""
    x0, y0 = origin
    ax.plot(x0 + c / C_MAX * (C_RIGHT - x0), y0 + y * unit, color="#222222", lw=1.6, zorder=4)
    ax.plot([x0, C_RIGHT], [y0 + asymptote * unit] * 2, color=AXIS_COLOR, lw=1.1,
            dashes=(4, 3), zorder=3)


def draw_fisher_plot(ax):
    draw_plot_axes(ax, FISHER_ORIGIN, FISHER_TIPS, y_label=rf"${FISHER_RATE}$",
                   x_label="$c$")
    c = np.linspace(0, C_MAX, 200)
    plot_curve(ax, FISHER_ORIGIN, c, binding_rate(c), FISHER_UNIT, K_UNBIND)


def draw_decision_plot(ax):
    draw_plot_axes(ax, DECISION_ORIGIN, DECISION_TIPS, y_label=r"$\langle T_{\rm dec}\rangle$",
                   x_label="$c$")
    c = np.linspace(K_UNBIND / (K_BIND * (DECISION_TOP - 1)), C_MAX, 200)
    plot_curve(ax, DECISION_ORIGIN, c, 1 + K_UNBIND / (c * K_BIND), DECISION_UNIT, 1.0)


def draw_slide(name, fisher, burst, decision):
    fig, ax = new_slide()
    draw_drift(ax)
    if fisher:
        equation(ax, (EQUATIONS_X, FISHER_Y), FISHER)
    if burst:
        draw_burst(ax)
        draw_fisher_plot(ax)
        for text, y in zip(NOTES, NOTES_Y):
            draw_label(ax, (10.0, y), text)
    if decision:
        equation(ax, (EQUATIONS_X, DECISION_Y), DECISION)
        draw_decision_plot(ax)
    save(fig, "scene14", name)


if __name__ == "__main__":
    for name, spec in SLIDES.items():
        draw_slide(name, **spec)
