"""
Scene 26 -- "We can also compute the Fisher information rate gain from this
nonlinear dependence of local concentration on background concentration,
highlighted by the term in red. A single promoter can produce a strong peak
in Fisher rate. [...] If the sphere is attractive enough, it can even dampen
the signal across all concentrations, which is unproductive."

The cartoon of scene 25 on the left. Right:

    scene26_1  the Fisher rate of a free probe, and of the probe on the anchor
               (SI eq. B5): the free rate at the local concentration, times
               the surface gain (d log phi / d log c)^2, in red
    scene26_2  Fig. 4G: that Fisher rate against c/c_sat, varying h0, with the
               free probe's (dashed)

The plot is re-drawn in the slides' font from the solver output behind the
paper's figure, as figure_production/fig4/fig4.ipynb draws it (the critical
h* run, shown in Fig. 4E, is left out there too).

Outputs: slides/figures/scene26/scene26_<k>.{pdf,svg}, full 128 x 96 mm slides.

Usage (from slides/):
    python3 scripts/scene26/scene26.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # slides/scripts, for `common`

from common.drawing import TEXT_COLOR, draw_label, draw_underbrace, new_slide, save
from scene25.scene25 import (CURVE_PT, EQUATIONS_X, LEGEND_PT, PANELS, draw_cartoon, draw_jumps,
                             draw_title, equation, fisher_rate, load_runs, plot_axes, run_label)

SURFACE_RED = "#D93A45"  # the red of scenes 2-3

# ---------------------------------------------------------------------------
# Equations
# ---------------------------------------------------------------------------

FREE = (r"\dot{ℱ}_{\rm free}(c) = \dfrac{1}{\tau_r(c)}", 11.0)
# B5, in pieces so that the surface term can be red and each term get its brace.
FISHER_LHS = r"\dot{ℱ}(c) ="
PROBE_TERM = r"\dfrac{1}{\tau_r(\phi)}"
SURFACE_TERM = r"\left(\dfrac{d\,\log\phi}{d\,\log c}\right)^2"
ANCHORED_Y = 7.7
TERM_GAP = 0.12   # between the pieces (slide units)
BRACE_DROP = 0.42  # from the baseline to the top of the braces
HEADINGS = [("Free probe", 12.45), ("Probe on the anchor", 9.2)]
SYMBOLS = [(r"$\tau_r$", "mean cycle time of the probe"),
           (r"$\phi$", "TF volume fraction at the probe")]
SYMBOLS_Y, SYMBOLS_STEP, SYMBOLS_INDENT = 4.6, 0.6, 0.85


def draw_equations(ax):
    for text, y in HEADINGS:
        draw_label(ax, (EQUATIONS_X, y), text, ha="left", va="baseline")
    equation(ax, (EQUATIONS_X, FREE[1]), FREE[0])

    def extent(text):
        """The x range of a piece's ink, in slide units."""
        box = text.get_window_extent().transformed(ax.transData.inverted())
        return box.x0, box.x1

    # Placed by ink: each piece starts just after the one before.
    lhs = equation(ax, (EQUATIONS_X, ANCHORED_Y), FISHER_LHS)
    probe = equation(ax, (extent(lhs)[1] + TERM_GAP, ANCHORED_Y), PROBE_TERM)
    surface = equation(ax, (extent(probe)[1] + TERM_GAP, ANCHORED_Y), SURFACE_TERM)
    surface.set_color(SURFACE_RED)
    draw_underbrace(ax, extent(probe), ANCHORED_Y - BRACE_DROP, "probe")
    draw_underbrace(ax, extent(surface), ANCHORED_Y - BRACE_DROP, "surface", color=SURFACE_RED)
    for k, (symbol, meaning) in enumerate(SYMBOLS):
        y = SYMBOLS_Y - k * SYMBOLS_STEP
        draw_label(ax, (EQUATIONS_X, y), symbol, ha="left", va="baseline")
        draw_label(ax, (EQUATIONS_X + SYMBOLS_INDENT, y), meaning, ha="left", va="baseline")


# ---------------------------------------------------------------------------
# Fig. 4G
# ---------------------------------------------------------------------------

FISHER_LIMITS = (0.1, 1e4)


def draw_fisher_plot(fig, ax_slide):
    spec, runs = PANELS["E"], load_runs("E")
    ax = plot_axes(fig, r"Fisher info. rate $\dot{ℱ}(c)$ (bits/s)", log_y=True,
                   y_limits=FISHER_LIMITS)
    for k, (run, color) in enumerate(zip(runs, spec["colors"])):
        if k == spec["critical"]:
            continue
        rate, free = fisher_rate(run, spec["converged_only"])
        ax.plot(run["x"], rate, color=color, lw=CURVE_PT, zorder=4, label=run_label("E", run, k))
        draw_jumps(ax, run, rate, color)
    ax.plot(run["x"], free, "--", color="black", lw=1.2, zorder=3,
            label=r"$\dot{ℱ}_{\rm free}(c)$")
    ax.legend(loc="upper left", fontsize=LEGEND_PT, frameon=False, handlelength=1.6,
              title=rf"$R = {runs[0]['R_nm']:.0f}$ nm",
              title_fontsize=LEGEND_PT, alignment="left", labelcolor=TEXT_COLOR)
    draw_title(ax_slide, r"varying $h_0$")


SLIDES = {"scene26_1": draw_equations, "scene26_2": None}


if __name__ == "__main__":
    for name, draw_right in SLIDES.items():
        fig, ax = new_slide()
        draw_cartoon(fig, ax)
        if draw_right:
            draw_right(ax)
        else:
            draw_fisher_plot(fig, ax)
        save(fig, "scene26", name)
