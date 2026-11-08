"""
Scene 25 -- "[...] Like good physicists, suppose we model a promoter as a
sphere. The standard way to describe mixtures of liquids at equilibrium is
Flory-Huggins theory. [...] Adding an attractive surface can create a local
phase separation while the rest of the mixture is homogeneous. [...] Here is
local concentration at the surface, as a function of concentration of a
protein in the mixture. Here is how it depends on the radius of the sphere,
and here dependence on how strongly it attracts transcription factors."

Left, on every slide (after Fig. 4A of the droplet paper): the anchor, a
stretch of chromatin folded into a sphere of radius R, the layer of TFs it
condenses around itself, and the probe on its surface, switching ON when a TF
binds it at the local concentration. Right:

    scene25_1  the Flory-Huggins potential (eq. 9, in units of kBT/v_sol), the
               Cahn-Hilliard equation with the anchor's boundary condition
               (eqs. 10-12), and the symbols
    scene25_2  Fig. 4E: surface volume fraction against c/c_sat, varying h0
    scene25_3  Fig. 4D: the same, varying R

The plots are re-drawn in the slides' font from the solver output behind the
paper's figure (reduced to fig4_curves.json by prepare_fig4_data.py), the way
figure_production/fig4/fig4.ipynb draws them: stable branch solid, metastable
branch dashed, open circles at the jump.

Outputs: slides/figures/scene25/scene25_<k>.{pdf,svg}, full 128 x 96 mm slides.

Usage (from slides/):
    python3 scripts/scene25/scene25.py
"""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # slides/scripts, for `common`

import numpy as np
from matplotlib.patches import Circle

from common.drawing import (AXIS_COLOR, LABEL_PT, SLIDE_UNITS, TEXT_COLOR, draw_arrow,
                            draw_label, new_slide, save, smooth_open, zoomed_axes)
from common.molecules import (DNA_BOND, DNA_LINE_PT, NUCLEOSOME_EXTENT, TF, TF_RADIUS, Polymer,
                              depth_scale, draw_dna, draw_legend, draw_tf, legend_area,
                              nucleosomes, scatter_proteins, tf_alpha)
from common.telegraph import draw_state_graph

# ---------------------------------------------------------------------------
# The paper's curves
# ---------------------------------------------------------------------------

CURVES = Path(__file__).with_name("fig4_curves.json")
V_TF = 67.0                 # nm^3, volume of a TF (Klf4)
K_PLUS = 4 * 3e-3 * 7.4     # um^3/s: 4 u D, probe size u = 3 nm, D = 7.4 um^2/s
K_D_FRACTION = 0.1          # K_D = 0.1 c_sat


def per_um3(phi):
    """Volume fraction -> TFs per um^3."""
    return np.asarray(phi) / V_TF * 1e9


def load_runs(panel):
    """The runs of a panel ("D" or "E"), as arrays."""
    runs = json.loads(CURVES.read_text())[panel]
    for run in runs:
        for key in ("phi_inf", "phi_surf", "success"):
            run[key] = np.array(run[key])
        run["branch"] = np.array(run["branch"])
        run["phi_other"] = np.array([np.nan if p is None else p for p in run["phi_other"]])
        run["x"] = run["phi_inf"] / run["binodal"]
        # The stable branch jumps from thin to thick film here (first order).
        run["jumps"] = np.flatnonzero((run["branch"][:-1] == "thin")
                                      & (run["branch"][1:] == "thick")) + 1
    return runs


def segments(run, converged_only):
    """Index ranges of the stable branch between jumps (D keeps only converged points)."""
    parts = np.split(np.arange(len(run["x"])), run["jumps"])
    return [s[run["success"][s]] for s in parts] if converged_only else parts


def fisher_rate(run, converged_only=False):
    """Fisher rate of the probe at the surface, (d log phi / d log c)^2 / tau_r(phi),
    each continuous branch differentiated separately; and that of a free probe."""
    k_minus = K_PLUS * K_D_FRACTION * per_um3(run["binodal"])
    tau = 1 / (K_PLUS * per_um3(run["phi_surf"])) + 1 / k_minus
    tau_free = 1 / (K_PLUS * per_um3(run["phi_inf"])) + 1 / k_minus
    rate = np.full(len(tau), np.nan)
    for s in segments(run, converged_only):
        gain = (np.gradient(run["phi_surf"][s], run["phi_inf"][s])
                * run["phi_inf"][s] / run["phi_surf"][s]) ** 2
        rate[s] = gain / tau[s]
    return rate, 1 / tau_free


# ---------------------------------------------------------------------------
# Cartoon (left half)
# ---------------------------------------------------------------------------

SPHERE_CENTER, SPHERE_RADIUS = np.array([4.6, 8.6]), 1.9
SPHERE_FILL = "#F3EFE6"
# The condensed layer: translucent discs stacked out to HALO_WIDTH beyond the
# surface, so the green deepens towards it.
HALO_WIDTH, HALO_LAYERS, HALO_ALPHA = 1.6, 14, 0.021

# Chromatin entering at the lower left, two loops folded inside the sphere,
# leaving at the top (traced in slide units, then smoothed), drawn at DNA_ZOOM
# of its usual size so its helix stays fine at this scale.
DNA_TRACE = [(0.5, 4.3), (1.4, 5.4), (2.4, 6.6), (3.2, 7.45), (3.6, 7.6), (3.75, 7.15),
             (4.3, 7.1), (4.65, 7.6), (4.45, 8.2), (3.9, 8.3), (3.7, 7.95), (4.2, 8.55),
             (4.9, 8.75), (5.65, 8.7), (6.0, 9.2), (5.75, 9.75), (5.1, 9.7), (4.95, 9.25),
             (5.25, 9.05), (5.55, 9.6), (5.6, 10.45), (6.4, 11.35), (7.8, 11.8), (9.6, 11.9)]

DNA_ZOOM = 0.6

PROBE_ANGLE = np.radians(-38)
PROBE_RADIUS = 0.11
R_ANGLE = np.radians(140)       # direction of the radius arrow

# The probe's two states, below the sphere, with a leader line from the probe.
MODEL_STATES = {"OFF": (4.3, 3.0), "ON": (7.0, 3.0)}
MODEL_LABELS = {("OFF", "ON"): r"$k^+c_{\rm local}$", ("ON", "OFF"): r"$k^-$"}
SITE_RADIUS = 0.13              # the probe inside each state: empty, or a TF bound
LEAD_END = (-0.4, 0.27)         # where the probe's leader line meets OFF (from its centre),
                                # its upper left, clear of the rate label

# TFs: a dense scatter thinned with distance r from the centre, keeping a
# fraction BACKGROUND + (1 - BACKGROUND) exp(-(r - R) / LAYER_DECAY).
TF_REGION = (0.3, 9.5, 1.0, 14.6)
TF_SCALE = 0.62                 # drawn smaller than in earlier scenes: there are many
TF_CANDIDATES, TF_SEED = 520, 2501
LAYER_DECAY, BACKGROUND = 0.7, 0.07
LEGEND, LEGEND_AT = ["DNA", "TF"], (1.0, 14.1)


def dna():
    """The chain, in the coordinates of the zoomed axes it is drawn on."""
    curve, _ = smooth_open(DNA_TRACE, n=600)
    return Polymer.from_trace([("coil", curve / DNA_ZOOM)], DNA_BOND)


def on_surface(angle, inset=0.0):
    return SPHERE_CENTER + (SPHERE_RADIUS - inset) * np.array([np.cos(angle), np.sin(angle)])


def scatter_tfs(chain):
    probe = on_surface(PROBE_ANGLE)
    lead = np.linspace(probe, np.array(MODEL_STATES["OFF"]) + LEAD_END, 30)
    model = np.array(list(MODEL_STATES.values()))
    model_box = np.array([(x, y) for x in np.linspace(model[0, 0] - 0.8, model[1, 0] + 0.8, 8)
                          for y in (model[0, 1] - 0.9, model[0, 1], model[0, 1] + 0.7)])
    avoid = [(SPHERE_CENTER, SPHERE_RADIUS + 0.05), (chain.nodes * DNA_ZOOM, 0.15),
             *[(site["center"] * DNA_ZOOM, NUCLEOSOME_EXTENT * DNA_ZOOM)
               for site in nucleosomes(chain.smoothed())],
             (lead, 0.08), (model_box, 0.4), legend_area(LEGEND, LEGEND_AT),
             (probe + np.array([[0.3, -0.15], [0.9, -0.15], [1.4, -0.15]]), 0.25)]  # "Probe"
    rng = np.random.default_rng(TF_SEED)
    candidates = scatter_proteins(TF_CANDIDATES, TF_REGION, rng, avoid,
                                  radius=TF_RADIUS * TF_SCALE)
    r = np.array([np.linalg.norm(q["center"] - SPHERE_CENTER) for q in candidates])
    keep = BACKGROUND + (1 - BACKGROUND) * np.exp(-(r - SPHERE_RADIUS) / LAYER_DECAY)
    return [q for q, p in zip(candidates, keep) if rng.uniform() < p]


def draw_cartoon(fig, ax):
    for k in range(HALO_LAYERS):
        ax.add_patch(Circle(SPHERE_CENTER, SPHERE_RADIUS + HALO_WIDTH * (1 - k / HALO_LAYERS),
                            fc=TF["fill"][:3], ec="none", alpha=HALO_ALPHA, zorder=0))
    ax.add_patch(Circle(SPHERE_CENTER, SPHERE_RADIUS, fc=SPHERE_FILL, ec=AXIS_COLOR, lw=1.3,
                        ls=(0, (4, 3)), zorder=1))
    chain = dna()
    draw_dna(zoomed_axes(fig, (0, 0), DNA_ZOOM), chain, z=1,
             line_pt=tuple(DNA_ZOOM * w for w in DNA_LINE_PT))
    for q in scatter_tfs(chain):
        draw_tf(ax, q["center"], scale=TF_SCALE * depth_scale(q["depth"]), angle=q["angle"],
                alpha=tf_alpha(q["depth"]), z=1.5 if q["depth"] < 0.5 else 6)

    # Radius.
    tip = on_surface(R_ANGLE, inset=0.05)
    draw_arrow(ax, SPHERE_CENTER, tip, AXIS_COLOR, z=4, lw=1.2, head=7)
    ax.add_patch(Circle(SPHERE_CENTER, 0.05, fc=AXIS_COLOR, ec="none", zorder=4))
    normal = np.array([np.cos(R_ANGLE - np.pi / 2), np.sin(R_ANGLE - np.pi / 2)])
    draw_label(ax, (SPHERE_CENTER + tip) / 2 + 0.28 * normal, "$R$")

    # Probe, and its two states.
    probe = on_surface(PROBE_ANGLE)
    off = np.array(MODEL_STATES["OFF"])
    ax.plot(*np.transpose([probe, off + LEAD_END]), color=AXIS_COLOR, lw=0.8, zorder=3)
    ax.add_patch(Circle(probe, PROBE_RADIUS, fc="black", ec="none", zorder=7))
    draw_label(ax, probe + (0.25, -0.15), "Probe", ha="left", va="top")
    draw_state_graph(ax, MODEL_STATES, list(MODEL_LABELS), labels=MODEL_LABELS)
    for name, center in MODEL_STATES.items():
        ax.add_patch(Circle(center, SITE_RADIUS, lw=1.0, zorder=4, ec=AXIS_COLOR, fc="white"))
        if name == "ON":
            draw_tf(ax, center, z=5)
        draw_label(ax, (center[0], center[1] - 0.68), name)
    draw_legend(ax, LEGEND, LEGEND_AT)


# ---------------------------------------------------------------------------
# Right half: equations, or a plot
# ---------------------------------------------------------------------------

EQUATION_PT = 10
EQUATIONS_X = 9.95
HEADINGS = [("Flory–Huggins potential", 13.55), ("Cahn–Hilliard equation", 10.1)]
POTENTIAL = [(r"f(\phi) = \dfrac{\phi}{\nu}\ln\phi + (1-\phi)\ln(1-\phi)", 12.45),
             (r"-\,\chi\phi^2 - \mu\phi", 11.2)]   # second line, after the "="
CAHN_HILLIARD = [(r"0 = \nabla^2\left[f'(\phi) - \kappa\nabla^2\phi\right]", 9.0),
                 (r"\kappa\,\partial_r\phi\,|_{r=R} = -h_0", 7.85)]
SYMBOLS = [(r"$\phi$", "TF volume fraction"),
           (r"$\nu$", "TF / solvent volume ratio"),
           (r"$\chi$", "TF–TF attraction"),
           (r"$\mu$", "chemical potential"),
           (r"$h_0$", "surface attraction"),
           (r"$\kappa$", "interface stiffness")]
SYMBOLS_Y, SYMBOLS_STEP, SYMBOLS_INDENT = 6.5, 0.6, 0.85

# Plots: axes box (left, bottom, width, height, slide units), title above.
PLOT_BOX = (11.2, 2.5, 8.2, 8.2)
PLOT_TITLE_XY = (15.3, 12.55)
TITLE_PT, LEGEND_PT, TICK_PT = 11, 8, 8.4
CURVE_PT = 1.6
PANELS = {
    "E": dict(colors=["orange", "red", "#2563EB", "gray", "green"], critical=3,
              title=r"varying $h_0$", converged_only=False),
    "D": dict(colors=["#7FFF5F", "magenta", "cyan"], critical=None,
              title=r"varying $R$", converged_only=True),
}
SURFACE_LIMITS = (0.0, 0.2)


def equation(ax, xy, tex):
    return ax.text(*xy, rf"${tex}$", ha="left", va="baseline", fontsize=EQUATION_PT,
                   color=TEXT_COLOR)


def draw_equations(ax):
    for text, y in HEADINGS:
        draw_label(ax, (EQUATIONS_X, y), text, ha="left", va="baseline")
    (first, y1), (second, y2) = POTENTIAL
    equation(ax, (EQUATIONS_X, y1), first)
    # The second line starts under the first term, after "f(phi) =" (placed by ink).
    lhs = equation(ax, (EQUATIONS_X, 0), r"f(\phi) =")
    indent = lhs.get_window_extent().transformed(ax.transData.inverted()).x1 + 0.2
    lhs.remove()
    equation(ax, (indent, y2), second)
    for tex, y in CAHN_HILLIARD:
        equation(ax, (EQUATIONS_X, y), tex)
    for k, (symbol, meaning) in enumerate(SYMBOLS):
        y = SYMBOLS_Y - k * SYMBOLS_STEP
        draw_label(ax, (EQUATIONS_X, y), symbol, ha="left", va="baseline")
        draw_label(ax, (EQUATIONS_X + SYMBOLS_INDENT, y), meaning, ha="left", va="baseline")


def plot_axes(fig, y_label, log_y=False, y_limits=SURFACE_LIMITS):
    """Axes on the right half, in the slides' font, the y label on top."""
    left, bottom, width, height = PLOT_BOX
    ax = fig.add_axes([left / SLIDE_UNITS[0], bottom / SLIDE_UNITS[1],
                       width / SLIDE_UNITS[0], height / SLIDE_UNITS[1]])
    if log_y:
        ax.set_yscale("log")
    ax.set_xlim(0, 1)
    ax.set_ylim(*y_limits)
    ax.set_xticks(np.linspace(0, 1, 6))
    ax.tick_params(labelsize=TICK_PT, colors=TEXT_COLOR)
    ax.set_xlabel(r"Concentration $c/c_{\rm sat}$", fontsize=LABEL_PT, color=TEXT_COLOR)
    ax.text(0, 1.03, y_label, transform=ax.transAxes, ha="left", va="bottom",
            fontsize=LABEL_PT, color=TEXT_COLOR)
    return ax


def run_label(panel, run, k):
    if panel == "D":
        return rf"$R = {run['R_nm']:.0f}$ nm"
    name = "h_*" if k == PANELS["E"]["critical"] else "h_0"
    units = r"\ (k_BT/{\rm nm}^2)" if k == 0 else ""  # on the first entry, as in the paper
    return rf"${name} = {run['h0']:.2f}{units}$"


def legend_title(panel, runs):
    if panel == "D":
        return rf"$h_0 = {runs[0]['h0']:.2f}\,k_BT/{{\rm nm}}^2$"
    return rf"$R = {runs[0]['R_nm']:.0f}$ nm"


def draw_jumps(ax, run, y, color):
    ax.plot(run["x"][np.r_[run["jumps"] - 1, run["jumps"]]],
            y[np.r_[run["jumps"] - 1, run["jumps"]]], "o", ms=5, mfc="white", mec=color,
            mew=1.2, zorder=5)


def draw_title(ax_slide, text):
    draw_label(ax_slide, PLOT_TITLE_XY, text, size=TITLE_PT)


def draw_surface_plot(fig, ax_slide, panel):
    """Fig. 4D or E: surface volume fraction against c / c_sat."""
    spec, runs = PANELS[panel], load_runs(panel)
    ax = plot_axes(fig, r"Surface vol. fraction $\phi_{\rm o}$")
    for k, (run, color) in enumerate(zip(runs, spec["colors"])):
        x, phi = run["x"], run["phi_surf"]
        dashed = k == spec["critical"]
        style, meta_style = ("--", ":") if dashed else ("-", "--")
        for j, s in enumerate(segments(run, spec["converged_only"])):
            ax.plot(x[s], phi[s], style, color=color, lw=CURVE_PT, zorder=4,
                    label=run_label(panel, run, k) if j == 0 else None)
        # Metastable continuations, anchored at the edges of the stable branch.
        ok = run["success"] if spec["converged_only"] else np.ones(len(x), bool)
        meta_thick = np.flatnonzero((run["branch"] == "thin") & ok)
        meta_thin = np.flatnonzero((run["branch"] == "thick") & ok)
        for J in run["jumps"]:
            ax.plot(x[np.r_[meta_thick, J]], np.r_[run["phi_other"][meta_thick], phi[J]],
                    meta_style, color=color, lw=0.9, zorder=3)
            ax.plot(x[np.r_[J - 1, meta_thin]], np.r_[phi[J - 1], run["phi_other"][meta_thin]],
                    meta_style, color=color, lw=0.9, zorder=3)
        draw_jumps(ax, run, phi, color)
    ax.legend(loc="upper left", fontsize=LEGEND_PT, frameon=False, handlelength=1.6,
              title=legend_title(panel, runs), title_fontsize=LEGEND_PT,
              alignment="left", labelcolor=TEXT_COLOR)
    draw_title(ax_slide, spec["title"])


SLIDES = {"scene25_1": None, "scene25_2": "E", "scene25_3": "D"}


if __name__ == "__main__":
    for name, panel in SLIDES.items():
        fig, ax = new_slide()
        draw_cartoon(fig, ax)
        if panel:
            draw_surface_plot(fig, ax, panel)
        else:
            draw_equations(ax)
        save(fig, "scene25", name)
