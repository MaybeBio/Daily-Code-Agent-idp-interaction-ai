"""
Scene 17 -- "This has an interesting consequence once we consider multiple
binding sites [...] at high concentration, cooperative binding is beneficial
but at low concentration, when binding is diffusion-limited, the best
promoters are not cooperative."

The cartoon of scene 11 with a second binding site on the promoter; the
six-state model (0, 1 or 2 TFs bound x transcription OFF / ON); and, from
17_5 on, Figure 3a of the SPRT paper (two binding sites), simplified as in
scene 16, with the optimal circuit at the concentration marked.

    scene17_1  nothing bound; top-left state highlighted
    scene17_2  one TF bound; top middle
    scene17_3  two TFs bound; top right
    scene17_4  two TFs bound, transcribing; bottom right
    scene17_5  the plot appears, low c marked: the optimal circuit there (the
               paper's left-hand inset) highlighted, the double-occupancy
               states greyed out; one TF bound, transcribing
    scene17_6  high c: the cooperative circuit (right-hand inset); two TFs
               bound, transcribing, more TFs around

The highlighted circuits are the arrows drawn black or red in the paper's
insets. When the script runs, it checks them against the data: each carries
at least CIRCUIT_FLUX of the largest steady-state flux of the optimal circuit
at that concentration, and the other arrows carrying more than
REPORTED_FLUX are printed.

Data: fig3a_full.json, the "full" curve of Figure 3a: ConstrainedPromoters/
data/fig3/a1 (commit dea9418), merged into one file. One entry per c0: the
11 log10 rates of the optimal circuit (in the order of PAPER_EDGES), its
Fisher rate, ..., c0 (entry 16).

Outputs: slides/figures/scene17/scene17_<k>.{pdf,svg}, full 128 x 96 mm slides.

Usage (from slides/):
    python3 scripts/scene17/scene17.py
"""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # slides/scripts, for `common`

import numpy as np
from matplotlib.patches import Circle

from common.drawing import TEXT_COLOR, draw_label, new_slide, save
from common.telegraph import STATE_LABEL_X, STATE_LABELS, draw_state_graph
from scene11.scene11 import LOCUS, draw_cartoon, molecules
from scene16.scene16 import FOUR_AD, LIMIT_STYLE, Y_LABEL, label_along, paper_axes

# ---------------------------------------------------------------------------
# The paper's data
# ---------------------------------------------------------------------------

DATA = Path(__file__).with_name("fig3a_full.json")
# The paper's two-site model ("2bs_no_konu" in ConstrainedPromoters/spec/
# circuits.json): its states, as named here, and its edges, the binding ones
# (rate proportional to c) marked.
PAPER_STATES = {5: "TL", 1: "TM", 3: "TR", 6: "BL", 2: "BM", 4: "BR"}
PAPER_EDGES = [(5, 1), (1, 5), (2, 6), (6, 2), (1, 2), (1, 3), (2, 4), (3, 1), (3, 4), (4, 2),
               (6, 5)]
PAPER_BINDING = {(5, 1), (1, 3), (6, 2), (2, 4)}
FISHER_ENTRY, C0_ENTRY = 11, 16
ACTIVATION_LIMIT = 0.32 * 5  # as plotted in the paper (0.3 kmax)
DIFFUSION_LIMIT = 2 * FOUR_AD  # two sites


def paper_rows():
    rows = json.loads(DATA.read_text())
    return [rows[str(i)] for i in range(1, len(rows) + 1)]


def paper_curve():
    rows = paper_rows()
    return np.array([r[C0_ENTRY] for r in rows]), np.array([r[FISHER_ENTRY] for r in rows])


def edge_fluxes(row):
    """Steady-state probability flux along each edge of the optimal circuit."""
    c = row[C0_ENTRY]
    rates = [10 ** w * (c if e in PAPER_BINDING else 1) for e, w in zip(PAPER_EDGES, row)]
    generator = np.zeros((6, 6))
    for (u, v), rate in zip(PAPER_EDGES, rates):
        generator[u - 1, v - 1] += rate
    generator -= np.diag(generator.sum(axis=1))
    values, vectors = np.linalg.eig(generator.T)
    occupancy = np.real(vectors[:, np.argmin(np.abs(values))])
    occupancy /= occupancy.sum()
    return {(PAPER_STATES[u], PAPER_STATES[v]): occupancy[u - 1] * rate
            for (u, v), rate in zip(PAPER_EDGES, rates)}


# ---------------------------------------------------------------------------
# Slides
# ---------------------------------------------------------------------------

# Rows of the data (from 1) marked as low and high concentration.
LOW_ROW, HIGH_ROW = 9, 26  # c0 = 1.08 and 168

# The circuits of the paper's insets (arrows drawn black or red there).
CIRCUITS = {
    "low": [("TL", "TM"), ("TM", "BM"), ("BM", "BL"), ("BL", "TL")],
    "high": [("TL", "TM"), ("TM", "TL"), ("TM", "TR"), ("TR", "BR"), ("BR", "BM"),
             ("BM", "BL"), ("BL", "TL")],
}
CIRCUIT_FLUX, REPORTED_FLUX = 0.4, 0.1  # relative to the largest flux

# Cartoon: TFs bound (first, then second site), transcription; highlighted
# state; or, from 17_5, the circuit and the data row marked.
SLIDES = {
    "scene17_1": dict(bound=0, transcribing=False, state="TL"),
    "scene17_2": dict(bound=1, transcribing=False, state="TM"),
    "scene17_3": dict(bound=2, transcribing=False, state="TR"),
    "scene17_4": dict(bound=2, transcribing=True, state="BR"),
    "scene17_5": dict(bound=1, transcribing=True, circuit="low", row=LOW_ROW),
    "scene17_6": dict(bound=2, transcribing=True, circuit="high", row=HIGH_ROW, tfs=2.0),
}

# ---------------------------------------------------------------------------
# Layout (slide units)
# ---------------------------------------------------------------------------

BOUND_SITES = [(-6.5, -1), (-11.5, -1)]  # scene 11's site, then one further upstream

# Six-state graph: the two-state graph's rows, three states each.
COLUMNS = (11.5, 13.5, 15.5)
ROWS = {"T": 13.22, "B": 10.10}
STATES = {row + col: (x, y) for row, y in ROWS.items() for col, x in zip("LMR", COLUMNS)}
NEIGHBOURS = [("TL", "TM"), ("TM", "TR"), ("BL", "BM"), ("BM", "BR"),
              ("TL", "BL"), ("TM", "BM"), ("TR", "BR")]
ARROWS = [arrow for a, b in NEIGHBOURS for arrow in ((a, b), (b, a))]
RATE_LABELS = {("TL", "TM"): r"$c\,k_b$", ("TM", "TR"): r"$c\,k_b$",
               ("BR", "BM"): r"$k_u$", ("BM", "BL"): r"$k_u$",
               ("TR", "BR"): r"$k_{\rm on}$", ("BR", "TR"): r"$k_{\rm off}$"}
STATE_LABELS_DX = 1.5  # state labels, moved right of the two-state graph's
DOUBLE_OCCUPANCY = ("TR", "BR")

# Figure 3a, below the graph.
PLOT_BOX = (12.6, 2.0, 6.6, 5.4)
X_LIMITS, Y_LIMITS = (1e-1, 1e3), (1e-2, 3e1)  # as in the paper
PLOT_LABEL_PT, TICK_PT = 8.4, 7.6
# Both labelled on the left, the diffusion limit low enough to end below the
# activation limit.
DIFFUSION_LABEL_C, ACTIVATION_LABEL_C = 0.15, 0.13
MARK_LABELS = {LOW_ROW: dict(text="low $c$", offset=(6, -3), ha="left"),
               HIGH_ROW: dict(text="high $c$", offset=(0, -7), ha="center")}


def check_circuits():
    rows = paper_rows()
    for name, row in (("low", LOW_ROW), ("high", HIGH_ROW)):
        flux = edge_fluxes(rows[row - 1])
        top = max(flux.values())
        print(f"{name} c = {rows[row - 1][C0_ENTRY]:.3g}:")
        for edge in CIRCUITS[name]:
            print(f"   {edge[0]}->{edge[1]}  flux {flux[edge] / top:.2f}")
            assert flux[edge] >= CIRCUIT_FLUX * top
        for edge, f in flux.items():
            if edge not in CIRCUITS[name] and f > REPORTED_FLUX * top:
                print(f"   also {edge[0]}->{edge[1]}  flux {f / top:.2f} (grey in the paper's inset)")


def draw_six_state_graph(ax, state=None, circuit=None):
    styles, faded = {}, ()
    if circuit:
        styles = {a: "highlight" if a in CIRCUITS[circuit] else "faint" for a in ARROWS}
        faded = DOUBLE_OCCUPANCY if circuit == "low" else ()
    draw_state_graph(ax, STATES, ARROWS, state, faded, RATE_LABELS, styles)
    for labels in STATE_LABELS.values():
        for text, y in labels:
            draw_label(ax, (STATE_LABEL_X + STATE_LABELS_DX, y), text)


def draw_paper_plot(fig, row):
    c, fisher = paper_curve()
    ax = paper_axes(fig, PLOT_BOX, X_LIMITS, Y_LIMITS, Y_LABEL, tick_pt=TICK_PT)
    ax.plot(c, fisher, color="black", lw=1.6, zorder=4)
    ax.plot(c, DIFFUSION_LIMIT * c, **LIMIT_STYLE)
    ax.axhline(ACTIVATION_LIMIT, linestyle="dotted", **LIMIT_STYLE)
    label_along(ax, DIFFUSION_LABEL_C, 1.6 * DIFFUSION_LIMIT * DIFFUSION_LABEL_C,
                "Diffusion limit", slope=1, size=PLOT_LABEL_PT)
    label_along(ax, ACTIVATION_LABEL_C, 1.2 * ACTIVATION_LIMIT, "Activation limit", slope=0,
                size=PLOT_LABEL_PT)
    mark = (c[row - 1], fisher[row - 1])
    ax.plot(*mark, "o", color="black", ms=5, zorder=5)
    label = MARK_LABELS[row]
    ax.annotate(label["text"], mark, xytext=label["offset"], textcoords="offset points",
                ha=label["ha"], va="top", fontsize=PLOT_LABEL_PT, color=TEXT_COLOR)


def draw_slide(name, mol, bound, transcribing, state=None, circuit=None, row=None, tfs=1.0):
    fig, ax = new_slide()
    sites = LOCUS.bound_tfs(BOUND_SITES)
    draw_cartoon(ax, mol, tfs, transcribing=transcribing, bound_tfs=sites[:bound])
    draw_six_state_graph(ax, state, circuit)
    if row:
        draw_paper_plot(fig, row)
    save(fig, "scene17", name)


if __name__ == "__main__":
    check_circuits()
    mol = molecules()
    for name, spec in SLIDES.items():
        draw_slide(name, mol, **spec)
