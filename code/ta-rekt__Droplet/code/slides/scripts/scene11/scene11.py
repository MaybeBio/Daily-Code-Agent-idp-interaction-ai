"""
Scene 11 -- "Here it is in action. Let's consider for now that transcription
gets activated as soon as a single green molecule is bound [...] How long would
you need to observe this to detect a 1.5x change in concentration? 1.2x?"

The locus with its promoter, TFs and free RNA polymerases (pink) scattered in
3D; the two-state model of transcription (OFF -> ON at rate c kb,
TF binding; ON -> OFF at rate ku, TF unbinding); a time series of transcription
simulated from that model; and, from 11_3 on, a concentration slider.

    scene11_1  OFF: nothing on the promoter or gene; the OFF state and every OFF
               stretch of the time series highlighted
    scene11_2  ON: a TF on the promoter, two polymerases elongating; the ON
               state and every ON stretch highlighted
    scene11_3  slider at 1/2x: half the TFs gone; OFF times visibly longer
    scene11_4  2x: twice the TFs; OFF times shorter, more bursts
    scene11_5  two time series over ZOOM_OUT x the time, in plots of the same
               shape: 1.5x on top, 1x below (the trace of 11_1-2 continued),
               each labelled on the right; the slider and TFs show 1.5x
    scene11_6  the same, re-sampled: 1.2x on top, a new 1x trace below, and no
               labels: which is which? (the slider and TFs show 1.2x)

Time series are exact simulations of the two-state Markov process (see
common/telegraph.py); each is the first random seed (from a fixed starting
seed) meeting the criteria noted in make_traces(), so they are reproducible
and easy to re-draw.

Outputs: slides/figures/scene11/scene11_<k>.{pdf,svg}, full 128 x 96 mm slides.

Usage (from slides/):
    python3 scripts/scene11/scene11.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # slides/scripts, for `common`

import numpy as np

from common.drawing import new_slide, save
from common.locus import Locus
from common.molecules import (TF_RADIUS, circle_points, depth_alpha, depth_scale, draw_legend,
                              draw_rnap_at, draw_scattered_tfs, scatter_proteins)
from common.telegraph import (TIME_SERIES_AXES, draw_slider, draw_time_series,
                              draw_two_state_graph, first_trace, stats, truncated, typical)

# ---------------------------------------------------------------------------
# Simulated time series (the model is in common/telegraph.py)
# ---------------------------------------------------------------------------

WINDOW = 12.0   # time shown on the plot (same units throughout; ZOOM_OUT x in 11_5-6)
ZOOM_OUT = 4


def make_traces():
    # 1x over ZOOM_OUT x the window, typical over that long run; its first window
    # has at least 5 complete ON bursts and is ON about half the time.
    base_long = first_trace(
        1.0, ZOOM_OUT * WINDOW,
        lambda tr: (typical(1.0)(tr) and stats(truncated(tr, WINDOW))["complete"] >= 5
                    and typical(1.0, 0.05)(truncated(tr, WINDOW))), seed=1)
    base = truncated(base_long, WINDOW)
    b = stats(base)
    return {
        "base": base,            # 11_1-2: the first window of base_long
        "base_long": base_long,  # 11_5: the same, over ZOOM_OUT x the time
        # 1/2x: OFF times visibly longer than at 1x
        "half": first_trace(0.5, WINDOW, lambda tr: stats(tr)["mean_off"] >= 1.5 * b["mean_off"],
                            seed=300),
        # 2x: OFF times shorter, more bursts
        "double": first_trace(2.0, WINDOW,
                              lambda tr: (stats(tr)["mean_off"] <= 0.6 * b["mean_off"]
                                          and stats(tr)["bursts"] >= b["bursts"] + 2), seed=400),
        # 1.5x and 1.2x over ZOOM_OUT x the window, typical over that run
        "high": first_trace(1.5, ZOOM_OUT * WINDOW, typical(1.5), seed=600),
        "raised": first_trace(1.2, ZOOM_OUT * WINDOW, typical(1.2), seed=700),
        # 11_6: 1x again, a new trace (unlike 11_5's, which the audience has seen)
        "base_new": first_trace(1.0, ZOOM_OUT * WINDOW, typical(1.0), seed=800),
    }


# ---------------------------------------------------------------------------
# Layout (slide units, traced from the storyboard)
# ---------------------------------------------------------------------------

# Where the sketch copies it, moved 0.4 left to leave room before the plots.
LOCUS = Locus(offset=(-0.43, -0.57))
POLYMERASES = [0.95, 0.35]  # two polymerases, at different stages of elongation
TRANSCRIPT_SEED = 11
BOUND_SITES = [(-6.5, -1)]  # the TF on the promoter (nodes from the gene start, side)

N_TFS_1X = 24
N_FREE_POLYMERASES = 9
RNAP_RADIUS = 0.37
SCATTER_SEED = 1100
REGION = (3.0, 10.4, 4.2, 14.2)  # right of the slider, clear of the plots

# Stacked time series are copies of the top one, STACK_STEP lower each.
STACK_STEP = 2.9
TRACE_PT = {1: 1.3, ZOOM_OUT: 0.55}  # line width, by window (in units of WINDOW)

# Each slide: concentration (slider and TFs), time series from top to bottom
# as (trace, label on the right or None), molecular state, highlight, slider,
# time window (in units of WINDOW).
SLIDES = {
    "scene11_1": dict(c=1.0, plots=[("base", None)], on=False, highlight=True, slider=False),
    "scene11_2": dict(c=1.0, plots=[("base", None)], on=True, highlight=True, slider=False),
    "scene11_3": dict(c=0.5, plots=[("half", None)]),
    "scene11_4": dict(c=2.0, plots=[("double", None)]),
    "scene11_5": dict(c=1.5, plots=[("high", r"$1.5\times$"), ("base_long", r"$1\times$")],
                      window=ZOOM_OUT),
    "scene11_6": dict(c=1.2, plots=[("raised", None), ("base_new", None)], window=ZOOM_OUT),
}


# ---------------------------------------------------------------------------
# Molecules: one population of TFs, thinned or extended with concentration
# ---------------------------------------------------------------------------

def molecules():
    tfs_bound = LOCUS.bound_tfs(BOUND_SITES)
    states = LOCUS.transcribe(POLYMERASES, TRANSCRIPT_SEED,
                              [circle_points(c, TF_RADIUS) for c in tfs_bound])
    # Everything scattered keeps clear of the ON state, so nothing moves when
    # the polymerases and the bound TF come and go.
    avoid = LOCUS.obstacles(states, tfs_bound)
    rng = np.random.default_rng(SCATTER_SEED)
    rnaps = scatter_proteins(N_FREE_POLYMERASES, REGION, rng, avoid, radius=RNAP_RADIUS)
    avoid += [(q["center"], RNAP_RADIUS * depth_scale(q["depth"])) for q in rnaps]

    base = scatter_proteins(N_TFS_1X, REGION, rng, avoid)
    keep = sorted(rng.choice(len(base), len(base) // 2, replace=False))
    tfs = {1.0: base, 0.5: [base[i] for i in keep]}
    for c in (2.0, 1.5, 1.2):
        tfs[c] = scatter_proteins(round(c * N_TFS_1X), REGION, rng, avoid, existing=base)
    return dict(bound=tfs_bound, states=states, rnaps=rnaps, tfs=tfs)


def draw_cartoon(ax, mol, c, transcribing, bound_tfs):
    """The locus, with its polymerases if `transcribing` and TFs bound at
    `bound_tfs` (centres), among the free polymerases and the TFs at
    concentration c."""
    LOCUS.draw(ax, mol["states"] if transcribing else [], bound_tfs)
    for q in mol["rnaps"]:
        draw_rnap_at(ax, q["center"], q["angle"], scale=depth_scale(q["depth"]),
                     alpha=depth_alpha(q["depth"]), z=0.4 if q["depth"] < 0.5 else 5.5)
    draw_scattered_tfs(ax, mol["tfs"][c], z_back=0.5, z_front=6)
    draw_legend(ax, ["DNA", "Promoter", "Gene", "TF", "RNAP"] + ["mRNA"] * transcribing)


def draw_slide(name, spec, traces, mol):
    c, on = spec["c"], spec.get("on", True)
    highlight = on if spec.get("highlight") else None
    fig, ax = new_slide()

    draw_cartoon(ax, mol, c, transcribing=on, bound_tfs=mol["bound"] if on else [])

    if spec.get("slider", True):
        draw_slider(ax, c)
    draw_two_state_graph(ax, highlight={True: "ON", False: "OFF", None: None}[highlight])
    plots, window = spec["plots"], spec.get("window", 1)
    for row, (trace, label) in enumerate(plots):
        draw_time_series(ax, traces[trace], window * WINDOW,
                         TIME_SERIES_AXES.moved(-row * STACK_STEP), lw=TRACE_PT[window],
                         highlight=highlight, title=row == 0,
                         time_label=row == len(plots) - 1, label=label)
    save(fig, "scene11", name)


if __name__ == "__main__":
    traces = make_traces()
    for key, tr in traces.items():
        s = stats(tr)
        print(f"{key:12s} bursts {s['bursts']:3d}  complete {s['complete']:3d}  "
              f"mean OFF {s['mean_off']:.2f}  ON fraction {s['on_fraction']:.2f}")
    mol = molecules()
    for name, spec in SLIDES.items():
        draw_slide(name, spec, traces, mol)
