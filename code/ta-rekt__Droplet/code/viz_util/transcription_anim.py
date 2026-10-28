#!/usr/bin/env python3
"""
Interactive telegraph-process animation.

Top subplot   : oscilloscope-style ON/OFF signal. A fixed-width time window that
                scrolls with the current time (newest data enters at the right
                edge). The window width doubles as the zoom / flow-speed knob:
                a live slider in the interactive view, or the WINDOW_WIDTH
                constant when rendering a GIF. Wider window  = zoomed out =
                signal appears to flow more slowly across the screen.

Bottom subplot: cumulative mRNA, produced at a constant rate while the gene is
                ON (flat while OFF). Its axes start at a fixed width and rescale
                outward whenever the curve reaches an edge, so the whole signal
                from t=0 always stays in view.

ON  durations ~ Exp(rate = K_OFF)   (mean = 1 / K_OFF)
OFF durations ~ Exp(rate = K_ON)    (mean = 1 / K_ON)

Run:
    python transcription_anim.py            # interactive window (default)
    python transcription_anim.py gif        # render a GIF to GIF_PATH
    python transcription_anim.py interactive
"""

import sys
import numpy as np
import matplotlib

# ---------------------------------------------------------------------------
# Config  (edit these)
# ---------------------------------------------------------------------------
RENDER        = "interactive"   # "interactive" or "gif"; a CLI arg overrides
GIF_PATH      = "telegraph.gif"

SEED          = 0

K_ON          = float(sys.argv[1]) if len(sys.argv) > 1 else 1.0
K_OFF         = 1.0     # ON  -> OFF rate   (mean ON  duration = 1 / K_OFF)
PROD_RATE     = 1.0     # mRNA produced per unit time while ON

T_TOTAL       = 400.0    # total simulated time
SIM_DT        = 0.2     # simulated time advanced per animation frame
FPS           = 25      # animation frames per second

WINDOW_WIDTH  = 30.0    # top-plot visible time span (zoom / flow speed)
WINDOW_MIN    = 2.0     # slider range (interactive mode)
WINDOW_MAX    = T_TOTAL

BOTTOM_XWIN0  = 300.0    # bottom-plot initial (fixed) x-axis width
BOTTOM_YWIN0  = 150.0     # bottom-plot initial y-axis height
RESCALE_GROW  = 1.6     # bottom-plot axes grow by this factor on hitting an edge

# colors
C_SIGNAL = "#264653"    # ON/OFF trace
C_ON     = "#2a9d8f"    # ON-region fill
C_MRNA   = "#e76f51"    # mRNA curve

# ---------------------------------------------------------------------------

# if len(sys.argv) > 1 and sys.argv[1] in ("gif", "interactive"):
#     RENDER = sys.argv[1]
# if RENDER == "gif":
#     matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
if RENDER == "interactive":
    from matplotlib.widgets import Slider


# ---------------------------------------------------------------------------
# Simulate the telegraph process and the cumulative mRNA
# ---------------------------------------------------------------------------
def simulate_telegraph(seed):
    """Return (starts, states): segment-start times and the state (0=OFF,1=ON)
    held on [starts[i], starts[i+1]). `starts` has a closing vertex at T_TOTAL,
    so len(starts) == len(states) + 1."""
    rng = np.random.default_rng(seed)
    starts = [0.0]
    states = [0]                                  # start OFF
    t = 0.0
    while t < T_TOTAL:
        rate = K_ON if states[-1] == 0 else K_OFF
        t += rng.exponential(1.0 / rate)
        if t >= T_TOTAL:
            break
        starts.append(t)
        states.append(1 - states[-1])
    starts.append(T_TOTAL)                        # closing vertex
    return np.asarray(starts), np.asarray(states)


def mrna_at_vertices(starts, states):
    """Cumulative mRNA at each vertex time (piecewise linear between them)."""
    m = np.zeros(len(starts))
    for i in range(len(starts) - 1):
        rate = PROD_RATE if states[i] == 1 else 0.0
        m[i + 1] = m[i] + rate * (starts[i + 1] - starts[i])
    return m


starts, states = simulate_telegraph(SEED)
mrna = mrna_at_vertices(starts, states)
y_step = np.append(states, states[-1])           # aligned with `starts` for steps-post
on_mask = (y_step == 1)


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------
fig, (ax_top, ax_bot) = plt.subplots(
    2, 1, figsize=(9, 6.2), gridspec_kw={"height_ratios": [1, 1.2]}
)
fig.subplots_adjust(hspace=0.35, bottom=0.16 if RENDER == "interactive" else 0.1)

# --- top: oscilloscope ON/OFF ---
ax_top.step(starts, y_step, where="post", color=C_SIGNAL, lw=2, zorder=3)
ax_top.fill_between(starts, 0, 1, where=on_mask, step="post",
                    color=C_ON, alpha=0.25, zorder=1)
edge_dot, = ax_top.plot([], [], "o", color=C_SIGNAL, ms=7, zorder=4)
time_text = ax_top.text(0.015, 0.86, "", transform=ax_top.transAxes,
                        fontsize=10, color="#555")
ax_top.set_ylim(-0.18, 1.18)
ax_top.set_yticks([0, 1])
ax_top.set_yticklabels(["OFF", "ON"])
ax_top.set_xlabel("time")
ax_top.set_title(f"mRNA transcription", fontsize=11)
ax_top.grid(True, axis="x", alpha=0.2)

# --- bottom: cumulative mRNA (auto-rescaling axes) ---
mrna_line, = ax_bot.plot([], [], color=C_MRNA, lw=2)
mrna_dot,  = ax_bot.plot([], [], "o", color=C_MRNA, ms=6)
ax_bot.set_xlim(0, BOTTOM_XWIN0)
ax_bot.set_ylim(0, BOTTOM_YWIN0)
ax_bot.set_xlabel("time")
ax_bot.set_ylabel("cumulative mRNA")
ax_bot.set_title("accumulated mRNA", fontsize=11)
# ax_bot.grid(True, alpha=0.2)

# --- speed/zoom slider (interactive only) ---
if RENDER == "interactive":
    slider_ax = fig.add_axes([0.15, 0.04, 0.7, 0.03])
    win_slider = Slider(slider_ax, "window\n(zoom)", WINDOW_MIN, WINDOW_MAX,
                        valinit=WINDOW_WIDTH)

view = {"xmax": BOTTOM_XWIN0, "ymax": BOTTOM_YWIN0}
n_frames = int(T_TOTAL / SIM_DT) + 1


def current_state(t_now):
    i = np.searchsorted(starts, t_now, side="right") - 1
    return states[min(max(i, 0), len(states) - 1)]


def update(frame):
    t_now = min(frame * SIM_DT, T_TOTAL)

    # top: scrolling window of width W (right edge = current time)
    W = win_slider.val if RENDER == "interactive" else WINDOW_WIDTH
    ax_top.set_xlim(t_now - W, t_now)
    edge_dot.set_data([t_now], [current_state(t_now)])
    time_text.set_text(f"t = {t_now:5.1f}")

    # bottom: draw the whole curve up to t_now, rescale axes at the edges
    mask = starts <= t_now
    xs = np.append(starts[mask], t_now)
    m_now = float(np.interp(t_now, starts, mrna))
    ys = np.append(mrna[mask], m_now)
    mrna_line.set_data(xs, ys)
    mrna_dot.set_data([t_now], [m_now])

    if t_now >= view["xmax"]:
        while view["xmax"] <= t_now:
            view["xmax"] *= RESCALE_GROW
        ax_bot.set_xlim(0, view["xmax"])
    if m_now >= view["ymax"]:
        while view["ymax"] <= m_now:
            view["ymax"] *= RESCALE_GROW
        ax_bot.set_ylim(0, view["ymax"])

    return mrna_line, mrna_dot, edge_dot, time_text


anim = FuncAnimation(fig, update, frames=n_frames,
                     interval=1000.0 / FPS, blit=False,
                     repeat=(RENDER == "interactive"))

if RENDER == "gif":
    anim.save(GIF_PATH, writer=PillowWriter(fps=FPS))
    print(f"wrote {GIF_PATH}  ({n_frames} frames)")
else:
    plt.show()
