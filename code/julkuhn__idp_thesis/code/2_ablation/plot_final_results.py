#!/usr/bin/env python3
"""Overview figure for the IDP binder design campaign."""
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch

OUT = "/capstor/store/cscs/swissai/ab011/jkuhn/idp_bindcraft/exp_panel_ablation/analysis/fig_FINAL_results_overview.png"
BLUE, BASEC, ORANGE = "#2a78d6", "#34475e", "#eb6834"
INK, SECOND, MUTED, GRID = "#0b0b0b", "#52514e", "#8a8a85", "#d9d7d1"
plt.rcParams.update({"figure.facecolor": "white", "axes.facecolor": "white", "font.family": "DejaVu Sans",
    "font.size": 10, "axes.edgecolor": SECOND, "axes.linewidth": 0.8, "axes.spines.top": False,
    "axes.spines.right": False, "xtick.color": SECOND, "ytick.color": SECOND})

fig = plt.figure(figsize=(14, 8.4))
gs = fig.add_gridspec(2, 2, height_ratios=[1.35, 1], width_ratios=[1.15, 1],
                      hspace=0.42, wspace=0.24, left=0.07, right=0.975, top=0.88, bottom=0.07)

# ---------- Panel A: acc/800 across scales ----------
axA = fig.add_subplot(gs[0, 0])
scales = ["Panel-5\n(n=5)", "16-target\n(n=16)", "scale-100\n(n=100, paired)"]
full = [141, 153, 117]; actif = [368, 341, 265]
x = np.arange(len(scales)); w = 0.38
axA.bar(x - w/2, full, w, color=BASEC, edgecolor='white', lw=0.8, label='IDP-BC (full)', zorder=3)
axA.bar(x + w/2, actif, w, color=BLUE, edgecolor='white', lw=0.8, label='actif λ=1', zorder=3)
for xi, (f, a) in zip(x, zip(full, actif)):
    axA.text(xi - w/2, f + 6, str(f), ha='center', fontsize=8, color=MUTED)
    axA.text(xi + w/2, a + 6, str(a), ha='center', fontsize=8, color=INK, fontweight='bold')
    axA.text(xi, max(f, a) + 34, f"{a/f:.1f}×", ha='center', fontsize=10, color=BLUE, fontweight='bold')
axA.set_xticks(x); axA.set_xticklabels(scales, fontsize=8.5)
axA.set_ylabel("accepted designs / 800 trajectories"); axA.set_ylim(0, 430)
axA.set_title("A   Throughput: actif ≈ 2× across every scale", loc='left', fontweight='bold', color=INK, fontsize=11)
axA.legend(frameon=False, fontsize=9, loc='upper right')
axA.grid(axis='y', ls=':', color=GRID, lw=0.7); axA.set_axisbelow(True); axA.tick_params(length=0)

# ---------- Panel B: scale-100 paired metrics ----------
axB = fig.add_subplot(gs[0, 1])
metrics = ["acc / 800", "designs /\nnode-h", "i_pTM_iface\n(×100, ↑)", "100 − i_pAE_iface\n(×100, ↑)"]
fv = [117, 7.5, 38.0, 73.2]; av = [265, 17.8, 41.3, 77.5]   # i_pAE flipped so higher=better for visual
x = np.arange(len(metrics))
axB.bar(x - w/2, fv, w, color=BASEC, edgecolor='white', lw=0.8, zorder=3)
axB.bar(x + w/2, av, w, color=BLUE, edgecolor='white', lw=0.8, zorder=3)
pv = ["p<10⁻⁷", "≈2.4×", "p=7×10⁻³", "p=6×10⁻¹²"]
for xi, (f, a, p) in zip(x, zip(fv, av, pv)):
    axB.text(xi + w/2, a + max(av)*0.015, f"{a:g}", ha='center', fontsize=7.5, color=INK, fontweight='bold')
    axB.text(xi - w/2, f + max(av)*0.015, f"{f:g}", ha='center', fontsize=7.5, color=MUTED)
    axB.text(xi, max(f, a) + max(av)*0.07, p, ha='center', fontsize=7.4, color=BLUE)
axB.set_xticks(x); axB.set_xticklabels(metrics, fontsize=7.8)
axB.set_title("B   scale-100 paired test: actif wins on every axis", loc='left', fontweight='bold', color=INK, fontsize=11)
axB.set_ylim(0, max(av)*1.2); axB.grid(axis='y', ls=':', color=GRID, lw=0.7); axB.set_axisbelow(True); axB.tick_params(length=0)
axB.set_ylabel("value")

# ---------- Panel C: coverage + binders in hand ----------
axC = fig.add_subplot(gs[1, 0])
axC.axis('off')
rows = [("", "IDP-BC", "actif"),
        ("binders (100 targets)", "729", "1467"),
        ("mean / target", "7.3", "14.7"),
        ("genuine 0-design targets", "10", "2"),
        ("designs per node-hour", "7.5", "17.8")]
yt = np.linspace(0.9, 0.18, len(rows))
for (lab, f, a), y in zip(rows, yt):
    hd = lab == ""
    axC.text(0.0, y, lab, fontsize=10, color=INK if hd else SECOND, fontweight='bold' if hd else 'normal', transform=axC.transAxes)
    axC.text(0.62, y, f, fontsize=10, ha='center', color=BASEC, fontweight='bold', transform=axC.transAxes)
    axC.text(0.86, y, a, fontsize=10, ha='center', color=BLUE, fontweight='bold', transform=axC.transAxes)
    if hd: axC.plot([0, 0.95], [y-0.05, y-0.05], color=GRID, lw=0.8, transform=axC.transAxes)
axC.text(0.0, 0.0, "scale-100 (BC1-fork) · actif tail-crop top-up still finishing", fontsize=7.5, color=MUTED, transform=axC.transAxes)
axC.set_title("C   Coverage & yield (scale-100)", loc='left', fontweight='bold', color=INK, fontsize=11)

# ---------- Panel D: conclusions ----------
axD = fig.add_subplot(gs[1, 1]); axD.axis('off')
axD.set_title("D   Conclusions", loc='left', fontweight='bold', color=INK, fontsize=11)
lines = [
 ("● Production method:", INK, True),
 ("   BC1-fork + actif λ=1  (ipTM loss → actif-pTM-energy)", BLUE, False),
 ("   ≈2× throughput + better interface quality, same coverage.", SECOND, False),
 ("● One controlled variable:", INK, True),
 ("   only the confidence loss changes; IDP stack + filters identical.", SECOND, False),
 ("● BC2 engine:", INK, True),
 ("   dead end for actif (degenerate no-contact minimiser);", SECOND, False),
 ("   best BC2 = V8b, but can't engage G3BP1/MAPT. Uncapped rerun in progress.", SECOND, False),
 ("● Open test:", INK, True),
 ("   actif-loss + ipSAE filter (running) — quality vs throughput.", SECOND, False),
]
y = 0.9
for txt, c, b in lines:
    axD.text(0.0, y, txt, fontsize=8.6 if not b else 9.2, color=c, fontweight='bold' if b else 'normal', transform=axD.transAxes)
    y -= 0.105

fig.suptitle("IDP binder design — final results:  actif λ=1 beats the IDP-BC baseline",
             x=0.07, ha='left', fontsize=14.5, fontweight='bold', y=0.955)
plt.savefig(OUT, dpi=200); print("wrote", OUT)
