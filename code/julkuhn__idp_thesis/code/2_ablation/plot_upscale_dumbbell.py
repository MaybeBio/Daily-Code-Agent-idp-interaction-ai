#!/usr/bin/env python3
"""Per-crop paired comparison of the 100-crop runs."""
import glob, os, csv, statistics as st
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

OUT = "/capstor/store/cscs/swissai/ab011/jkuhn/idp_bindcraft/exp_panel_ablation/analysis/fig_upscale_dumbbell.png"
OUTDIR = "/capstor/store/cscs/swissai/ab011/jkuhn/idp_bindcraft/upscale_500/output"
BLUE, BASEC, WINC, LOSEC = "#2a78d6", "#34475e", "#1baf7a", "#d9d7d1"
INK, SECOND, MUTED, GRID = "#0b0b0b", "#52514e", "#8a8a85", "#d9d7d1"
plt.rcParams.update({"figure.facecolor": "white", "axes.facecolor": "white", "font.family": "DejaVu Sans",
    "font.size": 10, "axes.edgecolor": SECOND, "axes.linewidth": 0.8, "axes.spines.top": False,
    "axes.spines.right": False, "xtick.color": SECOND, "ytick.color": SECOND})

COLS = {'i_pTM': ('Average_i_pTM', True), 'i_pAE': ('Average_i_pAE', False),
        'i_pTM_iface': ('Average_i_pTM_iface', True), 'i_pAE_iface': ('Average_i_pAE_iface', False)}

def crop_means(d):
    fd = f"{d}/final_design_stats.csv"
    if not os.path.isfile(fd): return None
    rows = list(csv.DictReader(open(fd)))
    out = {}
    for m, (c, _) in COLS.items():
        v = []
        for r in rows:
            try: v.append(float(r[c]))
            except: pass
        out[m] = st.mean(v) if v else None
    return out

full = {os.path.basename(d)[:-len("__full")]: d for d in glob.glob(f"{OUTDIR}/*__full")}
data = {}   # crop -> {'actif':{...}, 'full':{...}}
for c in sorted(full):
    a = crop_means(f"{OUTDIR}/{c}__actif"); f = crop_means(f"{OUTDIR}/{c}__full")
    if a and f: data[c] = {'actif': a, 'full': f}
crops = list(data)
n = len(crops)

fig, axes = plt.subplots(2, 2, figsize=(14, 8.5))
for ax, (m, (col, hi)) in zip(axes.ravel(), COLS.items()):
    # keep crops where both arms have this metric, sort by IDP(full) value
    pts = [(c, data[c]['full'][m], data[c]['actif'][m]) for c in crops
           if data[c]['full'][m] is not None and data[c]['actif'][m] is not None]
    pts.sort(key=lambda t: t[1])
    xs = np.arange(len(pts))
    fv = [p[1] for p in pts]; av = [p[2] for p in pts]
    wins = 0
    for x, f, a in zip(xs, fv, av):
        better = (a > f) if hi else (a < f)
        wins += better
        ax.plot([x, x], [f, a], color=(WINC if better else LOSEC), lw=1.0, zorder=2, alpha=0.9)
    ax.scatter(xs, fv, s=14, color=BASEC, zorder=3, label='IDP (full)', edgecolor='white', lw=0.3)
    ax.scatter(xs, av, s=14, color=BLUE, zorder=4, label='actif', edgecolor='white', lw=0.3)
    # mean reference lines
    ax.axhline(st.mean(fv), color=BASEC, ls=(0, (4, 3)), lw=1.0, zorder=1)
    ax.axhline(st.mean(av), color=BLUE, ls=(0, (4, 3)), lw=1.0, zorder=1)
    arrow = '↑ higher=better' if hi else '↓ lower=better'
    ax.set_title(f"{m}   ({arrow})", loc='left', fontweight='bold', color=INK, fontsize=11, pad=6)
    ax.set_ylabel(m, fontsize=9.5)
    ax.set_xlabel(f"crop (n={len(pts)}, sorted by IDP value)", fontsize=8.5)
    ax.text(0.015, 0.97, f"actif better on {wins}/{len(pts)} crops\nmeans: actif {st.mean(av):.3f}  ·  IDP {st.mean(fv):.3f}",
            transform=ax.transAxes, va='top', fontsize=8.5, color=SECOND,
            bbox=dict(boxstyle='round,pad=0.3', fc='white', ec=GRID, lw=0.6))
    ax.grid(axis='y', ls=':', color=GRID, lw=0.7); ax.set_axisbelow(True); ax.tick_params(length=0)
    ax.set_xticks([])

handles = [Line2D([0],[0], marker='o', color='none', markerfacecolor=BASEC, markersize=7, label='IDP (full)'),
           Line2D([0],[0], marker='o', color='none', markerfacecolor=BLUE, markersize=7, label='actif'),
           Line2D([0],[0], color=WINC, lw=2.5, label='actif better (per crop)'),
           Line2D([0],[0], color=LOSEC, lw=2.5, label='IDP better (per crop)')]
fig.legend(handles=handles, loc='upper right', frameon=False, fontsize=9, ncol=4,
           bbox_to_anchor=(0.995, 0.945))
fig.suptitle("Upscale production runs — per-crop paired quality (100 crops run both ways)",
             x=0.008, ha='left', fontsize=13.5, fontweight='bold', y=0.985)
fig.text(0.008, 0.005, f"n={len(crops)} paired crops; each vertical pair = one crop's per-crop mean over its accepted designs. "
         "Dashed lines = arm means. Crops sorted by the IDP value per panel, so the actif dot's offset shows the per-crop shift. "
         "Runs are live (trajectory-matched).", fontsize=7.4, color=MUTED)
plt.subplots_adjust(top=0.88, bottom=0.07, left=0.05, right=0.99, hspace=0.26, wspace=0.14)
plt.savefig(OUT, dpi=200); print("wrote", OUT, "| crops:", n)
