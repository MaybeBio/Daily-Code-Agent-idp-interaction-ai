#!/usr/bin/env python3
"""Paired comparison of actif lam1 and idp ref on the 100 crops run with both."""
import glob, os, csv, math, statistics as st
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

OUT = "/capstor/store/cscs/swissai/ab011/jkuhn/idp_bindcraft/exp_panel_ablation/analysis/fig_upscale_actif_vs_full.png"
OUTDIR = "/capstor/store/cscs/swissai/ab011/jkuhn/idp_bindcraft/upscale_500/output"
BLUE, BASEC = "#2a78d6", "#34475e"          # actif = blue, full = slate (baseline)
INK, SECOND, MUTED, GRID = "#0b0b0b", "#52514e", "#8a8a85", "#d9d7d1"
plt.rcParams.update({"figure.facecolor": "white", "axes.facecolor": "white", "font.family": "DejaVu Sans",
    "font.size": 10, "axes.edgecolor": SECOND, "axes.linewidth": 0.8, "axes.spines.top": False,
    "axes.spines.right": False, "xtick.color": SECOND, "ytick.color": SECOND})
def se(v): return st.stdev(v)/math.sqrt(len(v)) if len(v) > 1 else 0.0

COLS = {'i_pTM': 'Average_i_pTM', 'i_pAE': 'Average_i_pAE',
        'i_pTM_iface': 'Average_i_pTM_iface', 'i_pAE_iface': 'Average_i_pAE_iface'}

def crop_stats(d):
    """per-crop: acc/800 and mean of each score over that crop's accepted designs."""
    fd, ts = f"{d}/final_design_stats.csv", f"{d}/trajectory_stats.csv"
    if not os.path.isfile(fd): return None
    rows = list(csv.DictReader(open(fd)))
    ntr = (sum(1 for _ in open(ts)) - 1) if os.path.isfile(ts) else 0
    a8 = len(rows)*800/ntr if ntr else None
    sc = {}
    for k, c in COLS.items():
        v = [float(r[c]) for r in rows if _f(r.get(c))]
        sc[k] = st.mean(v) if v else None
    return a8, sc

def _f(x):
    try: float(x); return True
    except: return False

full = {os.path.basename(d)[:-len("__full")]: d for d in glob.glob(f"{OUTDIR}/*__full")}
pairs = sorted(full)
# collect paired per-crop values
ACC = {'actif': [], 'full': []}
SCO = {m: {'actif': [], 'full': []} for m in COLS}
npair = 0
for c in pairs:
    sa = crop_stats(f"{OUTDIR}/{c}__actif"); sf = crop_stats(f"{OUTDIR}/{c}__full")
    if not sa or not sf: continue
    npair += 1
    for arm, s in (('actif', sa), ('full', sf)):
        if s[0] is not None: ACC[arm].append(s[0])
        for m in COLS:
            if s[1][m] is not None: SCO[m][arm].append(s[1][m])

def mv(vals): return (st.mean(vals) if vals else float('nan'), se(vals) if vals else 0.0)

PANELS = [  # (key, ylabel, higher_better, ymax)
    ('acc', "accepted / 800 traj", True, None),
    ('i_pTM', "i_pTM  (global ↑)", True, 0.95),
    ('i_pAE', "i_pAE  (global ↓)", False, None),
    ('i_pTM_iface', "i_pTM_iface  (↑)", True, 0.6),
    ('i_pAE_iface', "i_pAE_iface  (↓)", False, None),
]
arms = ['full', 'actif']; cols = {'full': BASEC, 'actif': BLUE}
labels = {'full': 'IDP (full)', 'actif': 'actif'}

fig, axes = plt.subplots(1, 5, figsize=(14, 4.4))
for ax, (key, ylab, hi, ymax) in zip(axes, PANELS):
    data = ACC if key == 'acc' else SCO[key]
    means = {a: mv(data[a]) for a in arms}
    x = np.arange(len(arms))
    bars = ax.bar(x, [means[a][0] for a in arms], color=[cols[a] for a in arms],
                  width=0.62, edgecolor='white', lw=0.8, zorder=3)
    ax.errorbar(x, [means[a][0] for a in arms], yerr=[means[a][1] for a in arms],
                fmt='none', ecolor=SECOND, elinewidth=0.9, capsize=3, zorder=4)
    top = max(means[a][0]+means[a][1] for a in arms)
    ax.set_ylim(0, ymax if ymax else top*1.22)
    for xi, a in zip(x, arms):
        m, e = means[a]
        fmt = f"{m:.0f}" if key == 'acc' else f"{m:.3f}"
        ax.text(xi, m + e + ax.get_ylim()[1]*0.02, fmt, ha='center', va='bottom',
                fontsize=8.5, color=INK, fontweight='bold')
    ax.set_xticks(x); ax.set_xticklabels([labels[a] for a in arms], fontsize=9)
    ax.set_ylabel(ylab, fontsize=9.5)
    ax.grid(axis='y', ls=':', color=GRID, lw=0.7); ax.set_axisbelow(True); ax.tick_params(length=0)
    # ratio / delta annotation
    ma, mf = means['actif'][0], means['full'][0]
    if key == 'acc':
        ax.set_title(f"{ma/mf:.1f}×", loc='right', color=BLUE, fontsize=11, fontweight='bold', pad=4)
    else:
        better = (ma > mf) if hi else (ma < mf)
        ax.set_title(("actif better" if better else "full better"), loc='right',
                     color=(BLUE if better else BASEC), fontsize=8.5, pad=4)

fig.suptitle("Upscale production runs — actif vs IDP(full), 100 crops run both ways (paired)",
             x=0.008, ha='left', fontsize=13, fontweight='bold', y=0.99)
leg = [Patch(fc=BASEC, label='IDP (full) — baseline'), Patch(fc=BLUE, label='actif')]
axes[0].legend(handles=leg, loc='upper left', frameon=False, fontsize=8.5)
fig.text(0.008, 0.01, f"n={npair} paired crops (same crop run both ways). Macro over crops, error = SE; "
         "acc/800 includes zero-yield crops, score panels average only crops with ≥1 accepted design. "
         "Runs are live (trajectory-matched: ~184 vs ~187 traj/crop), so counts still rising. "
         "↑ higher=better, ↓ lower=better.", fontsize=7.4, color=MUTED)
plt.subplots_adjust(top=0.86, bottom=0.17, left=0.055, right=0.99, wspace=0.42)
plt.savefig(OUT, dpi=200); print("wrote", OUT, "| paired crops:", npair)
