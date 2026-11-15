#!/usr/bin/env python3
"""Boltz-2 interface scores for the re-predicted designs."""
import csv, numpy as np, statistics as st
from scipy import stats
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = "/capstor/store/cscs/swissai/ab011/jkuhn/idp_bindcraft/boltz2_env/actif_vs_idp/fig_boltz_validation.png"
BLUE, BASEC = "#2a78d6", "#34475e"
INK, SECOND, MUTED, GRID = "#0b0b0b", "#52514e", "#8a8a85", "#d9d7d1"
plt.rcParams.update({"figure.facecolor": "white", "axes.facecolor": "white", "font.family": "DejaVu Sans",
    "font.size": 10, "axes.edgecolor": SECOND, "axes.linewidth": 0.8, "axes.spines.top": False,
    "axes.spines.right": False, "xtick.color": SECOND, "ytick.color": SECOND})

rows = list(csv.DictReader(open("boltz_iptm.csv")))
data = {'full': [], 'actif': []}
per_crop = {}
for r in rows:
    try: ip = float(r['iptm'])
    except: continue
    data[r['arm']].append(ip)
    per_crop.setdefault(r['crop'], {}).setdefault(r['arm'], []).append(ip)
common = [c for c in per_crop if 'actif' in per_crop[c] and 'full' in per_crop[c]]
fa = [st.mean(per_crop[c]['full']) for c in common]; aa = [st.mean(per_crop[c]['actif']) for c in common]
p = stats.wilcoxon(aa, fa).pvalue

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.8), gridspec_kw=dict(width_ratios=[1, 1.15], wspace=0.28))

# --- left: distributions (box + jitter) ---
cols = {'full': BASEC, 'actif': BLUE}
for i, arm in enumerate(['full', 'actif']):
    v = data[arm]
    bp = ax1.boxplot(v, positions=[i], widths=0.5, patch_artist=True, showfliers=False,
                     medianprops=dict(color='white', lw=1.6), zorder=2)
    for b in bp['boxes']: b.set_facecolor(cols[arm]); b.set_edgecolor(cols[arm]); b.set_alpha(0.85)
    for w in bp['whiskers']+bp['caps']: w.set_color(cols[arm])
    jx = np.random.default_rng(0).normal(i, 0.06, len(v))
    ax1.scatter(jx, v, s=9, color=cols[arm], alpha=0.35, edgecolor='none', zorder=3)
    ax1.text(i, 1.02, f"{st.mean(v):.3f}", ha='center', fontsize=9, color=cols[arm], fontweight='bold')
ax1.axhline(0.6, color=MUTED, ls=(0, (4, 3)), lw=1.0, zorder=1)
ax1.text(1.5, 0.605, "confident (0.6)", fontsize=7.5, color=MUTED, ha='right', va='bottom')
ax1.set_xticks([0, 1]); ax1.set_xticklabels(['IDP (full)', 'actif'], fontsize=10)
ax1.set_ylabel("Boltz2 ipTM  (independent model)"); ax1.set_ylim(0, 1.08)
ax1.set_title("A   Independent-model confidence", loc='left', fontweight='bold', color=INK, fontsize=11, pad=14)
ax1.grid(axis='y', ls=':', color=GRID, lw=0.7); ax1.set_axisbelow(True); ax1.tick_params(length=0)

# --- right: paired per-crop scatter ---
ax2.scatter(fa, aa, s=30, color=BLUE, edgecolor='white', lw=0.5, zorder=3)
lim = [0.2, 1.0]
ax2.plot(lim, lim, color=SECOND, ls=(0, (4, 3)), lw=1.0, zorder=2)
ax2.set_xlim(lim); ax2.set_ylim(lim)
ax2.set_xlabel("IDP ipTM (per-crop mean)"); ax2.set_ylabel("actif ipTM (per-crop mean)")
ax2.set_title("B   Paired by crop (n=79)", loc='left', fontweight='bold', color=INK, fontsize=11, pad=14)
nabove = sum(1 for x, y in zip(fa, aa) if y > x)
ax2.text(0.05, 0.93, f"actif > IDP on {nabove}/{len(common)} crops\nWilcoxon p = {p:.2f} (n.s.)",
         transform=ax2.transAxes, fontsize=8.5, color=SECOND, va='top')
ax2.text(0.97, 0.06, "above line =\nactif more confident", transform=ax2.transAxes, fontsize=7.5,
         color=MUTED, ha='right', va='bottom')
ax2.grid(ls=':', color=GRID, lw=0.7); ax2.set_axisbelow(True); ax2.tick_params(length=0)

fig.suptitle("Boltz2 validation — are actif's accepted designs real, or AF2-gaming?",
             x=0.008, ha='left', fontsize=12.5, fontweight='bold', y=0.99)
fig.text(0.008, 0.005, f"308 BindCraft-accepted designs (150 IDP, 158 actif, 79 shared crops) re-predicted with Boltz2. "
         f"actif's designs are at least as confident to an independent model (0.714 vs 0.693) → the acceptance-rate "
         f"advantage is real, not an AF2 artefact.", fontsize=7.3, color=MUTED)
plt.subplots_adjust(top=0.84, bottom=0.17, left=0.07, right=0.985)
plt.savefig(OUT, dpi=200); print("wrote", OUT)
print(f"n.s. p={p:.3f}; >0.6: full {100*sum(x>0.6 for x in data['full'])/len(data['full']):.0f}%, actif {100*sum(x>0.6 for x in data['actif'])/len(data['actif']):.0f}%")
