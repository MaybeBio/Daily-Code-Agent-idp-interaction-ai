#!/usr/bin/env python3
"""Sectioned overview figure for the five-target ablation."""
import os, csv, json, glob, math, statistics as st
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

EXP = "/capstor/store/cscs/swissai/ab011/jkuhn/idp_bindcraft/exp_panel_ablation"
OUTPNG = f"{EXP}/analysis/fig_ablation_overview_errorbars.png"

BLUE, ORANGE = "#2a78d6", "#eb6834"
INK, SECOND, MUTED = "#0b0b0b", "#52514e", "#8a8a85"
GREY_BAR, BASE_C = "#c9c7c1", "#34475e"
plt.rcParams.update({
    "figure.facecolor": "white", "axes.facecolor": "white", "font.family": "DejaVu Sans",
    "font.size": 10, "axes.edgecolor": SECOND, "axes.linewidth": 0.8,
    "axes.spines.top": False, "axes.spines.right": False,
    "xtick.color": SECOND, "ytick.color": SECOND, "axes.titlesize": 11, "axes.labelsize": 10,
})

# ---- macro data ----
def dir_metrics(d):
    fd, tr = f"{d}/final_design_stats.csv", f"{d}/trajectory_stats.csv"
    ntr = (sum(1 for _ in open(tr)) - 1) if os.path.isfile(tr) else 0
    rows = list(csv.DictReader(open(fd))) if os.path.isfile(fd) else []
    ifc = []
    for r in rows:
        try: ifc.append(float(r["Average_i_pTM_iface"]))
        except: pass
    return len(rows), ntr, ifc

ARMDATA, SEC = {}, {}
for j in open(f"{EXP}/jobs.txt").read().splitlines():
    if not j.strip(): continue
    sp = j.split('\t')[0]
    sec, arm, crop = os.path.basename(sp)[:-5].split('__', 2)
    gene = crop.split('_')[0]
    d = json.load(open(sp))["design_path"].rstrip('/')
    ac, tr, ifc = dir_metrics(d)
    ARMDATA.setdefault(arm, {})[gene] = (ac, tr, ifc)
    SEC[arm] = sec

def macro(arm):
    a8, ifcm = [], []
    for gene, (ac, tr, ifc) in ARMDATA[arm].items():
        if tr: a8.append(ac*800/tr)
        if ifc: ifcm.append(st.mean(ifc))
    def ms(v): return (st.mean(v), st.stdev(v)/math.sqrt(len(v)) if len(v) > 1 else 0.0)
    return ms(a8), ms(ifcm)

SECTIONS = [("1  Loss term", ['no_losses', 'no_pae_loss', 'no_plddt_loss', 'with_template', 'vanilla_filters']),
            ("2  ipSAE / BECraft", ['becraft_w1_lam0', 'ipsae', 'ipsae_replace']),
            ("3  Energy", ['actif_w1_lam0', 'energy_w100_lam1', 'energy_w100_lam0', 'energy_w015_lam1']),
            ("4  Preprocessing", ['random'])]
PRETTY = {'actif_w1_lam0': 'actif', 'becraft_w1_lam0': 'becraft', 'energy_w100_lam1': 'energy·λ1',
          'energy_w100_lam0': 'energy·λ0', 'energy_w015_lam1': 'energy·w.15', 'ipsae_replace': 'ipSAE-rep',
          'ipsae': 'ipSAE', 'no_plddt_loss': 'no pLDDT', 'no_pae_loss': 'no PAE', 'no_losses': 'no losses',
          'with_template': 'w/ templ', 'vanilla_filters': 'vanilla', 'random': 'random', 'full': 'full'}
CTRL = {'random', 'vanilla_filters'}
def barcol(arm):
    if arm == 'full': return BASE_C
    if 'actif' in arm: return BLUE
    if arm in CTRL: return ORANGE
    return GREY_BAR

base_a8, base_ifc = macro('full')

fig, axes = plt.subplots(2, 4, figsize=(15, 7.4),
                         gridspec_kw=dict(width_ratios=[6, 4, 5, 2.4], hspace=0.30, wspace=0.18))
ACC_MAX = 430
for col, (title, arms) in enumerate(SECTIONS):
    full_first = ['full'] + arms
    for row, (which, ymax, ylab) in enumerate([('a8', ACC_MAX, 'accepted / 800 trajectories'),
                                               ('ifc', 0.55, r'i_pTM$_{iface}$')]):
        ax = axes[row][col]
        x = np.arange(len(full_first))
        vals, ses, cols = [], [], []
        for arm in full_first:
            (a8m, a8e), (ifm, ife) = macro(arm)
            v, e = (a8m, a8e) if which == 'a8' else (ifm, ife)
            vals.append(v); ses.append(e); cols.append(barcol(arm))
        ax.bar(x, vals, color=cols, edgecolor='white', linewidth=0.8, width=0.74, zorder=3)
        ax.errorbar(x, vals, yerr=ses, fmt='none', ecolor=SECOND, elinewidth=0.9, capsize=2.5, zorder=4)
        bval = base_a8[0] if which == 'a8' else base_ifc[0]
        ax.axhline(bval, color=BASE_C, ls=(0, (4, 3)), lw=1.0, zorder=2)
        ax.set_ylim(0, ymax)
        ax.set_xticks(x)
        ax.set_xticklabels([PRETTY.get(a, a) for a in full_first], rotation=38, ha='right', fontsize=8)
        for xi, (a, v, e) in zip(x, zip(full_first, vals, ses)):
            lab = f"{v:.0f}" if which == 'a8' else f"{v:.2f}"
            ax.text(xi, v + e + ymax*0.015, lab, ha='center', va='bottom', fontsize=7.6, color=MUTED)
        if col == 0:
            ax.set_ylabel(f"{ylab}\n(macro ± SE)", fontsize=9.5)
        else:
            ax.tick_params(labelleft=False)
        if row == 0:
            ax.set_title(title, loc='left', fontweight='bold', color=INK, fontsize=10.5, pad=6)
        ax.grid(axis='y', ls=':', color="#d9d7d1", lw=0.7); ax.set_axisbelow(True); ax.tick_params(length=0)

leg = [Patch(fc=BLUE, label='actif (standout)'), Patch(fc=BASE_C, label='IDP-BC baseline (full)'),
       Patch(fc=GREY_BAR, label='candidate arm'), Patch(fc=ORANGE, label='control')]
fig.legend(handles=leg, loc='upper right', frameon=False, fontsize=8.5, ncol=4, bbox_to_anchor=(0.995, 0.995))
fig.suptitle("Panel-5 ablation — throughput & interface quality by section",
             x=0.012, ha='left', fontsize=13, fontweight='bold', y=0.995)
fig.text(0.012, 0.015,
         "Macro-averaged over 5 targets (per-target mean, then mean); dashed line = IDP-BC baseline.  "
         "No significance stars: at n=5 the paired test cannot reach p<0.05 (min p = 2/2⁵ = 0.063) — "
         "powered tests come from the n=100 scale-up.", fontsize=8, color=MUTED)
plt.subplots_adjust(top=0.90, bottom=0.14, left=0.055, right=0.995)
plt.savefig(OUTPNG, dpi=200)
print("updated", OUTPNG)
