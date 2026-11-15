#!/usr/bin/env python3
"""Ablation figures for accepted-design yield and confidence scores."""
import os, csv, json, glob, math, statistics as st
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

OUT = "/capstor/store/cscs/swissai/ab011/jkuhn/idp_bindcraft/exp_panel_ablation/analysis/updated_2026-10-22"
os.makedirs(OUT, exist_ok=True)

# ---- validated palette (dataviz skill, canonical slots) ----
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, SECOND, MUTED = "#0b0b0b", "#52514e", "#8a8a85"
GREY_BAR, GREY_CTRL = "#c9c7c1", "#e4e2dc"
BASE_C = "#34475e"   # baseline slate

plt.rcParams.update({
    "figure.facecolor": "white", "axes.facecolor": "white",
    "font.family": "DejaVu Sans", "font.size": 10,
    "axes.edgecolor": SECOND, "axes.linewidth": 0.8,
    "axes.spines.top": False, "axes.spines.right": False,
    "xtick.color": SECOND, "ytick.color": SECOND,
    "xtick.labelsize": 9, "ytick.labelsize": 9,
    "axes.titlesize": 12, "axes.labelsize": 10.5,
})

# ============ data loading (macro) ============
def dir_metrics(d):
    fd, tr = f"{d}/final_design_stats.csv", f"{d}/trajectory_stats.csv"
    ntr = (sum(1 for _ in open(tr)) - 1) if os.path.isfile(tr) else 0
    rows = list(csv.DictReader(open(fd))) if os.path.isfile(fd) else []
    def col(c):
        v = []
        for r in rows:
            try: v.append(float(r[c]))
            except: pass
        return v
    return dict(ac=len(rows), tr=ntr, ifc=col("Average_i_pTM_iface"), ipaei=col("Average_i_pAE_iface"))

def macro(arm_dirs):
    out = {}
    for arm, gdirs in arm_dirs.items():
        a8, ifc, ipaei, zero = [], [], [], 0
        for g, d in gdirs.items():
            m = dir_metrics(d)
            if m['tr']: a8.append(m['ac']*800/m['tr'])
            if m['ac'] == 0: zero += 1
            if m['ifc']: ifc.append(st.mean(m['ifc']))
            if m['ipaei']: ipaei.append(st.mean(m['ipaei']))
        n = len(a8)
        out[arm] = dict(acc800=st.mean(a8) if a8 else 0.0,
                        acc800_se=(st.stdev(a8)/math.sqrt(n)) if n > 1 else 0.0,
                        ifc=st.mean(ifc) if ifc else float('nan'),
                        ifc_se=(st.stdev(ifc)/math.sqrt(len(ifc))) if len(ifc) > 1 else 0.0,
                        ipaei=st.mean(ipaei) if ipaei else float('nan'),
                        zero=zero, n=n)
    return out

EXP5 = "/capstor/store/cscs/swissai/ab011/jkuhn/idp_bindcraft/exp_panel_ablation"
p5_dirs = {}
for j in open(f"{EXP5}/jobs.txt").read().splitlines():
    if not j.strip(): continue
    sp = j.split('\t')[0]
    sec, arm, crop = os.path.basename(sp)[:-5].split('__', 2)
    d = json.load(open(sp))["design_path"].rstrip('/')
    p5_dirs.setdefault(arm, {})[crop.split('_')[0]] = d
m5 = macro(p5_dirs)

EXP16 = "/capstor/store/cscs/swissai/ab011/jkuhn/idp_bindcraft/exp_ablation_v2/output"
t16_dirs = {}
for d in glob.glob(f"{EXP16}/*"):
    b = os.path.basename(d)
    if '__' not in b: continue
    crop, arm = b.split('__', 1)
    t16_dirs.setdefault(arm, {})[crop] = d
m16 = macro(t16_dirs)

def bc2(base, n=5):
    tg = ['TARDBP_Q13148','TP53_P04637','APP_P05067','G3BP1_Q13283','MAPT_P10636']
    a8, zero = [], 0
    for t in tg:
        acc = len(glob.glob(f"{base}/{t}/3_Ranked/*.cif")); a8.append(acc/50*800)
        if acc == 0: zero += 1
    return dict(acc800=st.mean(a8), acc800_se=st.stdev(a8)/math.sqrt(len(a8)),
                ifc=float('nan'), ipaei=float('nan'), zero=zero, n=n)
bc2_arms = {
 'BC2 V7':  bc2("/capstor/store/cscs/swissai/ab011/jkuhn/idp_testing/V7_bc2v104_idphead/output"),
 'BC2 V8a': bc2("/capstor/store/cscs/swissai/ab011/jkuhn/idp_testing/V8_bc2v104_idpfull/output_a"),
 'BC2 V8b': bc2("/capstor/store/cscs/swissai/ab011/jkuhn/idp_testing/V8_bc2v104_idpfull/output_b"),
}

CTRL = {'random', 'vanilla_filters'}
def classify(arm):
    if arm == 'full': return 'base'
    if 'actif' in arm: return 'actif'
    if arm.startswith('BC2'): return 'bc2'
    if arm in CTRL: return 'ctrl'
    return 'bc1'
COL = {'base': BASE_C, 'actif': BLUE, 'bc2': ORANGE, 'ctrl': GREY_CTRL, 'bc1': GREY_BAR}
PRETTY = {'actif_w1_lam0': 'actif (λ0)', 'actif_w1_lam1': 'actif (λ1)',
          'becraft_w1_lam0': 'becraft', 'energy_w100_lam1': 'energy w1 λ1',
          'energy_w100_lam0': 'energy w1 λ0', 'energy_w015_lam1': 'energy w.15 λ1',
          'ipsae_replace': 'ipSAE-replace', 'ipsae': 'ipSAE',
          'no_plddt_loss': 'no pLDDT-loss', 'no_pae_loss': 'no PAE-loss',
          'no_losses': 'no losses', 'with_template': 'with template',
          'vanilla_filters': 'vanilla filters', 'random': 'random', 'full': 'IDP-BC (full)'}
def lab(a): return PRETTY.get(a, a)

def acc_panel(ax, items, title, base_val):
    items = sorted(items, key=lambda kv: kv[1]['acc800'])
    y = np.arange(len(items))
    vals = [v['acc800'] for _, v in items]
    ses = [v['acc800_se'] for _, v in items]
    cols = [COL[classify(k)] for k, _ in items]
    bars = ax.barh(y, vals, color=cols, edgecolor="white", linewidth=0.8, height=0.72, zorder=3)
    ax.errorbar(vals, y, xerr=ses, fmt='none', ecolor=SECOND, elinewidth=0.9, capsize=2.5, zorder=4)
    # baseline reference line
    if base_val is not None:
        ax.axvline(base_val, color=BASE_C, ls=(0, (4, 3)), lw=1.1, zorder=2)
    ax.set_yticks(y)
    ax.set_yticklabels([lab(k) for k, _ in items], fontsize=8.5)
    for tick, (k, _) in zip(ax.get_yticklabels(), items):
        tick.set_color(INK if classify(k) in ('actif', 'base') else SECOND)
    xmax = max(v + s for v, s in zip(vals, ses)) * 1.17
    ax.set_xlim(0, xmax)
    for yi, (k, v) in zip(y, items):
        t = f"{v['acc800']:.0f}"
        if v.get('zero', 0): t += f"  {v['zero']}/{v['n']}✗"
        ax.text(v['acc800'] + v['acc800_se'] + xmax*0.012, yi, t, va='center',
                fontsize=7.8, color=MUTED)
    ax.set_xlabel("accepted designs per 800 trajectories  (macro ± SE)")
    ax.set_title(title, loc='left', fontweight='bold', color=INK, pad=8)
    ax.grid(axis='x', ls=':', color="#d9d7d1", lw=0.7, zorder=0)
    ax.set_axisbelow(True)
    ax.tick_params(length=0)

# ===== FIG 1: Panel-5 acc/800 =====
fig, ax = plt.subplots(figsize=(7.6, 6.4))
acc_panel(ax, list(m5.items()) + list(bc2_arms.items()),
          "Panel-5 ablation — throughput", base_val=m5['full']['acc800'])
leg = [Patch(fc=BLUE, label='actif (best loss)'), Patch(fc=BASE_C, label='IDP-BC baseline'),
       Patch(fc=GREY_BAR, label='other BC1 loss arms'), Patch(fc=ORANGE, label='BC2 engine'),
       Patch(fc=GREY_CTRL, label='controls')]
ax.legend(handles=leg, loc='lower right', frameon=False, fontsize=8, handlelength=1.1,
          title="✗ = targets with 0 designs", title_fontsize=7.5)
fig.text(0.013, 0.012, "dashed line = IDP-BC baseline · n=5 targets", fontsize=7, color=MUTED)
plt.tight_layout(); plt.savefig(f"{OUT}/fig_panel5_acc800.png", dpi=200, bbox_inches='tight'); plt.close()

# ===== FIG 2: 16-target acc/800 =====
fig, ax = plt.subplots(figsize=(7.2, 3.4))
acc_panel(ax, list(m16.items()), "16-target ablation — throughput  (full = running, partial)",
          base_val=m16.get('full', {}).get('acc800'))
ax.legend(handles=leg[:3], loc='lower right', frameon=False, fontsize=8, handlelength=1.1)
fig.text(0.013, 0.02, "dashed line = IDP-BC baseline (partial) · n=16 targets", fontsize=7, color=MUTED)
plt.tight_layout(); plt.savefig(f"{OUT}/fig_16target_acc800.png", dpi=200, bbox_inches='tight'); plt.close()

# ===== FIG 3: quality vs quantity (true 2-axis scatter) =====
LABEL_ONLY = {'actif_w1_lam0', 'actif_w1_lam1', 'full', 'ipsae_replace', 'becraft_w1_lam0', 'vanilla_filters'}
# quality band of the non-control arms (context)
qvals = [v['ifc'] for m in (m5, m16) for a, v in m.items() if v['ifc'] == v['ifc'] and a not in CTRL]
fig, ax = plt.subplots(figsize=(7.4, 5.8))
ax.axhspan(min(qvals), max(qvals), color="#eef1f5", zorder=0)
for m, mk, setlab in [(m5, 'o', 'Panel-5'), (m16, 's', '16-target')]:
    for arm, v in m.items():
        if v['ifc'] != v['ifc']: continue
        cls = classify(arm); c = COL[cls] if cls in ('actif', 'base') else SECOND
        ax.scatter(v['acc800'], v['ifc'], marker=mk, s=66, facecolor=c, edgecolor='white',
                   linewidth=0.9, zorder=3)
        # label actif in both sets; everything else only once (Panel-5) to avoid duplicate text
        do_label = arm in LABEL_ONLY and (cls == 'actif' or setlab == 'Panel-5')
        if do_label:
            dy = -14 if arm == 'becraft_w1_lam0' else 9
            tag = lab(arm) + (f" [{setlab[:3]}]" if cls == 'actif' else "")
            ax.annotate(tag, (v['acc800'], v['ifc']), fontsize=7.4,
                        color=(INK if cls in ('actif', 'base') else SECOND),
                        xytext=(0, dy), textcoords='offset points', ha='center')
ax.text(0.015, 0.5*(min(qvals)+max(qvals)), "quality band\nof all loss arms", transform=ax.get_yaxis_transform(),
        fontsize=7, color=MUTED, va='center', ha='left')
ax.set_xlabel("throughput — accepted / 800 trajectories (macro)")
ax.set_ylabel("interface quality — i_pTM_iface (macro)")
ax.set_title("Quality vs. throughput", loc='left', fontweight='bold', color=INK, pad=8)
ax.grid(ls=':', color="#d9d7d1", lw=0.7); ax.set_axisbelow(True); ax.tick_params(length=0)
from matplotlib.lines import Line2D
sh = [Line2D([],[],marker='o',ls='',mfc=BLUE,mec='white',ms=8,label='actif'),
      Line2D([],[],marker='o',ls='',mfc=BASE_C,mec='white',ms=8,label='IDP-BC'),
      Line2D([],[],marker='o',ls='',mfc=SECOND,mec='white',ms=8,label='other BC1'),
      Line2D([],[],marker='o',ls='',mfc='none',mec=SECOND,ms=8,label='○ Panel-5'),
      Line2D([],[],marker='s',ls='',mfc='none',mec=SECOND,ms=8,label='□ 16-target')]
ax.legend(handles=sh, loc='lower right', frameon=False, fontsize=8, handlelength=1.0, ncol=1)
plt.tight_layout(); plt.savefig(f"{OUT}/fig_quality_vs_throughput.png", dpi=200, bbox_inches='tight'); plt.close()

print("wrote:", OUT)
for f in sorted(glob.glob(f"{OUT}/*.png")): print("  ", os.path.basename(f))
