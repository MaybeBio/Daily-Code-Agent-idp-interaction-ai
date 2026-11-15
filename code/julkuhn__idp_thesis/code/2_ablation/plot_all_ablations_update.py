#!/usr/bin/env python3
"""Five-target ablation figure: accepted-design yield and the four confidence scores."""
import os, csv, glob, math, statistics as st
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

OUT = "/capstor/store/cscs/swissai/ab011/jkuhn/idp_bindcraft/exp_panel_ablation/analysis/fig_all_ablations_errorbars.png"
BLUE, GREEN, RED, PURPLE, ORANGE = "#2a78d6", "#1baf7a", "#e34948", "#6250d6", "#eb6834"
INK, SECOND, MUTED, GRID, BASEC = "#0b0b0b", "#52514e", "#8a8a85", "#d9d7d1", "#34475e"
plt.rcParams.update({"figure.facecolor": "white", "axes.facecolor": "white", "font.family": "DejaVu Sans",
    "font.size": 10, "axes.edgecolor": SECOND, "axes.linewidth": 0.8, "axes.spines.top": False,
    "axes.spines.right": False, "xtick.color": SECOND, "ytick.color": SECOND})
def se(v): return st.stdev(v)/math.sqrt(len(v)) if len(v) > 1 else 0.0
NAN = float('nan')

EXP5 = "/capstor/store/cscs/swissai/ab011/jkuhn/idp_bindcraft/exp_panel_ablation"
IDT = "/capstor/store/cscs/swissai/ab011/jkuhn/idp_testing"
import json
GENE_OF = {}

def bc1_arm_dirs():
    """arm -> {gene: design_dir} from the ablation jobs.txt."""
    out = {}
    for j in open(f"{EXP5}/jobs.txt").read().splitlines():
        if not j.strip(): continue
        sp = j.split('\t')[0]; sec, arm, crop = os.path.basename(sp)[:-5].split('__', 2)
        d = json.load(open(sp))["design_path"].rstrip('/')
        out.setdefault(arm, {})[crop.split('_')[0]] = d
    return out
P5 = bc1_arm_dirs()

def bc1_scores(arm):
    cols = {'i_pTM': 'Average_i_pTM', 'i_pAE': 'Average_i_pAE',
            'i_pTM_iface': 'Average_i_pTM_iface', 'i_pAE_iface': 'Average_i_pAE_iface'}
    acc = {'a8': []}; sc = {k: [] for k in cols}
    for gene, d in P5[arm].items():
        fd, ts = f"{d}/final_design_stats.csv", f"{d}/trajectory_stats.csv"
        ntr = (sum(1 for _ in open(ts)) - 1) if os.path.isfile(ts) else 0
        rows = list(csv.DictReader(open(fd))) if os.path.isfile(fd) else []
        if ntr: acc['a8'].append(len(rows)*800/ntr)
        for k, c in cols.items():
            v = []
            for r in rows:
                try: v.append(float(r[c]))
                except: pass
            if v: sc[k].append(st.mean(v))
    return acc, sc

def bc2_scores(outpath):
    """BC2 UNCAPPED runs. acc/800 = per-target acc/traj*800 (macro); traj from campaign_state;
    i_pTM/i_pDAE from the live 3_Ranked/!_Ranked.csv (summary.csv absent: runs were stopped mid-leg)."""
    import json
    tg = ['TARDBP_Q13148','TP53_P04637','APP_P05067','G3BP1_Q13283','MAPT_P10636']
    a8, iptm, pdae, zero = [], [], [], 0
    for t in tg:
        d = f"{outpath}/{t}"
        acc = len(glob.glob(f"{d}/3_Ranked/*.cif"))
        if acc == 0: zero += 1
        tr = 0; sf = f"{d}/.campaign_state.json"
        if os.path.exists(sf):
            try: tr = json.load(open(sf)).get('trajectories', 0) or 0
            except: pass
        if tr: a8.append(acc/tr*800)
        rc = f"{d}/3_Ranked/!_Ranked.csv"
        if os.path.exists(rc):
            for r in csv.DictReader(open(rc)):
                try: iptm.append(float(r['i_pTM'])); pdae.append(float(r['i_pDAE']))
                except: pass
    return a8, iptm, pdae, zero
BC2MAP = {'BC2_V7': f"{IDT}/V7_bc2v104_idphead/output_V7_uncapped",
          'BC2_V8a': f"{IDT}/V8_bc2v104_idpfull/output_V8a_uncapped",
          'BC2_V8b': f"{IDT}/V8_bc2v104_idpfull/output_V8b_uncapped"}

ORDER = [('full','1_loss'),('no_losses','1_loss'),('no_pae_loss','1_loss'),('no_plddt_loss','1_loss'),
         ('with_template','1_loss'),('vanilla_filters','1_loss'),
         ('becraft_w1_lam0','2_ipsae'),('ipsae','2_ipsae'),('ipsae_replace','2_ipsae'),
         ('actif_w1_lam0','3_energy'),('energy_w100_lam1','3_energy'),('energy_w100_lam0','3_energy'),
         ('energy_w015_lam1','3_energy'),('random','4_preproc'),
         ('BC2_V7','BC2'),('BC2_V8a','BC2'),('BC2_V8b','BC2')]
SECCOL = {'1_loss': BLUE, '2_ipsae': GREEN, '3_energy': RED, '4_preproc': PURPLE, 'BC2': ORANGE}
PRETTY = {'full':'IDP','no_losses':'no_losses','no_pae_loss':'no_pae','no_plddt_loss':'no_plddt',
          'with_template':'w/_templ','vanilla_filters':'vanilla','becraft_w1_lam0':'becraft','ipsae':'ipsae',
          'ipsae_replace':'ipsae_rep','actif_w1_lam0':'actif','energy_w100_lam1':'energy·λ1',
          'energy_w100_lam0':'energy·λ0','energy_w015_lam1':'energy·w.15','random':'random',
          'BC2_V7':'BC2\nnative','BC2_V8a':'BC2\n+loss','BC2_V8b':'BC2\n+IDP gates'}

# metric dicts: value + SE, plus hatch/na flags
MET = {k: {} for k in ['acc','i_pTM','i_pAE','i_pTM_iface','i_pAE_iface']}
MSE = {k: {} for k in MET}; HATCH = {}
for arm, sec in ORDER:
    HATCH[arm] = (sec == 'BC2')
    if sec == 'BC2':
        a8, iptm, pdae, _zero = bc2_scores(BC2MAP[arm])
        MET['acc'][arm] = st.mean(a8) if a8 else 0; MSE['acc'][arm] = se(a8) if a8 else 0
        MET['i_pTM'][arm] = st.mean(iptm) if iptm else NAN; MSE['i_pTM'][arm] = se(iptm) if iptm else 0
        MET['i_pAE'][arm] = NAN; MSE['i_pAE'][arm] = 0
        MET['i_pTM_iface'][arm] = st.mean(pdae) if pdae else NAN; MSE['i_pTM_iface'][arm] = se(pdae) if pdae else 0  # i_pDAE
        MET['i_pAE_iface'][arm] = NAN; MSE['i_pAE_iface'][arm] = 0
    else:
        acc, sc = bc1_scores(arm)
        MET['acc'][arm] = st.mean(acc['a8']) if acc['a8'] else 0; MSE['acc'][arm] = se(acc['a8']) if acc['a8'] else 0
        for k in ['i_pTM','i_pAE','i_pTM_iface','i_pAE_iface']:
            MET[k][arm] = st.mean(sc[k]) if sc[k] else NAN; MSE[k][arm] = se(sc[k]) if sc[k] else 0

arms = [a for a, _ in ORDER]; x = np.arange(len(arms))
cols = [BASEC if a == 'full' else SECCOL[s] for a, s in ORDER]

PANELS = [  # (metric key, ylabel, ymax, higher_better)
    ('acc', "accepted / 800 traj\n(macro ± SE)", None, True),
    ('i_pTM', "i_pTM  (global ↑)", 0.95, True),
    ('i_pAE', "i_pAE  (global ↓)", None, False),
    ('i_pTM_iface', "i_pTM_iface ↑\n(BC2: i_pDAE, hatched)", 0.9, True),
    ('i_pAE_iface', "i_pAE_iface  (↓)", None, False),
]
fig, axes = plt.subplots(len(PANELS), 1, figsize=(14, 15), sharex=True)
for ax, (key, ylab, ymax, hi) in zip(axes, PANELS):
    V, E = MET[key], MSE[key]
    vals = [V[a] for a in arms]; ers = [E[a] for a in arms]
    plotv = [0 if v != v else v for v in vals]
    bars = ax.bar(x, plotv, color=cols, edgecolor='white', linewidth=0.8, width=0.76, zorder=3)
    for b, a in zip(bars, arms):
        if HATCH[a]: b.set_hatch('////'); b.set_edgecolor('white')
    ax.errorbar(x, plotv, yerr=[0 if e != e else e for e in ers], fmt='none', ecolor=SECOND,
                elinewidth=0.9, capsize=2.5, zorder=4)
    base = V['full']
    if base == base: ax.axhline(base, color=BASEC, ls=(0, (4, 3)), lw=1.0, zorder=2)
    top = max([v for v in plotv] + [ (V[a]+E[a]) for a in arms if V[a]==V[a] ]) * 1.18
    ax.set_ylim(0, ymax if ymax else top)
    for xi, a in zip(x, arms):
        v, e = V[a], E[a]
        if v != v:
            ax.text(xi, ax.get_ylim()[1]*0.03, 'n/a', ha='center', fontsize=6.5, color=MUTED, style='italic'); continue
        fmt = f"{v:.0f}" if key == 'acc' else f"{v:.2f}"
        ax.text(xi, v + e + ax.get_ylim()[1]*0.015, fmt, ha='center', va='bottom', fontsize=7, color=MUTED)
    ax.set_ylabel(ylab, fontsize=9.5)
    ax.grid(axis='y', ls=':', color=GRID, lw=0.7); ax.set_axisbelow(True); ax.tick_params(length=0)
    for xi in [5.5, 8.5, 12.5, 13.5]: ax.axvline(xi, color="#e4e2dc", lw=1.0, zorder=0)

axes[0].set_title("Panel-5 loss ablation + BC2 engine — same 5 IDR targets · all four confidence scores",
                  loc='left', fontweight='bold', color=INK, fontsize=13, pad=10)
leg = [Patch(fc=BASEC, label='IDP (baseline)'), Patch(fc=BLUE, label='loss-term'),
       Patch(fc=GREEN, label='ipSAE/BECraft'), Patch(fc=RED, label='energy'), Patch(fc=PURPLE, label='preproc'),
       Patch(fc=ORANGE, hatch='////', ec='white', label='BC2 engine')]
axes[0].legend(handles=leg, loc='upper right', frameon=False, fontsize=8.5, ncol=3)
axes[-1].set_xticks(x); axes[-1].set_xticklabels([PRETTY[a] for a in arms], rotation=50, ha='right', fontsize=8)
fig.text(0.008, 0.006, "5 IDR targets; macro over targets, dashed line = IDP baseline. "
         "↑ higher=better (i_pTM, i_pTM_iface) · ↓ lower=better (i_pAE, i_pAE_iface). "
         "BC2 (separate codebase) has no i_pTM_iface/i_pAE_iface: interface bar = i_pDAE (hatched), i_pAE/i_pAE_iface = n/a; "
         "BC2 = uncapped reruns (~700 traj/target, all 5 targets now covered). No significance stars (n=5 caps paired p at 0.063).",
         fontsize=7.2, color=MUTED)
plt.subplots_adjust(top=0.955, bottom=0.075, left=0.06, right=0.99, hspace=0.14)
plt.savefig(OUT, dpi=190); print("wrote", OUT)
