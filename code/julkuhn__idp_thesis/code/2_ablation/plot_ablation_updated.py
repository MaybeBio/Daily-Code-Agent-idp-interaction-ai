#!/usr/bin/env python3
"""Ablation figures for the five-target panel and the BindCraft 2 arms."""
import os, csv, json, glob, statistics as st
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = "/capstor/store/cscs/swissai/ab011/jkuhn/idp_bindcraft/exp_panel_ablation/analysis/updated_2026-10-22"
os.makedirs(OUT, exist_ok=True)

def dir_metrics(d):
    """accepted count, trajectory count, per-design score means from one output dir."""
    fd, tr = f"{d}/final_design_stats.csv", f"{d}/trajectory_stats.csv"
    ntr = (sum(1 for _ in open(tr)) - 1) if os.path.isfile(tr) else 0
    rows = list(csv.DictReader(open(fd))) if os.path.isfile(fd) else []
    def col(c):
        v = []
        for r in rows:
            try: v.append(float(r[c]))
            except: pass
        return v
    return dict(ac=len(rows), tr=ntr,
                ifc=col("Average_i_pTM_iface"), iptm=col("Average_i_pTM"),
                ipae=col("Average_i_pAE"), ipaei=col("Average_i_pAE_iface"))

def macro(arm_dirs):
    """arm -> {gene:dir}. Return per-arm macro dict."""
    out = {}
    for arm, gdirs in arm_dirs.items():
        a8, ifc, ipae_i, zero, tgt = [], [], [], 0, 0
        for g, d in gdirs.items():
            m = dir_metrics(d); tgt += 1
            if m['tr']: a8.append(m['ac']*800/m['tr'])
            if m['ac'] == 0: zero += 1
            if m['ifc']: ifc.append(st.mean(m['ifc']))
            if m['ipaei']: ipae_i.append(st.mean(m['ipaei']))
        out[arm] = dict(acc800=st.mean(a8) if a8 else 0.0,
                        acc800_sd=st.stdev(a8) if len(a8) > 1 else 0.0,
                        ifc=st.mean(ifc) if ifc else float('nan'),
                        ipaei=st.mean(ipae_i) if ipae_i else float('nan'),
                        zero=zero, ntgt=tgt)
    return out

# ---- Panel-5 (BC1) ----
EXP5 = "/capstor/store/cscs/swissai/ab011/jkuhn/idp_bindcraft/exp_panel_ablation"
p5_dirs = {}
for j in open(f"{EXP5}/jobs.txt").read().splitlines():
    if not j.strip(): continue
    sp = j.split('\t')[0]
    sec, arm, crop = os.path.basename(sp)[:-5].split('__', 2)
    d = json.load(open(sp))["design_path"].rstrip('/')
    p5_dirs.setdefault(arm, {})[crop.split('_')[0]] = d
m5 = macro(p5_dirs)

# ---- 16-target (BC1) ----
EXP16 = "/capstor/store/cscs/swissai/ab011/jkuhn/idp_bindcraft/exp_ablation_v2/output"
t16_dirs = {}
for d in glob.glob(f"{EXP16}/*"):
    b = os.path.basename(d)
    if '__' not in b: continue
    crop, arm = b.split('__', 1)
    t16_dirs.setdefault(arm, {})[crop] = d
m16 = macro(t16_dirs)

# ---- BC2 (Panel-5) ----
def bc2(base):
    tg = ['TARDBP_Q13148','TP53_P04637','APP_P05067','G3BP1_Q13283','MAPT_P10636']
    a8, zero = [], 0
    for t in tg:
        acc = len(glob.glob(f"{base}/{t}/3_Ranked/*.cif"))
        a8.append(acc/50*800)
        if acc == 0: zero += 1
    return dict(acc800=st.mean(a8), zero=zero, ntgt=5)
bc2_arms = {
 'BC2 V7 (native)':  bc2("/capstor/store/cscs/swissai/ab011/jkuhn/idp_testing/V7_bc2v104_idphead/output"),
 'BC2 V8a (+loss)':  bc2("/capstor/store/cscs/swissai/ab011/jkuhn/idp_testing/V8_bc2v104_idpfull/output_a"),
 'BC2 V8b (+gateswap)': bc2("/capstor/store/cscs/swissai/ab011/jkuhn/idp_testing/V8_bc2v104_idpfull/output_b"),
}

# ================= FIG 1: Panel-5 acc/800 + i_pTM_iface =================
def barplot(ax, items, title, baseline='full', show_ifc=True):
    items = sorted(items, key=lambda kv: kv[1]['acc800'])
    names = [k for k, _ in items]
    vals = [v['acc800'] for _, v in items]
    cols = []
    for k, v in items:
        if k == baseline: cols.append('#d62728')       # baseline rot
        elif 'actif' in k: cols.append('#2ca02c')       # actif grün
        elif k.startswith('BC2'): cols.append('#9467bd')# BC2 lila
        else: cols.append('#1f77b4')
    y = np.arange(len(names))
    ax.barh(y, vals, color=cols, edgecolor='k', lw=0.4)
    ax.set_yticks(y); ax.set_yticklabels(names, fontsize=8)
    ax.set_xlabel('acc / 800 Trajektorien (macro)')
    ax.set_title(title, fontsize=10, fontweight='bold')
    for yi, (k, v) in zip(y, items):
        lab = f"{v['acc800']:.0f}"
        if v.get('zero', 0): lab += f"  ({v['zero']}/{v['ntgt']}✗)"
        ax.text(v['acc800']+max(vals)*0.01, yi, lab, va='center', fontsize=7)
    if show_ifc:
        ax2 = ax.twiny()
        xi = [v.get('ifc', float('nan')) for _, v in items]
        ax2.plot(xi, y, 'D', color='k', ms=4, label='i_pTM_iface')
        ax2.set_xlim(0.30, 0.45); ax2.set_xlabel('i_pTM_iface (◆)', fontsize=8)
    ax.grid(axis='x', ls=':', alpha=0.5)

fig, ax = plt.subplots(figsize=(8.5, 6))
barplot(ax, list(m5.items()) + list(bc2_arms.items()), 'Panel-5 Ablation (BC1-Loss-Arme + BC2)')
plt.tight_layout(); plt.savefig(f"{OUT}/fig_panel5_ablation.png", dpi=150); plt.close()

# ================= FIG 2: 16-target acc/800 =================
fig, ax = plt.subplots(figsize=(8, 4))
barplot(ax, list(m16.items()), '16-Target Ablation (BC1)  — full = läuft noch (partiell)')
plt.tight_layout(); plt.savefig(f"{OUT}/fig_16target_ablation.png", dpi=150); plt.close()

# ================= FIG 3: Quantität vs Qualität (scatter) =================
fig, ax = plt.subplots(figsize=(8, 6))
for label, m, mk in [('Panel-5', m5, 'o'), ('16-Target', m16, 's')]:
    for arm, v in m.items():
        if v['ifc'] != v['ifc']: continue
        c = '#d62728' if arm == 'full' else '#2ca02c' if 'actif' in arm else '#1f77b4'
        ax.scatter(v['acc800'], v['ifc'], marker=mk, s=70, color=c, edgecolor='k', lw=0.5, zorder=3)
        ax.annotate(f"{arm}\n[{label[:3]}]", (v['acc800'], v['ifc']), fontsize=6,
                    xytext=(4, 3), textcoords='offset points')
ax.set_xlabel('acc / 800 (macro)  →  Durchsatz')
ax.set_ylabel('i_pTM_iface (macro)  →  Interface-Qualität')
ax.set_title('Quantität vs. Qualität  (rot=IDP-BC Baseline, grün=actif)', fontsize=10, fontweight='bold')
ax.grid(ls=':', alpha=0.5)
plt.tight_layout(); plt.savefig(f"{OUT}/fig_quality_vs_quantity.png", dpi=150); plt.close()

# ---- CSV-Dump der geplotteten Werte ----
with open(f"{OUT}/plotted_values.csv", 'w', newline='') as fh:
    w = csv.writer(fh); w.writerow(['set','arm','acc800_macro','acc800_sd','i_pTM_iface','i_pAE_iface','zero_targets','n_targets'])
    for s, m in [('panel5', m5), ('16target', m16)]:
        for arm, v in sorted(m.items(), key=lambda kv: -kv[1]['acc800']):
            w.writerow([s, arm, f"{v['acc800']:.1f}", f"{v['acc800_sd']:.1f}",
                        f"{v['ifc']:.3f}" if v['ifc']==v['ifc'] else '', f"{v['ipaei']:.3f}" if v['ipaei']==v['ipaei'] else '',
                        v['zero'], v['ntgt']])
    for arm, v in bc2_arms.items():
        w.writerow(['bc2_panel5', arm, f"{v['acc800']:.1f}", '', '', '', v['zero'], v['ntgt']])

print("geschrieben nach", OUT)
for f in sorted(glob.glob(f"{OUT}/*")): print("  ", os.path.basename(f))
