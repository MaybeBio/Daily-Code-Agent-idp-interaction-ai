#!/usr/bin/env python3
"""Per-target ablation statistics and paired Wilcoxon tests."""
import os, csv, json, math, sys, statistics as st
from collections import defaultdict
import numpy as np
from scipy import stats

EXP = "/capstor/store/cscs/swissai/ab011/jkuhn/idp_bindcraft/exp_panel_ablation"
OUT = f"{EXP}/analysis"
GENES = ["TARDBP", "APP", "TP53", "G3BP1", "MAPT"]

def load():
    PER = defaultdict(dict)
    for j in open(f"{EXP}/jobs.txt").read().splitlines():
        if not j.strip(): continue
        sp = j.split('\t')[0]
        sec, arm, crop = os.path.basename(sp)[:-5].split('__', 2)
        gene = crop.split('_')[0]
        d = json.load(open(sp))["design_path"].rstrip('/')
        trf, fdf = f"{d}/trajectory_stats.csv", f"{d}/final_design_stats.csv"
        n_tr = (sum(1 for _ in open(trf)) - 1) if os.path.isfile(trf) else 0
        R = ([r for r in csv.DictReader(open(fdf)) if r.get("Average_i_pTM_iface")]
             if os.path.isfile(fdf) else [])
        PER[arm][gene] = dict(
            sec=sec, tr=n_tr, ac=len(R),
            ifc=[float(r['Average_i_pTM_iface']) for r in R],
            iptm=[float(r['Average_i_pTM']) for r in R],
            ipae=[float(r['Average_i_pAE']) for r in R])
    return PER

PER = load()

def per_target_vectors(arm, key):
    """Return (genes, values) of the per-target mean for `key`, genes with data only."""
    g, v = [], []
    for gene in GENES:
        c = PER[arm].get(gene)
        if not c: continue
        if key == 'acc':
            if c['tr'] == 0: continue
            g.append(gene); v.append(c['ac'] * 800 / c['tr'])
        else:
            if not c[key]: continue
            g.append(gene); v.append(st.mean(c[key]))
    return g, v

def paired_wilcoxon(arm, base, key):
    """Paired test across targets present in BOTH arms."""
    ga, va = per_target_vectors(arm, key)
    gb, vb = per_target_vectors(base, key)
    da, db = dict(zip(ga, va)), dict(zip(gb, vb))
    common = [g for g in GENES if g in da and g in db]
    if len(common) < 3: return float('nan'), len(common)
    a = [da[g] for g in common]; b = [db[g] for g in common]
    if all(abs(x - y) < 1e-12 for x, y in zip(a, b)): return 1.0, len(common)
    try:
        return stats.wilcoxon(a, b).pvalue, len(common)
    except Exception:
        return float('nan'), len(common)

def two_prop_z(k1, n1, k2, n2):
    if min(n1, n2) == 0: return 1.0
    p1, p2, p = k1/n1, k2/n2, (k1+k2)/(n1+n2)
    se = math.sqrt(p*(1-p)*(1/n1 + 1/n2))
    return 2*(1 - stats.norm.cdf(abs((p1-p2)/se))) if se > 0 else 1.0

def star(p):
    if p != p: return "n/a"
    return '***' if p < .001 else '**' if p < .01 else '*' if p < .05 else 'ns'

BASE = 'full'
rows = []
for arm in PER:
    cells = PER[arm]
    tr = sum(c['tr'] for c in cells.values()); ac = sum(c['ac'] for c in cells.values())
    _, acc_v = per_target_vectors(arm, 'acc')
    _, ifc_v = per_target_vectors(arm, 'ifc')
    _, ipt_v = per_target_vectors(arm, 'iptm')
    _, ipa_v = per_target_vectors(arm, 'ipae')
    n_contrib = sum(1 for c in cells.values() if c['ac'] > 0)
    complete = all(c['tr'] >= 48 for c in cells.values())
    p_acc_w, n_pair = paired_wilcoxon(arm, BASE, 'acc')
    bt = sum(c['tr'] for c in PER[BASE].values()); ba = sum(c['ac'] for c in PER[BASE].values())
    rows.append(dict(
        arm=arm, section=list(cells.values())[0]['sec'],
        n_targets_with_designs=n_contrib, budget_complete=complete,
        traj_total=tr, acc_total=ac,
        acc_per_800_macro=st.mean(acc_v) if acc_v else float('nan'),
        acc_per_800_macro_sd=st.stdev(acc_v) if len(acc_v) > 1 else 0.0,
        acc_per_800_micro=ac*800/tr if tr else float('nan'),
        i_pTM_iface_macro=st.mean(ifc_v) if ifc_v else float('nan'),
        i_pTM_iface_macro_sd=st.stdev(ifc_v) if len(ifc_v) > 1 else 0.0,
        i_pTM_macro=st.mean(ipt_v) if ipt_v else float('nan'),
        i_pAE_macro=st.mean(ipa_v) if ipa_v else float('nan'),
        p_acc_paired=p_acc_w, n_paired_targets=n_pair,
        p_acc_pooled_z=1.0 if arm == BASE else two_prop_z(ac, tr, ba, bt),
        p_ifc_paired=paired_wilcoxon(arm, BASE, 'ifc')[0],
    ))

rows.sort(key=lambda r: -(r['acc_per_800_macro'] if r['acc_per_800_macro'] == r['acc_per_800_macro'] else -1))
with open(f"{OUT}/macro_summary.csv", 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)

# per-target breakdown
with open(f"{OUT}/per_target_breakdown.csv", 'w', newline='') as fh:
    w = csv.writer(fh)
    w.writerow(['arm','target','traj','accepted','acc_per_800','i_pTM_iface_mean','i_pTM_mean','i_pAE_mean'])
    for arm in PER:
        for gene in GENES:
            c = PER[arm].get(gene)
            if not c: continue
            w.writerow([arm, gene, c['tr'], c['ac'],
                        round(c['ac']*800/c['tr'],1) if c['tr'] else '',
                        round(st.mean(c['ifc']),3) if c['ifc'] else '',
                        round(st.mean(c['iptm']),3) if c['iptm'] else '',
                        round(st.mean(c['ipae']),3) if c['ipae'] else ''])

B = [r for r in rows if r['arm'] == BASE][0]
print(f"{'arm':<20}{'acc/800':>9}{'±sd':>7}{'p(pair)':>9}{'p(pool)':>10}  {'ifc':>7}{'±sd':>7}{'p':>8}  {'tgt':>4} {'budget':>7}")
print("-"*104)
for r in rows:
    flag = ""
    if r['n_targets_with_designs'] < 5: flag += f" [{r['n_targets_with_designs']}/5 tgt]"
    if not r['budget_complete']: flag += " [budget incomplete]"
    pp = "—" if r['arm'] == BASE else f"{r['p_acc_paired']:.3f}{star(r['p_acc_paired'])}"
    pz = "—" if r['arm'] == BASE else f"{r['p_acc_pooled_z']:.1e}{star(r['p_acc_pooled_z'])}"
    pi = "—" if r['arm'] == BASE else f"{r['p_ifc_paired']:.3f}"
    print(f"{r['arm']:<20}{r['acc_per_800_macro']:>9.1f}{r['acc_per_800_macro_sd']:>7.0f}"
          f"{pp:>12}{pz:>13}  {r['i_pTM_iface_macro']:>7.3f}{r['i_pTM_iface_macro_sd']:>7.3f}{pi:>8}"
          f"  {r['n_targets_with_designs']:>4} {'ok' if r['budget_complete'] else 'INCOMPL':>7}{flag}")
print(f"\nwrote {OUT}/macro_summary.csv and {OUT}/per_target_breakdown.csv")
