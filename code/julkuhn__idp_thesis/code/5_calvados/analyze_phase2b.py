#!/usr/bin/env python3
"""CALVADOS phase 2b: the phase 2 metrics for runs initialised from the bound pose."""
import os, glob, re, json, csv
import numpy as np
import MDAnalysis as mda

ROOT = "/capstor/store/cscs/swissai/ab011/jkuhn/idp_bindcraft"
P2   = f"{ROOT}/calvados_sanity/phase2_complexes"
SUF  = "_prebound"   # <<< einziger Unterschied zu 2a wo es "_complex" war

def iface_from_bc(name):
    crop, arm = name.rsplit('__', 1)
    base = f"{ROOT}/exp_ablation/output/{crop}__full" if arm == 'full' \
           else f"{ROOT}/exp_ablation_v2/output/{crop}__{arm}"
    csvp = f"{base}/final_design_stats.csv"
    if not os.path.isfile(csvp): return set()
    iface = set()
    for row in csv.DictReader(open(csvp)):
        for tok in re.findall(r'[A-Za-z]?(\d+)', row.get('InterfaceResidues','')):
            iface.add(int(tok))
    return iface

def analyze_one(name):
    sim = f"{P2}/sims_prebound/{name}/{name}{SUF}"
    data = f"{sim}/data"
    base = f"{name}{SUF}"
    tgt_rg = f"{data}/{base}_{name}_target_rg.npy"
    bnd_rg = f"{data}/{base}_{name}_binder_rg.npy"
    cmap   = f"{data}/{base}_{name}_target_{name}_binder_cmap.npy"
    com_top = f"{data}/{base}_com_top.pdb"
    com_dcd = f"{data}/{base}_com_traj.dcd"
    if not all(os.path.isfile(p) for p in (tgt_rg, bnd_rg, cmap, com_top, com_dcd)):
        return None
    rg_t = np.load(tgt_rg); rg_b = np.load(bnd_rg); cm = np.load(cmap)
    u = mda.Universe(com_top, com_dcd)
    coms = []
    for ts in u.trajectory:
        if u.atoms.n_atoms >= 2:
            coms.append(u.atoms.positions[:2].copy())
    coms = np.array(coms)
    com_dist = np.linalg.norm(coms[:,0,:] - coms[:,1,:], axis=1) if len(coms) else np.array([])
    bc_iface = iface_from_bc(name)
    calv_iface_counts = cm.sum(axis=1)
    thr = np.quantile(calv_iface_counts, 0.8) if len(calv_iface_counts)>0 else 0
    calv_iface = set(int(i+1) for i in np.where(calv_iface_counts > thr)[0])
    overlap = bc_iface & calv_iface
    jaccard = len(overlap) / max(len(bc_iface | calv_iface), 1)
    return dict(
        name=name,
        n_frames=len(rg_t),
        rg_target_med=float(np.median(rg_t)) if len(rg_t) else 0,
        rg_binder_med=float(np.median(rg_b)) if len(rg_b) else 0,
        com_dist_med=float(np.median(com_dist)) if len(com_dist) else 0,
        com_dist_min=float(com_dist.min()) if len(com_dist) else 0,
        com_dist_max=float(com_dist.max()) if len(com_dist) else 0,
        bound_fraction=float(np.mean(com_dist < 30.0)) if len(com_dist) else 0,
        bc_iface_n=len(bc_iface),
        calv_iface_n=len(calv_iface),
        iface_overlap=len(overlap),
        iface_jaccard=jaccard,
    )

results = []
for d in sorted(glob.glob(f"{P2}/sims_prebound/*/")):
    name = os.path.basename(d.rstrip('/'))
    if name == 'logs': continue
    r = analyze_one(name)
    if r:
        r['_crop'] = name.rsplit('__', 1)[0]
        r['_arm']  = name.rsplit('__', 1)[1]
        results.append(r)

print(f"  {len(results)} Komplexe analysiert")
from collections import defaultdict
by_arm = defaultdict(list)
for r in results: by_arm[r['_arm']].append(r)
print()
print(f"  {'Arm':16s} n   com_dist_med   bound_frac   iface_jaccard   calv_iface_med")
for arm in ('full','becraft_w1_lam0','ipsae_replace'):
    rs = by_arm.get(arm, [])
    if not rs: continue
    cmed = np.median([r['com_dist_med'] for r in rs])
    bf   = np.median([r['bound_fraction'] for r in rs])
    jc   = np.median([r['iface_jaccard'] for r in rs])
    civ  = np.median([r['calv_iface_n'] for r in rs])
    print(f"  {arm:16s} {len(rs):3d}  {cmed:8.1f}       {bf:5.2f}         {jc:5.2f}           {civ:5.1f}")

out = f"{P2}/phase2b_analysis.csv"
if results:
    keys = ['_crop','_arm','name','n_frames','rg_target_med','rg_binder_med',
            'com_dist_med','com_dist_min','com_dist_max','bound_fraction',
            'bc_iface_n','calv_iface_n','iface_overlap','iface_jaccard']
    with open(out,'w',newline='') as f:
        w=csv.DictWriter(f, fieldnames=keys); w.writeheader(); w.writerows(results)
    print(f"\n  -> {out}")
