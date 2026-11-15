#!/usr/bin/env python3
"""CALVADOS phase 3: target-binder proximity for the experimentally determined controls."""
import os, glob
import numpy as np
import MDAnalysis as mda
BASE="/capstor/store/cscs/swissai/ab011/jkuhn/idp_bindcraft/calvados_sanity/phase3_pos_control/sims"
rows=[]
for ctrl in sorted(os.path.basename(d) for d in glob.glob(f"{BASE}/*") if os.path.isdir(d) and os.path.basename(d)!="logs"):
    data=f"{BASE}/{ctrl}/data"
    top=f"{data}/{ctrl}_prebound_com_top.pdb"; dcd=f"{data}/{ctrl}_prebound_com_traj.dcd"
    trg=f"{data}/{ctrl}_prebound_{ctrl}_target_rg.npy"; bnd=f"{data}/{ctrl}_prebound_{ctrl}_binder_rg.npy"
    if not (os.path.isfile(top) and os.path.isfile(dcd)):
        print(f"{ctrl}: missing com traj"); continue
    u=mda.Universe(top,dcd)
    coms=np.array([u.atoms.positions[:2].copy() for _ in u.trajectory])
    cd=np.linalg.norm(coms[:,0,:]-coms[:,1,:],axis=1)
    rt=np.load(trg) if os.path.isfile(trg) else np.array([]); rb=np.load(bnd) if os.path.isfile(bnd) else np.array([])
    rows.append(dict(ctrl=ctrl,n=len(cd),
        com_med=float(np.median(cd)),com_min=float(cd.min()),com_max=float(cd.max()),
        bound_frac=float(np.mean(cd<30.0)),
        rg_t=float(np.median(rt)) if len(rt) else None, rg_b=float(np.median(rb)) if len(rb) else None))
print(f"{'control':8} {'n':>4} {'com_med(A)':>11} {'com_min':>8} {'bound_frac(<30A)':>16}")
for r in rows:
    print(f"{r['ctrl']:8} {r['n']:>4} {r['com_med']:>11.1f} {r['com_min']:>8.1f} {r['bound_frac']:>16.2f}")
import json; json.dump(rows,open(os.path.dirname(__file__)+"/phase3_results.json","w"),indent=1)
print("\nReference: phase2/2b designed complexes had bound_fraction ~0.00, com_dist_med ~190-224 A.")
print("wrote phase3_results.json")
