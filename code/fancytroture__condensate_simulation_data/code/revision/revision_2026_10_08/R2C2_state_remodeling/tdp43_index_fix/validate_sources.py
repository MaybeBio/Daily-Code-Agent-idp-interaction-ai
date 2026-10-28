#!/usr/bin/env python3
"""Check residue-pair ordering against trajectory coordinates and fingerprint inputs."""
import hashlib
import json
from pathlib import Path
import sys
import mdtraj as md
import numpy as np

OUT=Path(__file__).resolve().parent
ROOT=OUT.parents[2]
sys.path.insert(0,str(ROOT/'project_figure/scripts/supplementary_figures'))
import shared as s


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    results=[]
    manifest=[]
    for family in s.FAMILIES:
        for role in s.ROLES:
            d=s.model_dir('TDP43',family,role)
            temp=s.TEMPS['TDP43'][role]
            top=d/'starts/seed_000.pdb'
            if family=='free':traj=d/f'outputs_{temp}_cut_stride1/tdp43_{temp}_s00_r00.dcd'
            else:traj=d/'extracted_500ns_bulk_chains/chain000_bulk_500ns.dcd'
            t=md.load_frame(str(traj),0,top=str(top))
            assert t.n_atoms==148
            i,j=np.triu_indices(148,3)
            fpath=d/'intermediate_tdp43/combined_features.npy'
            x=np.load(fpath,mmap_mode='r')
            # build_npz.py uses periodic=True unless --no-pbc was supplied.
            dist=md.compute_distances(t,np.column_stack((i,j)),periodic=True)[0]
            np.testing.assert_allclose(x[0],dist,atol=1e-5,rtol=0)
            results.append(dict(family=family,role=role,atoms=t.n_atoms,features=list(x.shape),
                                first_frame_distance_max_error=float(abs(x[0]-dist).max()),periodic=True,
                                first_pairs=np.column_stack((i,j))[:5].tolist()))
            small=[d/'avg_contacts.npy',d/'avg_mu.npy',d/'plot_fig3b.py',d/'scripts/build_npz.py',top]
            small+=sorted(d.glob('state_map_att*_raw_to_final.npy'))
            for p in small:
                manifest.append(dict(path=str(p.relative_to(ROOT)),bytes=p.stat().st_size,
                                     mtime_ns=p.stat().st_mtime_ns,sha256=sha(p)))
            for p in [traj,fpath]:
                manifest.append(dict(path=str(p.relative_to(ROOT)),bytes=p.stat().st_size,
                                     mtime_ns=p.stat().st_mtime_ns,sha256=None))
    previous=OUT/'source_manifest_before.json'
    if previous.exists():
        assert json.loads(previous.read_text())==manifest, 'Original source file changed during this repair'
    else:previous.write_text(json.dumps(manifest,indent=2)+'\n')
    (OUT/'trajectory_pair_validation.json').write_text(json.dumps(results,indent=2)+'\n')
    figures=[]
    for name in ['project_figure','project_figure_set_Tl_Th']:
        for p in sorted((ROOT/name).rglob('*')):
            if p.is_file() and p.suffix.lower() in ['.png','.pdf','.svg','.eps']:
                figures.append(dict(path=str(p.relative_to(ROOT)),bytes=p.stat().st_size,
                                    mtime_ns=p.stat().st_mtime_ns,sha256=sha(p)))
    guard=OUT/'original_figures_before.json'
    if guard.exists():
        assert json.loads(guard.read_text())==figures, 'An original figure changed; investigate without reverting user edits'
    else:guard.write_text(json.dumps(figures,indent=2)+'\n')
    print('Original figure guard:',len(figures),'files')
    print(json.dumps(results,indent=2))


if __name__=='__main__':main()
