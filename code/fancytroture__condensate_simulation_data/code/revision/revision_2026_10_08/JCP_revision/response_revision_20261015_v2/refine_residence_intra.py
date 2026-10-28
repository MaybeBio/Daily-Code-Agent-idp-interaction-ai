"""Resolve fast Nintra decay at the original intrachain 0.1 ns cadence."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
from pathlib import Path
import json,time
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from MDAnalysis.coordinates.DCD import DCDReader
from MDAnalysis.lib.nsgrid import FastNS
from complete_residence import ROOT,BASE,OUT,PLAN,STAGES

def process(item):
    case,bid=item['case'],item['block'];target=OUT/f'{case}_block{bid}_intra_0p1ns.npz'
    if target.exists():return str(target)
    stage=next(s for s in STAGES if s['case']==case)
    meta=json.loads(Path(stage['source_metadata']).read_text());oldids=stage['inter_parent_chain_ids_reconstructed']
    with np.load(BASE/'raw'/f'{case}_block{bid}.npz') as z:
        frames=z['frames'];short=z['bulk'].reshape(5,500,200).all(axis=1)
    selected=[]
    for sub in range(5):
        start=int(frames[sub*500]);end=start+9999;old=np.zeros(200,bool)
        for s in meta:
            if s['chain_id'] in oldids and start>=s['start'] and end<=s['end']:old[s['chain_id']]=True
        selected.append(np.flatnonzero(short[sub]|old) if old.any() else np.array([],int))
    with np.load(OUT/f'{case}_block{bid}_intramolecular.npz') as z:coarse=z['intra']
    values=np.full((5000,200),np.nan);start=time.monotonic()
    for sub,part in enumerate(item['parts']):
        ids=selected[sub]
        if not len(ids):continue
        path=ROOT/part['path']
        with DCDReader(str(path)) as reader:
            cl=reader.n_atoms//200
            for j in range(0,10000,10):
                dest=sub*1000+j//10
                if j%50==0:
                    values[dest,ids]=coarse[sub*200+j//50,ids];continue
                ts=reader[j];box=ts.dimensions[:3].astype(float)*.1
                pos=(ts.positions.reshape(200,cl,3)[ids].reshape(-1,3).astype(float)*.1)%box
                pairs=FastNS(.6601,pos.astype(np.float32),np.array([*box,90,90,90],np.float32),True).self_search().get_pairs()
                pairs=pairs[(pairs[:,0]//cl==pairs[:,1]//cl)&(np.abs(pairs[:,0]-pairs[:,1])>2)]
                delta=pos[pairs[:,0]]-pos[pairs[:,1]];delta-=np.rint(delta/box)*box
                pairs=pairs[np.einsum('ij,ij->i',delta,delta)<.66**2]
                values[dest,ids]=np.bincount(pairs[:,0]//cl,minlength=len(ids))
        print(case,bid,'sub',sub,'chains',len(ids),'seconds',round(time.monotonic()-start,1),flush=True)
    np.savez_compressed(target,intra=values,dt_ns=.1)
    return str(target)

if __name__=='__main__':
    with ProcessPoolExecutor(max_workers=6) as pool:list(pool.map(process,PLAN))
