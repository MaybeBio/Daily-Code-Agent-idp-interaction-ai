"""Add Rg and Nintra to the existing matched-residence comparison."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
from pathlib import Path
import json, csv, sys, time
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from MDAnalysis.coordinates.DCD import DCDReader
from MDAnalysis.lib.nsgrid import FastNS
from scipy.spatial.distance import pdist

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
BASE=ROOT/'revision_2026_10_08/R2C5_chain_selection/central60_2026_10_10'
OUT=HERE/'residence'
OUT.mkdir(exist_ok=True)
PLAN=json.loads((BASE/'plan.json').read_text())
STAGES=json.loads((ROOT/'revision_2026_09_30/selection_density/stage_data.json').read_text())['stages']
LAGS=np.arange(61)*.5

def extract(item):
    case,bid=item['case'],item['block']
    target=OUT/f'{case}_block{bid}_intramolecular.npz'
    if target.exists():return str(target)
    rg,ni,frames=[],[],[]
    start=time.monotonic()
    for part in item['parts']:
        path=ROOT/part['path'];assert path.stat().st_size==part['bytes']
        number=int(path.stem.split('part')[-1])
        with DCDReader(str(path)) as reader:
            cl=reader.n_atoms//200
            for j in range(0,10000,50):
                ts=reader[j];box=ts.dimensions[:3].astype(float)*.1
                xyz=ts.positions.astype(float)*.1
                raw=xyz.reshape(200,cl,3)
                delta=np.diff(raw,axis=1);delta-=np.rint(delta/box)*box
                whole=np.concatenate([raw[:,:1],raw[:,:1]+np.cumsum(delta,axis=1)],axis=1)
                centered=whole-whole.mean(axis=1,keepdims=True)
                rg.append(np.sqrt(np.mean(np.sum(centered**2,axis=2),axis=1)))
                pos=xyz%box
                pairs=FastNS(.6601,pos.astype(np.float32),np.array([*box,90,90,90],np.float32),True).self_search().get_pairs()
                keep=(pairs[:,0]//cl==pairs[:,1]//cl)&(np.abs(pairs[:,0]-pairs[:,1])>2)
                pairs=pairs[keep]
                disp=pos[pairs[:,0]]-pos[pairs[:,1]];disp-=np.rint(disp/box)*box
                pairs=pairs[np.einsum('ij,ij->i',disp,disp)<.66**2]
                counts=np.bincount(pairs[:,0]//cl,minlength=200)
                if j==0:
                    ii,jj=np.triu_indices(cl,k=3)
                    for chain in [0,17,99,199]:
                        d=raw[chain,ii]-raw[chain,jj];d-=np.rint(d/box)*box
                        assert np.count_nonzero(np.sum(d*d,axis=1)<.66**2)==counts[chain]
                ni.append(counts);frames.append((number-1)*10000+j)
        print(case,bid,path.name,'elapsed',round(time.monotonic()-start,1),flush=True)
    with np.load(BASE/'raw'/f'{case}_block{bid}.npz') as z:
        np.testing.assert_array_equal(frames,z['contact_frames'])
    np.savez_compressed(target,rg=np.array(rg),intra=np.array(ni),frames=np.array(frames))
    return str(target)

def acf(x,maxlag=60):
    y=x-x.mean(axis=0);v=np.mean(y*y,axis=0)
    assert np.all(v>0)
    return np.array([np.mean(y[:len(y)-k]*y[k:],axis=0)/v for k in range(maxlag+1)]).T

def cross(y,threshold=.2,dt=.5):
    hits=np.flatnonzero(y<=threshold)
    if not len(hits):return np.inf
    k=hits[0]
    return 0. if k==0 else (k-1)*dt+(threshold-y[k-1])/(y[k]-y[k-1])*dt

def summarize():
    blocks={}
    for item in PLAN:
        case,bid=item['case'],item['block']
        stage=next(s for s in STAGES if s['case']==case)
        metadata=json.loads(Path(stage['source_metadata']).read_text())
        ids=stage['inter_parent_chain_ids_reconstructed']
        with np.load(BASE/'raw'/f'{case}_block{bid}.npz') as z:
            short=z['bulk'].reshape(5,500,200).all(axis=1);frames=z['frames']
        with np.load(OUT/f'{case}_block{bid}_intramolecular.npz') as z:
            series={m:z[m] for m in ['rg','intra']}
        with np.load(OUT/f'{case}_block{bid}_intra_0p1ns.npz') as z:
            series['intra']=z['intra']
        groupcurves={m:{g:[] for g in ['original_selection','short_matched_original']} for m in series}
        for sub in range(5):
            start=int(frames[sub*500]);end=start+9999
            old=np.zeros(200,bool)
            for s in metadata:
                if s['chain_id'] in ids and start>=s['start'] and end<=s['end']:old[s['chain_id']]=True
            if not old.any():continue
            assert short[sub].any()
            for m,values in series.items():
                cadence=1000 if m=='intra' else 200
                cohort_ids=np.flatnonzero(old|short[sub])
                curves=acf(values[sub*cadence:(sub+1)*cadence,cohort_ids],300 if m=='intra' else 60)
                for g,mask in [('original_selection',old),('short_matched_original',short[sub])]:
                    groupcurves[m][g].append(curves[mask[cohort_ids]].mean(axis=0))
        blocks.setdefault(case,{})[bid]={m:{g:np.mean(v,axis=0) for g,v in groups.items() if v} for m,groups in groupcurves.items()}
    storage={};payload={}
    for case,bb in sorted(blocks.items()):
        for m in ['rg','intra']:
            for g in ['original_selection','short_matched_original']:
                bids=[b for b in sorted(bb) if g in bb[b][m]]
                arr=np.array([bb[b][m][g] for b in bids]);assert len(arr)>=3
                ix=np.random.default_rng(319).integers(len(arr),size=(4000,len(arr)))
                boot=arr[ix].mean(axis=1)
                storage[case,m,g]=(arr.mean(axis=0),boot,len(arr))
                payload[f'{case}_{m}_{g}']=arr
    rows=[]
    for system,tl,th in [('FUS',280,325),('A1',230,260),('TDP43',260,300)]:
        for m in ['rg','intra']:
            for g in ['original_selection','short_matched_original']:
                ml,bl,nl=storage[f'{system}_{tl}',m,g];mh,bh,nh=storage[f'{system}_{th}',m,g]
                bh=bh[np.random.default_rng(1337).permutation(4000)]
                dt=.1 if m=='intra' else .5
                l,h=cross(ml,dt=dt),cross(mh,dt=dt);ls=np.array([cross(x,dt=dt) for x in bl]);hs=np.array([cross(x,dt=dt) for x in bh])
                lower=np.where(np.isinf(hs),30.,hs)/ls
                upper=hs/np.where(np.isinf(ls),30.,ls)
                ci=[np.quantile(lower,.025,method='inverted_cdf'),np.quantile(upper,.975,method='inverted_cdf')]
                rows.append(dict(system=system,metric=m,group=g,quantity='tau',ratio=h/l,ci_low=ci[0],ci_high=ci[1],low_ns=l,high_ns=h,low_blocks=nl,high_blocks=nh,low_censored=int(np.isinf(ls).sum()),high_censored=int(np.isinf(hs).sum()),dt_ns=dt))
    np.savez_compressed(OUT/'block_curves.npz',**payload)
    with (OUT/'new_ratios.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    print(json.dumps(rows,indent=2),flush=True)

if __name__=='__main__':
    with ProcessPoolExecutor(max_workers=3) as pool:
        list(pool.map(extract,PLAN))
    summarize()
