"""Matched definitions and explicit intervals for all four dynamical observables."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
from pathlib import Path
import numpy as np
import json,csv
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
OLD=ROOT/'revision_2026_09_30/interchain_recompute/analysis'
N=5000

def cross(t,y,q):
    y=np.atleast_2d(y);ok=y<=q;idx=ok.argmax(axis=1);found=ok.any(axis=1)
    assert not np.any(idx[found]==0)
    rows=np.arange(len(y));before=np.maximum(idx-1,0)
    with np.errstate(divide='ignore',invalid='ignore'):
        out=t[before]+(q-y[rows,before])/(y[rows,idx]-y[rows,before])*(t[idx]-t[before])
    return np.where(found,out,np.inf)

def load(system,temp,metric):
    if metric in ['rg','intra']:
        with np.load(HERE/'intramolecular'/f'{system}_{temp}.npz') as z:
            t=z[metric+'_time'];units=z[metric+'_units']
        mean=units.mean(axis=0)
        seed=20261015+temp+{'FUS':1000,'A1':2000,'TDP43':3000}[system]
        weights=np.random.default_rng(seed).multinomial(len(units),np.ones(len(units))/len(units),size=N)
        boot=weights@units/len(units)
        method='selected_trajectory_bootstrap'
    else:
        with np.load(OLD/f'{system}_{temp}_{metric}_curves.npz') as z:
            t=z['time_ns'];mean=z['mean'];boot=z['block_250_boot']
        assert len(boot)==N
        method='synchronized_250ns_origin_block_bootstrap'
    return t,mean,boot,method

def main():
    rows=[]
    for system,tl,th in [('FUS',280,325),('A1',230,260),('TDP43',260,300)]:
        for metric in ['rg','intra','inter','partner']:
            t,l,bl,method=load(system,tl,metric);ht,h,bh,hm=load(system,th,metric)
            np.testing.assert_array_equal(t,ht);assert method==hm
            for kind,values in [('threshold',[.4,.45,.5,.55,.6] if metric=='partner' else [.1,.15,.2,.25,.3]),('integral',[10.,20.,30.])]:
                for value in values:
                    if kind=='threshold':
                        lv=cross(t,l,value)[0];hv=cross(t,h,value)[0]
                        ls=cross(t,bl,value);hs=cross(t,bh,value)
                        horizon=t[-1]
                        lower=np.where(np.isinf(hs),horizon,hs)/ls
                        upper=hs/np.where(np.isinf(ls),horizon,ls)
                    else:
                        k=np.searchsorted(t,value)+1
                        lv=np.trapz(l[:k],t[:k]);hv=np.trapz(h[:k],t[:k])
                        ls=np.trapz(bl[:,:k],t[:k],axis=1);hs=np.trapz(bh[:,:k],t[:k],axis=1)
                        assert np.all(ls>0)
                        lower=upper=hs/ls;horizon=value
                    assert np.isfinite(lv) and np.isfinite(hv)
                    ci=[np.quantile(lower,.025,method='inverted_cdf'),np.quantile(upper,.975,method='inverted_cdf')]
                    assert np.isfinite(ci).all(),(system,metric,kind,value,ci)
                    rows.append(dict(system=system,metric=metric,kind=kind,value=value,ratio=hv/lv,ci_low=ci[0],ci_high=ci[1],low_ns=lv,high_ns=hv,low_censored=int(np.isinf(ls).sum()),high_censored=int(np.isinf(hs).sum()),horizon_ns=horizon,bootstrap=N,method=method))
            print(system,metric,'complete',flush=True)
    with (HERE/'sensitivity_ratios.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    (HERE/'sensitivity_ratios.json').write_text(json.dumps(rows,indent=2)+'\n')

if __name__=='__main__':main()
