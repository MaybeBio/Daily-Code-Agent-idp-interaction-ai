#!/usr/bin/env python3
"""Matched-length residence sensitivity; chains are not bootstrap replicas."""
import csv
import json
import os
from pathlib import Path
import sys
from concurrent.futures import ProcessPoolExecutor

os.environ.setdefault('MPLCONFIGDIR','/tmp/mpl-r2c5')
import numpy as np
from extract import HERE, ROOT, AUDIT, CONFIG, dump

DT=.5
LAGS=np.arange(61)*DT
NBOOT=4000
POPCOUNT=np.array([i.bit_count() for i in range(256)],dtype=np.uint8)


def table(name,rows):
    if not rows:
        return
    keys=list(dict.fromkeys(k for r in rows for k in r))
    with (HERE/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=keys)
        w.writeheader(); w.writerows(rows)


def eligible(flags):
    assert flags.shape==(2500,200)
    long=flags.all(axis=0)
    short=flags.reshape(5,500,200).all(axis=1)
    assert np.all(short[:,long])
    return long,short


def unwrap(com,box):
    inc=np.diff(com.astype(float),axis=0)
    inc-=np.rint(inc/box)*box
    inc-=inc.mean(axis=1,keepdims=True)
    return np.concatenate([np.zeros_like(com[:1]),np.cumsum(inc,axis=0)],axis=0)


def all_inside(flags,length):
    c=np.concatenate([np.zeros((1,flags.shape[1]),int),np.cumsum(~flags,axis=0)])
    return c[length:]-c[:-length]==0


def curves(counts):
    n=counts.sum(axis=2).astype(float)
    fluct=n-n.mean(axis=0)
    var=(fluct**2).mean(axis=0)
    present=counts>0
    packed=np.packbits(present,axis=2)
    npar=present.sum(axis=2)
    inter,partner=[],[]
    for lag in range(len(LAGS)):
        end=len(n)-lag
        cov=(fluct[:end]*fluct[lag:]).mean(axis=0)
        inter.append(np.divide(cov,var,out=np.full_like(var,np.nan),where=var>1e-12))
        overlap=POPCOUNT[packed[:end]&packed[lag:]].sum(axis=2)
        values=np.divide(overlap,npar[:end],out=np.full(overlap.shape,np.nan),where=npar[:end]>0)
        valid=np.isfinite(values)
        partner.append(np.divide(np.nansum(values,axis=0),valid.sum(axis=0),
            out=np.full(values.shape[1],np.nan),where=valid.sum(axis=0)>0))
    return np.array(inter).T,np.array(partner).T


def crossing(y,threshold):
    hit=np.flatnonzero(np.asarray(y)<=threshold)
    if not len(hit):
        return np.inf
    j=hit[0]
    if j==0:
        return 0.
    return float(LAGS[j-1]+(threshold-y[j-1])/(y[j]-y[j-1])*DT)


def per_block(path):
    case=path.stem.rsplit('_block',1)[0]
    bid=int(path.stem.rsplit('_block',1)[1])
    stage=next(x for x in json.loads((AUDIT/'stage_data.json').read_text())['stages'] if x['case']==case)
    metadata=json.loads(Path(stage['source_metadata']).read_text())
    ids=np.asarray(stage['inter_parent_chain_ids_reconstructed'],int)
    with np.load(path) as z:
        flags=z['bulk']; uniform=z['uniform_bulk']; frames=z['frames']
        unwrapped=unwrap(z['com_nm'],z['box_nm'])
        pc=z['partner_contacts']
        time=z['time_ns']
    old=np.zeros(200,bool);old[ids]=True
    mobility=[]
    for lag in [1,10,50]:
        k=int(round(lag/.2))
        delta=unwrapped[k:]-unwrapped[:-k]
        inside=all_inside(flags,k+1)
        active=np.zeros((len(delta),200),bool)
        for s in metadata:
            if s['chain_id'] in ids:
                active[:,s['chain_id']]=(frames[:-k]>=s['start'])&(frames[k:]<=s['end'])
        for dims,label in [(2,'xy'),(3,'xyz')]:
            msd=(delta[:,:,:dims]**2).sum(axis=2)
            perchain=msd.mean(axis=0)
            for name,mask in [('all',np.ones(200,bool)),('original_ids',old),('not_selected',~old)]:
                mobility.append(dict(case=case,block=bid,lag_ns=lag,coordinates=label,group=name,
                    msd_nm2=float(perchain[mask].mean()),n_chains=int(mask.sum())))
            for name,mask in [('all_bulk_origins',inside),('original_ids_bulk_origins',inside&old[None])]:
                denom=mask.sum(axis=0); good=denom>0
                means=(msd*mask).sum(axis=0)[good]/denom[good]
                mobility.append(dict(case=case,block=bid,lag_ns=lag,coordinates=label,group=name,
                    msd_nm2=float(means.mean()),n_chains=int(good.sum())))
            if active.any():
                mobility.append(dict(case=case,block=bid,lag_ns=lag,coordinates=label,group='original_window_matched',
                    msd_nm2=float(msd[active].mean()),n_chains=int(active.any(axis=0).sum())))
                mobility.append(dict(case=case,block=bid,lag_ns=lag,coordinates=label,group='all_window_matched',
                    msd_nm2=float((msd.mean(axis=1)*active.sum(axis=1)).sum()/active.sum()),n_chains=200))
    occupancy=[]
    for s in metadata:
        ix=(frames>=s['start'])&(frames<=s['end'])
        if ix.any():
            c=s['chain_id']
            occupancy.append(dict(case=case,block=bid,chain_id=c,metadata_id=s['id'],
                observed_frames=int(ix.sum()),observed_span_ns=float(time[ix][-1]-time[ix][0]),
                full_original_window=bool(ix.sum()==2500),
                source_rule_bulk_fraction=float(flags[ix,c].mean()),
                uniform_rule_bulk_fraction=float(uniform[ix,c].mean()),
                used_in_inter=bool(c in ids)))
    long,short=eligible(flags)
    counts=[]
    historical_groups=['original_selection','short_matched_original']
    blocks={m:{g:[] for g in ['residence_500','residence_100','short_only']+historical_groups} for m in ['inter','partner']}
    for sub in range(5):
        keep=short[sub]
        historical=np.zeros(200,bool)
        start=int(frames[sub*500]);end=start+9999
        for s in metadata:
            if s['chain_id'] in ids and start>=s['start'] and end<=s['end']:
                historical[s['chain_id']]=True
        subpc=pc[sub*200:(sub+1)*200]
        counts.append(dict(case=case,block=bid,subblock=sub,first_time_ns=float(time[sub*500]),
            n_500=int(long.sum()),n_100=int(keep.sum()),n_new=int((keep&~long).sum()),
            n_original=int(historical.sum()),
            mean_N_inter_500=float(subpc[:,long].sum(axis=2).mean()) if long.any() else '',
            mean_N_inter_100=float(subpc[:,keep].sum(axis=2).mean()) if keep.any() else '',
            n_uniform_500=int(uniform.all(axis=0).sum()),
            n_uniform_100=int(uniform[sub*500:(sub+1)*500].all(axis=0).sum())))
        # Compute only eligible tagged chains; partner identities still span all 200 chains.
        ix=np.flatnonzero(keep|historical)
        if not len(ix):
            raise RuntimeError(f'No 100 ns resident chains: {case} block{bid} sub{sub}')
        cc=pc[sub*200:(sub+1)*200,ix,:]
        ic,ps=curves(cc)
        for metric,arr in [('inter',ic),('partner',ps)]:
            groups=[('residence_500',long[ix]),('residence_100',keep[ix]),('short_only',keep[ix]&~long[ix])]
            if historical.any():
                groups += [('original_selection',historical[ix]),('short_matched_original',keep[ix])]
            for group,mask in groups:
                if mask.any():
                    if not np.isfinite(arr[mask]).all():
                        raise RuntimeError('Nonfinite curve for resident chain')
                    blocks[metric][group].append(arr[mask].mean(axis=0))
    result={m:{g:np.mean(v,axis=0) for g,v in gg.items() if len(v)==5 or (g in historical_groups and len(v)>0)} for m,gg in blocks.items()}
    print('ANALYZED',case,bid,'long',int(long.sum()),'short',short.sum(axis=1).tolist(),flush=True)
    return case,bid,result,mobility,occupancy,counts


def bootstrap(arr,seed):
    rng=np.random.default_rng(seed)
    idx=rng.integers(len(arr),size=(NBOOT,len(arr)))
    return arr[idx].mean(axis=1)


def main():
    paths=sorted((HERE/'raw').glob('*_block*.npz'))
    assert len(paths)==24, f'Need all 24 blocks, got {len(paths)}'
    results={}; mobility=[]; occupancy=[]; counts=[]
    with ProcessPoolExecutor(max_workers=4) as pool:
        for case,bid,r,m,o,c in pool.map(per_block,paths):
            results.setdefault(case,{})[bid]=r
            mobility+=m;occupancy+=o;counts+=c
    table('mobility_blocks.csv',mobility)
    table('vamp_sampled_residence.csv',occupancy)
    table('residence_counts.csv',counts)
    curve_rows=[];stats=[];ratio_rows=[];storage={}
    for case,bb in sorted(results.items()):
        for metric,th in [('inter',.2),('partner',.5)]:
            matched=[b for b in sorted(bb) if 'residence_500' in bb[b][metric]]
            for group in ['residence_500','residence_100','short_only','residence_100_all_blocks','original_selection','short_matched_original']:
                source='residence_100' if group=='residence_100_all_blocks' else group
                candidates=sorted(bb) if group in ['residence_100_all_blocks','original_selection','short_matched_original'] else matched
                bids=[b for b in candidates if source in bb[b][metric]]
                if not bids:
                    continue
                arr=np.array([bb[b][metric][source] for b in bids])
                boot=bootstrap(arr,319)
                mean=arr.mean(axis=0)
                storage[(case,metric,group)]=(mean,boot)
                low,high=np.quantile(boot,[.025,.975],axis=0)
                for j,t in enumerate(LAGS):
                    curve_rows.append(dict(case=case,metric=metric,group=group,lag_ns=t,
                        mean=float(mean[j]),ci_low=float(low[j]),ci_high=float(high[j]),n_blocks=len(bids)))
                for name,val,bs in [('tau',crossing(mean,th),np.array([crossing(y,th) for y in boot])),
                    ('I30',np.trapz(mean,LAGS),np.trapz(boot,LAGS,axis=1))]:
                    ci=np.quantile(bs,[.025,.975],method='inverted_cdf')
                    stats.append(dict(case=case,metric=metric,group=group,quantity=name,threshold=th if name=='tau' else '',
                        estimate_ns=float(val),ci_low_ns=float(ci[0]),ci_high_ns=float(ci[1]),
                        noncrossing_boot_fraction=float(np.isinf(bs).mean()),n_blocks=len(bids)))
            if len(matched)<4:
                print('WARNING missing 500ns resident controls',case,metric,matched,flush=True)
    for system,lo,hi in [('FUS','280','325'),('A1','230','260'),('TDP43','260','300')]:
        for metric,th in [('inter',.2),('partner',.5)]:
            for group in ['residence_500','residence_100','short_only','residence_100_all_blocks','original_selection','short_matched_original']:
                kl=(f'{system}_{lo}',metric,group); kh=(f'{system}_{hi}',metric,group)
                if kl not in storage or kh not in storage:
                    continue
                ml,bl=storage[kl];mh,bh=storage[kh]
                # Independent resampling across temperatures, paired across selection rules.
                bh=bh[np.random.default_rng(1337).permutation(NBOOT)]
                for name in ['tau','I30']:
                    if name=='tau':
                        l=crossing(ml,th);h=crossing(mh,th)
                        ls=np.array([crossing(y,th) for y in bl]);hs=np.array([crossing(y,th) for y in bh])
                    else:
                        l=np.trapz(ml,LAGS);h=np.trapz(mh,LAGS)
                        ls=np.trapz(bl,LAGS,axis=1);hs=np.trapz(bh,LAGS,axis=1)
                    if np.isfinite(ls).all() and np.isfinite(hs).all() and np.all(ls>0):
                        ci=np.quantile(hs/ls,[.025,.975])
                        status='finite'
                    else:
                        ci=[np.nan,np.nan];status='censored_or_nonpositive_denominator'
                    ratio_rows.append(dict(system=system,metric=metric,group=group,quantity=name,
                        ratio=float(h/l) if np.isfinite(l) and np.isfinite(h) and l>0 else '',
                        ci_low=float(ci[0]),ci_high=float(ci[1]),status=status))
    table('curves.csv',curve_rows);table('timescales.csv',stats);table('temperature_ratios.csv',ratio_rows)
    changes=[]
    for case in CONFIG:
        for metric,th in [('inter',.2),('partner',.5)]:
            longkey=(case,metric,'residence_500');shortkey=(case,metric,'residence_100')
            if longkey not in storage or shortkey not in storage:continue
            lm,lb=storage[longkey];sm,sb=storage[shortkey]
            for quantity in ['tau','I30']:
                if quantity=='tau':
                    lv=crossing(lm,th);sv=crossing(sm,th)
                    ls=np.array([crossing(x,th) for x in lb]);ss=np.array([crossing(x,th) for x in sb])
                else:
                    lv=np.trapz(lm,LAGS);sv=np.trapz(sm,LAGS)
                    ls=np.trapz(lb,LAGS,axis=1);ss=np.trapz(sb,LAGS,axis=1)
                ci=np.quantile(ss/ls,[.025,.975]) if np.isfinite(ls).all() and np.isfinite(ss).all() and np.all(ls>0) else [np.nan,np.nan]
                changes.append(dict(case=case,metric=metric,quantity=quantity,
                    short_over_long=float(sv/lv) if np.isfinite(lv) and lv>0 else '',ci_low=float(ci[0]),ci_high=float(ci[1])))
    table('paired_residence_effects.csv',changes)
    msummary=[]
    for case in CONFIG:
        for lag in [1,10,50]:
            for dims in ['xy','xyz']:
                for selected,reference in [('original_ids','all'),('original_ids_bulk_origins','all_bulk_origins'),
                    ('original_window_matched','all_window_matched')]:
                    a={r['block']:r['msd_nm2'] for r in mobility if r['case']==case and r['lag_ns']==lag and
                        r['coordinates']==dims and r['group']==selected}
                    b={r['block']:r['msd_nm2'] for r in mobility if r['case']==case and r['lag_ns']==lag and
                        r['coordinates']==dims and r['group']==reference}
                    bids=sorted(set(a)&set(b))
                    if not bids:continue
                    arr=np.array([[a[x],b[x]] for x in bids]);boot=bootstrap(arr,731)
                    ratio=arr[:,0].mean()/arr[:,1].mean();ci=np.quantile(boot[:,0]/boot[:,1],[.025,.975])
                    msummary.append(dict(case=case,lag_ns=lag,coordinates=dims,selected=selected,reference=reference,
                        selected_msd_nm2=float(arr[:,0].mean()),all_msd_nm2=float(arr[:,1].mean()),ratio=float(ratio),
                        ci_low=float(ci[0]),ci_high=float(ci[1]),n_blocks=len(bids)))
    table('mobility_summary.csv',msummary)
    dump(HERE/'analysis_completed.json',dict(status='complete',blocks=24,conditions=6,chains_per_frame=200,
        bootstrap_unit='500 ns temporal block, all chains and five subblocks kept together',bootstrap_draws=NBOOT,
        uncertainty_limitation='Four blocks from one slab per condition; not independent simulation replicas.',
        contact_cutoff_nm=.66,acf_length_ns=100,max_lag_ns=30,residence_ns=[500,100],
        populations_or_models_retrained=False,original_files_modified=False))


if __name__=='__main__':
    main()
