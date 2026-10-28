#!/usr/bin/env python3
"""Reconstruct per-seed contacts and repair residue indexing, without training."""
import ast
import hashlib
import json
import os
from pathlib import Path
import sys
import time

os.environ.setdefault('OPENBLAS_NUM_THREADS', '4')
os.environ.setdefault('OMP_NUM_THREADS', '4')
import h5py
import numpy as np
from scipy.optimize import linear_sum_assignment
from scipy.sparse import csr_matrix

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[2]
sys.path.insert(0, str(ROOT/'project_figure/scripts/supplementary_figures'))
import shared as s
SAM = ROOT/'revision_2026_09_30/state_sampling'
BLOCK = 4096


def log(*args):
    print(time.strftime('%H:%M:%S'), *args, flush=True)


def make_maps(v, nres):
    gap = 0 if nres == 146 else 2
    i,j = np.triu_indices(nres, gap+1)
    assert v.shape[-1] == len(i) == 10585
    a = np.zeros(v.shape[:-1]+(nres,nres), dtype=np.float32)
    a[...,i,j] = v
    a[...,j,i] = v
    return a


def alignment_function(path):
    # Reuse only the original pure alignment functions, without GPU/training imports.
    wanted = {'_flat_contact_for_dist','cost_ssd_between_sets','_pick_central_attempt','covamp_within_system_align'}
    tree = ast.parse(path.read_text())
    funcs = [n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name in wanted]
    assert len(funcs) == len(wanted)
    namespace = {'np':np, 'linear_sum_assignment':linear_sum_assignment}
    exec(compile(ast.Module(body=funcs, type_ignores=[]), str(path), 'exec'), namespace)
    return namespace['covamp_within_system_align']


def raw_mus(d, attempts, q):
    mus=[]
    for a,att in enumerate(attempts):
        path=d/f'intermediate_tdp43/data_combined_att{att}.hdf5'
        mu=q[:,a*5:(a+1)*5].sum(0,dtype=np.float64)
        if path.exists():
            with h5py.File(path,'r') as f:
                key=f'combined/{att}/5/50/mu'
                if key in f:mu=np.abs(f[key][:].astype(float))
        mus.append(mu/mu.sum())
    return np.array(mus)


def saved_average_error(v, old_perms, saved):
    maps=make_maps(v, saved.shape[-1])
    mean=np.mean([x[p] for x,p in zip(maps,old_perms)],axis=0)
    return float(abs(mean-saved).max()), float(np.sqrt(np.mean((mean-saved)**2)))


def robust_prototypes(binary,q,proto):
    nf,npairs=binary.shape
    nm=q.shape[1]//5
    labels=q.reshape(nf,nm,5).argmax(-1)
    norms=(proto*proto).sum(-1)
    distances=np.empty((nf,nm),dtype=np.float32)
    for start in range(0,nf,BLOCK):
        end=min(nf,start+BLOCK)
        c=csr_matrix(binary[start:end])
        dots=c@proto.T
        columns=np.arange(nm)[None,:]*5+labels[start:end]
        distances[start:end]=np.asarray(c.sum(-1),dtype=np.float32)-2*np.take_along_axis(dots,columns,axis=1)+norms[columns]
        if start//BLOCK%100==0:log('Robust distances',end,'/',nf)
    thresholds=np.full((nm,5),np.inf,dtype=np.float32)
    for m in range(nm):
        for k in range(5):
            vals=distances[labels[:,m]==k,m]
            if len(vals)>=500:thresholds[m,k]=np.percentile(vals,90.)
    sums=np.zeros((nm*5,npairs),dtype=np.float64)
    denom=np.zeros(nm*5,dtype=np.float64)
    for start in range(0,nf,BLOCK):
        end=min(nf,start+BLOCK)
        lab=labels[start:end]
        cut=thresholds[np.arange(nm)[None,:],lab]
        use=(distances[start:end]<=cut)&np.isfinite(cut)
        weights=q[start:end].reshape(-1,nm,5).copy()
        weights*=((np.arange(5)[None,None,:]==lab[:,:,None])&use[:,:,None])
        weights=weights.reshape(-1,nm*5)
        sums+=(csr_matrix(binary[start:end]).T@weights.astype(np.float64)).T
        denom+=weights.sum(0,dtype=np.float64)
        if start//BLOCK%100==0:log('Robust averages',end,'/',nf)
    result=proto.copy()
    valid=denom>1e-12
    result[valid]=(sums[valid]/denom[valid,None]).astype(np.float32)
    return result.reshape(nm,5,npairs),thresholds


def condition(family,role):
    d=s.model_dir('TDP43',family,role)
    tag=f'TDP43_{family}_{s.TEMPS["TDP43"][role]}'
    dest=OUT/tag
    dest.mkdir(exist_ok=True)
    src=SAM/f'{tag}_sufficient_statistics.npz'
    with np.load(src) as z:
        attempts=z['attempts'].astype(int)
        lengths=z['lengths'].astype(int)
    old_raw_to_final=np.stack([np.load(d/f'state_map_att{att}_raw_to_final.npy') for att in attempts])
    old_perms=np.argsort(old_raw_to_final,axis=1)
    saved=np.load(d/'avg_contacts.npy')
    if role=='Th':
        # No index change at high temperature: retain the saved mapping exactly.
        assert saved.shape==(5,148,148)
        for gap in [0,1,2]:
            idx=np.arange(148-gap)
            assert np.all(saved[:,idx,idx+gap]==0)
        avg_mu=np.load(d/'avg_mu.npy')
        identity=np.tile(np.arange(5),(len(attempts),1))
        info=dict(reference_attempt=int((d/'ref_att.txt').read_text()),changed_models=0,
                  new_to_old_by_attempt=identity.tolist(),population=avg_mu.tolist())
        report=dict(condition=tag,mode='unchanged_high_temperature_source',saved_shape=list(saved.shape),
                    rebuild_performed=False,legacy=info,fixed=info,
                    note='Correct 148-residue indexing already present; no alignment change requested here')
        for label in ['legacy','fixed']:
            np.savez_compressed(dest/f'{label}_alignment.npz',avg_contacts=saved,avg_mu=avg_mu,
                                new_to_old=identity,raw_to_final=old_raw_to_final,attempts=attempts,
                                ref=info['reference_attempt'])
        (dest/'validation.json').write_text(json.dumps(report,indent=2)+'\n')
        log(tag,'retained unchanged: correct 148-residue source')
        return report
    cache=dest/'raw_per_seed_contacts.npz'
    if cache.exists():
        with np.load(cache) as z:
            raw=z['raw'];mu=z['mu'];mode=str(z['mode'].item())
        log(tag,'using completed per-seed cache',mode)
    else:
        features=np.load(d/'intermediate_tdp43/combined_features.npy',mmap_mode='r')
        nf,npairs=features.shape
        assert nf==lengths.sum() and npairs==10585
        q=np.empty((nf,len(attempts)*5),dtype=np.float32)
        for a,att in enumerate(attempts):
            start=0
            qd=d/f'intermediate_tdp43/qcache_combined_n5_att{att}_lag50'
            for ti,length in enumerate(lengths):
                arr=np.load(qd/f'traj{ti}.npy',mmap_mode='r')
                assert arr.shape==(length,5)
                q[start:start+length,a*5:(a+1)*5]=arr
                start+=length
        assert np.isfinite(q).all() and q.min()>=-1e-7
        np.testing.assert_allclose(q.reshape(nf,-1,5).sum(-1),1.,atol=1e-5)
        mu=raw_mus(d,attempts,q)
        binary=np.empty(features.shape,dtype=bool)
        sums=np.zeros((len(attempts)*5,npairs),dtype=np.float64)
        for start in range(0,nf,BLOCK):
            end=min(nf,start+BLOCK)
            f=features[start:end]
            assert np.isfinite(f).all() and np.all(f>=0)
            binary[start:end]=f<.66
            sums+=(csr_matrix(binary[start:end]).T@q[start:end].astype(np.float64)).T
            if start//BLOCK%100==0:log(tag,'Raw contact averages',end,'/',nf)
        proto=(sums/q.sum(0,dtype=np.float64)[:,None]).astype(np.float32)
        unfiltered=proto.reshape(len(attempts),5,npairs)
        err=saved_average_error(unfiltered,old_perms,saved)
        log(tag,'unfiltered saved-average max/rms error',err)
        np.savez_compressed(dest/'unfiltered_per_seed_contacts.npz',raw=unfiltered,mu=mu,attempts=attempts)
        if err[0]<2e-5:
            raw=unfiltered;mode='none'
        else:
            raw,thresholds=robust_prototypes(binary,q,proto)
            mode='robust_keep90'
            err=saved_average_error(raw,old_perms,saved)
            log(tag,'robust saved-average max/rms error',err)
            np.save(dest/'robust_thresholds.npy',thresholds)
            assert err[0]<2e-4, ('Cannot reproduce old maps; do not silently change filtering',tag,err)
        np.savez_compressed(cache,raw=raw,mu=mu,mode=mode,attempts=attempts)
        del q,binary,features,sums
    align=alignment_function(d/'plot_fig3b.py')
    report=dict(condition=tag,mode=mode,saved_shape=list(saved.shape),rebuild_performed=True,
                old_average_error=saved_average_error(raw,old_perms,saved))
    for label,nres in [('legacy',saved.shape[-1]),('fixed',148)]:
        contacts=make_maps(raw,nres)
        valid=[(int(att),x,p) for att,x,p in zip(attempts,contacts,mu)]
        perm,ref=align(valid,5,mask_gap=4,verbose=False)
        order=np.argsort(np.mean([p[perm[int(att)]] for att,p in zip(attempts,mu)],axis=0))[::-1]
        perms=np.array([perm[int(att)][order] for att in attempts])
        avg=np.mean([x[p] for x,p in zip(contacts,perms)],axis=0)
        avg_mu=np.mean([x[p] for x,p in zip(mu,perms)],axis=0)
        # new final -> raw -> old final, used to permute saved sufficient statistics.
        new_to_old=np.take_along_axis(old_raw_to_final,perms,axis=1)
        report[label]=dict(reference_attempt=int(ref),changed_models=int(np.any(new_to_old!=np.arange(5),axis=1).sum()),
                           new_to_old_by_attempt=new_to_old.tolist(),population=avg_mu.tolist())
        if label=='legacy':
            report[label]['saved_average_max_error']=float(abs(avg-saved).max())
            assert np.array_equal(new_to_old,np.tile(np.arange(5),(len(attempts),1))), ('Legacy alignment not reproduced',tag,report[label])
        np.savez_compressed(dest/f'{label}_alignment.npz',avg_contacts=avg,avg_mu=avg_mu,
                            new_to_old=new_to_old,raw_to_final=np.argsort(perms,axis=1),attempts=attempts,ref=ref)
    (dest/'validation.json').write_text(json.dumps(report,indent=2)+'\n')
    log(tag,'DONE. Changed seed alignments:',report['fixed']['changed_models'])
    return report


def main():
    reports=[]
    for family in ['free','condensate']:
        for role in ['Tl','Th']:
            reports.append(condition(family,role))
    (OUT/'contact_rebuild_validation.json').write_text(json.dumps(reports,indent=2)+'\n')
    log('ALL CONTACT REBUILDS COMPLETE')


if __name__=='__main__':
    main()
