"""Recover per-trajectory intramolecular curves from the existing DCDs."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
from pathlib import Path
import json,time,csv
from concurrent.futures import ProcessPoolExecutor
import numpy as np
import mdtraj as md
from mdtraj.formats import DCDTrajectoryFile

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
OUT=HERE/'intramolecular';OUT.mkdir(exist_ok=True)
BASES={'FUS':'fus/fus_protein_single_simulation','A1':'A1/A1_protein_single_simulation','TDP43':'TDP43/tdp43_protein_single_simulation'}
RGDIRS={'FUS':'paper_thermo_figures_2temps_280_325_vampnet','A1':'paper_thermo_figures_a1_4p17_relaxed_mu10','TDP43':'paper_thermo_figures_tdp43_4p17_relaxed_mu10'}

def acf(x,maxlag,pair):
    x=np.asarray(x,float);x-=x.mean();n=len(x)
    f=np.fft.rfft(x,n=2*n);c=np.fft.irfft(f*f.conjugate(),n=2*n)[:maxlag+1]
    c/=c[0]
    if pair:c*=n/(n-np.arange(maxlag+1))
    return c

def compute(task):
    system,temp=task;target=OUT/f'{system}_{temp}.npz'
    if target.exists():return str(target)
    base=ROOT/BASES[system];d=base/'condensate'/f'tica_100dimension_10frames_{temp}K'
    paths=sorted((d/'extracted_500ns_bulk_chains').glob('chain*_bulk_500ns.dcd'))
    lengths=np.load(d/f'intermediate_{system.lower()}'/'combined_tica_lengths.npy').astype(int)
    features=np.load(d/f'intermediate_{system.lower()}'/'combined_features.npy',mmap_mode='r')
    assert len(paths)==len(lengths)
    top=md.load(str(next((d/'starts').glob('seed_000.pdb')))).topology
    ii,jj=np.triu_indices(top.n_atoms,k=3);pairs=np.array([ii,jj]).T
    rgcurves=[];incurves=[];series={};checks=[];offset=0;start=time.monotonic()
    for unit,(path,length) in enumerate(zip(paths,lengths)):
        rg=[];counts=[]
        with DCDTrajectoryFile(str(path)) as f:
            assert len(f)==length
            done=0
            while done<length:
                xyz,cell,angle=f.read(n_frames=min(5000,length-done))
                tr=md.Trajectory(xyz*.1,top,unitcell_lengths=cell*.1,unitcell_angles=angle)
                rg.extend(md.compute_rg(tr).astype(float))
                distances=md.compute_distances(tr[::10],pairs,periodic=True)
                cc=np.count_nonzero(distances<.66,axis=1)
                counts.extend(cc)
                for k in [0,len(cc)-1]:
                    index=offset+done+10*k
                    expected=np.count_nonzero(features[index]<.66)
                    if expected!=cc[k]:
                        raise ValueError(('Feature/DCD contact mismatch',system,temp,path.name,index,int(expected),int(cc[k])))
                done+=len(xyz)
        rg=np.asarray(rg);counts=np.asarray(counts)
        series[f'intra_{unit}']=counts
        if length==50000:
            rgcurves.append(acf(rg,6000,False))
        incurves.append(acf(counts,600,True))
        checks.append({'file':str(path),'frames':int(length),'feature_checks':'two frames per 5000-frame chunk agree exactly'})
        offset+=length
        if unit%5==4:print(system,temp,unit+1,'elapsed',round(time.monotonic()-start,1),flush=True)
    rgcurves=np.array(rgcurves);incurves=np.array(incurves)
    with np.load(base/RGDIRS[system]/'figure1_rg_ree_relaxation'/f'cache_condensate_{temp}K.npz') as old:
        rg_error=float(np.max(np.abs(rgcurves.mean(axis=0)-np.interp(np.arange(6001)*.01,old['time_axis'],old['acf_rg']))))
    with np.load(ROOT/f'project_figure_set_Tl_Th/source_data/intrachain_relaxation/intrachain_{system}_condensate_{temp}K.npz') as old:
        intra_error=float(np.max(np.abs(incurves[:,:301].mean(axis=0)-old['acf_mean'])));
        assert len(incurves)==int(old['n_trajectories'])
    assert rg_error<2e-5,('Rg mismatch',system,temp,rg_error)
    assert intra_error<2e-5,('Intra mismatch',system,temp,intra_error)
    np.savez_compressed(target,rg_units=rgcurves,intra_units=incurves,rg_time=np.arange(6001)*.01,intra_time=np.arange(601)*.1,**series)
    target.with_suffix('.json').write_text(json.dumps({'system':system,'temperature':temp,'rg_max_abs_error':rg_error,'intra_max_abs_error':intra_error,'sources':checks},indent=2))
    print('DONE',system,temp,rg_error,intra_error,flush=True)
    return str(target)

if __name__=='__main__':
    tasks=[('FUS',280),('FUS',325),('A1',230),('A1',260),('TDP43',260),('TDP43',300)]
    with ProcessPoolExecutor(max_workers=2) as pool:list(pool.map(compute,tasks))
