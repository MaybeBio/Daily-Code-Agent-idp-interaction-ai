#!/usr/bin/env python3
"""Read-only, physical-time extraction for the R2C5 residence control."""
import argparse
from concurrent.futures import ProcessPoolExecutor
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import time
import warnings

os.environ.setdefault('MPLCONFIGDIR', '/tmp/mpl-r2c5')
import numpy as np
from MDAnalysis.coordinates.DCD import DCDReader
from MDAnalysis.lib.nsgrid import FastNS
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
AUDIT = ROOT / 'revision_2026_09_30/selection_density'
SELECTOR = ROOT.parent / 'preprocessing/extract_bulk_chain_trajectories.py'
spec = importlib.util.spec_from_file_location('selector', SELECTOR)
selector = importlib.util.module_from_spec(spec)
spec.loader.exec_module(selector)
warnings.filterwarnings('ignore', category=DeprecationWarning)
CONFIG = {'FUS_280': (4.,20,80), 'FUS_325': (5.,10,90),
          'A1_230': (4.,20,80), 'A1_260': (4.5,20,80),
          'TDP43_260': (4.,20,80), 'TDP43_300': (4.,20,80)}


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def dump(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix('.tmp.json')
    tmp.write_text(json.dumps(obj, indent=2, allow_nan=False)+'\n')
    tmp.replace(path)


def bulk_flags(com, box, params):
    cutoff, low, high = params
    pairs = FastNS(cutoff, np.asarray(com, np.float32),
                   np.array([*box,90,90,90], np.float32), True).self_search().get_pairs()
    ds = selector.DSU(len(com))
    for a,b in pairs:
        ds.union(int(a),int(b))
    largest = np.asarray(ds.largest_component())
    flags = np.zeros(len(com),bool)
    if len(largest)>=30:
        split, lo, hi = selector.solve_pbc_z_range(com[largest,2],box[2],low,high)
        zz = (com[largest,2]-split)%box[2]
        flags[largest] = (zz>=lo)&(zz<=hi)
    return flags


def whole_com(xyz, box, mass):
    raw = xyz.reshape(200,-1,3)
    delta = np.diff(raw,axis=1)
    delta -= np.rint(delta/box)*box
    whole = np.concatenate([raw[:,:1],raw[:,:1]+np.cumsum(delta,axis=1)],axis=1)
    return ((whole*mass[None,:,None]).sum(axis=1)/mass.sum())%box


def contacts(xyz,box,chain_len,cutoff=.66):
    pos = np.asarray(xyz,dtype=float)%box
    # Expand the float32 neighbor search, then apply the cutoff in float64.
    pairs = FastNS(cutoff+1e-4,pos.astype(np.float32),
                   np.array([*box,90,90,90],np.float32),True).self_search().get_pairs()
    inter = (pairs[:,0]//chain_len)!=(pairs[:,1]//chain_len)
    pairs = pairs[inter]
    delta = pos[pairs[:,0]]-pos[pairs[:,1]]
    delta -= np.rint(delta/box)*box
    pairs = pairs[np.einsum('ij,ij->i',delta,delta)<=cutoff**2]
    ci = pairs//chain_len
    counts = np.bincount(ci[:,0]*200+ci[:,1],minlength=40000).reshape(200,200)
    counts = counts+counts.T
    assert not counts.diagonal().any()
    assert counts.sum()==2*len(pairs)
    return counts.astype(np.int32)


def check_counts(xyz,box,chain_len,result):
    pos = xyz%box
    tree = cKDTree(pos,boxsize=box)
    for c in [0,17,99,199]:
        neigh = tree.query_ball_point(pos[c*chain_len:(c+1)*chain_len],.66)
        counts = np.bincount(np.concatenate(neigh).astype(int)//chain_len,minlength=200)
        counts[c]=0
        np.testing.assert_array_equal(counts,result[c])


def plan():
    inventory = json.loads((AUDIT/'inventory.json').read_text())
    result = []
    for entry in inventory:
        parts = entry['parts']
        last_full = max(i+1 for i,p in enumerate(parts) if p['frames']==10000 and 'read_error' not in p)
        starts = np.rint(np.linspace(51,last_full-4,4)).astype(int)
        for bid,first in enumerate(starts):
            ps = parts[first-1:first+4]
            assert len(ps)==5 and all(p['frames']==10000 for p in ps)
            result.append(dict(case=entry['case'],block=bid,first_part=int(first),
                simulation=entry['simulation'],parts=ps,selection_params=CONFIG[entry['case']],
                selection_dt_ns=.2,contact_dt_ns=.5,raw_dt_ns=.01,duration_ns=500.,
                cutoff_nm=.66,source_selection_sha256=digest(SELECTOR),script_sha256=digest(__file__)))
    dump(HERE/'plan.json',result)
    return result


def process(item):
    from openabc.forcefields.parsers.hps_parser import _hps_amino_acid_mass_dict
    case, bid = item['case'],item['block']
    out = HERE/'raw'/f'{case}_block{bid}.npz'
    ident = hashlib.sha256(json.dumps(item,sort_keys=True).encode()).hexdigest()
    if out.exists():
        with np.load(out) as z:
            if str(z['identity'])!=ident:
                raise RuntimeError(f'Incompatible cache {out}')
        return str(out)
    sim = Path(item['simulation'])
    pdb = next(sim.glob('*_slab_vampnet_start.pdb'))
    names = [line[17:20].strip() for line in pdb.read_text().splitlines() if line.startswith(('ATOM  ','HETATM'))]
    masses = np.array([_hps_amino_acid_mass_dict[x] for x in names]).reshape(200,-1)
    assert np.all(masses==masses[0])
    cl = masses.shape[1]
    frames, times, coms, flags, uniform = [],[],[],[],[]
    cf,ct,pc = [],[],[]
    started = time.monotonic()
    checks=[]
    for part in item['parts']:
        path = ROOT/part['path']
        assert path.stat().st_size==part['bytes']
        index = int(path.stem.split('part')[-1])
        with DCDReader(str(path)) as reader:
            assert reader.n_atoms==200*cl and len(reader)==10000
            assert np.isclose(reader.dt,10.,atol=1e-4)
            initial=float(reader[0].time/1000)
            for j in sorted(set(range(0,10000,20))|set(range(0,10000,50))):
                ts = reader[j]
                assert np.allclose(ts.dimensions,[200,200,2800,90,90,90])
                box=ts.dimensions[:3].astype(float)*.1
                xyz=ts.positions.astype(float)*.1
                t=float(ts.time/1000)
                if not np.isclose(t,initial+j*.01,atol=.001,rtol=0):
                    raise RuntimeError('DCD time mismatch')
                fr=(index-1)*10000+j
                if j%20==0:
                    original = selector.get_com_nm(ts.positions,200)
                    frames.append(fr); times.append(t)
                    coms.append(whole_com(xyz,box,masses[0]))
                    flags.append(bulk_flags(original,box,item['selection_params']))
                    uniform.append(bulk_flags(original,box,(4.,20,80)))
                if j%50==0:
                    counts=contacts(xyz,box,cl)
                    if j==0:
                        check_counts(xyz,box,cl,counts)
                        checks.append(dict(part=index,frame=fr,reference_chains=4,exact_match=True))
                    cf.append(fr); ct.append(t); pc.append(counts)
        print(f'{case} block{bid} part{index}: {len(pc)} contact frames, {time.monotonic()-started:.1f}s',flush=True)
    assert np.all(np.diff(frames)==20) and np.all(np.diff(cf)==50)
    assert np.allclose(np.diff(times),.2,atol=.001,rtol=0)
    assert np.allclose(np.diff(ct),.5,atol=.001,rtol=0)
    out.parent.mkdir(parents=True,exist_ok=True)
    tmp=out.with_suffix('.tmp.npz')
    np.savez_compressed(tmp,identity=np.asarray(ident),frames=np.asarray(frames),time_ns=times,
        com_nm=np.asarray(coms,np.float32),bulk=np.asarray(flags),uniform_bulk=np.asarray(uniform),
        contact_frames=np.asarray(cf),contact_time_ns=ct,partner_contacts=np.asarray(pc,np.int32),box_nm=box)
    tmp.replace(out)
    dump(out.with_suffix('.json'),dict(case=case,block=bid,status='complete',identity=ident,
        first_time_ns=times[0],last_time_ns=times[-1],contact_frames=len(pc),selection_frames=len(frames),
        elapsed_seconds=time.monotonic()-started,operator_checks=checks,topology_sha256=digest(pdb)))
    return str(out)


if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('--workers',type=int,default=4)
    ap.add_argument('--pilot',action='store_true')
    args=ap.parse_args()
    items=plan()
    if args.pilot:
        items=items[:1]
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for value in pool.map(process,items):
            print('COMPLETE',value,flush=True)
