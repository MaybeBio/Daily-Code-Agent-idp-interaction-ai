#!/usr/bin/env python3
"""Independent raw-Angstrom DCD reader check, with unit-conversion precision recorded."""
import argparse
import json
import mdtraj as md
import numpy as np
from compute import HERE, ROOT, pilot, dump
from investigate_precision import pairs


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cases',required=True)
    args=parser.parse_args()
    legacy,_=pilot.load_legacy()
    for case in args.cases.split(','):
        c=next(c for c in legacy.CONDITIONS if f'{c.system}_{c.temp[:-1]}'==case)
        out=HERE/'raw'/case
        source=json.loads((out/'inputs.json').read_text())
        with np.load(out/'timeseries.npz') as z:
            frame,chain,counts=z['frame'],z['chain'],z['partner_contacts']
            length=int(z['chain_len'])
        checks=[]
        for row in np.unique(np.linspace(0,len(frame)-1,9,dtype=int)):
            f,cid=int(frame[row]),int(chain[row])
            part,local=f//10000+1,f%10000
            path=ROOT/source['parts'][part-1]['path']
            with md.formats.DCDTrajectoryFile(str(path)) as reader:
                reader.seek(local)
                pos,cell,angle=reader.read(n_frames=1)
            assert np.allclose(angle,90)
            x,box=pos[0].astype(float)*.1,cell[0].astype(float)*.1
            ref=legacy.count_contacts_for_chains(x,box,length,[cid],.66,np.repeat(np.arange(200),length))[cid]
            high=md.load_frame(str(path),local,top=str(c.topology))
            highref=legacy.count_contacts_for_chains(high.xyz[0],high.unitcell_lengths[0],length,
                     [cid],.66,np.repeat(np.arange(200),length))[cid]
            differences=[]
            if highref!=ref:
                rawpairs=pairs(x,box,cid,length)
                highpairs=pairs(high.xyz[0],high.unitcell_lengths[0].astype(float),cid,length)
                for pair in sorted(rawpairs.symmetric_difference(highpairs)):
                    delta=x[pair[0]]-x[pair[1]]
                    delta-=box*np.round(delta/box)
                    distance=float(np.linalg.norm(delta))
                    differences.append(dict(atoms=list(map(int,pair)),distance_raw_nm=distance))
                    if abs(distance-.66)>1e-4:
                        raise AssertionError('Reader discrepancy is not confined to cutoff rounding')
            got=int(counts[row].sum())
            checks.append(dict(case=case,frame=f,part=part,local=local,chain=cid,
                independent_raw_float64_count=ref,production_count=got,match=ref==got,
                mdtraj_highlevel_float32_count=highref,highlevel_match=highref==got,
                highlevel_precision_difference_pairs=differences))
        dump(HERE/f'validation_rawunits_{case}.json',checks)
        if not all(r['match'] for r in checks):
            raise AssertionError(f'{case}: raw-units independent contact mismatch')
        print(f'{case}: {len(checks)} exact raw-unit checks passed; '
              f'{sum(not r["highlevel_match"] for r in checks)} high-level rounding differences recorded',flush=True)


if __name__=='__main__':
    main()
