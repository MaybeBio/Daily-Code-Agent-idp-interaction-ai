#!/usr/bin/env python3
"""Independent mdtraj coordinate/contact check of completed MDAnalysis extraction."""
import argparse
import json
import numpy as np
import mdtraj as md
from compute import HERE, ROOT, pilot, dump


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cases',required=True)
    args=parser.parse_args()
    legacy,_=pilot.load_legacy()
    results=[]
    for case in args.cases.split(','):
        c=next(c for c in legacy.CONDITIONS if f'{c.system}_{c.temp[:-1]}'==case)
        path=HERE/'raw'/case
        manifest=json.loads((path/'inputs.json').read_text())
        with np.load(path/'timeseries.npz') as z:
            f,chain,counts=z['frame'],z['chain'],z['partner_contacts']
            chain_len=int(z['chain_len'])
        rows=np.unique(np.linspace(0,len(f)-1,9,dtype=int))
        for row in rows:
            frame=int(f[row]);part=frame//10000+1;local=frame%10000
            source=ROOT/manifest['parts'][part-1]['path']
            tr=md.load_frame(str(source),local,top=str(c.topology))
            cid=int(chain[row])
            ref=legacy.count_contacts_for_chains(tr.xyz[0],tr.unitcell_lengths[0],chain_len,
                                                [cid],.66,np.repeat(np.arange(200),chain_len))[cid]
            got=int(counts[row].sum())
            results.append(dict(case=case,frame=frame,part=part,local=local,chain=cid,
                                mdtraj_reference=ref,MDAnalysis_result=got,match=ref==got))
        if not all(r['match'] for r in results if r['case']==case):
            dump(HERE/f'validation_{case}.json',results)
            raise AssertionError(f'{case}: independent reader contact mismatch; investigate precision/cutoff')
        dump(HERE/f'validation_{case}.json',[r for r in results if r['case']==case])
        print(f'{case}: {len(rows)} independent reader/contact checks passed',flush=True)


if __name__=='__main__':
    main()
