#!/usr/bin/env python3
"""Cross-validate saved identities, counts and coordinates without editing inputs."""
import argparse
import json
from pathlib import Path
import numpy as np
from compute import HERE, ROOT, sha, dump


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cases',required=True)
    args=parser.parse_args()
    reports=[]
    for case in args.cases.split(','):
        out=HERE/'raw'/case
        manifest=json.loads((out/'inputs.json').read_text())
        done=json.loads((out/'completed.json').read_text())
        assert manifest['script_sha256']==sha(HERE/'compute.py')
        for p in manifest['parts']:
            stat=(ROOT/p['path']).stat()
            assert stat.st_size==p['bytes'] and stat.st_mtime_ns==p['mtime_ns']
        with np.load(out/'timeseries.npz') as z:
            frame,seg,chain,t,c=[z[k] for k in ['frame','segment','chain','time_ns','partner_contacts']]
            assert str(z['identity_hash'])==done['identity_hash']
        assert len(frame)==done['chain_frames']
        assert len(np.unique(frame))==done['source_frames']
        assert len(np.unique(seg))==done['sampled_segments']
        assert np.all(np.isfinite(t)) and np.all(c>=0)
        assert np.all(c[np.arange(len(chain)),chain]==0)
        assert np.all(frame%50==0)
        assert len(np.unique(np.column_stack([frame,seg]),axis=0))==len(frame)
        pair_checks=0
        changes=np.r_[0,np.flatnonzero(np.diff(frame)!=0)+1,len(frame)]
        for a,b in zip(changes[:-1],changes[1:]):
            sub=c[a:b][:,chain[a:b]]
            np.testing.assert_array_equal(sub,sub.T)
            pair_checks+=(b-a)*(b-a-1)//2
        conditions=dict(case=case,status='passed',source_inventory_unchanged=True,
                        chain_frames=len(frame),cross_chain_reciprocity_checks=int(pair_checks),
                        no_duplicate_samples=True,no_self_chain_contacts=True,
                        minimum_interchain_contacts=int(c.sum(axis=1).min()),
                        maximum_interchain_contacts=int(c.sum(axis=1).max()))
        reports.append(conditions)
        dump(out/'output_audit.json',conditions)
        print(conditions,flush=True)
    dump(HERE/'output_audits.json',reports)


if __name__=='__main__':
    main()
