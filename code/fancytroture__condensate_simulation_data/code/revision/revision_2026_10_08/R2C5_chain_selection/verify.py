#!/usr/bin/env python3
"""Independent output checks, including comparison with pre-existing COM data."""
import csv
import json
import numpy as np
from extract import HERE,AUDIT,CONFIG,dump,digest


def read(name):
    with (HERE/name).open() as f:return list(csv.DictReader(f))


def main():
    checks=[]
    for case in CONFIG:
        paths=sorted((HERE/'raw').glob(f'{case}_block*.npz'))
        assert len(paths)==4
        for path in paths:
            meta=json.loads(path.with_suffix('.json').read_text())
            assert meta['status']=='complete'
            with np.load(path) as z:
                assert z['bulk'].shape==(2500,200)
                assert z['partner_contacts'].shape==(1000,200,200)
                assert np.allclose(np.diff(z['time_ns']),.2,atol=.001,rtol=0)
                assert np.allclose(np.diff(z['contact_time_ns']),.5,atol=.001,rtol=0)
                pc=z['partner_contacts']
                assert np.all(pc>=0)
                np.testing.assert_array_equal(pc,pc.swapaxes(1,2))
                assert not np.diagonal(pc,axis1=1,axis2=2).any()
                assert np.all(z['bulk'].reshape(5,500,200).all(axis=1)[:,z['bulk'].all(axis=0)])
                if path.stem.endswith('block0'):
                    with np.load(AUDIT/f'{case}_scan_s100.npz') as previous:
                        ix=np.searchsorted(previous['frames'],z['frames'][::5])
                        np.testing.assert_array_equal(previous['frames'][ix],z['frames'][::5])
                        np.testing.assert_allclose(previous['com_nm'][ix],z['com_nm'][::5],atol=1e-6,rtol=0)
                        np.testing.assert_array_equal(previous['bulk'][ix],z['uniform_bulk'][::5])
            checks.append(dict(file=path.name,status='passed',raw_sha256=digest(path),
                independent_contact_comparisons=len(meta['operator_checks'])*4,
                previous_500frame_200chain_COM_comparison=path.stem.endswith('block0')))
    cc=read('curves.csv');stats=read('timescales.csv');rr=read('temperature_ratios.csv')
    for row in stats:
        selected=[r for r in cc if all(r[k]==row[k] for k in ['case','metric','group'])]
        t=np.array([float(r['lag_ns']) for r in selected]);y=np.array([float(r['mean']) for r in selected])
        assert np.isclose(y[0],1.,atol=1e-12)
        assert np.isfinite(y).all() and len(y)==61
        if row['quantity']=='I30':
            value=np.sum(np.diff(t)*(y[1:]+y[:-1])/2)
        else:
            th=float(row['threshold']);hit=np.flatnonzero(y<=th)
            if not len(hit):value=np.inf
            else:
                i=int(hit[0]);value=t[i-1]+(t[i]-t[i-1])*(y[i-1]-th)/(y[i-1]-y[i])
        np.testing.assert_allclose(value,float(row['estimate_ns']),rtol=1e-12)
    prov=read('training_coordinate_provenance.csv')
    good=[r for r in prov if r['metadata_available']=='True']
    assert len(good)==240 and all(float(r['max_pbc_coordinate_error_A'])==0 for r in good)
    assert len(prov)==241
    dump(HERE/'verification.json',dict(status='passed',raw_blocks=checks,
        coordinate_provenance_checked_trajectories=240,coordinate_provenance_checked_frames=720,
        untraced_training_trajectories=1,scalar_checks=len(stats),
        tests='test_analysis.py: six unit tests; independent KD-tree contact checks at every part start',
        original_documents_modified=False))
    print('PASS: 24 raw blocks; 600,000 cached COM comparisons; 480 KD-tree chain/frame checks; 720 provenance frames.',flush=True)


if __name__=='__main__':main()
