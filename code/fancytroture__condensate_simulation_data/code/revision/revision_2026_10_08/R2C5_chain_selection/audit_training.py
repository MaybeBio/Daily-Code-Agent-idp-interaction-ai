#!/usr/bin/env python3
"""Verify current VAMP chain coordinates against their raw parent windows."""
import ast
import csv
import json
from pathlib import Path
import re
import warnings
import numpy as np
from MDAnalysis.coordinates.DCD import DCDReader
from extract import ROOT,HERE,AUDIT,digest,dump,CONFIG
from analyze import table

warnings.filterwarnings('ignore',category=DeprecationWarning)


def constants(path):
    wanted={'CUTOFF_NM','PERCENTILE_LOW','PERCENTILE_HIGH','STRIDE','SKIP_US','WINDOW_US','N_TARGET','TIME_BINS'}
    out={}
    for node in ast.parse(path.read_text()).body:
        if isinstance(node,ast.Assign):
            for name in node.targets:
                if isinstance(name,ast.Name) and name.id in wanted:
                    out[name.id]=ast.literal_eval(node.value)
    return out


def main():
    inv=json.loads((AUDIT/'inventory.json').read_text())
    stages={s['case']:s for s in json.loads((AUDIT/'stage_data.json').read_text())['stages']}
    sources=[];rows=[];summary=[]
    for entry in inv:
        case=entry['case'];stage=stages[case]
        metapath=Path(stage['source_metadata']);meta=json.loads(metapath.read_text())
        bychain={s['chain_id']:s for s in meta}
        parts={int(re.search(r'part(\d+)',p['path'])[1]):p for p in entry['parts']}
        sim=Path(entry['simulation']);train=Path(entry['training'])
        extract=sim/'extract_20chains_1us_dcd.py'
        if extract.exists():
            cfg=constants(extract)
            if cfg:
                assert (cfg['CUTOFF_NM'],cfg['PERCENTILE_LOW'],cfg['PERCENTILE_HIGH'])==CONFIG[case]
                sources.append(dict(case=case,path=str(extract),sha256=digest(extract),parameters=cfg))
        wrapper=train/'run_prepare_tica.sh'
        if case.startswith('TDP43'):
            txt=wrapper.read_text()
            assert 'EXTRACT_MIN_BULK_FRACTION:-0.8' in txt
            sources.append(dict(case=case,path=str(wrapper),sha256=digest(wrapper),
                wrapper_min_bulk_fraction_default=.8,python_default=1.,
                actual_historical_environment_override='not archived'))
        good=0;unknown=[]
        for file in entry['selections'][1]['dcd']:
            path=ROOT/file['path'];cid=int(re.search(r'chain(\d+)',path.name)[1])
            with DCDReader(str(path)) as selected:
                if cid not in bychain:
                    unknown.append(dict(chain=cid,path=str(path),frames=len(selected)))
                    rows.append(dict(case=case,chain_id=cid,trajectory=str(path),frames=len(selected),
                        metadata_available=False,checked_frames=0,max_pbc_coordinate_error_A='',source_start_frame=''))
                    continue
                segment=bychain[cid]
                assert len(selected)==segment['end']-segment['start']+1==50000
                errors=[]
                for j in [0,len(selected)//2,len(selected)-1]:
                    target=selected[j].positions.copy()
                    frame=segment['start']+j
                    part=parts[frame//10000+1]
                    with DCDReader(str(ROOT/part['path'])) as source:
                        ts=source[frame%10000]
                        xyz=ts.positions.reshape(200,-1,3)[cid]
                        delta=target.astype(float)-xyz
                        box=ts.dimensions[:3]
                        delta-=np.rint(delta/box)*box
                        errors.append(float(np.abs(delta).max()))
                assert max(errors)<1e-5,(case,cid,errors)
                good+=1
                rows.append(dict(case=case,chain_id=cid,trajectory=str(path),frames=len(selected),
                    metadata_available=True,checked_frames=3,max_pbc_coordinate_error_A=max(errors),
                    source_start_frame=segment['start']))
        summary.append(dict(case=case,metadata_trajectories=len(meta),coordinate_verified_trajectories=good,
            inter_used_chains=len(stage['inter_parent_chain_ids_reconstructed']),
            training_trajectories=good+len(unknown),untraced_training_trajectories=unknown,
            metadata_sha256=digest(metapath),
            statement='The 40 traced VAMP trajectories share parent selected windows with the interchain analysis; not all 200 chains.'))
        print('VERIFIED',case,good,'untraced',len(unknown),flush=True)
    table('training_coordinate_provenance.csv',rows)
    dump(HERE/'training_provenance.json',summary)
    dump(HERE/'selection_source_parameters.json',sources)


if __name__=='__main__':main()
