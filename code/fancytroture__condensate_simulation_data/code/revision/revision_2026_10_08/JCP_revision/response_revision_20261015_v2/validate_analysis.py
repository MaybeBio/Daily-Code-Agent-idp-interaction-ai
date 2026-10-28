"""Independent checks of sampling, selection support, and plotted statistics."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
from pathlib import Path
import csv,json
import numpy as np
from MDAnalysis.coordinates.DCD import DCDReader
from complete_residence import ROOT,BASE,OUT,PLAN,STAGES
HERE=Path(__file__).resolve().parent

def read(path):
    with path.open() as f:return list(csv.DictReader(f))

inventory={(r['case'],int(r['block']),int(r['subblock'])):r for r in read(BASE/'residence_counts.csv')}
checked_counts=0;raw_checks=0
for item in PLAN:
    case,bid=item['case'],item['block']
    stage=next(s for s in STAGES if s['case']==case)
    meta=json.loads(Path(stage['source_metadata']).read_text());original_ids=stage['inter_parent_chain_ids_reconstructed']
    with np.load(BASE/'raw'/f'{case}_block{bid}.npz') as z:
        short=z['bulk'].reshape(5,500,200).all(axis=1);frames=z['frames']
    with np.load(OUT/f'{case}_block{bid}_intra_0p1ns.npz') as z:fine=z['intra']
    with np.load(OUT/f'{case}_block{bid}_intramolecular.npz') as z:coarse=z['intra']
    valid=np.isfinite(fine[::5]);np.testing.assert_array_equal(fine[::5][valid],coarse[valid])
    first_supported=None
    for sub in range(5):
        first=int(frames[sub*500]);last=first+9999
        ids={s['chain_id'] for s in meta if s['chain_id'] in original_ids and first>=s['start'] and last<=s['end']}
        r=inventory[case,bid,sub]
        assert len(ids)==int(r['n_original'])
        assert short[sub].sum()==int(r['n_100'])
        checked_counts+=1
        expected=set(np.flatnonzero(short[sub]))|ids if ids else set()
        present=set(np.flatnonzero(np.isfinite(fine[sub*1000:(sub+1)*1000]).all(axis=0)))
        assert expected==present
        if ids and first_supported is None:first_supported=sub
    if first_supported is not None:
        sub=first_supported
        with DCDReader(str(ROOT/item['parts'][sub]['path'])) as reader:
            ts=reader[10];cl=reader.n_atoms//200;box=ts.dimensions[:3].astype(float)*.1
            xyz=ts.positions.astype(float).reshape(200,cl,3)*.1
        ids=np.flatnonzero(np.isfinite(fine[sub*1000+1]));ii,jj=np.triu_indices(cl,k=3)
        for chain in [ids[0],ids[-1]]:
            delta=xyz[chain,ii]-xyz[chain,jj];delta-=np.rint(delta/box)*box
            expected=np.count_nonzero(np.sum(delta*delta,axis=1)<.66**2)
            assert expected==fine[sub*1000+1,chain]
            raw_checks+=1

stats=read(OUT/'new_ratios.csv')
previous=(HERE/'residence_run.log').read_text();previous=json.loads(previous[previous.index('[\n  {'):])
original_rg={(r['system'],r['group']):r for r in previous if r['metric']=='rg'}
old_times=read(BASE/'timescales.csv')
temps={'FUS':(280,325),'A1':(230,260),'TDP43':(260,300)}
for r in stats:
    assert int(r['low_censored'])==int(r['high_censored'])==0
    assert float(r['low_ns'])>float(r['dt_ns']) and float(r['high_ns'])>float(r['dt_ns'])
    for role,temp in zip(['low','high'],temps[r['system']]):
        old=next(t for t in old_times if t['case']==f"{r['system']}_{temp}" and t['metric']=='inter' and t['group']==r['group'] and t['quantity']=='tau')
        assert int(r[role+'_blocks'])==int(old['n_blocks'])
    if r['metric']=='rg':
        old=original_rg[r['system'],r['group']]
        for field in ['low_ns','high_ns','ratio','ci_low','ci_high']:
            np.testing.assert_allclose(float(r[field]),float(old[field]),atol=1e-12,rtol=1e-12)
    else:assert float(r['dt_ns'])==.1

sen=json.loads((HERE/'sensitivity_ratios.json').read_text())
assert len(sen)==96
assert all(r['bootstrap']==5000 and r['low_censored']==r['high_censored']==0 for r in sen)
result=dict(selection_subwindows_verified=checked_counts,independent_raw_contact_checks=raw_checks,
            original_Rg_residence_results_preserved=True,intra_sampling_ns=.1,
            no_residence_crossing_at_first_sample=True,all_residence_bootstrap_draws_uncensored=True,
            sensitivity_ratios_with_intervals=len(sen),all_sensitivity_bootstrap_draws_uncensored=True)
(HERE/'analysis_validation.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
