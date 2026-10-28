#!/usr/bin/env python3
"""Apply existing trajectory bootstrap draws to the unchanged C3 display support."""
import hashlib
import json
import numpy as np
from recalculate import OUT, ROOT, old, read, write, distance

cache = ROOT/'revision_2026_10_08/R2C2_state_remodeling/followup/bootstrap_draws.npz'
references = read(OUT/'distance_comparison.csv')
results, contrasts, checks = [], [], []
values = {}
with np.load(cache) as z:
    for system in ('FUS','A1','TDP43'):
        states = old.original.transition_states(system)
        index = {s:i for i,s in enumerate(states)}
        for family in ('free','condensate'):
            original = [old.load_condition(system,family,r,states) for r in ('Tl','Th')]
            for source in ('model_seeds','sampling','sampling_10ns','sampling_100ns'):
                mapped = []
                for role,d in zip(('Tl','Th'),original):
                    source_key = 'sampling' if family=='free' and source.startswith('sampling') else source
                    prefix = f'{system}_{family}_{role}_{source_key}'
                    ts, ps = z[prefix+'_t'],z[prefix+'_mu']
                    assert ts.shape==(4001,5,5) and ps.shape==(4001,5)
                    np.testing.assert_allclose(ts.sum(-1),1,atol=1e-12)
                    t = np.zeros((4001,len(states),len(states)))
                    p = np.zeros((4001,len(states)))
                    for local_i,display_i in d['mapping'].items():
                        i=index[display_i]
                        p[:,i]=ps[:,local_i]
                        for local_j,display_j in d['mapping'].items():
                            t[:,i,index[display_j]]=ts[:,local_i,local_j]
                    mapped.append((t,p))
                for masked,definition in ((True,'submitted_edge_mask'),(False,'all_entries_original_display_union')):
                    ts = [item[0]*d['mask'] if masked else item[0] for item,d in zip(mapped,original)]
                    ds = distance(ts[0],ts[1],mapped[0][1],mapped[1][1])
                    assert np.isfinite(ds).all()
                    lo,hi = np.quantile(ds[1:],[.025,.975])
                    ref=next(r for r in references if (r['system'],r['family'],r['definition'])==(system,family,definition))
                    if source=='model_seeds':
                        np.testing.assert_allclose(ds[0],float(ref['distance']),atol=2e-8,rtol=0)
                        checks.append(dict(system=system,family=family,definition=definition,full_point_error=float(abs(ds[0]-float(ref['distance'])))))
                    results.append(dict(system=system,family=family,definition=definition,uncertainty=source,
                                        estimate=float(ds[0]),ci_low=float(lo),ci_high=float(hi),n_bootstrap=4000))
                    values[system,family,definition,source]=ds
        for definition in ('submitted_edge_mask','all_entries_original_display_union'):
            for source in ('model_seeds','sampling','sampling_10ns','sampling_100ns'):
                ds=values[system,'condensate',definition,source]-values[system,'free',definition,source]
                lo,hi=np.quantile(ds[1:],[.025,.975])
                contrasts.append(dict(system=system,definition=definition,uncertainty=source,
                                      estimate=float(ds[0]),ci_low=float(lo),ci_high=float(hi),n_bootstrap=4000))
write('trajectory_and_model_intervals.csv',results)
write('trajectory_and_model_contrasts.csv',contrasts)
(OUT/'sampling_validation.json').write_text(json.dumps(dict(checks=checks,cache=str(cache),
    cache_sha256=hashlib.sha256(cache.read_bytes()).hexdigest(),n_bootstrap=4000,
    no_new_resampling_or_simulation=True,original_display_mapping_and_zero_padding_retained=True,
    population_normalization='Reuse normalized full-five-state masses from existing bootstrap cache; no retained-support renormalization; full-point rounding agreement checked'),indent=2)+'\n')
for row in results:
    if row['uncertainty']=='sampling' and row['definition']=='all_entries_original_display_union':
        print(row)
