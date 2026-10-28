#!/usr/bin/env python3
"""Sampling-grid and unit-level cross-checks without rereading raw DCDs."""
import argparse
import json
import numpy as np
from analyze import HERE, SYSTEMS, runs, crossing, write_csv


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cases',default=','.join(f'{s}_{v}' for s,l,h in SYSTEMS for v in (l,h)))
    args=parser.parse_args()
    rows, checks=[],[]
    for case in args.cases.split(','):
        raw=HERE/'raw'/case
        with np.load(raw/'timeseries.npz') as z:
            frame,seg,time,contacts=[z[k] for k in ['frame','segment','time_ns','partner_contacts']]
        metadata=json.loads((raw/'segments.json').read_text())
        manifest=json.loads((raw/'inputs.json').read_text())
        parts=manifest['parts']
        lo=(parts[-1]['part']-70)*10000
        hi=(parts[-1]['part']-1)*10000+parts[-1]['frames']-1
        # Independent interval arithmetic, without importing the schedule builder.
        for s,entry in enumerate(metadata):
            a=max(lo,entry['start']); b=min(hi,entry['end'])
            expected=np.arange(((a+49)//50)*50,b+1,50,dtype=np.int64)
            np.testing.assert_array_equal(frame[seg==s],expected)
        assert np.all(contacts>=0)
        all_units=[]
        for s in np.unique(seg):
            ix=np.flatnonzero(seg==s)
            for sl in runs(frame[ix],time[ix]):
                if sl.stop-sl.start>=8:
                    all_units.append(ix[sl])
        for stride in [1,2,4]:
            for metric,horizon,threshold in [('inter',30.,.2),('partner',50.,.5)]:
                dt=.5*stride
                t=np.arange(int(round(horizon/dt))+1)*dt
                curves=[]
                for ix in all_units:
                    # Absolute grid, not segment-relative resampling.
                    ix=ix[frame[ix]%(50*stride)==0]
                    if len(ix)<8:
                        continue
                    c=contacts[ix]
                    curve=np.full(len(t),np.nan)
                    if metric=='inter':
                        x=c.sum(1).astype(float)
                        x-=x.mean()
                        var=np.mean(x*x)
                        if var<=1e-12:
                            continue
                        ac=np.correlate(x,x,mode='full')[len(x)-1:]
                        n=min(len(curve),len(ac))
                        curve[:n]=ac[:n]/(var*np.arange(len(x),len(x)-n,-1))
                    else:
                        p=c>0
                        for lag in range(min(len(p),len(t))):
                            x=p[:len(p)-lag]; y=p[lag:]
                            den=x.sum(1); ok=den>0
                            if ok.any():
                                curve[lag]=np.mean((x[ok]&y[ok]).sum(1)/den[ok])
                    curves.append(curve)
                mean=np.nanmean(curves,axis=0)
                if stride==1:
                    with np.load(HERE/'analysis'/f'{case}_{metric}_curves.npz') as z:
                        err=float(np.max(np.abs(mean-z['mean'])))
                        np.testing.assert_allclose(mean,z['mean'],atol=1e-12)
                    checks.append(dict(case=case,metric=metric,independent_mean_max_error=err,
                                       independent_window_schedule_match=True,n_units=len(curves)))
                keep=t<=30
                rows.append(dict(case=case,metric=metric,dt_ns=dt,n_units=len(curves),
                                 first_crossing_ns=crossing(t,mean,threshold),threshold=threshold,
                                 I_30ns=float(np.trapz(mean[keep],t[keep])),
                                 first_interval_crossing=bool(crossing(t,mean,threshold)<=dt)))
    write_csv(HERE/'analysis'/'sampling_grid_sensitivity.csv',rows)
    write_csv(HERE/'analysis'/'independent_curve_checks.csv',checks)
    ratio=[]
    for system,low,high in SYSTEMS:
        for metric in ['inter','partner']:
            for dt in [.5,1.,2.]:
                l=next((r for r in rows if r['case']==f'{system}_{low}' and r['metric']==metric and r['dt_ns']==dt),None)
                h=next((r for r in rows if r['case']==f'{system}_{high}' and r['metric']==metric and r['dt_ns']==dt),None)
                if l is None or h is None:
                    continue
                ratio.append(dict(system=system,metric=metric,dt_ns=dt,
                                  crossing_ratio=h['first_crossing_ns']/l['first_crossing_ns'],
                                  I_30ns_ratio=h['I_30ns']/l['I_30ns']))
    write_csv(HERE/'analysis'/'sampling_grid_ratios.csv',ratio)
    print(f'{len(checks)} independent curve/window checks passed',flush=True)


if __name__=='__main__':
    main()
