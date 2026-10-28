"""Extend partner-retention lags using the existing paired 20 ns trajectories."""
import ast
import csv
import hashlib
import json
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
REV = HERE.parent.parent
SOURCE = REV/'paired_slab_temperature_jump/analysis_full20ns_2026_10_09'
SCRIPT = REV/'FUS_slab_temperature_jump/analyze_pilot_acf.py'
ANALYSIS = REV/'paired_slab_temperature_jump/analyze_full20ns.py'
SYSTEMS = {'FUS': (280,325), 'A1': (230,260), 'TDP43': (260,300)}
HORIZONS = [9.95, 12.5, 15., 17.5, 19.]
OUT = HERE/'extended_partner'
OUT.mkdir(exist_ok=True)


def original_function(path, name):
    node = next(n for n in ast.parse(path.read_text()).body if isinstance(n, ast.FunctionDef) and n.name==name)
    scope = {'np': np}
    exec(compile(ast.Module(body=[node], type_ignores=[]), str(path), 'exec'), scope)
    return scope[name]


def branch(job):
    system, temperature, snapshot = job
    path = SOURCE/'observables'/f'{system}_snapshot{snapshot:02d}_{temperature}K.npz'
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    metadata = json.loads(path.with_suffix('.json').read_text())
    assert digest == metadata['output_sha256']
    with np.load(path) as source:
        time = source['time_ns']
        selected = (time > .1+1e-8) & (time <= 20.+1e-8)
        assert selected.sum()==1990
        contacts = source['partner_contacts'][selected]
        assert (contacts.sum(axis=2)>0).all()
        curve = original_function(SCRIPT, 'partner_curve')(contacts, 1900)
    with np.load(SOURCE/'per_snapshot_curves.npz') as source:
        baseline = source[f'{system}__{temperature}__full__partner__original'][snapshot-1]
    np.testing.assert_allclose(curve[:996], baseline, rtol=0, atol=2e-12)
    return job, curve, digest


def main():
    jobs = [(s,t,n) for s,temps in SYSTEMS.items() for t in temps for n in range(1,9)]
    curves = {}; hashes = {}
    with ProcessPoolExecutor(max_workers=4) as pool:
        for n,(job,curve,digest) in enumerate(pool.map(branch,jobs),1):
            curves[job] = curve
            hashes['%s_%s_%02d'%job] = digest
            if n%8==0: print(f'Checked and extended {n}/48 branches', flush=True)
    weights = np.load(SOURCE/'bootstrap_pair_weights.npz')['weights']
    assert weights.shape==(4000,8)
    cross = original_function(ANALYSIS,'crossing')
    saved = {}; records = []; condition_records = []
    for system,temps in SYSTEMS.items():
        values = [np.stack([curves[system,t,n] for n in range(1,9)]) for t in temps]
        for t, v in zip(temps,values):saved[f'{system}_{t}'] = v
        means = [v.mean(0) for v in values]
        boots = [weights@v for v in values]
        for horizon in HORIZONS:
            end = int(round(horizon*100))+1
            point = [cross(v[:end],.5)[0] for v in means]
            tau = [cross(v[:,:end],.5) for v in boots]
            unit = [cross(v[:,:end],.5) for v in values]
            finite = np.isfinite(tau[0]) & np.isfinite(tau[1])
            ratio = point[1]/point[0] if np.isfinite(point).all() else None
            ci = np.quantile(tau[1]/tau[0],[.025,.975],method='inverted_cdf') if finite.all() else [None,None]
            records.append(dict(system=system,horizon_ns=horizon,low_crossing_ns=float(point[0]) if np.isfinite(point[0]) else None,
                                high_crossing_ns=float(point[1]) if np.isfinite(point[1]) else None,ratio=ratio,
                                ci_low=ci[0],ci_high=ci[1],n_bootstrap=4000,both_uncensored=int(finite.sum()),
                                min_origins_per_chain=1990-(end-1),low_individual_censored=int((~np.isfinite(unit[0])).sum()),
                                high_individual_censored=int((~np.isfinite(unit[1])).sum())))
            for t,b,mean in zip(temps,tau,point):
                condition_records.append(dict(system=system,temperature=t,horizon_ns=horizon,
                                              mean_crossing_ns=float(mean) if np.isfinite(mean) else None,
                                              bootstrap_censored=int((~np.isfinite(b)).sum())))
    np.savez_compressed(OUT/'per_snapshot_partner_curves.npz',**saved)
    for name,data in [('ratios',records),('conditions',condition_records)]:
        (OUT/f'{name}.json').write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')
        with (OUT/f'{name}.csv').open('w',newline='') as stream:
            writer=csv.DictWriter(stream,fieldnames=list(data[0]));writer.writeheader();writer.writerows(data)
    (OUT/'provenance.json').write_text(json.dumps({'source':str(SOURCE),'max_lag_ns':19.,'analysis_window_ns':[.1,20.],
        'source_hashes':hashes,'original_lags_reproduced':True,'all_4000_original_paired_weights_retained':True,
        'partner_definition':'Endpoint overlap normalized by each origin partner count; returns allowed; chains then snapshots averaged.'},indent=2)+'\n')
    print(json.dumps(records,indent=2),flush=True)


if __name__=='__main__':main()
