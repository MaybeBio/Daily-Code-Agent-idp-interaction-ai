#!/usr/bin/env python3
"""Gap-safe ACFs, partner overlaps and synchronized origin-block uncertainty."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
DT = .5
BLOCKS = (100., 250., 500.)
SYSTEMS = [('FUS', '280', '325'), ('A1', '230', '260'), ('TDP43', '260', '300')]


def write_csv(path, rows):
    if not rows:
        return
    fields = list(dict.fromkeys(k for row in rows for k in row))
    with path.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def runs(frame, time):
    """Return slices without joining missing raw frames or DCD time jumps."""
    if len(frame) == 0:
        return []
    boundaries = np.flatnonzero((np.diff(frame) != 50) |
                               (~np.isclose(np.diff(time), DT, atol=.001, rtol=0)))+1
    ends = np.r_[0, boundaries, len(frame)]
    return [slice(int(a), int(b)) for a,b in zip(ends[:-1], ends[1:])]


def correlation_values(values, lag):
    """Per-origin connected covariance, normalized by full-run population variance."""
    x = np.asarray(values, dtype=float)
    fluct = x-x.mean()
    var = np.mean(fluct**2)
    if len(x) <= lag or var <= 1e-12:
        return np.empty(0)
    return fluct[:len(x)-lag]*fluct[lag:] / var


def partner_values(contacts, lag):
    present = contacts > 0
    if len(present) <= lag:
        return np.empty(0)
    origin = present[:len(present)-lag]
    n = origin.sum(axis=1)
    overlap = np.sum(origin & present[lag:], axis=1)
    return np.divide(overlap, n, out=np.full(len(n), np.nan), where=n>0)


def crossing(t, y, threshold):
    hit = np.flatnonzero(np.asarray(y) <= threshold)
    if not len(hit):
        return np.inf
    i = int(hit[0])
    if i == 0:
        return float(t[0])
    return float(t[i-1] + (threshold-y[i-1])*(t[i]-t[i-1])/(y[i]-y[i-1]))


def interval(values):
    # Order-statistic quantiles retain right-censored (+inf) draws without
    # interpolation such as inf-inf, and never discard noncrossing draws.
    x = np.asarray(values)
    if np.isnan(x).any():
        raise ValueError('NaN bootstrap statistic')
    return tuple(float(v) for v in np.quantile(x, [.025, .975], method='inverted_cdf'))


def ratio_bounds(num, den, horizon):
    """Conservative ratios when a first crossing is only known to exceed horizon."""
    num, den = np.asarray(num), np.asarray(den)
    nlo = np.where(np.isinf(num), horizon, num)
    dlo = np.where(np.isinf(den), horizon, den)
    if np.any(dlo <= 0):
        raise ValueError('Nonpositive denominator time')
    lower = np.divide(nlo, den, out=np.zeros_like(nlo), where=np.isfinite(den))
    upper = num / dlo
    return lower, upper


def analyze_case(case, out, n_boot, seed):
    raw = HERE/'raw'/case
    complete = json.loads((raw/'completed.json').read_text())
    if complete['status'] != 'complete' or complete['failed_rows']:
        raise ValueError(f'Incomplete raw computation {case}')
    with np.load(raw/'timeseries.npz') as z:
        frame, seg, chains, times, contacts = [z[k] for k in
                            ['frame', 'segment', 'chain', 'time_ns', 'partner_contacts']]
        raw_identity = str(z['identity_hash'])
    units, unit_info = [], []
    excluded = []
    for s in np.unique(seg):
        ids = np.flatnonzero(seg == s)
        order = np.argsort(frame[ids])
        ids = ids[order]
        for run_id, sl in enumerate(runs(frame[ids], times[ids])):
            ix = ids[sl]
            info = dict(case=case, segment=int(s), chain=int(chains[ix[0]]), run=run_id,
                        n_samples=len(ix), span_ns=float(times[ix[-1]]-times[ix[0]]),
                        start_ns=float(times[ix[0]]), end_ns=float(times[ix[-1]]),
                        start_frame=int(frame[ix[0]]), end_frame=int(frame[ix[-1]]))
            if len(ix) < 8:
                excluded.append(dict(info, reason='fewer_than_8_contiguous_samples'))
                continue
            info.update(mean_N_inter=float(contacts[ix].sum(axis=1).mean()),
                        mean_N_partner=float((contacts[ix]>0).sum(axis=1).mean()))
            units.append(ix)
            unit_info.append(info)
    write_csv(out/f'{case}_units.csv', unit_info)
    write_csv(out/f'{case}_excluded.csv', excluded)
    result = {}
    for metric, maxlag, thresholds in [('inter', 60, (.1,.2,.3)), ('partner', 100, (.3,.5,.7))]:
        t = np.arange(maxlag+1)*DT
        unit_curves = np.full((len(units), maxlag+1), np.nan)
        block_stats = {}
        for block in BLOCKS:
            start = int(np.floor(times.min()/block))
            stop = int(np.floor(times.max()/block))
            block_stats[block] = (start, np.zeros((stop-start+1, maxlag+1)),
                                 np.zeros((stop-start+1, maxlag+1)))
        for u, ix in enumerate(units):
            values = contacts[ix].sum(axis=1).astype(float)
            for lag in range(maxlag+1):
                origin_values = correlation_values(values, lag) if metric=='inter' else partner_values(contacts[ix], lag)
                if not len(origin_values):
                    continue
                valid = np.isfinite(origin_values)
                if not valid.any():
                    continue
                y = origin_values[valid]
                unit_curves[u, lag] = y.mean()
                origin_t = times[ix[:len(origin_values)]][valid]
                for block, (start, numer, denom) in block_stats.items():
                    bi = np.floor(origin_t/block).astype(int)-start
                    # Each run contributes total weight one at a given lag.
                    numer[:, lag] += np.bincount(bi, weights=y/len(y), minlength=len(numer))
                    denom[:, lag] += np.bincount(bi, weights=np.full(len(y),1/len(y)), minlength=len(numer))
        mean = np.nanmean(unit_curves, axis=0)
        if not np.all(np.isfinite(mean)) or not np.isclose(mean[0],1):
            raise ValueError(f'Invalid {case} {metric} ACF')
        boots = {}
        for block, (start, num, den) in block_stats.items():
            np.testing.assert_allclose(num.sum(axis=0)/den.sum(axis=0), mean, atol=2e-13)
            rng = np.random.default_rng(seed + int(block))
            weights = rng.multinomial(len(num), np.full(len(num),1/len(num)), size=n_boot)
            boot_den = weights@den
            if np.any(boot_den <= 0):
                raise ValueError('Bootstrap replicate has no supported pairs; revise blocks explicitly')
            boots[block] = (weights@num)/boot_den
        rng = np.random.default_rng(seed+1)
        unit_weights = rng.multinomial(len(units), np.full(len(units), 1/len(units)), size=n_boot)
        finite = np.isfinite(unit_curves)
        segment_den = unit_weights@finite.astype(float)
        segment_boot = (unit_weights@np.nan_to_num(unit_curves))/segment_den
        payload = dict(time_ns=t, mean=mean, unit_curves=unit_curves, segment_boot=segment_boot,
                       raw_identity=np.asarray(raw_identity))
        for b in BLOCKS:
            payload[f'block_{int(b)}_boot'] = boots[b]
            payload[f'block_{int(b)}_numerator'] = block_stats[b][1]
            payload[f'block_{int(b)}_denominator'] = block_stats[b][2]
        np.savez_compressed(out/f'{case}_{metric}_curves.npz', **payload)
        rows = []
        for block, boot in boots.items():
            lo, hi = np.quantile(boot,[.025,.975],axis=0)
            for j in range(len(t)):
                rows.append(dict(case=case,metric=metric,block_ns=block,lag_ns=t[j],mean=mean[j],
                                 ci_low=lo[j],ci_high=hi[j],n_units=int(finite[:,j].sum())))
        write_csv(out/f'{case}_{metric}_curves.csv', rows)
        result[metric] = dict(t=t,mean=mean,units=unit_curves,boots=boots,
                              thresholds=thresholds,segment_boot=segment_boot)
    return result


def summarize(results, out):
    rows, ratios, censoring = [], [], []
    stats = {}
    for case, metrics in results.items():
        for metric, r in metrics.items():
            t, mean = r['t'], r['mean']
            horizon = float(t[-1])
            for threshold in r['thresholds']:
                ut = np.array([crossing(t,y,threshold) for y in r['units']])
                censoring.append(dict(case=case,metric=metric,threshold=threshold,horizon_ns=horizon,
                                 n_units=len(ut),right_censored_units=int(np.isinf(ut).sum())))
            for method, boot in [(f'origin_block_{int(b)}ns',r['boots'][b]) for b in BLOCKS]+[('segment_bootstrap',r['segment_boot'])]:
                quantities = []
                for end in [10.,20.,30.] + ([50.] if metric=='partner' else []):
                    keep = t <= end
                    quantities.append((f'I_{int(end)}ns',float(np.trapz(mean[keep],t[keep])),
                                       np.trapz(boot[:,keep],t[keep],axis=1),end,None))
                for th in r['thresholds']:
                    quantities.append((f'crossing_{th}',crossing(t,mean,th),
                                       np.array([crossing(t,y,th) for y in boot]),horizon,th))
                for name, val, samples, end, th in quantities:
                    lo, hi = interval(samples)
                    stats[(case,metric,method,name)] = (val,samples,end,th)
                    rows.append(dict(case=case,metric=metric,method=method,quantity=name,
                                 estimate_ns=val,ci_low_ns=lo,ci_high_ns=hi,window_ns=end,
                                 status='right_censored' if np.isinf(val) else 'observed',
                                 censoring_fraction=float(np.isinf(samples).mean()),
                                 first_interval_crossing=bool(th is not None and val<=DT)))
    for system, low, high in SYSTEMS:
        low_case, high_case = f'{system}_{low}', f'{system}_{high}'
        if low_case not in results or high_case not in results:
            continue
        for key, (l, ls, end, th) in stats.items():
            case, metric, method, name = key
            if case != low_case:
                continue
            h, hs, _, _ = stats[(high_case,metric,method,name)]
            if th is None:
                ratio = h/l if l>0 else np.nan
                if np.any(ls<=0):
                    lo, hi = np.nan,np.nan
                    bound_status='nonpositive_integral_denominator'
                else:
                    lo, hi = interval(hs/ls)
                    bound_status='finite'
                lower, upper = ratio,ratio
            else:
                lb, ub = ratio_bounds(hs,ls,end)
                lo = float(np.quantile(lb,.025,method='inverted_cdf'))
                hi = float(np.quantile(ub,.975,method='inverted_cdf'))
                lb0, ub0 = ratio_bounds(np.array([h]),np.array([l]),end)
                lower,upper=float(lb0[0]),float(ub0[0])
                ratio = h/l if np.isfinite(h) and np.isfinite(l) else np.nan
                bound_status='censoring_bounds' if np.isinf(hs).any() or np.isinf(ls).any() else 'finite'
            ratios.append(dict(system=system,metric=metric,method=method,quantity=name,
                               low_ns=l,high_ns=h,ratio=ratio,ratio_lower_bound=lower,
                               ratio_upper_bound=upper,ci_low=lo,ci_high=hi,
                               interval_status=bound_status,
                               low_boot_censored_fraction=float(np.isinf(ls).mean()),
                               high_boot_censored_fraction=float(np.isinf(hs).mean())))
    write_csv(out/'timescale_summary.csv',rows)
    write_csv(out/'temperature_ratios.csv',ratios)
    write_csv(out/'unit_censoring.csv',censoring)
    return rows,ratios


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cases', default=','.join(f'{s}_{v}' for s,l,h in SYSTEMS for v in (l,h)))
    parser.add_argument('--bootstrap', type=int, default=5000)
    args=parser.parse_args()
    out=HERE/'analysis'
    out.mkdir(exist_ok=True)
    results={}
    for i, case in enumerate(args.cases.split(',')):
        print(f'Analyzing {case}',flush=True)
        # Per-case seed must not change when cases are run separately.
        seed=20260930+int.from_bytes(hashlib.sha256(case.encode()).digest()[:3], 'little')
        results[case]=analyze_case(case,out,args.bootstrap,seed)
    summarize(results,out)
    (out/'methods.json').write_text(json.dumps(dict(bootstrap=args.bootstrap,dt_ns=DT,
        origin_blocks_ns=BLOCKS,script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        cases=list(results),primary_uncertainty='Synchronized absolute-time origin-block bootstrap; 250 ns primary, 100/500 ns sensitivity',
        centering='Per-contiguous-run sample mean and population variance; fixed during origin bootstrap',
        averaging='Equal supported-run weighting at each lag; bootstrap reweights original origins only, no splicing',
        partner='Two-time partner-set overlap divided by origin set size; returns allowed; not connected ACF',
        integral='Signed trapezoidal finite-window integral, no tail extrapolation or negative clipping',
        censoring='Noncrossing threshold draws retained as right-censored; interval bounds retain all draws',
        limitations='Historical selected chains; one slab per condition; fixed centering; no independent-slab or reselection uncertainty',
        segment_bootstrap='Secondary comparison only; chains from the same slab are not independent simulations'),indent=2)+'\n')


if __name__=='__main__':
    main()
