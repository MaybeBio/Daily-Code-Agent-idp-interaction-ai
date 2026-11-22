"""
rg_dynamics.py — single-chain dynamics from the stored 300 K frames (deck wetware row 'minute perturbations': single-
molecule switching is testable in our trajectories).

Per chain and model, from the per-frame Rg series of the FDT npz (1,000 evenly spaced production frames):
  tau_int    integrated autocorrelation time of Rg (frames; Sokal window, c = 5), and n_eff = frames / (2 tau_int)
  switching  two-state hysteresis on Rg: compact below the 25th percentile, expanded above the 75th; transitions
             compact <-> expanded per 1,000 frames and the mean dwell (frames) in each state
  bimodality Sarle's coefficient of the Rg distribution (> 0.555 suggests two populations)
Times in model frames; CALVADOS-2 also in ns of simulation time (run length 100 ns up to 150 residues, x (N/150)^1.5
above; the production half holds the 1,000 frames). Coarse-grained Langevin times are not physical times: compare chains
within a model, not with experiment.
Then: do faster-fluctuating or more frequently switching chains respond more strongly to heat? Spearman of tau_int and
switching with chi and the calibrated total response, raw and with nu and log N held fixed (rank partial).
    python rg_dynamics.py [--workers N] [--out DIR]
Output: Analysis_FDT/total_response/dynamics/ (per_chain.csv, summary.txt)
"""
import sys

sys.dont_write_bytecode = True

import argparse  # noqa: E402
import os  # noqa: E402
from concurrent.futures import ThreadPoolExecutor  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

HERE = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, HERE)
import total_response_evidence as E  # noqa: E402

TLA = E.TLA
OUT = f'{TLA}/Analysis_FDT/total_response/dynamics'
NPZ = {'C': f'{TLA}/Analysis_FDT/npz/300K', 'M': f'{TLA}/Analysis_FDT_MPIPI/npz_300K', 'R': f'{TLA}/Analysis_FDT_RECHARGED/npz_300K'}


def tau_int(x, c=5.0):
    x = np.asarray(x, float) - np.mean(x)
    n = len(x)
    f = np.fft.rfft(x, 2 * n)
    acf = np.fft.irfft(f * np.conj(f))[:n]
    acf /= acf[0]
    tau = 0.5
    for w in range(1, n):
        tau += acf[w]
        if w >= c * tau:
            break
    return tau


def switching(x):
    lo, hi = np.percentile(x, [25, 75])
    trans, dwell, cur, run = 0, {-1: [], 1: []}, 0, 0
    for v in x:
        s = -1 if v < lo else (1 if v > hi else 0)
        if s != 0 and s != cur:
            if cur != 0:
                trans += 1
                dwell[cur].append(run)
            cur, run = s, 0
        run += 1
    if cur != 0:
        dwell[cur].append(run)
    return trans, np.mean(dwell[-1]) if dwell[-1] else np.nan, np.mean(dwell[1]) if dwell[1] else np.nan


def one(args):
    k, uid = args
    try:
        p = os.path.join(NPZ[k], f'{uid}.npz')
        if k == 'M' and not os.path.exists(p):   # the 250 compact chains
            p = os.path.join(TLA, 'Analysis_FDT_MPIPI', 'compact', 'npz_300K', f'{uid}.npz')
        z = np.load(p)
        rg = z['rg'].astype(float)
        rg = rg[np.isfinite(rg)]
        t = tau_int(rg)
        tr, dc, de = switching(rg)
        g = stats.skew(rg) ** 2 + 1
        kur = stats.kurtosis(rg) + 3 * (len(rg) - 1) ** 2 / ((len(rg) - 2) * (len(rg) - 3))
        return {'model': E.MODELS[k], 'UniProt_ID': uid, 'frames': len(rg), 'tau_int_frames': t, 'n_eff': len(rg) / (2 * t),
                'switches_per_1000': 1000 * tr / len(rg), 'dwell_compact': dc, 'dwell_expanded': de,
                'bimodality': g / kur, 'rg_cv': rg.std() / rg.mean()}
    except Exception as ex:  # noqa: BLE001
        return {'model': E.MODELS[k], 'UniProt_ID': uid, 'error': repr(ex)}


def rn(v):
    return stats.norm.ppf((stats.rankdata(v) - 0.5) / len(v))


def partial(y, x, Z):
    Zr = np.column_stack([np.ones(len(y))] + [rn(z) for z in Z])
    ry = rn(y) - Zr @ np.linalg.lstsq(Zr, rn(y), rcond=None)[0]
    rx = rn(x) - Zr @ np.linalg.lstsq(Zr, rn(x), rcond=None)[0]
    return np.corrcoef(ry, rx)[0, 1]


def main(workers, out):
    os.makedirs(out, exist_ok=True)
    D, _ = E.load()
    tasks = [(k, u) for k in ('C', 'M', 'R') for u in D[k].index]
    with ThreadPoolExecutor(workers) as ex:
        P = pd.DataFrame(list(ex.map(one, tasks)))
    P = P.dropna(subset=['tau_int_frames'])
    P['N'] = [D[{v: k for k, v in E.MODELS.items()}[m]].N.get(u) for m, u in zip(P.model, P.UniProt_ID)]
    run_ns = 100.0 * np.maximum(1.0, (P.N / 150.0) ** 1.5)
    P['tau_int_ns_C'] = np.where(P.model == 'CALVADOS-2', P.tau_int_frames * (run_ns / 2) / 1000, np.nan)
    rows = []
    for k, m in E.MODELS.items():
        x = P[P.model == m].set_index('UniProt_ID').join(D[k][['chi', 'total', 'nu']], how='inner')
        x = x.dropna(subset=['chi', 'total', 'nu', 'N', 'tau_int_frames', 'switches_per_1000', 'bimodality', 'rg_cv'])
        logN = np.log(x.N)
        for dyn in ('tau_int_frames', 'switches_per_1000', 'bimodality', 'rg_cv'):
            for resp in ('chi', 'total'):
                rows.append({'model': m, 'dynamics': dyn, 'response': resp, 'n': len(x),
                             'spearman': stats.spearmanr(x[dyn], x[resp])[0],
                             'partial | nu, logN': partial(x[resp].to_numpy(), x[dyn].to_numpy(), [x.nu.to_numpy(), logN.to_numpy()])})
    C = pd.DataFrame(rows)
    P.to_csv(os.path.join(out, 'per_chain.csv'), index=False)
    C.to_csv(os.path.join(out, 'correlations.csv'), index=False)
    S = P.groupby('model').agg(chains=('UniProt_ID', 'size'), tau_frames_median=('tau_int_frames', 'median'),
                               tau_frames_p95=('tau_int_frames', lambda s: s.quantile(0.95)), n_eff_median=('n_eff', 'median'),
                               n_eff_p5=('n_eff', lambda s: s.quantile(0.05)), switches_median=('switches_per_1000', 'median'),
                               frac_bimodal=('bimodality', lambda s: (s > 0.555).mean()), tau_ns_median=('tau_int_ns_C', 'median'))
    pd.set_option('display.width', 200)
    lines = ['== per model', S.round(3).to_string(), '', '== dynamics vs response (Spearman; rank partial with nu and log N fixed)',
             C.round(3).to_string(index=False)]
    open(os.path.join(out, 'summary.txt'), 'w').write('\n'.join(lines) + '\n')
    print('\n'.join(lines))


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--workers', type=int, default=8)
    ap.add_argument('--out', default=OUT)
    a = ap.parse_args()
    main(a.workers, a.out)
