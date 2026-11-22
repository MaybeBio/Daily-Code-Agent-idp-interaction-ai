"""
charge_dependence_check.py — how much of the charge dependence of the total response under the calibrated model follows
from the hydrophilic forcing acting on charged residues, and how much from the ensembles (2026-10-25).

Exact split of the hydrophilic term of each chain (as total_response_evidence.factors):
    H / <Rg> = -r(Rg, D) x CV(Rg) x gamma sd(D) / RT       (coupling r x breathing x size of the forcing fluctuations)
Check: the coupling r is a correlation, so |r| <= 1 for every chain (this tests the stored SDs against the stored
covariance; the product itself reproduces H by construction).
Only the last factor scales by definition with the number of charged residues. Reported, for CALVADOS-2 and
Mpipi-Recharged, as rank partial correlations at fixed nu and log N (total_response_evidence.partial):
  each factor, H and the total against the fraction of D, E, K and R, and against |NCPR| at fixed D+E+K+R as well
  the share of the variance of log|H| beyond nu and log N carried by each factor (covariance with log|H| / its variance;
  the shares add to 1 because log|H| is the sum of the factors' log magnitudes)
  controls: stronger contraction (-chi, -H, -total) against the charged fraction and against the polar + Gly fraction
  (N, Q, S, T, H, G) for the plain force fields (chi) and for the linear form (A0), whose forcing also acts on polar
  residues and Gly
    BENDER_REANALYSIS_DIR=... python charge_dependence_check.py      (openmpipi env)
Output: Traj Level Analysis/Analysis_FDT/total_response/charge_dependence.txt
"""
import sys

sys.dont_write_bytecode = True

import os  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

HERE = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, os.environ.get('BENDER_REANALYSIS_DIR', HERE))
sys.path.insert(0, HERE)
import total_response_evidence as E  # noqa: E402

OUTF = os.path.join(E.OUT, 'charge_dependence.txt')
FDT = {'C': os.path.join(E.TLA, 'Analysis_FDT', 'fdt_results_300K.csv'),
       'R': os.path.join(E.TLA, 'Analysis_FDT_RECHARGED', 'recharged_fdt_300K.csv')}
RT = 8.314462618e-3 * 300.0


def composition(bio):
    b = bio[['f_DEK', 'f_R', 'NCPR', 'f_N', 'f_Q', 'f_S', 'f_T', 'f_G', 'f_H']].copy()
    b['fDEKR'] = b.f_DEK + b.f_R
    b['fpolG'] = b[['f_N', 'f_Q', 'f_S', 'f_T', 'f_G', 'f_H']].sum(axis=1)
    b['absNCPR'] = b.NCPR.abs()
    return b[['fDEKR', 'fpolG', 'absNCPR']]


def resid_on(x, v):
    A = np.column_stack([np.ones(len(x)), E.rank_normal(x.nu.values), np.log(x.N.values)])
    return v - A @ np.linalg.lstsq(A, v, rcond=None)[0]


def main():
    bio = composition(E.bio_table())
    L = ['Charge dependence of the total response, calibrated model (rank partial correlations at fixed nu and log N)', '']
    E.FORCING = 'capped'
    D, g = E.load()
    for k, lab in (('C', 'CALVADOS-2'), ('R', 'Mpipi-Recharged')):
        d = D[k]
        f = pd.read_csv(FDT[k])
        f = f[f.metric == 'rg'].set_index('UniProt_ID').reindex(d.index)
        x = pd.DataFrame({'coupling': -d.dchi_rg_hyd * RT / (f.A_sd * d.D_sd_hyd), 'CV': f.A_sd / f.A_mean,
                          'fluct': g[k] * d.D_sd_hyd / RT, 'H': d.H, 'total': d.total, 'nu': d.nu, 'N': d.N}).join(bio).dropna()
        dev = float((-x.coupling * x.CV * x.fluct * 100 - x.H).abs().max())
        rmax = float(x.coupling.abs().max())
        if dev > 1e-6 or not rmax <= 1 + 1e-9:
            raise SystemExit(f'{lab}: factor split inconsistent (H to {dev:.1e}, max |r| {rmax:.3f})')
        Z = np.column_stack([x.nu, np.log(x.N)])
        Zq = np.column_stack([x.nu, np.log(x.N), x.fDEKR])
        L.append(f'== {lab}: {len(x)} chains; split reproduces H to {dev:.1e}; max |r| {rmax:.3f}; median coupling r(Rg, D) '
                 f'{x.coupling.median():.2f}')
        for nm in ('coupling', 'CV', 'fluct', 'H', 'total'):
            L.append(f'  {nm:8s} vs D+E+K+R {E.partial(x[nm].values, x.fDEKR.values, Z):+.2f}; '
                     f'vs |NCPR| at fixed D+E+K+R {E.partial(x[nm].values, x.absNCPR.values, Zq):+.2f}')
        lh = resid_on(x, np.log(np.abs(x.H)).values)            # log|H| = sum of the factors' log magnitudes (exact)
        L.append('  share of the variance of log|H| beyond nu and log N: ' + ', '.join(
            f'{nm} {np.cov(resid_on(x, np.log(np.abs(x[nm])).values), lh)[0, 1] / lh.var():.2f}' for nm in ('coupling', 'CV', 'fluct')))
        L.append('')
    L.append('Controls: partial rank correlation with stronger contraction (-chi, -H, -total) | nu, log N')
    for forcing, name in (('A0', 'linear form'), ('capped', 'calibrated model')):
        E.FORCING = forcing
        D, _ = E.load()
        for k, lab in (('C', 'CALVADOS-2'), ('R', 'Mpipi-Recharged')):
            x = D[k][['chi', 'H', 'total', 'nu', 'N']].join(bio).dropna()
            Z = np.column_stack([x.nu, np.log(x.N)])
            r = lambda col, v: E.partial(-x[col].values, x[v].values, Z)  # noqa: E731
            L.append(f'  {name:16s} {lab:16s} charged fraction: chi {r("chi", "fDEKR"):+.2f}, H {r("H", "fDEKR"):+.2f}, '
                     f'total {r("total", "fDEKR"):+.2f} | polar + Gly fraction: H {r("H", "fpolG"):+.2f}, total {r("total", "fpolG"):+.2f}')
    E.FORCING = 'capped'
    open(OUTF, 'w').write('\n'.join(L) + '\n')
    print('\n'.join(L))


if __name__ == '__main__':
    main()
