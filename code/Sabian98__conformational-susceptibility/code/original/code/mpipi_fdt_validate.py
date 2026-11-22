"""
mpipi_fdt_validate.py — integration check of the Mpipi-GG FDT across the five temperatures, the counterpart of
`TODO_T2_trajectory_loader.py validate` for CALVADOS-2.

Per protein and observable (runs with status ok at all five temperatures; a run whose run.in temperature does not
match its folder is not a run at that temperature): chi_T at each T, the local slope d<A>/dT (np.gradient), and
integral of chi_T dT (trapezoid, 278 -> 320 K) vs <A>_320 - <A>_278, with z from the SEs in quadrature
(TODO_T2_fdt_computation.global_integration_check). Per adjacent temperature pair, mean trapezoid vs mean
observed difference.

Inputs: Traj Level Analysis/Analysis_FDT_MPIPI/mpipi_<T>K.csv and mpipi_<T>K_fdt.csv (mpipi_temperature_hpc.py).
Outputs (same folder): mpipi_fdt_validation.csv, mpipi_fdt_validation_summary.csv, mpipi_fdt_validation_segments.csv.

Usage: python3 mpipi_fdt_validate.py
"""
import sys

sys.dont_write_bytecode = True  # keep the NAS code folder free of __pycache__

import os  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

HERE = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.environ.get('BENDER_REANALYSIS_DIR', '/Users/taseefrahman/Downloads/Scripts/bender_reanalysis'))
import TODO_T2_fdt_computation as fdt  # noqa: E402
import bender_common as bc  # noqa: E402

OUT = '/Volumes/User Homes/Nu Project/BENDER-MT/Traj Level Analysis/Analysis_FDT_MPIPI'
TEMPS = [278, 288, 300, 310, 320]


def main(out=OUT):
    temps = np.array(TEMPS, dtype=float)
    ok = None
    for T in TEMPS:
        s = pd.read_csv(os.path.join(out, f'mpipi_{T}K.csv'), usecols=['UniProt_ID', 'status'])
        ids = set(s.loc[s['status'] == 'ok', 'UniProt_ID'])
        ok = ids if ok is None else ok & ids
    df = pd.concat([pd.read_csv(os.path.join(out, f'mpipi_{T}K_fdt.csv'),
                                usecols=['UniProt_ID', 'T', 'metric', 'chi_total', 'chi_total_se', 'A_mean', 'A_mean_se'])
                    for T in TEMPS])
    df['T'] = df['T'].round().astype(int)
    df = df[df['UniProt_ID'].isin(ok)]
    vals = ['chi_total', 'chi_total_se', 'A_mean', 'A_mean_se']
    wide = df.pivot_table(index=['UniProt_ID', 'metric'], columns='T', values=vals)
    wide = wide.dropna(subset=[(v, T) for v in vals for T in TEMPS])

    rows = []
    for (uid, metric), r in wide.iterrows():
        chi, chi_se, A, A_se = (np.array([r[(v, T)] for T in TEMPS]) for v in vals)
        row = {'UniProt_ID': uid, 'metric': metric}
        for T, c, s in zip(TEMPS, chi, fdt.local_slopes(A, temps)):
            row[f'chi_{T}K'] = c
            row[f'slope_{T}K'] = s
        for T in TEMPS:
            row[f'A_{T}K'] = r[('A_mean', T)]
        row.update(fdt.global_integration_check(chi, chi_se, A, A_se, temps))
        row['linear_slope'] = fdt.linear_slope(A, temps)
        rows.append(row)
    val = pd.DataFrame(rows)
    path = os.path.join(out, 'mpipi_fdt_validation.csv')
    val.to_csv(path, index=False)

    summ, seg = [], []
    for metric, g in val.groupby('metric'):
        for T in TEMPS:
            x, y = g[f'chi_{T}K'], g[f'slope_{T}K']
            summ.append({'metric': metric, 'T': T, 'n': len(g), 'mean_chi': x.mean(), 'mean_slope': y.mean(),
                         'pearson_r': x.corr(y), 'spearman_r': x.corr(y, method='spearman')})
        summ.append({'metric': metric, 'T': 'integral', 'n': len(g), 'mean_chi': g['integral'].mean(),
                     'mean_slope': g['observed_diff'].mean(),
                     'ratio_observed_to_integral': g['observed_diff'].mean() / g['integral'].mean(),
                     'pearson_r': g['integral'].corr(g['observed_diff']),
                     'spearman_r': g['integral'].corr(g['observed_diff'], method='spearman'),
                     'frac_abs_z_lt_2': (g['z'].abs() < 2).mean()})
        for T0, T1 in zip(TEMPS[:-1], TEMPS[1:]):
            trap = (T1 - T0) * (g[f'chi_{T0}K'] + g[f'chi_{T1}K']) / 2
            obs = g[f'A_{T1}K'] - g[f'A_{T0}K']
            seg.append({'metric': metric, 'segment': f'{T0}-{T1}', 'n': len(g), 'mean_trapezoid': trap.mean(),
                        'mean_observed': obs.mean(), 'ratio_observed_to_trapezoid': obs.mean() / trap.mean(),
                        'pearson_r': trap.corr(obs)})
    summ, seg = pd.DataFrame(summ), pd.DataFrame(seg)
    summ.to_csv(os.path.join(out, 'mpipi_fdt_validation_summary.csv'), index=False)
    seg.to_csv(os.path.join(out, 'mpipi_fdt_validation_segments.csv'), index=False)
    bc.write_provenance(path, __file__, {'temps': TEMPS, 'runs': 'status ok at all five temperatures',
                                         'slopes': 'np.gradient on the T grid',
                                         'integral': 'trapezoid, SEs in quadrature'},
                        [os.path.join(out, f'mpipi_{T}K{s}.csv') for T in TEMPS for s in ('', '_fdt')])
    with pd.option_context('display.width', 200, 'display.max_columns', 20):
        print(summ[summ['T'] == 'integral'].drop(columns='T').to_string(index=False, float_format=lambda v: f'{v:.3g}'))
        print(seg.pivot(index='metric', columns='segment', values='ratio_observed_to_trapezoid').round(3).to_string())
    print(f'saved {path} ({val["UniProt_ID"].nunique()} proteins, {len(val)} rows)')


if __name__ == '__main__':
    main()
