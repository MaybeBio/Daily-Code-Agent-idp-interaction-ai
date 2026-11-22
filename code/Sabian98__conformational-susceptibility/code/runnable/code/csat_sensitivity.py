"""
csat_sensitivity.py — how the headline results of the calibrated model depend on the saturation c_sat (2026-10-26).
For each c_sat in csat_scan.CSATS, the hydrophilic strength is the one refitted to the four calibration constructs at that
c_sat (saturating_forcing.py, variant S1, TEA strength held), Mpipi-GG strengths follow the SD-ratio rule of
total_response_evidence.load, and the per-chain terms come from csat_scan.py. Recomputed per c_sat:
  held-out mean absolute error over the five series to which no parameter is fitted (two alpha-synuclein loops, p53,
  sfAFP, Aap) and the PNt error, from saturating_<model>.csv
  share of chains that compact on heating, median total response, medians below nu 0.45 and above 0.63
  rank agreement of the total response between force fields (CALVADOS-2 against Mpipi-Recharged on all shared chains,
  against Mpipi-GG on the 1,107 chains of all three)
  rank correlation of heat sensitivity with nu, and the rank partial correlation of stronger compaction on heating with
  the Asp, Glu, Lys and Arg fraction at fixed nu and log N
Check: c_sat 2 must reproduce total_response_evidence.load() to 1e-9.
    BENDER_REANALYSIS_DIR=... python csat_sensitivity.py      (openmpipi env)
Output: Traj Level Analysis/Analysis_FDT/capped_300K/csat_scan/{sensitivity.csv, sensitivity.txt}
"""
import sys

sys.dont_write_bytecode = True

import os  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

HERE = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, os.environ.get('BENDER_REANALYSIS_DIR', HERE))
sys.path.insert(0, HERE)
import csat_scan as CS  # noqa: E402
import total_response_evidence as T  # noqa: E402

BM = (__import__('os').environ['CS_DATA'] + '/BENDER-MT')
SAT = os.path.join(BM, 'heldout_temperature', 'results', 'exploratory', 'saturating')
FILES = {'C': 'calvados', 'R': 'recharged', 'M': 'gg'}
NAME = {'C': 'CALVADOS-2', 'R': 'Mpipi-Recharged', 'M': 'Mpipi-GG'}
GT = {'C': 2.0, 'R': 0.4}


def heldout(cs):
    """(gamma_h, MAE over the five unfitted series, PNt MAE) per directly calibrated model at this c_sat."""
    out = {}
    for k, f in (('C', 'calvados'), ('R', 'recharged')):
        t = pd.read_csv(os.path.join(SAT, f'saturating_{f}.csv'))
        r = t[(t.variant == 'S1') & np.isclose(t.c_sat, cs)].iloc[0]
        out[k] = (float(r.gamma_h), (3 * r.MAE_aSyn_p53 + 2 * r.MAE_controls) / 5, float(r.MAE_PNt))
    return out


def totals(cs, base, ho):
    D = {}
    for k in ('C', 'R', 'M'):
        t = pd.read_csv(os.path.join(CS.OUT, f'csat_dchi_{FILES[k]}.csv')).set_index('UniProt_ID')
        d = base[k][['chi', 'nu', 'N']].join(t[[f'dchi_rg_rel_chg_cap{cs:g}', f'dchi_rg_rel_tea_cap{cs:g}']].rename(
            columns={f'dchi_rg_rel_chg_cap{cs:g}': 'hyd', f'dchi_rg_rel_tea_cap{cs:g}': 'tea'}), how='inner')
        D[k] = d
    g = {k: ho[k][0] for k in ('C', 'R')}
    g['M'] = np.mean([(g[k] * D[k].hyd).std() / D[k].chi.std() for k in ('C', 'R')]) * D['M'].chi.std() / D['M'].hyd.std()
    gt = dict(GT)
    gt['M'] = np.mean([(gt[k] * D[k].tea).std() / D[k].chi.std() for k in ('C', 'R')]) * D['M'].chi.std() / D['M'].tea.std()
    for k in D:
        D[k]['total'] = D[k].chi + gt[k] * D[k].tea + g[k] * D[k].hyd
    return D, g


def main():
    base, gbase = T.load()
    bio = T.bio_table()
    rows, L = [], ['Sensitivity of the calibrated model to the saturation c_sat (strength refitted at each value)', '']
    for cs in CS.CSATS:
        ho = heldout(cs)
        D, g = totals(cs, base, ho)
        if np.isclose(cs, 2.0):            # with the paper's strengths, c_sat 2 must give the stored totals exactly
            ho2 = {k: (gbase[k],) + ho[k][1:] for k in ('C', 'R')}
            D2, _ = totals(cs, base, ho2)
            for k in D2:
                dev = (D2[k].total - base[k].total.reindex(D2[k].index)).abs().max()
                if dev > 1e-9:
                    raise SystemExit(f'c_sat 2 does not reproduce the stored totals for {k} ({dev:.1e})')
            L.append('check: c_sat 2 with the paper strengths reproduces the stored totals (< 1e-9); the scan refit differs by '
                     + ', '.join(f'{NAME[k]} {ho[k][0] - gbase[k]:+.1e}' for k in ('C', 'R')) + ' in gamma_hyd')
        r = {'c_sat': cs, 'gamma_C': g['C'], 'gamma_R': g['R'], 'gamma_M': g['M'],
             'heldout_MAE_C': ho['C'][1], 'heldout_MAE_R': ho['R'][1], 'PNt_MAE_C': ho['C'][2], 'PNt_MAE_R': ho['R'][2]}
        for k in ('C', 'M', 'R'):
            d = D[k]
            r[f'compact_pct_{k}'] = 100 * (d.total < 0).mean()
            r[f'median_{k}'] = d.total.median()
            r[f'median_nu_lt045_{k}'] = d.total[(d.nu > 0.30) & (d.nu < 0.45)].median()   # the paper's most compact bin
            r[f'median_nu_gt063_{k}'] = d.total[d.nu > 0.63].median()
            r[f'rho_heat_nu_{k}'] = stats.spearmanr(-d.total, d.nu, nan_policy='omit').statistic
        trio = D['C'].index.intersection(D['M'].index).intersection(D['R'].index)
        cr = D['C'].index.intersection(D['R'].index)
        r['agree_CR'] = stats.spearmanr(D['C'].total[cr], D['R'].total[cr]).statistic
        r['agree_CM'] = stats.spearmanr(D['C'].total[trio], D['M'].total[trio]).statistic
        for k in ('C', 'R'):
            x = D[k].join(bio[['f_DEK', 'f_R']], how='inner').dropna(subset=['nu', 'N', 'f_DEK', 'f_R'])
            Z = np.column_stack([x.nu, np.log(x.N)])
            r[f'partial_charge_{k}'] = T.partial(-x.total.values, (x.f_DEK + x.f_R).values, Z)
        rows.append(r)
    S = pd.DataFrame(rows)
    os.makedirs(CS.OUT, exist_ok=True)
    S.to_csv(os.path.join(CS.OUT, 'sensitivity.csv'), index=False)
    show = S.set_index('c_sat').T
    L += [show.round(3).to_string(), '']
    open(os.path.join(CS.OUT, 'sensitivity.txt'), 'w').write('\n'.join(L) + '\n')
    print('\n'.join(L))


if __name__ == '__main__':
    main()
