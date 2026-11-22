"""
capped_vs_linear.py — what changes in the proteome results between the linear calibration (A0, the first manuscript
versions) and the capped main model (from 2026-10-24). Recombines stored per-chain terms through
total_response_evidence.load() with each forcing; no trajectories read. Same quantities as s1_proteome.py: share
contracting, medians overall and by compaction, rank agreement of the two totals, agreement between force fields
(the 1,107 chains of all three models, and every chain of both for CALVADOS-2 against Mpipi-Recharged), partial rank
correlations with composition at fixed compaction and length, kingdom medians, the chains that no longer contract.
    BENDER_REANALYSIS_DIR=... python capped_vs_linear.py       (openmpipi env)
Output: Traj Level Analysis/Analysis_FDT/total_response/capped_vs_linear/{summary.txt, by_nu.csv, agreement.csv,
        partials.csv, chains_<model>.csv}
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
import s1_proteome as S1P  # noqa: E402   (its partial() and bins; it pins the evidence module to A0 on import)
import total_response_evidence as T  # noqa: E402

OUT = os.path.join(T.TLA, 'Analysis_FDT', 'total_response', 'capped_vs_linear')
NAME = S1P.NAME
DESCS = ('f_DEK', 'f_pi', 'hydropathy_avg', 'f_K')


def totals(forcing):
    T.FORCING = forcing                  # load() and gamma_table() read the module global at call time
    D, g = T.load()
    return D, g


def main():
    os.makedirs(OUT, exist_ok=True)
    DL, gl = totals('A0')
    DC, gc = totals('capped')
    T.FORCING = 'capped'
    bio = T.bio_table()
    L = ['Linear calibration (A0) against the capped main model, proteome results (% per K)',
         'strengths: linear ' + ', '.join(f'{k} {v:.3f}' for k, v in gl.items()),
         '           capped ' + ', '.join(f'{k} {v:.3f}' for k, v in gc.items()), '']
    S, rows = {}, []
    for k in ('C', 'R', 'M'):
        d = DL[k][['nu', 'chi', 'total']].rename(columns={'total': 'total_linear'}).join(
            DC[k][['chi', 'total']].rename(columns={'chi': 'chi_c', 'total': 'total_capped'}), how='inner')
        S[k] = d.drop(columns='chi_c')
        S[k].to_csv(os.path.join(OUT, f'chains_{k}.csv'))
        rho = stats.spearmanr(d.total_linear, d.total_capped).statistic
        L.append(f'{NAME[k]} ({len(d)} chains): contracting {100 * (d.total_linear < 0).mean():.1f}% -> '
                 f'{100 * (d.total_capped < 0).mean():.1f}%; median {d.total_linear.median():+.3f} -> '
                 f'{d.total_capped.median():+.3f}; rank agreement of the two totals {rho:.2f}')
        for lab, x in d.groupby(pd.cut(d.nu, S1P.BINS), observed=True):
            rows.append({'model': NAME[k], 'nu_bin': str(lab), 'n': len(x), 'linear': x.total_linear.median(),
                         'capped': x.total_capped.median()})
    by = pd.DataFrame(rows)
    by.to_csv(os.path.join(OUT, 'by_nu.csv'), index=False)
    L += ['', 'median by compaction bin (linear -> capped)', by.round(3).to_string(index=False), '']
    agr, par = [], []
    trio = S['C'].index.intersection(S['M'].index).intersection(S['R'].index)
    for a_, b_ in (('C', 'R'), ('C', 'M'), ('M', 'R')):
        ids = S[a_].index.intersection(S[b_].index) if (a_, b_) == ('C', 'R') else trio
        r0 = stats.spearmanr(S[a_].total_linear[ids], S[b_].total_linear[ids]).statistic
        r1 = stats.spearmanr(S[a_].total_capped[ids], S[b_].total_capped[ids]).statistic
        L.append(f'agreement {NAME[a_]} vs {NAME[b_]} ({len(ids)} chains): linear {r0:.2f} -> capped {r1:.2f}')
        agr.append({'pair': f'{NAME[a_]} vs {NAME[b_]}', 'n': len(ids), 'linear': r0, 'capped': r1})
    pd.DataFrame(agr).to_csv(os.path.join(OUT, 'agreement.csv'), index=False)
    L.append('')
    for k in ('C', 'R'):
        d = S[k].join(bio[list(DESCS) + ['kingdom', 'length', 'f_R']], how='inner').dropna(subset=['nu', *DESCS])
        Z = np.column_stack([d.nu, np.log(d.length)])
        parts = []
        for desc in DESCS:
            p0 = S1P.partial(-d.total_linear.values, d[desc].values, Z)
            p1 = S1P.partial(-d.total_capped.values, d[desc].values, Z)
            parts.append(f'{desc} {p0:+.2f} -> {p1:+.2f}')
            par.append({'model': NAME[k], 'descriptor': desc, 'linear': p0, 'capped': p1})
        L.append(f'{NAME[k]} partial rank correlation with stronger contraction | nu, log N: ' + ', '.join(parts))
        L.append(f'{NAME[k]} kingdom medians linear {d.groupby("kingdom").total_linear.median().reindex(S1P.KINGDOMS).round(2).to_dict()}')
        L.append(f'{NAME[k]} kingdom medians capped {d.groupby("kingdom").total_capped.median().reindex(S1P.KINGDOMS).round(2).to_dict()}')
        x = d[d.total_capped >= 0]
        L.append(f'{NAME[k]} chains not contracting under the capped model: {len(x)}; median nu {x.nu.median():.3f} (all '
                 f'{d.nu.median():.3f}), D+E+K+R fraction {(x.f_DEK + x.f_R).median():.3f} (all {(d.f_DEK + d.f_R).median():.3f}), '
                 f'hydropathy {x.hydropathy_avg.median():.2f} (all {d.hydropathy_avg.median():.2f})')
        L.append('')
    pd.DataFrame(par).to_csv(os.path.join(OUT, 'partials.csv'), index=False)
    open(os.path.join(OUT, 'summary.txt'), 'w').write('\n'.join(L) + '\n')
    print('\n'.join(L))


if __name__ == '__main__':
    main()
