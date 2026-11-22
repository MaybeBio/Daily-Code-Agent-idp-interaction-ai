"""
s1_proteome.py — do the proteome conclusions survive the charged-group form of the hydrophilic forcing (S1, Amendment 3
of the held-out pre-registration)? Recombines stored per-chain terms, no trajectories read.

Frozen total (as total_response_evidence.load): chi + gamma_TEA TEA + gamma_hyd hyd, hyd on D E K R N Q S T H G.
S1 total: chi + gamma_TEA TEA + gamma_chg (hyd_DE + hyd_K + hyd_R), gamma_chg fitted on the Wuttke constructs alone
(0.469789 CALVADOS-2, 0.294234 Mpipi-Recharged, as registered); Mpipi-GG at the mean relative size of the other two, as
the frozen total. The group-split hydrophilic terms come from total_response/all_observables (tea_fdt / tea_crossff
--classes); they are first checked to add up to the hydrophilic term the frozen totals use.
Compared: share contracting, medians overall and by compaction, rank agreement of frozen and S1 totals, agreement between
force fields, partial rank correlations with composition at fixed compaction and length, kingdom medians, and the
composition of chains that no longer contract.
    BENDER_REANALYSIS_DIR=... python s1_proteome.py       (openmpipi env)
Output: Traj Level Analysis/Analysis_FDT/total_response/s1_sensitivity/{summary.txt, by_nu.csv, chains_<model>.csv}
"""
import sys

sys.dont_write_bytecode = True

import os  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

HERE = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, HERE)
os.environ['BENDER_FORCING'] = 'A0'          # this sensitivity compares S1 against the linear A0 model by design
import total_response_evidence as T  # noqa: E402

TR = os.path.join(T.TLA, 'Analysis_FDT', 'total_response')
AO = os.path.join(TR, 'all_observables')
OUT = os.path.join(TR, 's1_sensitivity')
FILES = {'C': 'tea_dchi_classes_300K.csv', 'M': 'tea_dchi_gg_classes_300K.csv', 'R': 'tea_dchi_recharged_classes_300K.csv'}
G_CHG = {'C': 0.469789, 'R': 0.294234}
NAME = {'C': 'CALVADOS-2', 'M': 'Mpipi-GG', 'R': 'Mpipi-Recharged'}
BINS = [0.30, 0.45, 0.50, 0.54, 0.58, 0.63, 1.0]
KINGDOMS = ['Bacteria', 'Plants', 'Fungi', 'Mammals', 'Protists', 'Viruses']


def rk(x):
    return stats.rankdata(x)


def partial(y, x, Z):
    """Spearman partial correlation of y and x given the columns of Z (rank-residualized least squares)."""
    Zr = np.column_stack([rk(z) for z in Z.T] + [np.ones(len(y))])
    ry = rk(y) - Zr @ np.linalg.lstsq(Zr, rk(y), rcond=None)[0]
    rx = rk(x) - Zr @ np.linalg.lstsq(Zr, rk(x), rcond=None)[0]
    return float(np.corrcoef(ry, rx)[0, 1])


def main():
    os.makedirs(OUT, exist_ok=True)
    T.FORCING = 'A0'                               # read by load() at call time, whatever the import order
    D, g = T.load()
    bio = T.bio_table()
    L = ['S1 sensitivity of the proteome results (hydrophilic forcing on D, E, K, R only; % per K)', '']
    S1 = {}
    for k in ('C', 'R', 'M'):
        a = pd.read_csv(os.path.join(AO, FILES[k])).set_index('UniProt_ID')
        grp = a[[f'dchi_rg_rel_hyd_{x}' for x in ('DE', 'K', 'R', 'pol', 'G')]].sum(axis=1)
        d = D[k].join(a[['dchi_rg_rel_hyd', 'dchi_rg_rel_hyd_DE', 'dchi_rg_rel_hyd_K', 'dchi_rg_rel_hyd_R']], how='inner')
        e1 = float((grp.reindex(d.index) - d.dchi_rg_rel_hyd).abs().max())
        e2 = float((d.dchi_rg_rel_hyd - d.hyd).abs().max())
        L.append(f'{NAME[k]}: {len(d)} chains; groups add up to the hydrophilic term to {e1:.1e}, and it equals the frozen '
                 f'totals\' term to {e2:.1e}')
        if max(e1, e2) > 1e-6:
            raise ValueError(f'{k}: hydrophilic terms disagree')
        d['chg'] = d.dchi_rg_rel_hyd_DE + d.dchi_rg_rel_hyd_K + d.dchi_rg_rel_hyd_R
        S1[k] = d
    # Mpipi-GG strength: mean relative size (SD of the forcing term over SD of chi) of the two calibrated models
    f = np.mean([(G_CHG[k] * S1[k].chg).std() / S1[k].chi.std() for k in ('C', 'R')])
    G_CHG['M'] = f * S1['M'].chi.std() / S1['M'].chg.std()
    L.append(f'Mpipi-GG gamma_chg by relative size: {G_CHG["M"]:.3f}')
    L.append('')
    rows = []
    for k, d in S1.items():
        d['total_S1'] = d.chi + d.TEA + G_CHG[k] * d.chg
        d.to_csv(os.path.join(OUT, f'chains_{k}.csv'), columns=['nu', 'chi', 'TEA', 'H', 'total', 'chg', 'total_S1'])
        rho = stats.spearmanr(d.total, d.total_S1).statistic
        L.append(f'{NAME[k]}: contracting {100 * (d.total < 0).mean():.1f}% (frozen) -> {100 * (d.total_S1 < 0).mean():.1f}% '
                 f'(S1); median {d.total.median():+.3f} -> {d.total_S1.median():+.3f}; rank agreement frozen vs S1 {rho:.2f}')
        for lab, x in d.groupby(pd.cut(d.nu, BINS), observed=True):
            rows.append({'model': NAME[k], 'nu_bin': str(lab), 'n': len(x), 'frozen': x.total.median(),
                         'S1': x.total_S1.median()})
    by = pd.DataFrame(rows)
    by.to_csv(os.path.join(OUT, 'by_nu.csv'), index=False)
    L += ['', 'median by compaction bin (frozen -> S1)', by.round(3).to_string(index=False), '']
    # agreement between force fields
    agr, par = [], []
    for a_, b_ in (('C', 'R'), ('C', 'M'), ('M', 'R')):
        ids = S1[a_].index.intersection(S1[b_].index)
        if a_ == 'C' and b_ == 'R':
            pass
        else:
            ids = ids.intersection(S1['C'].index).intersection(S1['M'].index).intersection(S1['R'].index)
        r0 = stats.spearmanr(S1[a_].total[ids], S1[b_].total[ids]).statistic
        r1 = stats.spearmanr(S1[a_].total_S1[ids], S1[b_].total_S1[ids]).statistic
        L.append(f'agreement {NAME[a_]} vs {NAME[b_]} ({len(ids)} chains): frozen {r0:.2f} -> S1 {r1:.2f}')
        agr.append({'pair': f'{NAME[a_]} vs {NAME[b_]}', 'n': len(ids), 'frozen': r0, 'S1': r1})
    pd.DataFrame(agr).to_csv(os.path.join(OUT, 'agreement.csv'), index=False)
    L.append('')
    # composition beyond compaction and length (partial Spearman), as in the Results
    for k in ('C', 'R'):
        d = S1[k].join(bio[['f_DEK', 'f_pi', 'hydropathy_avg', 'f_K', 'kingdom', 'length']], how='inner', rsuffix='_b')
        d = d.dropna(subset=['nu', 'f_DEK', 'f_pi', 'hydropathy_avg'])
        Z = np.column_stack([d.nu, np.log(d.length)])
        parts = []
        for desc in ('f_DEK', 'f_pi', 'hydropathy_avg', 'f_K'):
            p0 = partial(-d.total.values, d[desc].values, Z)
            p1 = partial(-d.total_S1.values, d[desc].values, Z)
            parts.append(f'{desc} {p0:+.2f} -> {p1:+.2f}')
            par.append({'model': NAME[k], 'descriptor': desc, 'frozen': p0, 'S1': p1})
        L.append(f'{NAME[k]} partial rank correlation with stronger contraction | nu, log N: ' + ', '.join(parts))
        km0 = d.groupby('kingdom').total.median().reindex(KINGDOMS).round(2).to_dict()
        km1 = d.groupby('kingdom').total_S1.median().reindex(KINGDOMS).round(2).to_dict()
        L.append(f'{NAME[k]} kingdom medians frozen {km0}')
        L.append(f'{NAME[k]} kingdom medians S1     {km1}')
        x = d[d.total_S1 >= 0]
        if len(x):
            L.append(f'{NAME[k]} chains not contracting under S1: {len(x)}; median nu {x.nu.median():.3f} (all '
                     f'{d.nu.median():.3f}), D+E+K+R fraction {(x.f_DEK + bio.f_R.reindex(x.index)).median():.3f} (all '
                     f'{(d.f_DEK + bio.f_R.reindex(d.index)).median():.3f}), hydropathy {x.hydropathy_avg.median():.2f} '
                     f'(all {d.hydropathy_avg.median():.2f})')
        L.append('')
    pd.DataFrame(par).to_csv(os.path.join(OUT, 'partials.csv'), index=False)
    open(os.path.join(OUT, 'summary.txt'), 'w').write('\n'.join(L) + '\n')
    print('\n'.join(L))


if __name__ == '__main__':
    main()
