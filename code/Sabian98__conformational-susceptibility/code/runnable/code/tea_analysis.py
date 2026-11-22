"""
tea_analysis.py — what TEA's hydration dH/dT (tea_fdt.py) does to CALVADOS-2 susceptibility at 300 K.

chi_TEA(gamma) = chi + gamma * dchi   (chi(Rg) as % of <Rg> per K; Arg kappa from CHARMM36m, the paper's choice;
                                       CGenFF Arg as a check)
  A  share of chains that contract on heating at gamma = 0..4; gamma* where each chain switches sign
  B  by CALVADOS-2 nu bin: chi, dchi, chi_TEA(2), share contracting
  C  ranking: Spearman of chi_TEA(gamma) with chi, and with Mpipi-Recharged / Mpipi-GG chi at 300 K; leftover
     (own nu spline + log N) correlations
  D  composition: Spearman of dchi and gamma* with hydrophobic fraction, pi fraction, FCR, hydropathy
  E  by kingdom: share contracting at gamma = 2, median gamma*
  F  expanded charged chains (nu > 0.63): chi, dchi, chi_TEA(2)
Inputs: Analysis_FDT/tea_300K/tea_dchi_300K.csv, Analysis_FDT/fdt_results_300K.csv, BENDER_BIO, and
Analysis_FDT_RECHARGED / Analysis_FDT_MPIPI at 300 K. Output: Analysis_FDT/tea_300K/*.csv and summary.txt.
"""
import sys

sys.dont_write_bytecode = True

import os  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402

TLA = (__import__('os').environ['CS_DATA'] + '/BENDER-MT/Traj Level Analysis')
BIO_CSV = (__import__('os').environ['CS_DATA'] + '/BENDER-BIO/BENDER_BIO.csv')
OUT = os.path.join(TLA, 'Analysis_FDT', 'tea_300K')
BINS = [0.30, 0.45, 0.50, 0.54, 0.58, 0.63, 1.0]
GAMMAS = [0, 1, 2, 3, 4]


def spline(x, knots=5):
    q = np.quantile(x, np.linspace(0, 1, knots + 2)[1:-1])
    return np.column_stack([x, x ** 2, x ** 3] + [np.clip(x - t, 0, None) ** 3 for t in q])


def leftover(y, nu, N):
    X = np.column_stack([np.ones(len(y)), spline(nu), np.log(N)])
    return y - X @ np.linalg.lstsq(X, y, rcond=None)[0]


def rel(path, extra=()):
    d = pd.read_csv(path)
    d = d[d.metric == 'rg'].set_index('UniProt_ID')
    return 100 * d.chi_total / d.A_mean


def load(out=OUT):
    t = pd.read_csv(os.path.join(out, 'tea_dchi_300K.csv')).set_index('UniProt_ID')
    c = pd.read_csv(os.path.join(TLA, 'Analysis_FDT', 'fdt_results_300K.csv'))
    c = c[c.metric == 'rg'].set_index('UniProt_ID')
    bio = pd.read_csv(BIO_CSV, usecols=['UniProt_ID', 'nu', 'nu_fit_r2', 'length', 'FCR', 'f_pi', 'f_hydrophobic',
                                        'hydropathy_avg', 'kingdom'], low_memory=False).set_index('UniProt_ID')
    d = pd.DataFrame({'chi': 100 * c.chi_total / c.A_mean,
                      'dchi': t.dchi_rg_rel_charmm36, 'dchi_cgenff': t.dchi_rg_rel_cgenff}).dropna().join(bio)
    g = np.where((d.chi > 0) & (d.dchi < 0), -d.chi / d.dchi, np.where(d.chi <= 0, 0.0, np.inf))
    d['gamma_star'] = g
    for k in GAMMAS:
        d[f'chi_g{k}'] = d.chi + k * d.dchi
    return d


def main():
    os.makedirs(OUT, exist_ok=True)
    d = load()
    T, lines = {}, [f'{len(d)} CALVADOS-2 proteins at 300 K; chi(Rg) and TEA dchi in % of <Rg> per K '
                    '(Arg kappa CHARMM36m; CGenFF Arg in the *_cgenff columns)']
    T['A_contract_by_gamma'] = pd.DataFrame([{
        'gamma': k, 'frac_contract': (d[f'chi_g{k}'] < 0).mean(), 'median_chi': d[f'chi_g{k}'].median(),
        'frac_contract_cgenff': ((d.chi + k * d.dchi_cgenff) < 0).mean()} for k in GAMMAS])
    gs = d.gamma_star.replace(np.inf, np.nan)
    T['A_gamma_star'] = pd.DataFrame([{'n': len(d), 'already_negative': (d.gamma_star == 0).mean(),
                                       'never_switches': np.isinf(d.gamma_star).mean(),
                                       'median_gamma_star': gs.median(), 'q25': gs.quantile(0.25), 'q75': gs.quantile(0.75),
                                       'frac_switch_below_2': (d.gamma_star < 2).mean(),
                                       'frac_switch_below_3': (d.gamma_star < 3).mean()}])
    b = pd.cut(d.nu, BINS)
    g = d.groupby(b, observed=True)
    T['B_by_nu'] = pd.DataFrame({'n': g.size(), 'chi': g.chi.median(), 'dchi': g.dchi.median(),
                                 'chi_TEA_g2': g.chi_g2.median(), 'frac_contract_g2': g.chi_g2.apply(lambda x: (x < 0).mean()),
                                 'median_gamma_star': g.gamma_star.median()}).reset_index().astype({'nu': str})
    # C ranking
    rows = []
    for k in GAMMAS[1:]:
        rows.append({'comparison': 'chi_TEA vs chi (CALVADOS-2)', 'gamma': k, 'n': len(d),
                     'spearman': spearmanr(d.chi, d[f'chi_g{k}'])[0]})
    other = {'Mpipi-Recharged': rel(os.path.join(TLA, 'Analysis_FDT_RECHARGED', 'recharged_fdt_300K.csv')),
             'Mpipi-GG': rel(os.path.join(TLA, 'Analysis_FDT_MPIPI', 'mpipi_300K_fdt.csv'))}
    for name, o in other.items():
        ids = d.index.intersection(o.index)
        for k in GAMMAS:
            rows.append({'comparison': f'chi_TEA vs {name} chi', 'gamma': k, 'n': len(ids),
                         'spearman': spearmanr(d.loc[ids, f'chi_g{k}'], o[ids])[0]})
    clean = d[d.nu_fit_r2 >= 0.99]
    L0 = leftover(clean.chi.to_numpy(), clean.nu.to_numpy(), clean.length.to_numpy())
    for k in GAMMAS[1:]:
        Lk = leftover(clean[f'chi_g{k}'].to_numpy(), clean.nu.to_numpy(), clean.length.to_numpy())
        rows.append({'comparison': 'leftover chi_TEA vs leftover chi', 'gamma': k, 'n': len(clean),
                     'spearman': spearmanr(L0, Lk)[0], 'pearson': np.corrcoef(L0, Lk)[0, 1]})
    T['C_ranking'] = pd.DataFrame(rows)
    # D composition
    rows = []
    for col in ('f_hydrophobic', 'f_pi', 'FCR', 'hydropathy_avg', 'nu', 'length'):
        x = d[col]
        rows.append({'descriptor': col, 'rho_dchi': spearmanr(x, d.dchi, nan_policy='omit')[0],
                     'rho_dchi_per_residue_hydrophobic': np.nan,
                     'rho_chi_TEA_g2': spearmanr(x, d.chi_g2, nan_policy='omit')[0],
                     'rho_chi': spearmanr(x, d.chi, nan_policy='omit')[0]})
    T['D_composition'] = pd.DataFrame(rows).drop(columns='rho_dchi_per_residue_hydrophobic')
    # E kingdom
    g = d.groupby('kingdom')
    T['E_by_kingdom'] = pd.DataFrame({'n': g.size(), 'chi': g.chi.median(), 'dchi': g.dchi.median(),
                                      'frac_contract_g0': g.chi.apply(lambda x: (x < 0).mean()),
                                      'frac_contract_g2': g.chi_g2.apply(lambda x: (x < 0).mean()),
                                      'median_gamma_star': g.gamma_star.median()}).reset_index()
    # F expanded
    e = d[d.nu > 0.63]
    T['F_expanded'] = pd.DataFrame([{'n': len(e), 'median_FCR': e.FCR.median(), 'chi': e.chi.median(),
                                     'dchi': e.dchi.median(), 'chi_TEA_g2': e.chi_g2.median(),
                                     'frac_contract_g0': (e.chi < 0).mean(), 'frac_contract_g2': (e.chi_g2 < 0).mean()}])
    d.to_csv(os.path.join(OUT, 'tea_per_protein.csv'))
    for name, t in T.items():
        t.to_csv(os.path.join(OUT, f'{name}.csv'), index=False)
        lines += ['', f'== {name}', t.to_string(index=False, float_format=lambda v: f'{v:.3f}')]
    open(os.path.join(OUT, 'summary.txt'), 'w').write('\n'.join(lines) + '\n')
    print('\n'.join(lines))


if __name__ == '__main__':
    main()
