"""
tea_crossff_compare.py — do the three force fields agree once the same TEA hydration forcing is added?

Per model m at 300 K: chi_m (FDT, % of <Rg> per K) and dchi_m (TEA term per unit of the model's own scale; CALVADOS-2
from tea_fdt.py with gamma in mol/kcal, Mpipi-GG / Mpipi-Recharged from tea_crossff.py with dimensionless gamma').
The scales are not calibrated for the Mpipi models, so the comparison uses quantities that do not depend on them:
  A  rank agreement of the TEA term itself (dchi) between models
  B  rank agreement of the switch point gamma* = -chi/dchi (chains that never switch ranked last)
  C  total response chi + s_m dchi with each model's scale s_m matched so that its median TEA term equals CALVADOS-2's
     at gamma = 2: rank agreement between models and share contracting; also with each s_m varied 0.5x-2x
  D  composition: Spearman of each model's dchi with hydrophobic fraction, pi fraction, FCR
Arg kappa: CHARMM36m (the paper's choice); CGenFF in the *_cgenff rows.
Output: Analysis_FDT/tea_300K/crossff/*.csv and summary.txt.
"""
import sys

sys.dont_write_bytecode = True

import os  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402

TLA = (__import__('os').environ['CS_DATA'] + '/BENDER-MT/Traj Level Analysis')
TEA_DIR = os.path.join(TLA, 'Analysis_FDT', 'tea_300K')
OUT = os.path.join(TEA_DIR, 'crossff')
BIO_CSV = (__import__('os').environ['CS_DATA'] + '/BENDER-BIO/BENDER_BIO.csv')
NAMES = {'C': 'CALVADOS-2', 'M': 'Mpipi-GG', 'R': 'Mpipi-Recharged'}
PAIRS = (('C', 'M'), ('C', 'R'), ('M', 'R'))


def chi_rel(path, ok=None):
    d = pd.read_csv(path)
    d = d[d.metric == 'rg'].set_index('UniProt_ID')
    s = 100 * d.chi_total / d.A_mean
    return s if ok is None else s[s.index.isin(ok)]


def load(arg='charmm36', tea_dir=TEA_DIR):
    ok_m = set(pd.read_csv(os.path.join(TLA, 'Analysis_FDT_MPIPI', 'mpipi_300K.csv')).query("status == 'ok'").UniProt_ID)
    rc = pd.read_csv(os.path.join(TLA, 'Analysis_FDT_RECHARGED', 'recharged_fdt_300K.csv'))
    ok_r = set(rc.loc[rc.max_bond_nm <= 1.0, 'UniProt_ID'])
    chi = {'C': chi_rel(os.path.join(TLA, 'Analysis_FDT', 'fdt_results_300K.csv')),
           'M': chi_rel(os.path.join(TLA, 'Analysis_FDT_MPIPI', 'mpipi_300K_fdt.csv'), ok_m),
           'R': chi_rel(os.path.join(TLA, 'Analysis_FDT_RECHARGED', 'recharged_fdt_300K.csv'), ok_r)}
    files = {'C': 'tea_dchi_300K.csv', 'M': 'tea_dchi_gg_300K.csv', 'R': 'tea_dchi_recharged_300K.csv'}
    dchi = {}
    for k, f in files.items():
        t = pd.read_csv(os.path.join(tea_dir, f)).set_index('UniProt_ID')
        dchi[k] = t[f'dchi_rg_rel_{arg}']
    return {k: pd.DataFrame({'chi': chi[k], 'dchi': dchi[k]}).dropna() for k in NAMES}


def gamma_star(d):
    return pd.Series(np.where((d.chi > 0) & (d.dchi < 0), -d.chi / d.dchi, np.where(d.chi <= 0, 0.0, np.inf)), index=d.index)


def rho(a, b):
    return spearmanr(a, b)[0]


def analyse(D, tag):
    rows_A, rows_B, rows_C = [], [], []
    trio = D['C'].index.intersection(D['M'].index).intersection(D['R'].index)
    s = {k: D['C'].dchi.median() * 2 / D[k].dchi.median() if k != 'C' else 2.0 for k in D}   # matched scales
    for a, b in PAIRS:
        for set_name, ids in (('pair', D[a].index.intersection(D[b].index)), ('three-model set', trio)):
            x, y = D[a].loc[ids], D[b].loc[ids]
            rows_A.append({'arg': tag, 'pair': f'{NAMES[a]} vs {NAMES[b]}', 'set': set_name, 'n': len(ids),
                           'rho_chi': rho(x.chi, y.chi), 'rho_dchi': rho(x.dchi, y.dchi)})
            ga, gb = gamma_star(x).replace(np.inf, 1e9), gamma_star(y).replace(np.inf, 1e9)
            rows_B.append({'arg': tag, 'pair': f'{NAMES[a]} vs {NAMES[b]}', 'set': set_name, 'n': len(ids),
                           'rho_gamma_star': rho(ga, gb)})
            for f in (0.5, 1.0, 2.0):
                ta, tb = x.chi + s[a] * x.dchi, y.chi + f * s[b] * y.dchi
                rows_C.append({'arg': tag, 'pair': f'{NAMES[a]} vs {NAMES[b]}', 'set': set_name, 'n': len(ids),
                               f'scale_{NAMES[b]}_x': f, 'rho_total': rho(ta, tb),
                               'contract_a': (ta < 0).mean(), 'contract_b': (tb < 0).mean(),
                               'same_sign': (np.sign(ta) == np.sign(tb)).mean()})
    scales = pd.DataFrame([{'arg': tag, 'model': NAMES[k], 'matched_scale': s[k], 'median_dchi': D[k].dchi.median(),
                            'median_chi': D[k].chi.median(), 'n': len(D[k])} for k in D])
    return pd.DataFrame(rows_A), pd.DataFrame(rows_B), pd.DataFrame(rows_C), scales


def composition(D):
    bio = pd.read_csv(BIO_CSV, usecols=['UniProt_ID', 'f_hydrophobic', 'f_pi', 'FCR', 'hydropathy_avg'],
                      low_memory=False).set_index('UniProt_ID')
    rows = []
    for k, d in D.items():
        b = bio.reindex(d.index)
        for c in b.columns:
            rows.append({'model': NAMES[k], 'descriptor': c, 'rho_dchi': spearmanr(b[c], d.dchi, nan_policy='omit')[0]})
    return pd.DataFrame(rows).pivot(index='descriptor', columns='model', values='rho_dchi').reset_index()


def main(tea_dir=TEA_DIR, out=OUT):
    os.makedirs(out, exist_ok=True)
    T, lines = {}, []
    parts = [analyse(load(arg, tea_dir), arg) for arg in ('charmm36', 'cgenff')]
    for k, name in enumerate(('A_rank_dchi', 'B_rank_gamma_star', 'C_total_matched_scale', 'scales')):
        T[name] = pd.concat([p[k] for p in parts])
    T['D_composition'] = composition(load('charmm36', tea_dir))
    for name, t in T.items():
        t.to_csv(os.path.join(out, f'{name}.csv'), index=False)
        lines += ['', f'== {name}', t.to_string(index=False, float_format=lambda v: f'{v:.3f}')]
    open(os.path.join(out, 'summary.txt'), 'w').write('\n'.join(lines) + '\n')
    print('\n'.join(lines))


if __name__ == '__main__':
    main(*sys.argv[1:3])
