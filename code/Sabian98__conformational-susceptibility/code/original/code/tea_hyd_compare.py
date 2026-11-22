"""
tea_hyd_compare.py — does the hydrophilic-solvation correction transfer between force fields?

Per model at 300 K (tea_fdt.py / tea_crossff.py --classes, column 'hyd'): chi(Rg) and the hydrophilic term dchi_hyd
per unit gamma_hyd (% of <Rg> per K). The Wuttke calibration (tea_calibrate.py fit --free hyd --tmax 320, on integrase,
CspM34 and both ProTa segments; leave-one-out folds within +/-7%) gives gamma_hyd = 0.358 +/- 0.026 for CALVADOS-2 and
0.209 +/- 0.018 for Mpipi-Recharged; Mpipi-GG has no Wuttke runs. Comparisons:
  A  rank agreement of the hydrophilic term between models
  B  total response chi + s_m dchi_hyd with s_m set so that SD(s_m dchi_hyd) = f SD(chi) in each model, for f from 0 to 2
     and for f_cal, the ratio the calibrated Mpipi-Recharged term has on its own proteome
  C  CALVADOS-2 and Mpipi-Recharged each at its own calibrated gamma_hyd, with the TEA term at the prior the calibration
     held it at (gamma_tea 2.0 / 0.4; corrected 2026-10-10, it was left out before): agreement of the total response,
     share of chains contracting on heating (by nu bin and kingdom); Mpipi-GG at the mean relative size of the other two
  D  composition: Spearman of each model's hydrophilic term with FCR, pi, hydropathy, length and nu
(nu bins in this script use CALVADOS-2 nu for all three models.)
Output: Analysis_FDT/tea_300K/hyd/*.csv, summary.txt.
"""
import sys

sys.dont_write_bytecode = True

import os  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402

HERE = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, HERE)
import tea_crossff_compare as X  # noqa: E402

TEA_DIR = X.TEA_DIR
OUT = os.path.join(TEA_DIR, 'hyd')
FILES = {'C': 'tea_dchi_classes_300K.csv', 'M': 'tea_dchi_gg_classes_300K.csv', 'R': 'tea_dchi_recharged_classes_300K.csv'}
GAMMA_HYD = {'C': 0.358, 'R': 0.209}     # Wuttke calibration, 278-320 K
GAMMA_TEA = {'C': 2.0, 'R': 0.4}         # TEA term held at this prior while gamma_hyd was fitted
TEA_CLASSES = ['aliph', 'pro', 'arom', 'argC36', 'rest']
BINS = [0.30, 0.45, 0.50, 0.54, 0.58, 0.63, 1.0]


def load(tea_dir=TEA_DIR):
    base = X.load('charmm36', tea_dir)
    out = {}
    for k, f in FILES.items():
        t = pd.read_csv(os.path.join(tea_dir, f)).set_index('UniProt_ID')
        out[k] = pd.DataFrame({'hyd': t['dchi_rg_rel_hyd'], 'tea': sum(t[f'dchi_rg_rel_{c}'] for c in TEA_CLASSES)}
                              ).join(base[k][['chi']], how='inner').dropna()
    return out


def main(tea_dir=TEA_DIR, out=OUT):
    os.makedirs(out, exist_ok=True)
    D = load(tea_dir)
    bio = pd.read_csv(X.BIO_CSV, usecols=['UniProt_ID', 'nu', 'FCR', 'f_pi', 'hydropathy_avg', 'length', 'kingdom'],
                      low_memory=False).set_index('UniProt_ID')
    trio = D['C'].index.intersection(D['M'].index).intersection(D['R'].index)
    T3 = {k: d.loc[trio] for k, d in D.items()}
    lines = [f'{len(trio)} chains in all three models; hydrophilic term per unit gamma_hyd, % of <Rg> per K']
    # A
    A = pd.DataFrame([{'pair': f'{X.NAMES[a]} vs {X.NAMES[b]}', 'n': len(trio), 'rho_chi': spearmanr(T3[a].chi, T3[b].chi)[0],
                       'rho_hyd': spearmanr(T3[a].hyd, T3[b].hyd)[0]} for a, b in X.PAIRS])
    # f_cal from the calibrated Mpipi-Recharged term on its own full set
    fcal = {X.NAMES[k]: float((g * D[k].hyd).std() / D[k].chi.std()) for k, g in GAMMA_HYD.items()}
    rows = []
    for f in [0.0, 0.25, 0.5, 1.0, 2.0] + sorted(fcal.values()):
        tt = {k: T3[k].chi + (f * T3[k].chi.std() / T3[k].hyd.std()) * T3[k].hyd for k in T3}
        r = {'hyd_spread_over_chi_spread': round(f, 3)}
        for a, b in X.PAIRS:
            r[f'{X.NAMES[a]} vs {X.NAMES[b]}'] = spearmanr(tt[a], tt[b])[0]
            r[f'same_sign {X.NAMES[a]}-{X.NAMES[b]}'] = (np.sign(tt[a]) == np.sign(tt[b])).mean()
        rows.append(r)
    B = pd.DataFrame(rows)
    # C: each calibrated model at its own gamma_hyd; Mpipi-GG at the mean relative size of the other two
    fM = float(np.mean(list(fcal.values())))
    g_eff = dict(GAMMA_HYD, M=fM * D['M'].chi.std() / D['M'].hyd.std())
    fT = float(np.mean([(GAMMA_TEA[k] * D[k].tea).std() / D[k].chi.std() for k in GAMMA_TEA]))
    g_tea = dict(GAMMA_TEA, M=fT * D['M'].chi.std() / D['M'].tea.std())
    tot = {k: D[k].chi + g_tea[k] * D[k].tea + g_eff[k] * D[k].hyd for k in D}
    rows = []
    for k in D:
        Tb = tot[k].to_frame('tot').join(D[k].chi).join(bio, how='left')
        rows.append({'model': X.NAMES[k], 'gamma_hyd': g_eff[k], 'gamma_tea': g_tea[k], 'calibrated': k in GAMMA_HYD, 'group': 'all', 'n': len(Tb),
                     'frac_contract': (Tb.tot < 0).mean(), 'frac_contract_chi_only': (Tb.chi < 0).mean(),
                     'median_total': Tb.tot.median(), 'median_chi': Tb.chi.median(), 'rho_total_vs_chi': spearmanr(Tb.tot, Tb.chi)[0]})
        for col, grp in (('nu_bin', pd.cut(Tb.nu, BINS).astype(str)), ('kingdom', Tb.kingdom)):
            for key, idx in Tb.groupby(grp).groups.items():
                rows.append({'model': X.NAMES[k], 'gamma_hyd': g_eff[k], 'calibrated': k in GAMMA_HYD,
                             'group': f'{col}={key}', 'n': len(idx), 'frac_contract': (Tb.tot[idx] < 0).mean(),
                             'frac_contract_chi_only': (Tb.chi[idx] < 0).mean(), 'median_total': Tb.tot[idx].median(),
                             'median_chi': Tb.chi[idx].median()})
    Cc = pd.DataFrame(rows)
    agree = []
    for a, b in X.PAIRS:
        for set_name, ids in (('pair', D[a].index.intersection(D[b].index)), ('three-model set', trio)):
            agree.append({'pair': f'{X.NAMES[a]} vs {X.NAMES[b]}', 'set': set_name, 'n': len(ids),
                          'rho_chi': spearmanr(D[a].chi[ids], D[b].chi[ids])[0],
                          'rho_total_calibrated': spearmanr(tot[a][ids], tot[b][ids])[0],
                          'same_sign_total': (np.sign(tot[a][ids]) == np.sign(tot[b][ids])).mean()})
    Ag = pd.DataFrame(agree)
    # D: composition
    rows = []
    for k, d in D.items():
        b = bio.reindex(d.index)
        for c in ['FCR', 'f_pi', 'hydropathy_avg', 'nu', 'length']:
            rows.append({'model': X.NAMES[k], 'descriptor': c, 'rho_hyd': spearmanr(b[c], d.hyd, nan_policy='omit')[0],
                         'rho_chi': spearmanr(b[c], d.chi, nan_policy='omit')[0]})
    Dd = pd.DataFrame(rows).pivot(index='descriptor', columns='model', values='rho_hyd').reset_index()
    T = {'A_rank_hyd': A, 'B_total_vs_relative_size': B, 'C_calibrated_by_model': Cc, 'C_calibrated_agreement': Ag,
         'D_composition': Dd}
    lines.append('calibrated terms: SD(gamma_hyd dchi_hyd)/SD(chi) = '
                 + ', '.join(f'{v:.2f} ({k})' for k, v in fcal.items()))
    for name, t in T.items():
        t.to_csv(os.path.join(out, f'{name}.csv'), index=False)
        lines += ['', f'== {name}', t.to_string(index=False, float_format=lambda v: f'{v:.3f}')]
    open(os.path.join(out, 'summary.txt'), 'w').write('\n'.join(lines) + '\n')
    print('\n'.join(lines))


if __name__ == '__main__':
    main(*sys.argv[1:3])
