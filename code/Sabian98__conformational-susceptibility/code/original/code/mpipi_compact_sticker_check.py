"""
mpipi_compact_sticker_check.py — are the chains that are compact in Mpipi-GG sticker-driven, and is that why their
susceptibility keeps rising below the theta point (CALVADOS-2 plateaus at ~0.135 % Rg/K below nu 0.50)?

Data at 300 K: Mpipi-GG = the 992-set (Analysis_FDT_MPIPI/mpipi_300K*.csv) plus the 250-protein compact set
(Analysis_FDT_MPIPI/compact/); CALVADOS-2 = Analysis_FDT/fdt_results_300K.csv (all 10,282). chi(Rg) as % of <Rg>
per K; nu = primary nu with power-law fit R^2 >= 0.99 in the model concerned. Sequence descriptors from Test 2
(test2_slope_vs_sequence.sequence_layers) plus the aromatic count.

  A  what sets compaction: Spearman of each model's nu, and of nu_Mpipi - nu_CALVADOS (paired), with descriptors
  B  who is compact: median descriptors of chains with own nu < 0.50 vs the rest, per model
  C  below theta (own nu < 0.50): partial Spearman of chi with descriptors, controlling for own nu (cubic spline)
     and log N; CV R^2 of the sticker vs charge groups on that residual (Test 2 CV code)
  D  energy-term split of chi by own-nu bin (Mpipi WF / DH / bond; CALVADOS AH / DH incl. dielectric / bond)

Output: Traj Level Analysis/Analysis_FDT_MPIPI/compact/sticker_check/ (A_compaction.csv, B_composition.csv,
C_below_theta.csv, C_below_theta_r2.csv, D_terms.csv)
Usage: python3 mpipi_compact_sticker_check.py
"""
import sys

sys.dont_write_bytecode = True  # keep the NAS code folder free of __pycache__

import os  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

HERE = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, HERE)
import test2_slope_vs_sequence as t2  # noqa: E402
import test2b_cross_model_residual as t2b  # noqa: E402

t1, L, bc = t2.t1, t2.L, t2.bc
MP = '/Volumes/User Homes/Nu Project/BENDER-MT/Traj Level Analysis/Analysis_FDT_MPIPI'
OUT_DIR = os.path.join(MP, 'compact', 'sticker_check')
R2_MIN, THETA = 0.99, 0.50
BINS = [0.30, 0.45, 0.50, 0.54, 0.58, 0.63, 1.0]
STICKERS = ['f_pi', 'f_pi_extended', 'lambda_avg', 'f_aro', 'n_aro', 'f_Y', 'f_F', 'f_W', 'f_R',
            'log1p_aromatic_pairs', 'log1p_cation_pi_pairs', 'SHD']
CHARGE = ['FCR', 'NCPR', 'uversky_R', 'f_E', 'f_D', 'f_K', 'kappa', 'SCD']
OTHER = ['hydropathy_avg', 'f_P', 'f_G', 'f_S', 'f_Q', 'f_N', 'logN']
DESC = STICKERS + CHARGE + OTHER


def mpipi_300() -> pd.DataFrame:
    rows = []
    for s, d in (('main', MP), ('compact', os.path.join(MP, 'compact'))):
        tag = 'mpipi_300K' if s == 'main' else 'mpipi_compact_300K'
        m = pd.read_csv(os.path.join(d, f'{tag}.csv')).set_index('UniProt_ID')
        f = pd.read_csv(os.path.join(d, f'{tag}_fdt.csv'))
        f = f[f['metric'] == 'rg'].set_index('UniProt_ID')
        m = m[m['status'] == 'ok']
        A = f['A_mean']
        rows.append(pd.DataFrame({'set': s, 'nu': m['nu'], 'nu_r2': m['nu_r2'], 'chi': 100 * f['chi_total'] / A,
                                  'WF': 100 * f['chi_WF'] / A, 'DH': 100 * f['chi_DH'] / A,
                                  'bond': 100 * f['chi_bond'] / A}).dropna(subset=['chi']))
    return pd.concat(rows)


def calvados_300(bio: pd.DataFrame) -> pd.DataFrame:
    f = pd.read_csv(os.path.join(L.OUT_ROOT, 'fdt_results_300K.csv'))
    f = f[f['metric'] == 'rg'].set_index('UniProt_ID')
    A = f['A_mean']
    return pd.DataFrame({'nu': bio['nu'], 'nu_r2': bio['nu_fit_r2'], 'chi': 100 * f['chi_total'] / A,
                         'AH': 100 * f['chi_AH'] / A, 'DH': 100 * f['chi_DH'] / A,
                         'dielectric': 100 * f['chi_dielectric'] / A, 'bond': 100 * f['chi_bond'] / A}).dropna(subset=['chi'])


def partial_spearman(y, x, Z):
    ry = t2b.residualise(stats.rankdata(y), Z)
    rx = t2b.residualise(stats.rankdata(x), Z)
    return np.corrcoef(ry, rx)[0, 1]


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    bio = pd.read_csv(bc.BIO_CSV, low_memory=False).set_index('UniProt_ID')
    seq, _ = t2.sequence_layers(bio)
    seq['n_aro'] = bio['f_aro'] * bio['length']
    seq['logN'] = np.log(bio['length'])
    models = {'Mpipi-GG': mpipi_300(), 'CALVADOS-2': calvados_300(bio)}
    for k, d in models.items():
        models[k] = d[d['nu_r2'] >= R2_MIN].join(seq[DESC], how='inner')
    M, C = models['Mpipi-GG'], models['CALVADOS-2']
    print(f"300 K, nu fit R2 >= {R2_MIN}: Mpipi-GG {len(M)} ({(M['set'] == 'compact').sum()} from the compact set), "
          f"CALVADOS-2 {len(C)}")

    # A: what sets compaction
    P = M[['nu']].join(C[['nu']], lsuffix='_M', rsuffix='_C', how='inner').join(seq[DESC])
    P['dnu'] = P['nu_M'] - P['nu_C']
    A = pd.DataFrame({'nu Mpipi-GG': {c: stats.spearmanr(M[c], M['nu'])[0] for c in DESC},
                      'nu CALVADOS-2 (same proteins)': {c: stats.spearmanr(P[c], P['nu_C'])[0] for c in DESC},
                      'nu Mpipi - nu CALVADOS': {c: stats.spearmanr(P[c], P['dnu'])[0] for c in DESC}})
    A.index.name = 'descriptor'
    A.to_csv(os.path.join(OUT_DIR, 'A_compaction.csv'))
    print(f"\nA  Spearman with compaction (paired n {len(P)}; negative = more compact). nu Mpipi vs CALVADOS "
          f"Spearman {stats.spearmanr(P['nu_M'], P['nu_C'])[0]:.2f}")
    print(A.reindex(A['nu Mpipi - nu CALVADOS'].abs().sort_values(ascending=False).index).round(3).to_string())

    # B: composition of the chains that are compact in each model
    rows = []
    for k, d in models.items():
        lo, hi = d[d['nu'] < THETA], d[d['nu'] >= THETA]
        for c in STICKERS + CHARGE + ['f_P', 'f_G', 'hydropathy_avg']:
            rows.append({'model': k, 'descriptor': c, 'n_compact': len(lo), 'median_compact': lo[c].median(),
                         'median_rest': hi[c].median(), 'spearman_nu_below_theta': stats.spearmanr(lo[c], lo['nu'])[0]})
    B = pd.DataFrame(rows)
    B.to_csv(os.path.join(OUT_DIR, 'B_composition.csv'), index=False)
    print(f'\nB  median descriptors, own nu < {THETA} vs the rest')
    print(B.pivot(index='descriptor', columns='model', values=['median_compact', 'median_rest']).round(3).to_string())

    # C: below theta, what sets chi at fixed nu
    rows, r2 = [], []
    for k, d in models.items():
        lo = d[d['nu'] < THETA]
        Z = np.column_stack([t2b.spline(lo['nu'].to_numpy(), knots=3), lo['logN']])
        for c in DESC:
            rows.append({'model': k, 'descriptor': c, 'n': len(lo), 'partial_spearman_chi': partial_spearman(lo['chi'], lo[c], Z),
                         'spearman_chi': stats.spearmanr(lo['chi'], lo[c])[0]})
        y = t2b.residualise(lo['chi'].to_numpy(), Z)
        rn = lambda cols: np.column_stack([t1.rank_normal(lo[c].to_numpy()) for c in cols])  # noqa: E731
        r2.append({'model': k, 'n': len(lo), 'R2 stickers': t2.cv_r2(rn(STICKERS), y, 'linear'),
                   'R2 charge': t2.cv_r2(rn(CHARGE), y, 'linear'),
                   'R2 stickers + charge': t2.cv_r2(rn(STICKERS + CHARGE), y, 'linear'),
                   'R2 all (quadratic)': t2.cv_r2(rn(DESC), y, 'quadratic'),
                   'spearman chi vs nu below theta': stats.spearmanr(lo['chi'], lo['nu'])[0]})
    Cp, Cr = pd.DataFrame(rows), pd.DataFrame(r2)
    Cp.to_csv(os.path.join(OUT_DIR, 'C_below_theta.csv'), index=False)
    Cr.to_csv(os.path.join(OUT_DIR, 'C_below_theta_r2.csv'), index=False)
    print(f'\nC  below theta (own nu < {THETA}): chi residual after own nu + log N')
    print(Cr.round(3).to_string(index=False))
    w = Cp.pivot(index='descriptor', columns='model', values='partial_spearman_chi')
    print(w.reindex(w.abs().max(axis=1).sort_values(ascending=False).index).head(14).round(3).to_string())

    # D: energy-term split by own-nu bin
    D = []
    for k, d in models.items():
        terms = ['WF', 'DH', 'bond'] if k == 'Mpipi-GG' else ['AH', 'DH', 'dielectric', 'bond']
        g = d.assign(b=pd.cut(d['nu'], BINS)).groupby('b', observed=True)
        t = g[['chi'] + terms].median()
        t.insert(0, 'n', g.size())
        t.insert(0, 'model', k)
        D.append(t.reset_index())
    D = pd.concat(D)
    D.to_csv(os.path.join(OUT_DIR, 'D_terms.csv'), index=False)
    print('\nD  median % Rg per K by own-nu bin (CALVADOS DH includes its dielectric part)')
    print(D.round(4).to_string(index=False))
    bc.write_provenance(os.path.join(OUT_DIR, 'A_compaction.csv'), __file__,
                        {'nu_fit_r2_min': R2_MIN, 'theta_nu': THETA, 'stickers': STICKERS, 'charge': CHARGE},
                        [bc.BIO_CSV, os.path.join(L.OUT_ROOT, 'fdt_results_300K.csv'),
                         os.path.join(MP, 'mpipi_300K_fdt.csv'), os.path.join(MP, 'compact', 'mpipi_compact_300K_fdt.csv')])
    print(f'\nsaved {OUT_DIR}')


if __name__ == '__main__':
    main()
