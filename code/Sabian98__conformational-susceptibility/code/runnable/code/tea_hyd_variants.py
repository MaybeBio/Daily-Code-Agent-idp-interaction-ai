"""
tea_hyd_variants.py — do the proteome-wide results for the calibrated hydrophilic term depend on which residues carry it?

For the chains all three models simulate, the hydrophilic term split by residue group (tea_fdt.py / tea_crossff.py
--classes: hyd_DE, hyd_K, hyd_R, hyd_pol (NQSTH), hyd_G) is recombined for each residue-set variant of
tea_calibrate.py sens, with each calibrated model at its own gamma_hyd for that variant (Wuttke fit, 278-320 K:
Wuttke2014/calibration/{calvados,recharged}/sensitivity.csv), the TEA term at the prior the calibration held it at
(gamma_tea 2.0 / 0.4; added 2026-10-10, left out before), and Mpipi-GG at the mean relative size of the other two.
Reported per variant: share contracting on heating, median slope, rank of the total vs chi, median by nu bin, and the
cross-model rank agreement of the calibrated total response.
    python tea_hyd_variants.py GROUPS_DIR OUT_DIR
GROUPS_DIR holds tea_dchi{,_gg,_recharged}_classes_300K.csv with the group columns.
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

CAL = (__import__('os').environ['CS_DATA'] + '/BENDER-MT/Traj Level Analysis/Wuttke2014/calibration')
VARIANTS = {'all (DEKRNQSTH+G)': ['DE', 'K', 'R', 'pol', 'G'], 'without G (backbone)': ['DE', 'K', 'R', 'pol'],
            'without R': ['DE', 'K', 'pol', 'G'], 'charged only (DEKR)': ['DE', 'K', 'R'],
            'polar + G only': ['pol', 'G'], 'acidic only (DE)': ['DE']}
FILES = {'C': 'tea_dchi_classes_300K.csv', 'M': 'tea_dchi_gg_classes_300K.csv', 'R': 'tea_dchi_recharged_classes_300K.csv'}
GAMMA_TEA = {'C': 2.0, 'R': 0.4}
TEA_CLASSES = ['aliph', 'pro', 'arom', 'argC36', 'rest']
BINS = [0.30, 0.45, 0.50, 0.54, 0.58, 0.63, 1.0]


def main(gdir, out):
    os.makedirs(out, exist_ok=True)
    base = X.load('charmm36')
    G = {k: pd.read_csv(os.path.join(gdir, f)).set_index('UniProt_ID') for k, f in FILES.items()}
    ids = G['C'].index.intersection(G['M'].index).intersection(G['R'].index)
    chi = {k: base[k].chi.reindex(ids) for k in G}
    tea = {k: sum(G[k].loc[ids, f'dchi_rg_rel_{c}'] for c in TEA_CLASSES) for k in G}
    # Mpipi-GG strengths at the others' relative size, scaled on each model's full set as total_response_evidence.load
    # (fixed 2026-10-10: scaling on the 1,106-chain three-model set gave gamma_tea 0.364 / gamma_hyd 0.1051 here and 0.300 /
    # 0.1059 elsewhere, hence two values, 0.87 and 0.86, for the same agreement)
    own = {k: G[k].index.intersection(base[k].chi.dropna().index) for k in G}
    chiF = {k: base[k].chi.reindex(own[k]) for k in G}
    teaF = {k: sum(G[k].loc[own[k], f'dchi_rg_rel_{c}'] for c in TEA_CLASSES) for k in G}
    fT = np.mean([(GAMMA_TEA[k] * teaF[k]).std() / chiF[k].std() for k in GAMMA_TEA])
    g_tea = dict(GAMMA_TEA, M=fT * chiF['M'].std() / teaF['M'].std())
    sens = {k: pd.read_csv(os.path.join(CAL, m, 'sensitivity.csv')) for k, m in (('C', 'calvados'), ('R', 'recharged'))}
    nu = pd.read_csv(X.BIO_CSV, usecols=['UniProt_ID', 'nu'], low_memory=False).set_index('UniProt_ID').nu.reindex(ids)
    rows, bins = [], []
    for var, groups in VARIANTS.items():
        hyd = {k: sum(G[k].loc[ids, f'dchi_rg_rel_hyd_{g}'] for g in groups) for k in G}
        gam = {}
        for k in ('C', 'R'):
            s = sens[k]
            s = s[(s['hyd residues'] == var) & (s.interp_deg == 3) & (s.tmax.astype(str).isin(['320.0', '320']))]
            gam[k] = float(s.gamma_hyd.iloc[0])
        hydF = {k: sum(G[k].loc[own[k], f'dchi_rg_rel_hyd_{g}'] for g in groups) for k in G}
        f = np.mean([(gam[k] * hydF[k]).std() / chiF[k].std() for k in ('C', 'R')])
        gam['M'] = f * chiF['M'].std() / hydF['M'].std()
        tot = {k: chi[k] + g_tea[k] * tea[k] + gam[k] * hyd[k] for k in G}
        r = {'variant': var, 'n': len(ids), **{f'gamma_{X.NAMES[k]}': gam[k] for k in G}}
        for k in G:
            r[f'contract_{X.NAMES[k]}'] = (tot[k] < 0).mean()
            r[f'median_{X.NAMES[k]}'] = tot[k].median()
            r[f'rho_vs_chi_{X.NAMES[k]}'] = spearmanr(tot[k], chi[k])[0]
        for a, b in X.PAIRS:
            r[f'agree {X.NAMES[a]}-{X.NAMES[b]}'] = spearmanr(tot[a], tot[b])[0]
        rows.append(r)
        for k in G:
            for lab, grp in tot[k].groupby(pd.cut(nu, BINS)):
                bins.append({'variant': var, 'model': X.NAMES[k], 'nu_bin': str(lab), 'n': len(grp), 'median': grp.median()})
    t = pd.DataFrame(rows)
    b = pd.DataFrame(bins).pivot_table(index=['variant', 'nu_bin'], columns='model', values='median').reset_index()
    t.to_csv(os.path.join(out, 'hyd_variants.csv'), index=False)
    b.to_csv(os.path.join(out, 'hyd_variants_by_nu.csv'), index=False)
    txt = (f'{len(ids)} chains in all three models; each variant with its own Wuttke-calibrated gamma_hyd (CALVADOS-2, '
           'Mpipi-Recharged; 278-320 K) and Mpipi-GG at their mean relative size\n\n' + t.round(3).T.to_string() +
           '\n\nmedian total slope (% per K) by nu bin\n' + b.round(3).to_string(index=False))
    open(os.path.join(out, 'hyd_variants.txt'), 'w').write(txt + '\n')
    print(txt)


if __name__ == '__main__':
    main(*sys.argv[1:3])
