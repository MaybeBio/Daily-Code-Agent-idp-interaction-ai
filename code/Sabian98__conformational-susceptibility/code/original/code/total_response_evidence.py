"""
total_response_evidence.py — the deck's evidence redone for the total heating response (2026-10-10).

Total response at 300 K, % of <Rg> per K, per model:
    total = chi (Boltzmann + the model's own electrostatic dH/dT; recharged_compare.load)
          + gamma_tea * dchi_tea (TEA hydrophobic term, CHARMM36m Arg; tea_fdt.py / tea_crossff.py, sum of the classes)
          + gamma_hyd * dchi_hyd (hydrophilic-solvation term; same scripts with --classes, column 'hyd')
This is the model the Wuttke calibration fitted (tea_calibrate.py fit/cv --free hyd, 278-320 K): TEA held at its prior
gamma_tea (CALVADOS-2 2.0 mol/kcal, TEA's value; Mpipi-Recharged 0.4) and gamma_hyd = 0.358 +/- 0.026 (CALVADOS-2),
0.209 +/- 0.018 (Mpipi-Recharged). Mpipi-GG (no Wuttke runs): each term at the mean relative size, SD(term)/SD(chi),
of the other two. (Before 2026-10-10 evening the proteome totals left the TEA term out; it is 5-14% of the spread of
the hydrophilic term and in the same direction.)
Parts: S sticker, Q charge (Boltzmann), dQ electrostatic dH/dT, B bond, TEA hydrophobic solvation, H hydrophilic
solvation (sum = total).

Sections, each mirroring deck slides (old slide numbers):
  A  by each model's own nu bin: chi, H, total (slide 12)
  B  cross-model agreement of the total, overall, by CALVADOS-2 nu bin (slide 13 table) and by kingdom
  C  expanded chains, CALVADOS-2 nu > 0.63: Boltzmann electrostatics, electrostatic dH/dT, H, total (slide 14)
  D  leftover after own nu (spline) + log N: share explained by compaction and length, cross-model r, variance share of
     each part (slides 15, 19, 21)
  E  composition of the leftover: Spearman and rank partial correlations, pi x (D+E+K) tercile grid, parts by decile
     (slides 16-18, 22)
  F  sequence predictability of the leftover and of the total: CV R^2 of descriptor groups, noise ceiling from
     total_response_noise.py (slides 20, 28)
  G  charge-held chains (Boltzmann definition, slide 24): their total response vs other chains at matched nu (slides 24-25)
  H  thirds by pi and by D+E+K fraction, median total by own nu bin (slides 17-18)
  I  factor split of the hydrophilic term: coupling r(Rg, E_hyd) x breathing CV(Rg) x fluctuation sd(E_hyd)/RT (slide 26)
  J  gamma sensitivity: B and A at gamma +/- 1 SE
Writes lineage_inputs/per_protein.csv (total for CALVADOS-2 as C, Mpipi-Recharged as M) for lineage_test.py (slide 32).
    python total_response_evidence.py [--noise DIR] [--out DIR]
Output: Traj Level Analysis/Analysis_FDT/total_response/*.csv, summary.txt
"""
import sys

sys.dont_write_bytecode = True

import argparse  # noqa: E402
import os  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

HERE = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, HERE)
import recharged_compare as RC  # noqa: E402

TLA = RC.TLA
TEA = os.path.join(TLA, 'Analysis_FDT', 'tea_300K')
OUT = os.path.join(TLA, 'Analysis_FDT', 'total_response')
MODELS = RC.MODELS
PAIRS = (('C', 'M'), ('C', 'R'), ('M', 'R'))
GAMMA = {'C': (0.358, 0.026), 'R': (0.209, 0.018)}
GAMMA_TEA = {'C': 2.0, 'R': 0.4}   # the calibration's prior (held fixed while gamma_hyd was fitted)
HYD_FILES = {'C': 'tea_dchi_classes_300K.csv', 'M': 'tea_dchi_gg_classes_300K.csv', 'R': 'tea_dchi_recharged_classes_300K.csv'}
# Forcing behind the totals: 'capped' (the main model from 2026-10-24: charged-group hydrophilic and TEA terms capped at
# 2 contacts per residue, capped_model.py, per-chain terms from capped_proteome.py) or 'A0' (the linear calibration of the
# first manuscript versions; its outputs are kept in total_response/A0_linear). Set with BENDER_FORCING.
FORCING = os.environ.get('BENDER_FORCING', 'capped')
if FORCING not in ('capped', 'A0'):
    raise ValueError(f'BENDER_FORCING must be capped or A0, not {FORCING!r}')
CAPPED = os.path.join(TLA, 'Analysis_FDT', 'capped_300K')
CAPPED_FILES = {'C': 'capped_dchi_calvados_300K.csv', 'M': 'capped_dchi_gg_300K.csv', 'R': 'capped_dchi_recharged_300K.csv'}
PARTS = ['S', 'Q', 'dQ', 'B', 'TEA', 'H']
TEA_CLASSES = ['aliph', 'pro', 'arom', 'argC36', 'rest']
BINS = RC.BINS
BINS13 = [0.30, 0.50, 0.54, 0.58, 0.63, 1.0]
R2_MIN = RC.R2_MIN
HYD_RES = 'DEKRNQSTHG'
GROUPS = {
    'stickers': ['f_F', 'f_W', 'f_Y', 'f_R', 'lambda_avg', 'f_aro', 'f_pi', 'f_pi_extended', 'log1p_aromatic_pairs',
                 'log1p_cation_pi_pairs'],
    'charge': ['FCR', 'NCPR', 'uversky_R', 'kappa', 'SCD', 'f_D', 'f_E', 'f_K'],
    'hydropathy': ['hydropathy_avg', 'uversky_H', 'SHD', 'f_hydrophobic', 'f_polar'],
    'hydrophilic residues': ['f_D', 'f_E', 'f_K', 'f_R', 'f_N', 'f_Q', 'f_S', 'f_T', 'f_H', 'f_G', 'f_hydres'],
}


def sp(a, b):
    return stats.spearmanr(a, b, nan_policy='omit')[0]


def rank_normal(x):
    return stats.norm.ppf((stats.rankdata(x) - 0.5) / len(x))


def partial(y, x, z):
    """Rank-based partial correlation of y and x with z (one or more columns) held fixed."""
    Z = np.atleast_2d(np.asarray(z, dtype=float).T).T
    ok = np.isfinite(np.asarray(y, dtype=float)) & np.isfinite(np.asarray(x, dtype=float)) & np.isfinite(Z).all(1)
    Z = np.column_stack([rank_normal(c) for c in Z[ok].T])
    ry = RC.residualise(rank_normal(np.asarray(y, dtype=float)[ok]), Z)
    rx = RC.residualise(rank_normal(np.asarray(x, dtype=float)[ok]), Z)
    return np.corrcoef(ry, rx)[0, 1]


def gamma_table():
    """(hydrophilic strength, SE) and TEA strength per calibrated model for the active forcing."""
    if FORCING == 'capped':
        import capped_model as CM
        return CM.params()
    return dict(GAMMA), dict(GAMMA_TEA)


def load(gscale=None):
    """Per model: recharged_compare parts + H; gscale = {'C': g, 'R': g} overrides the calibrated gamma."""
    D = RC.load()
    hyd = {}
    if FORCING == 'capped':                        # per-chain capped terms per unit strength (capped_proteome.py)
        for k, f in CAPPED_FILES.items():
            t = pd.read_csv(os.path.join(CAPPED, f)).set_index('UniProt_ID')
            hyd[k] = t[['dchi_rg_rel_chg_cap', 'dchi_rg_chg_cap', 'D_sd_chg_cap', 'dchi_rg_rel_tea_cap']].rename(
                columns={'dchi_rg_rel_chg_cap': 'dchi_rg_rel_hyd', 'dchi_rg_chg_cap': 'dchi_rg_hyd',
                         'D_sd_chg_cap': 'D_sd_hyd', 'dchi_rg_rel_tea_cap': 'tea'})
    else:
        for k, f in HYD_FILES.items():
            t = pd.read_csv(os.path.join(TEA, f)).set_index('UniProt_ID')
            t['tea'] = sum(t[f'dchi_rg_rel_{c}'] for c in TEA_CLASSES)
            hyd[k] = t[['dchi_rg_rel_hyd', 'dchi_rg_hyd', 'D_sd_hyd', 'tea']]
    gh, gtea = gamma_table()
    g = {k: (gscale or {}).get(k, v[0]) for k, v in gh.items()}
    for k in D:
        D[k] = D[k].join(hyd[k].rename(columns={'dchi_rg_rel_hyd': 'hyd'}), how='inner')
    # Mpipi-GG: SD(gamma hyd)/SD(chi) = mean of the calibrated models' ratios on their own full sets (tea_hyd_compare.py)
    f = np.mean([(g[k] * D[k].hyd).std() / D[k].chi.std() for k in ('C', 'R')])
    g['M'] = f * D['M'].chi.std() / D['M'].hyd.std()
    gt = dict(gtea)
    ft = np.mean([(gt[k] * D[k].tea).std() / D[k].chi.std() for k in ('C', 'R')])
    gt['M'] = ft * D['M'].chi.std() / D['M'].tea.std()
    for k in D:
        D[k]['H'] = g[k] * D[k].hyd
        D[k]['TEA'] = gt[k] * D[k].tea
        D[k]['total'] = D[k].chi + D[k].TEA + D[k].H
    g.update({f'tea_{k}': v for k, v in gt.items()})
    return D, g


def bio_table():
    bio = pd.read_csv(RC.BIO_CSV, low_memory=False).set_index('UniProt_ID')
    bio['f_DEK'] = bio.f_D + bio.f_E + bio.f_K
    bio['f_hydres'] = sum(bio[f'f_{a}'] for a in HYD_RES)
    bio['abs_NCPR'] = bio.NCPR.abs()
    return bio


def trio_ids(D):
    return D['C'].index.intersection(D['M'].index).intersection(D['R'].index)


# ---------------------------------------------------------------- A, B, C, J
def by_nu(D):
    rows = []
    for k, d in D.items():
        for lab, x in list(d.groupby(pd.cut(d.nu, BINS), observed=True)) + [('all', d)]:
            rows.append({'model': MODELS[k], 'own_nu_bin': str(lab), 'n': len(x), 'chi_median': x.chi.median(),
                         'TEA_median': x.TEA.median(), 'H_median': x.H.median(), 'total_median': x.total.median(),
                         'total_q25': x.total.quantile(0.25), 'total_q75': x.total.quantile(0.75),
                         'frac_contract': (x.total < 0).mean(), 'frac_contract_chi': (x.chi < 0).mean()})
    return pd.DataFrame(rows)


def agreement(D, bio):
    rows, trio = [], trio_ids(D)
    for a, b in PAIRS:
        for set_name, ids in (('pair', D[a].index.intersection(D[b].index)), ('three-model set', trio)):
            rows.append({'pair': f'{MODELS[a]} vs {MODELS[b]}', 'set': set_name, 'group': 'all', 'n': len(ids),
                         'rho_total': sp(D[a].total[ids], D[b].total[ids]), 'rho_chi': sp(D[a].chi[ids], D[b].chi[ids]),
                         'same_sign_total': (np.sign(D[a].total[ids]) == np.sign(D[b].total[ids])).mean()})
    # slide-13 table: three-model set by CALVADOS-2 nu, medians and rank agreement within bins
    tab = []
    nuC = D['C'].nu[trio]
    for lab, ids in pd.Series(trio, index=trio).groupby(pd.cut(nuC, BINS13), observed=True):
        ids = ids.index
        r = {'bin_nu_C': str(lab), 'n': len(ids)}
        for k in D:
            r[f'total_{MODELS[k]}'] = D[k].total[ids].median()
            r[f'chi_{MODELS[k]}'] = D[k].chi[ids].median()
        for a, b in PAIRS:
            r[f'rho_total {MODELS[a]}-{MODELS[b]}'] = sp(D[a].total[ids], D[b].total[ids])
        tab.append(r)
    # by kingdom: C-R on all paired chains, three-model set where > 20 chains
    king = bio.kingdom
    ids = D['C'].index.intersection(D['R'].index)
    for kg, x in pd.Series(ids, index=ids).groupby(king.reindex(ids)):
        i = x.index
        i3 = i.intersection(trio)
        r = {'pair': 'by kingdom', 'set': kg, 'group': kg, 'n': len(i), 'rho_total C-R': sp(D['C'].total[i], D['R'].total[i]),
             'rho_chi C-R': sp(D['C'].chi[i], D['R'].chi[i]), 'n_three_model': len(i3)}
        for a, b in (('C', 'M'), ('M', 'R')):
            r[f'rho_total {a}-{b}'] = sp(D[a].total[i3], D[b].total[i3]) if len(i3) > 20 else np.nan
        rows.append(r)
    return pd.DataFrame(rows), pd.DataFrame(tab)


def expanded(D):
    trio = trio_ids(D)
    ids = trio[D['C'].nu[trio] > 0.63]
    rows = []
    for k in D:
        d = D[k].loc[ids]
        rows.append({'model': MODELS[k], 'n': len(d), 'sticker_S': d.S.median(), 'charge_boltzmann_Q': d.Q.median(),
                     'electrostatic_dHdT': d.dQ.median(), 'electrostatic_net': (d.Q + d.dQ).median(),
                     'hydrophobic_TEA': d.TEA.median(), 'hydrophilic_H': d.H.median(), 'total': d.total.median(), 'chi': d.chi.median(),
                     'frac_contract_total': (d.total < 0).mean(), 'frac_contract_chi': (d.chi < 0).mean(),
                     'H_share_of_total_median': (d.H / d.total).median()})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------- D, E
def resid_frame(d, bio):
    X = np.column_stack([RC.spline(d.nu.to_numpy()), np.log(d.N.to_numpy())])
    out = pd.DataFrame({'L': RC.residualise(d.total.to_numpy(), X), 'L_chi': RC.residualise(d.chi.to_numpy(), X)}, index=d.index)
    for p in PARTS:
        out[p] = RC.residualise(d[p].to_numpy(), X)
    out.attrs['r2_compaction'] = 1 - out.L.var() / d.total.var()
    out.attrs['r2_compaction_chi'] = 1 - out.L_chi.var() / d.chi.var()
    return out


def leftovers(D, bio):
    full = {k: d[d.nu_r2 >= R2_MIN] for k, d in D.items()}
    trio = [u for u in trio_ids(D) if all(D[k].nu_r2.get(u, 0) >= R2_MIN for k in D)]
    res_full = {k: resid_frame(d, bio) for k, d in full.items()}
    res_trio = {k: resid_frame(D[k].loc[trio], bio) for k in D}
    z = {k: res_trio[k] / res_trio[k].L.std(ddof=1) for k in D}
    shared = sum(z.values()) / 3
    rows, corr = [], []

    def share(res, label, n, r2=np.nan, r2chi=np.nan):
        v = res.L.var(ddof=1)
        rows.append({'target': label, 'n': n, 'r2_compaction_length_total': r2, 'r2_compaction_length_chi': r2chi,
                     **{f'share_{p}': np.cov(res[p], res.L)[0, 1] / v for p in PARTS}})

    for k, r in res_full.items():
        share(r, f'{MODELS[k]}, own full set', len(r), r.attrs['r2_compaction'], r.attrs['r2_compaction_chi'])
    for k, r in res_trio.items():
        share(r, f'{MODELS[k]}, three-model set', len(r), r.attrs['r2_compaction'], r.attrs['r2_compaction_chi'])
    share(shared, 'shared leftover (three models)', len(trio))
    for a, b in PAIRS:
        corr.append({'pair': f'{MODELS[a]} vs {MODELS[b]}', 'set': 'three-model set', 'n': len(trio),
                     'r_leftover_total': np.corrcoef(res_trio[a].L, res_trio[b].L)[0, 1],
                     'r_leftover_chi': np.corrcoef(res_trio[a].L_chi, res_trio[b].L_chi)[0, 1]})
    pr = full['C'].index.intersection(full['R'].index)
    lc, lr = resid_frame(D['C'].loc[pr], bio), resid_frame(D['R'].loc[pr], bio)
    corr.append({'pair': 'CALVADOS-2 vs Mpipi-Recharged', 'set': 'all paired', 'n': len(pr),
                 'r_leftover_total': np.corrcoef(lc.L, lr.L)[0, 1], 'r_leftover_chi': np.corrcoef(lc.L_chi, lr.L_chi)[0, 1]})
    sharedCR = (lc / lc.L.std(ddof=1) + lr / lr.L.std(ddof=1)) / 2
    share(sharedCR, 'shared leftover (CALVADOS-2 + Mpipi-Recharged, all paired)', len(pr))
    return pd.DataFrame(rows), pd.DataFrame(corr), res_full, res_trio, shared, sharedCR


DESCRIPTORS = ['f_pi', 'f_DEK', 'FCR', 'abs_NCPR', 'NCPR', 'f_hydres', 'f_G', 'f_R', 'f_K', 'f_D', 'f_E', 'f_polar',
               'f_hydrophobic', 'hydropathy_avg', 'lambda_avg', 'kappa', 'SCD', 'length']


def composition(targets, bio):
    rows = []
    for name, y in targets.items():
        b = bio.reindex(y.index)
        r = {'target': name, 'n': len(y)}
        for c in DESCRIPTORS:
            r[f'rho {c}'] = sp(b[c], y)
        r['partial f_pi | f_DEK'] = partial(y, b.f_pi, b.f_DEK)
        r['partial f_DEK | f_pi'] = partial(y, b.f_DEK, b.f_pi)
        r['partial FCR | f_pi'] = partial(y, b.FCR, b.f_pi)
        r['partial f_hydres | f_pi'] = partial(y, b.f_hydres, b.f_pi)
        r['partial f_DEK | f_pi, f_hydres'] = partial(y, b.f_DEK, np.column_stack([b.f_pi, b.f_hydres]))
        r['partial f_R | FCR'] = partial(y, b.f_R, b.FCR)
        r['partial f_K | FCR'] = partial(y, b.f_K, b.FCR)
        rows.append(r)
    return pd.DataFrame(rows).set_index('target').T


def grid(y, bio):
    b = bio.reindex(y.index)
    tp = pd.cut(b.f_pi, [-1, 0.068, 0.156, 1], labels=['low', 'mid', 'high'])
    tc = pd.cut(b.f_DEK, [-1, 0.18, 0.34, 1], labels=['low', 'mid', 'high'])
    g = y.groupby([tp, tc], observed=False)
    return pd.DataFrame({'median': g.median(), 'n': g.size()}).reset_index().rename(columns={'f_pi': 'pi_third', 'f_DEK': 'DEK_third'})


def deciles(res, bio, label):
    rows = []
    b = bio.reindex(res.index)
    for c in ('f_DEK', 'f_pi', 'f_hydres'):
        q = pd.qcut(b[c].rank(method='first'), 10, labels=False)
        for dec, x in res.groupby(q):
            rows.append({'target': label, 'by': c, 'decile': int(dec) + 1, 'n': len(x),
                         **{f'{p}_mean': x[p].mean() for p in ['L'] + PARTS}, 'L_se': x.L.std() / np.sqrt(len(x))})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------- F
def predictability(targets, bio, noise):
    import test2_slope_vs_sequence as t2
    rows = []
    for name, (y, noise_var) in targets.items():
        b = bio.reindex(y.index)
        seq, layers = t2.sequence_layers(b)
        seq['f_hydres'] = b.f_hydres
        rn = lambda cols: np.column_stack([rank_normal(seq[c].to_numpy()) for c in cols])  # noqa: E731
        yy = y.to_numpy()
        r = {'target': name, 'n': len(yy), 'noise_ceiling': 1 - noise_var / yy.var(ddof=1) if noise_var is not None else np.nan}
        cols = []
        for lname, lc in layers:
            cols += lc
            r[f'R2 {lname} (cum.)'] = t2.cv_r2(rn(cols), yy, 'linear')
        for g, gc in GROUPS.items():
            r[f'R2 {g} only'] = t2.cv_r2(rn(gc), yy, 'linear')
        r['R2 quadratic all'] = t2.cv_r2(rn(cols), yy, 'quadratic')
        r['R2 best'] = max(r[f'R2 {layers[-1][0]} (cum.)'], r['R2 quadratic all'])
        r['best share of ceiling'] = r['R2 best'] / r['noise_ceiling']
        rows.append(r)
        print(f'  F {name}: ceiling {r["noise_ceiling"]:.3f}, best {r["R2 best"]:.3f}', flush=True)
    return pd.DataFrame(rows)


def noise_model(noise_dir, g=None):
    """Mean squared SE (relative units) of total and chi per model from total_response_noise.py / capped_noise.py
    subsets; with g (load()'s strengths), the noise files must have been computed at those strengths."""
    out = {}
    for k in MODELS:
        p = os.path.join(noise_dir, f'total_response_se_{k}.csv') if noise_dir else ''
        if os.path.exists(p):
            t = pd.read_csv(p)
            if g is not None and 'gamma_hyd' in t and abs(t.gamma_hyd.median() - g[k]) > 1e-3 * abs(g[k]):
                raise ValueError(f'{p}: noise computed at gamma_hyd {t.gamma_hyd.median():.4f}, the totals use {g[k]:.4f} '
                                 f'({FORCING} forcing); pass the matching --noise folder')
            out[k] = {'total': float((t.total_rel_se ** 2).mean()), 'chi': float((t.chi_rel_se ** 2).mean()), 'n': len(t),
                      'median_rel_se_total': float((t.total_rel_se / t.total_rel.abs()).median())}
    return out


# ---------------------------------------------------------------- G, H, I
def charge_held(D):
    rows = []
    for k, d in D.items():
        h = (d.Q > 0) & (d.Q > d.S)
        for lab, x in d.groupby(pd.cut(d.nu, BINS), observed=True):
            hh = h[x.index]
            rows.append({'model': MODELS[k], 'own_nu_bin': str(lab), 'n_held': int(hh.sum()), 'n_other': int((~hh).sum()),
                         'total_held': x.total[hh].median(), 'total_other': x.total[~hh].median(),
                         'H_held': x.H[hh].median(), 'H_other': x.H[~hh].median(),
                         'chi_held': x.chi[hh].median(), 'chi_other': x.chi[~hh].median()})
        z = (d.total - d.groupby(pd.cut(d.nu, BINS), observed=True).total.transform('median'))
        rows.append({'model': MODELS[k], 'own_nu_bin': 'all (total minus nu-bin median)', 'n_held': int(h.sum()),
                     'n_other': int((~h).sum()), 'total_held': z[h].median(), 'total_other': z[~h].median(),
                     'mannwhitney_p': stats.mannwhitneyu(z[h].dropna(), z[~h].dropna()).pvalue if h.sum() > 2 else np.nan})   # nu <= 0.30 chains are outside the bins (z NaN)
    return pd.DataFrame(rows)


def thirds(D, bio):
    rows = []
    for k, d in D.items():
        b = bio.reindex(d.index)
        for col, cuts in (('f_pi', [-1, 0.068, 0.156, 1]), ('f_DEK', [-1, 0.18, 0.34, 1])):
            t = pd.cut(b[col], cuts, labels=['low', 'mid', 'high'])
            for (lab, th), x in d.groupby([pd.cut(d.nu, BINS), t], observed=True):
                rows.append({'model': MODELS[k], 'by': col, 'own_nu_bin': str(lab), 'third': th, 'n': len(x),
                             'total_median': x.total.median() if len(x) >= 10 else np.nan,
                             'chi_median': x.chi.median() if len(x) >= 10 else np.nan})
    return pd.DataFrame(rows)


def factors(D):
    """H / <Rg> = -gamma Cov(Rg, D)/(RT <Rg>) = -r(Rg, D) x CV(Rg) x gamma sd(D)/RT (exact)."""
    rows = []
    fdt = {'C': pd.read_csv(os.path.join(TLA, 'Analysis_FDT', 'fdt_results_300K.csv')),
           'M': pd.read_csv(os.path.join(TLA, 'Analysis_FDT_MPIPI', 'mpipi_300K_fdt.csv')),
           'R': pd.read_csv(os.path.join(TLA, 'Analysis_FDT_RECHARGED', 'recharged_fdt_300K.csv'))}
    RT = {'C': 8.314462618e-3 * 300, 'R': 8.314462618e-3 * 300, 'M': 1.98720425864083e-3 * 300}
    for k, d in D.items():
        f = fdt[k]
        f = f[f.metric == 'rg'].set_index('UniProt_ID').reindex(d.index)
        cov = -d.dchi_rg_hyd * RT[k]                       # Cov(Rg, D) per unit gamma
        x = pd.DataFrame({'coupling': -cov / (f.A_sd * d.D_sd_hyd), 'breathing_CV': f.A_sd / f.A_mean,
                          'fluct_per_rootN': d.D_sd_hyd / RT[k] / np.sqrt(d.N), 'nu': d.nu, 'H': d.H})
        for lab, y in x.groupby(pd.cut(x.nu, BINS), observed=True):
            rows.append({'model': MODELS[k], 'own_nu_bin': str(lab), 'n': len(y), 'H_median': y.H.median(),
                         'coupling_median': y.coupling.median(), 'breathing_CV_median': y.breathing_CV.median(),
                         'fluct_per_rootN_median': y.fluct_per_rootN.median()})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------- main
def main(noise_dir, out):
    os.makedirs(out, exist_ok=True)
    pd.set_option('display.width', 250)
    pd.set_option('display.max_columns', 40)
    D, g = load()
    bio = bio_table()
    T = {}
    T['A_by_nu'] = by_nu(D)
    T['B_agreement'], T['B_table_slide13'] = agreement(D, bio)
    T['C_expanded'] = expanded(D)
    T['D_leftover_shares'], T['D_leftover_corr'], res_full, res_trio, shared, sharedCR = leftovers(D, bio)
    tg = {'shared (three models)': shared.L, 'shared (C + R, all paired)': sharedCR.L,
          **{f'{MODELS[k]} leftover (own full set)': r.L for k, r in res_full.items()},
          **{f'{MODELS[k]} total (own full set)': D[k].total for k in D}}
    T['E_composition'] = composition(tg, bio).reset_index().rename(columns={'index': 'statistic'})
    T['E_grid_shared'] = grid(shared.L, bio)
    T['E_grid_sharedCR'] = grid(sharedCR.L, bio)
    T['E_deciles'] = pd.concat([deciles(shared, bio, 'shared (three models)'), deciles(sharedCR, bio, 'shared (C + R)')])
    nm = noise_model(noise_dir, g)
    if nm:
        nz = {}
        for k, r in res_trio.items():
            nz[k] = nm[k]['total'] / r.L.std(ddof=1) ** 2 if k in nm else np.nan
        targets = {'shared leftover (three models)': (shared.L, np.nansum(list(nz.values())) / 9 if len(nz) == 3 else None)}
        for k in ('C', 'R'):
            r = res_full[k]
            targets[f'{MODELS[k]} leftover (own full set)'] = (r.L, nm[k]['total'] if k in nm else None)
            targets[f'{MODELS[k]} total (own full set)'] = (D[k].total.loc[r.index], nm[k]['total'] if k in nm else None)
        T['F_predictability'] = predictability(targets, bio, nm)
        T['F_noise'] = pd.DataFrame(nm).T.reset_index().rename(columns={'index': 'model'})
    T['G_charge_held'] = charge_held(D)
    T['H_thirds'] = thirds(D, bio)
    T['I_factors_hyd'] = factors(D)
    sens = []
    gh = gamma_table()[0]
    for lab, sc in (('gamma - 1 SE', {k: v[0] - v[1] for k, v in gh.items()}), ('gamma + 1 SE', {k: v[0] + v[1] for k, v in gh.items()})):
        Ds, gs = load(sc)
        a, tab = agreement(Ds, bio)
        a = a[a.group == 'all'].assign(variant=lab, **{f'gamma_{k}': v for k, v in gs.items()})
        sens.append(a)
        sens.append(by_nu(Ds).query("own_nu_bin == 'all'").assign(variant=lab))
    T['J_gamma_sensitivity'] = pd.concat(sens)
    # lineage_test.py inputs: CALVADOS-2 total as C, Mpipi-Recharged total as M (own nu, length)
    li = os.path.join(out, 'lineage_inputs')
    os.makedirs(li, exist_ok=True)
    per = pd.DataFrame({'chi_C': D['C'].total, 'nu_C': D['C'].nu, 'length_C': D['C'].N}).join(
        pd.DataFrame({'chi_M': D['R'].total, 'nu_M': D['R'].nu, 'length_M': D['R'].N}), how='outer')
    per.index.name = 'UniProt_ID'
    per.to_csv(os.path.join(li, 'per_protein.csv'))
    lines = ['gamma: ' + ', '.join(f'{k} {v:.4f}' for k, v in g.items()),
             *[f'{MODELS[k]}: {len(d)} chains' for k, d in D.items()], f'three-model set: {len(trio_ids(D))}']
    for name, t in T.items():
        t.to_csv(os.path.join(out, f'{name}.csv'), index=False)
        lines += ['', f'== {name}', t.to_string(index=False, float_format=lambda v: f'{v:.3f}')]
    open(os.path.join(out, 'summary.txt'), 'w').write('\n'.join(lines) + '\n')
    print('\n'.join(lines))


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--noise', default=None, help='directory with total_response_se_<model>.csv')
    ap.add_argument('--out', default=OUT)
    a = ap.parse_args()
    main(a.noise, a.out)
