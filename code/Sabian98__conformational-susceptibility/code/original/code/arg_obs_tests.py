"""
arg_obs_tests.py — four tests on the all-observables solvation pass (tea_fdt.py / tea_crossff.py --classes
--all-observables → Analysis_FDT/total_response/all_observables/; hydrophilic term split by residue group).

Calibrated total = chi + gamma_tea TEA + gamma_hyd hyd (as total_response_evidence.py). The Arg question: the
'without R' variant drops Arg from the hydrophilic forcing and uses its own Wuttke-calibrated strength
(Wuttke2014/calibration/<model>/sensitivity.csv, interpolation 3, 278-320 K; Mpipi-GG from
tea_300K/consolidation/hyd_variants.csv, at the other models' relative size).
  1  Arg ranking: per chain, total with Arg minus total without Arg (Rg); which chains differ most, what predicts it;
     a shortlist of tractable chains (<= 150 residues, CALVADOS-2 nu fit R^2 >= 0.99) where the two variants differ most
  2  Does the growth-temperature result depend on Arg? Total response and projected swelling at the OGT vs OGT (organism
     means, bacteria with phylum fixed effects; all kingdoms with kingdom fixed effects), with and without Arg
  3  Total response for all eight observables: median by nu bin, share of chains whose sign differs from chi, and
     CALVADOS-2 vs Mpipi-GG rank agreement of chi and of the total
  4  Residue-age result (deck slide 37) vs Arg: partial correlation of the age score with the leftover of the total,
     holding pi and D+E+K fixed (as lineage_test.py), then also Arg; and with an age score that leaves Arg out
Models: CALVADOS-2 and Mpipi-GG (complete); Mpipi-Recharged included when its pass is >= 99% complete.
    python arg_obs_tests.py [--out DIR]
Output: Analysis_FDT/total_response/arg_obs/ (A_arg_*.csv, B_ogt_arg.csv, C_observables*.csv, D_age_arg.csv, summary.txt)
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
import ogt_response as OR  # noqa: E402
import total_response_evidence as E  # noqa: E402

TLA = E.TLA
AO = f'{TLA}/Analysis_FDT/total_response/all_observables'
OUT = f'{TLA}/Analysis_FDT/total_response/arg_obs'
FILES = {'C': 'tea_dchi_classes_300K.csv', 'M': 'tea_dchi_gg_classes_300K.csv', 'R': 'tea_dchi_recharged_classes_300K.csv'}
FDT = {'C': f'{TLA}/Analysis_FDT/fdt_results_300K.csv', 'M': f'{TLA}/Analysis_FDT_MPIPI/mpipi_300K_fdt.csv',
       'R': f'{TLA}/Analysis_FDT_RECHARGED/recharged_fdt_300K.csv'}
N_EXPECTED = {'C': 10282, 'M': 1303, 'R': 10266}
TEA_CLASSES = ['aliph', 'pro', 'arom', 'argC36', 'rest']
OBS_C = {'rg': 'rg', 'ree': 'ree', 'delta': 'delta', 'fragmentation_index': 'graph_fragmentation_index',
         'global_efficiency': 'graph_global_efficiency', 'avg_clustering': 'graph_avg_clustering',
         'transitivity': 'graph_transitivity', 'degree_assortativity': 'graph_degree_assortativity'}
AGE = {**dict.fromkeys('GA', 1), **dict.fromkeys('VD', 2), 'P': 3, 'S': 4, **dict.fromkeys('EL', 5), 'T': 6, 'R': 7,
       'N': 8, 'K': 9, 'Q': 10, 'I': 11, 'C': 12, 'H': 13, 'F': 14, 'M': 15, 'Y': 16, 'W': 17}
BINS = E.BINS


def gammas():
    g = {}
    for k, m in (('C', 'calvados'), ('R', 'recharged')):
        s = pd.read_csv(f'{TLA}/Wuttke2014/calibration/{m}/sensitivity.csv')
        s = s[(s.interp_deg == 3) & (s.tmax.astype(str).isin(['320.0', '320']))].set_index('hyd residues').gamma_hyd
        g[k] = {'with Arg': float(s['all (DEKRNQSTH+G)']), 'without Arg': float(s['without R'])}
    v = pd.read_csv(f'{TLA}/Analysis_FDT/tea_300K/consolidation/hyd_variants.csv').set_index('variant')['gamma_Mpipi-GG']
    g['M'] = {'with Arg': float(v['all (DEKRNQSTH+G)']), 'without Arg': float(v['without R'])}
    return g


GAMMA_TEA = {'C': 2.0, 'R': 0.4, 'M': 0.3000}


def load(models):
    """Per model: one row per chain and observable: chi (relative), TEA, hyd, hyd_R (per unit strength)."""
    out = {}
    for k in models:
        t = pd.read_csv(f'{AO}/{FILES[k]}').set_index('UniProt_ID')
        f = pd.read_csv(FDT[k])
        if k == 'M':
            ok = pd.read_csv(f'{TLA}/Analysis_FDT_MPIPI/mpipi_300K.csv', usecols=['UniProt_ID', 'status'])
            f = f[f.UniProt_ID.isin(ok.loc[ok.status == 'ok', 'UniProt_ID'])]
        if k == 'R':
            f = f[f.max_bond_nm <= 1.0]
        rows = []
        for metric, col in OBS_C.items():
            c = col if k == 'C' else metric
            if f'dchi_{c}_rel_hyd' not in t:
                continue
            ff = f[f.metric == metric].set_index('UniProt_ID')
            d = pd.DataFrame({'chi': 100 * ff.chi_total / ff.A_mean, 'chi_abs': ff.chi_total, 'A_mean': ff.A_mean})
            d = d.join(pd.DataFrame({'tea': sum(t[f'dchi_{c}_rel_{x}'] for x in TEA_CLASSES), 'hyd': t[f'dchi_{c}_rel_hyd'],
                                     'hyd_R': t[f'dchi_{c}_rel_hyd_R'],
                                     'tea_abs': sum(t[f'dchi_{c}_{x}'] for x in TEA_CLASSES), 'hyd_abs': t[f'dchi_{c}_hyd'],
                                     'hydR_abs': t[f'dchi_{c}_hyd_R']}), how='inner')
            rows.append(d.assign(metric=metric))
        out[k] = pd.concat(rows).reset_index()
    return out


def totals(d, k, g):
    for var, gam in g[k].items():
        hyd = d.hyd if var == 'with Arg' else d.hyd - d.hyd_R
        hyd_abs = d.hyd_abs if var == 'with Arg' else d.hyd_abs - d.hydR_abs
        d[f'total {var}'] = d.chi + GAMMA_TEA[k] * d.tea + gam * hyd
        d[f'total_abs {var}'] = d.chi_abs + GAMMA_TEA[k] * d.tea_abs + gam * hyd_abs
    d['arg_effect'] = d['total with Arg'] - d['total without Arg']
    return d


def main(out):
    os.makedirs(out, exist_ok=True)
    pd.set_option('display.width', 220)
    pd.set_option('display.max_columns', 30)
    models = [k for k in ('C', 'M', 'R') if pd.read_csv(f'{AO}/{FILES[k]}', usecols=['UniProt_ID']).shape[0] >= 0.99 * N_EXPECTED[k]]
    g = gammas()
    L = load(models)
    for k in models:
        L[k] = totals(L[k], k, g)
    bio = E.bio_table()
    bio['age'] = bio.sequence.map(lambda q: np.mean([AGE[c] for c in q if c in AGE]))
    bio['age_noR'] = bio.sequence.map(lambda q: np.mean([AGE[c] for c in q if c in AGE and c != 'R']))
    nu = {'C': bio.nu.rename('nu'), 'M': pd.read_csv(f'{TLA}/Analysis_FDT_MPIPI/mpipi_300K.csv').set_index('UniProt_ID').nu,
          'R': pd.read_csv(f'{TLA}/Analysis_FDT_RECHARGED/recharged_300K_nu.csv').set_index('UniProt_ID').nu}
    lines = [f'models: {[E.MODELS[k] for k in models]}; strengths {g}; TEA {GAMMA_TEA}']
    # ---------------- 1 Arg ranking
    rk, cors = [], []
    for k in models:
        r = L[k][L[k].metric == 'rg'].set_index('UniProt_ID')
        r = r.join(bio[['length', 'f_R', 'f_K', 'f_DEK', 'f_pi', 'FCR', 'nu_fit_r2']]).join(nu[k].rename('nu_own'))
        for c in ('f_R', 'f_K', 'FCR', 'f_pi', 'nu_own', 'length'):
            cors.append({'model': E.MODELS[k], 'descriptor': c, 'spearman_with_arg_effect': stats.spearmanr(r[c], r.arg_effect, nan_policy='omit')[0]})
        r['model'] = E.MODELS[k]
        rk.append(r)
        lines += ['', f'== 1 Arg effect on the total response of Rg, {E.MODELS[k]} (% per K): with minus without Arg',
                  r.arg_effect.describe(percentiles=[0.05, 0.5, 0.95]).round(4).to_string(),
                  'by own nu bin, median arg effect: ' + r.groupby(pd.cut(r.nu_own, BINS), observed=True).arg_effect.median().round(3).to_string().replace('\n', '; ')]
    RK = pd.concat(rk)
    RK.to_csv(os.path.join(out, 'A_arg_per_chain.csv'))
    pd.DataFrame(cors).to_csv(os.path.join(out, 'A_arg_correlates.csv'), index=False)
    c = RK[RK.model == 'CALVADOS-2']
    short = c[(c.length <= 150) & (c.nu_fit_r2 >= 0.99)]
    if 'M' in models:
        mm = RK[RK.model == 'Mpipi-GG'].arg_effect
        short = short.join(mm.rename('arg_effect_GG'))
    short = short.sort_values('arg_effect').head(25)
    short.to_csv(os.path.join(out, 'A_arg_shortlist.csv'))
    lines += ['', 'Arg effect correlates:', pd.DataFrame(cors).round(3).to_string(index=False), '',
              'shortlist (CALVADOS-2, <= 150 residues, clean nu fit; most negative = Arg adds most contraction):',
              short[['length', 'f_R', 'f_K', 'nu_own', 'total with Arg', 'total without Arg', 'arg_effect'] +
                    (['arg_effect_GG'] if 'M' in models else [])].round(3).to_string()]
    # ---------------- 2 OGT with and without Arg
    bio['OGT'] = OR.ogt_table(bio)
    bio['phylum'] = OR.ogt_table(bio, 'phylum')
    rng = np.random.default_rng(5)
    rows = []
    for k in models:
        r = L[k][L[k].metric == 'rg'].set_index('UniProt_ID')
        rg = pd.read_csv(FDT[k]).query("metric == 'rg'").set_index('UniProt_ID').A_mean
        d = r.join(bio[['organism', 'kingdom', 'phylum', 'OGT', 'length']]).join(rg.rename('rg')).dropna(subset=['OGT'])
        dT = d.OGT + 273.15 - 300.0
        a300 = np.log(d.rg / np.sqrt(d.length))
        for var in ('with Arg', 'without Arg'):
            d[f'alphaOGT {var}'] = a300 + np.log1p(d[f'total {var}'] / 100 * dT)
        cols = ['total with Arg', 'total without Arg', 'alphaOGT with Arg', 'alphaOGT without Arg']
        for label, sub, fe in (('bacteria, phylum FE', d[d.kingdom == 'Bacteria'].dropna(subset=['phylum']), 'phylum'),
                               ('all kingdoms, kingdom FE', d, 'kingdom')):
            o = sub.groupby(['organism', fe]).agg(OGT=('OGT', 'first'), **{x: (x, 'mean') for x in cols}).reset_index()
            o = o[o[fe].isin(o[fe].value_counts()[lambda s: s >= 3].index)]
            for x in cols:
                def fit(dd):
                    X = np.column_stack([np.ones(len(dd)), dd.OGT, pd.get_dummies(dd[fe], drop_first=True).astype(float)])
                    return 10 * np.linalg.lstsq(X, dd[x], rcond=None)[0][1]
                b = [fit(o.iloc[rng.integers(0, len(o), len(o))]) for _ in range(2000)]
                rows.append({'model': E.MODELS[k], 'set': label, 'organisms': len(o), 'quantity': x, 'slope_per_10C': fit(o),
                             'ci_lo': np.percentile(b, 2.5), 'ci_hi': np.percentile(b, 97.5)})
    B = pd.DataFrame(rows)
    B.to_csv(os.path.join(out, 'B_ogt_arg.csv'), index=False)
    lines += ['', '== 2 growth temperature with and without Arg (per 10 C)', B.round(4).to_string(index=False)]
    # ---------------- 3 all observables
    rows, agree = [], []
    for k in models:
        d = L[k].merge(nu[k].rename('nu_own'), left_on='UniProt_ID', right_index=True, how='left')
        for metric, x in d.groupby('metric'):
            rel = metric != 'degree_assortativity'
            tot, chi = (x['total with Arg'], x.chi) if rel else (x['total_abs with Arg'], x.chi_abs)
            row = {'model': E.MODELS[k], 'observable': metric, 'units': '% per K' if rel else 'per K (absolute)',
                   'n': len(x), 'chi_median': chi.median(), 'total_median': tot.median(),
                   'frac_total_negative': (tot < 0).mean(), 'frac_chi_negative': (chi < 0).mean(),
                   'frac_sign_flip': (np.sign(tot) != np.sign(chi)).mean(), 'spearman_total_vs_chi': stats.spearmanr(tot, chi, nan_policy='omit')[0]}
            for lab, grp in x.assign(t=tot).groupby(pd.cut(x.nu_own, BINS), observed=True):
                row[f'total {lab}'] = grp.t.median()
            rows.append(row)
    C = pd.DataFrame(rows)
    if 'M' in models:
        for metric in OBS_C:
            a = L['C'][L['C'].metric == metric].set_index('UniProt_ID')
            b = L['M'][L['M'].metric == metric].set_index('UniProt_ID')
            ids = a.index.intersection(b.index)
            col = 'total with Arg' if metric != 'degree_assortativity' else 'total_abs with Arg'
            chic = 'chi' if metric != 'degree_assortativity' else 'chi_abs'
            agree.append({'observable': metric, 'n': len(ids), 'rho_chi C2-GG': stats.spearmanr(a.loc[ids, chic], b.loc[ids, chic], nan_policy='omit')[0],
                          'rho_total C2-GG': stats.spearmanr(a.loc[ids, col], b.loc[ids, col], nan_policy='omit')[0]})
    A3 = pd.DataFrame(agree)
    C.to_csv(os.path.join(out, 'C_observables.csv'), index=False)
    A3.to_csv(os.path.join(out, 'C_observables_agreement.csv'), index=False)
    lines += ['', '== 3 total response for every observable (with Arg)', C.round(4).to_string(index=False), '',
              A3.round(3).to_string(index=False)]
    # ---------------- 4 residue age vs Arg
    def rn(v):
        return stats.norm.ppf((stats.rankdata(v) - 0.5) / len(v))

    def partial(y, x, Z):
        Zr = np.column_stack([np.ones(len(y))] + [rn(z) for z in Z])
        ry = rn(y) - Zr @ np.linalg.lstsq(Zr, rn(y), rcond=None)[0]
        rx = rn(x) - Zr @ np.linalg.lstsq(Zr, rn(x), rcond=None)[0]
        return np.corrcoef(ry, rx)[0, 1]
    rows = []
    for k, f in (('C', f'{TLA}/Analysis_FDT/total_response/lineage/per_protein.csv'),
                 ('R', f'{TLA}/Analysis_FDT/total_response/lineage/recharged/per_protein.csv')):
        p = pd.read_csv(f).set_index('UniProt_ID').join(bio[['f_R', 'age_noR']]).dropna(subset=['L', 'age', 'pi', 'dek', 'f_R'])
        rows.append({'model': E.MODELS[k], 'n': len(p), 'raw age vs L': stats.spearmanr(p.age, p.L)[0],
                     'age | pi, D+E+K (rank-normal partial; lineage_test/deck: -0.10 / -0.14)': partial(p.L, p.age, [p.pi, p.dek]),
                     'age | pi, D+E+K, Arg': partial(p.L, p.age, [p.pi, p.dek, p.f_R]),
                     'age without Arg | pi, D+E+K': partial(p.L, p.age_noR, [p.pi, p.dek]),
                     'Arg | pi, D+E+K': partial(p.L, p.f_R, [p.pi, p.dek])})
    D4 = pd.DataFrame(rows)
    D4.to_csv(os.path.join(out, 'D_age_arg.csv'), index=False)
    lines += ['', '== 4 residue age vs leftover of the total (negative = more contraction), rank partial correlations',
              D4.round(3).to_string(index=False)]
    open(os.path.join(out, 'summary.txt'), 'w').write('\n'.join(lines) + '\n')
    print('\n'.join(lines))


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', default=OUT)
    main(ap.parse_args().out)
