"""
What holds a chain together, by energy rather than by composition.

Per protein and model, chi(Rg) (% of <Rg> per K, averaged over the five temperatures as in Test 2b) splits exactly
into energy-term parts:
  S    sticker part           Cov(Rg, S)/RT^2         CALVADOS-2 chi_AH, Mpipi-GG chi_WF
  Q    charge part            Cov(Rg, Q)/RT^2         CALVADOS-2 chi_DH - chi_dielectric, Mpipi-GG chi_DH
  dQ   charge dH/dT part      -Cov(Rg, dQ/dT)/RT      CALVADOS-2 chi_dielectric, Mpipi-GG 0
  B    bond part              Cov(Rg, B)/RT^2
  rest Mpipi-GG chi_total minus S + Q + B (xtc rounding of the recomputed terms), so the parts add up exactly
Q > 0: charges pull the chain in as it compacts (cohesive); Q < 0: charges push it apart (opposing).
Dynamic charge share phi = Q / (|S| + |Q|).

Static definition, 300 K, from the per-frame npz: mean sticker and charge energy per residue in kT;
cohesion c_S = -<S>/N, c_Q = -<Q>/N (c_Q > 0: charges net attractive); static charge share psi = c_Q / (|c_S| + |c_Q|).

A  classes: per nu bin, both models, full sets: how many chains have cohesive / opposing / dominant charge parts,
   and how many have net attractive charge energy
B  leftover decomposition (paired sets of Test 2b, compaction nu | rg | both): each part residualised on the same
   design as Test 2b (OLS is linear, so the part leftovers add up to the Test 2b leftover exactly); variance shares;
   each descriptor's correlation with the leftover split into the parts' contributions; CV R^2 of the sticker and
   charge groups on each part. Also per model on its full set.
C  matched-compaction comparison of mechanism classes (phi, psi): leftover and its parts by class; coupling per
   unit cohesion (S / c_S vs Q / c_Q) for chains with net attractive charges

Output: Traj Level Analysis/Analysis_FDT/energy_mechanism_split/
Usage: python3 energy_mechanism_split.py
"""
import sys

sys.dont_write_bytecode = True  # keep the NAS code folder free of __pycache__

import os  # noqa: E402
from concurrent.futures import ThreadPoolExecutor  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

HERE = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, HERE)
import test2_slope_vs_sequence as t2  # noqa: E402
import test2b_cross_model_residual as t2b  # noqa: E402

t1, L, bc = t2.t1, t2.L, t2.bc
MPIPI_DIR = t2b.MPIPI_DIR
OUT_DIR = os.path.join(L.OUT_ROOT, 'energy_mechanism_split')
CAL_NPZ = os.path.join(L.OUT_ROOT, 'npz', '300K')
MP_NPZ = os.path.join(MPIPI_DIR, 'npz_300K')
TEMPS = t2b.TEMPS
R_KJ, R_KCAL = 8.314462618e-3, 1.98720425864083e-3
BINS = [0.30, 0.45, 0.50, 0.54, 0.58, 0.63, 1.0]
PARTS = ['S', 'Q', 'dQ', 'B', 'rest']
PHI_CUT = 0.1                   # |phi| <= 0.1: charges neutral to the response
KEY_DESC = ['f_aro', 'f_pi', 'lambda_avg', 'f_R', 'f_K', 'f_E', 'f_D', 'FCR', 'NCPR', 'SCD', 'kappa', 'SHD',
            'hydropathy_avg']


def parts_cal() -> pd.DataFrame:
    cols = ['UniProt_ID', 'metric', 'chi_total', 'chi_AH', 'chi_DH', 'chi_dielectric', 'chi_bond', 'A_mean']
    d = pd.concat([pd.read_csv(os.path.join(L.OUT_ROOT, f'fdt_results_{T}K.csv'), usecols=cols).query("metric == 'rg'")
                   .assign(T=T) for T in TEMPS])
    d = d.assign(S=d.chi_AH, Q=d.chi_DH - d.chi_dielectric, dQ=d.chi_dielectric, B=d.chi_bond)
    d['rest'] = d.chi_total - d[['S', 'Q', 'dQ', 'B']].sum(axis=1)
    return d


def parts_mpipi() -> pd.DataFrame:
    cols = ['UniProt_ID', 'metric', 'n_residues', 'chi_total', 'chi_WF', 'chi_DH', 'chi_bond', 'A_mean']
    out, ok = [], None
    for T in TEMPS:
        s = pd.read_csv(os.path.join(MPIPI_DIR, f'mpipi_{T}K.csv'), usecols=['UniProt_ID', 'status'])
        ids = set(s.loc[s.status == 'ok', 'UniProt_ID'])
        ok = ids if ok is None else ok & ids
        out.append(pd.read_csv(os.path.join(MPIPI_DIR, f'mpipi_{T}K_fdt.csv'), usecols=cols).query("metric == 'rg'").assign(T=T))
    d = pd.concat(out)
    d = d[d.UniProt_ID.isin(ok)]
    d = d.assign(S=d.chi_WF, Q=d.chi_DH, dQ=0.0, B=d.chi_bond)
    d['rest'] = d.chi_total - d[['S', 'Q', 'B']].sum(axis=1)
    return d


def average(d: pd.DataFrame) -> pd.DataFrame:
    """Parts as % of <Rg> per K, averaged over the temperatures (same convention as Test 2b chi_rel)."""
    g = d.groupby('UniProt_ID')
    A = g['A_mean'].mean()
    out = pd.DataFrame({p: 100 * g[p].mean() / A for p in PARTS + ['chi_total']}).rename(columns={'chi_total': 'chi'})
    out['n_T'] = g.size()
    out['rg300'] = d[d['T'] == 300].set_index('UniProt_ID')['A_mean']
    out = out[out.n_T == len(TEMPS)]
    out['phi'] = out.Q / (out.S.abs() + out.Q.abs())
    return out


def static_energy(ids, model: str) -> pd.DataFrame:
    """Mean sticker and charge energy per residue at 300 K in kT, from the per-frame npz."""
    folder, kS, kQ, RT = ((CAL_NPZ, 'V_AH', 'V_DH', R_KJ * 300) if model == 'C' else (MP_NPZ, 'E_WF', 'E_DH', R_KCAL * 300))

    def one(u):
        p = os.path.join(folder, f'{u}.npz')
        if not os.path.exists(p):
            return u, np.nan, np.nan
        z = np.load(p)
        n = len(str(z['sequence']))
        return u, -float(z[kS].mean()) / RT / n, -float(z[kQ].mean()) / RT / n

    with ThreadPoolExecutor(16) as ex:
        rows = list(ex.map(one, ids))
    s = pd.DataFrame(rows, columns=['UniProt_ID', 'c_S', 'c_Q']).set_index('UniProt_ID')
    s['psi'] = s.c_Q / (s.c_S.abs() + s.c_Q.abs())
    return s


def classes(df: pd.DataFrame) -> pd.DataFrame:
    df = df.assign(bin=pd.cut(df.nu, BINS))
    g = df.groupby('bin', observed=True)
    out = pd.DataFrame({
        'n': g.size(),
        'chi_median': g.chi.median(),
        'S_median': g.S.median(), 'Q_median': g.Q.median(), 'dQ_median': g.dQ.median(), 'B_median': g.B.median(),
        'frac_Q_cohesive': g.phi.apply(lambda x: (x > PHI_CUT).mean()),
        'frac_Q_neutral': g.phi.apply(lambda x: (x.abs() <= PHI_CUT).mean()),
        'frac_Q_opposing': g.phi.apply(lambda x: (x < -PHI_CUT).mean()),
        'frac_Q_dominant': g.apply(lambda x: (x.Q.abs() > x.S.abs()).mean(), include_groups=False),
        'phi_median': g.phi.median(),
        'c_S_median': g.c_S.median(), 'c_Q_median': g.c_Q.median(),
        'frac_charge_attractive': g.c_Q.apply(lambda x: (x > 0).mean()),
        'frac_psi_ge_0.2': g.psi.apply(lambda x: (x >= 0.2).mean()),
        'psi_max': g.psi.max(),
    })
    out.index = out.index.astype(str)
    out.loc['all'] = [len(df), df.chi.median(), df.S.median(), df.Q.median(), df.dQ.median(), df.B.median(),
                      (df.phi > PHI_CUT).mean(), (df.phi.abs() <= PHI_CUT).mean(), (df.phi < -PHI_CUT).mean(),
                      (df.Q.abs() > df.S.abs()).mean(), df.phi.median(), df.c_S.median(), df.c_Q.median(),
                      (df.c_Q > 0).mean(), (df.psi >= 0.2).mean(), df.psi.max()]
    return out


def leftover_parts(j: pd.DataFrame, design: np.ndarray, tag: str) -> pd.DataFrame:
    """Residualise chi and every part on the same design; the part leftovers add up to the chi leftover."""
    out = pd.DataFrame(index=j.index)
    out[f'res_{tag}'] = t2b.residualise(j['chi'].to_numpy(), design)
    for p in PARTS:
        out[f'res_{tag}_{p}'] = t2b.residualise(j[p].to_numpy(), design)
    gap = np.abs(out[f'res_{tag}'] - out[[f'res_{tag}_{p}' for p in PARTS]].sum(axis=1)).max()
    assert gap < 1e-8 * max(1.0, out[f'res_{tag}'].abs().max()), f'part leftovers do not add up ({gap})'
    return out


def decompose(res: pd.Series, parts: pd.DataFrame, seq: pd.DataFrame, label: str, rn) -> tuple:
    """Variance shares, descriptor-correlation split and group CV R^2 for one leftover and its parts."""
    y = res.to_numpy()
    sd = y.std(ddof=1)
    var_rows = {'target': label, 'n': len(y)}
    for p in parts.columns:
        var_rows[f'share_{p}'] = np.cov(parts[p], y)[0, 1] / y.var(ddof=1)
    desc_rows = []
    for c in KEY_DESC:
        x = t1.rank_normal(seq[c].to_numpy())
        row = {'target': label, 'descriptor': c, 'r_total': np.corrcoef(x, y)[0, 1]}
        for p in parts.columns:
            row[f'r_from_{p}'] = np.cov(x, parts[p])[0, 1] / (x.std(ddof=1) * sd)
        desc_rows.append(row)
    r2_rows = []
    for g in ('stickers', 'charge'):
        X = rn(t2b.GROUPS[g])
        row = {'target': label, 'group': g, 'R2_total': t2.cv_r2(X, y, 'linear')}
        for p in parts.columns:
            yp = parts[p].to_numpy()
            row[f'R2_{p}'] = t2.cv_r2(X, yp, 'linear') if yp.std() > 1e-9 * sd else np.nan
        r2_rows.append(row)
    return var_rows, desc_rows, r2_rows


def class_compare(df: pd.DataFrame, tag: str) -> pd.DataFrame:
    """Leftover and its parts by mechanism class (compaction already removed by the residualisation)."""
    df = df.copy()
    df['class_phi'] = np.select([df.phi > PHI_CUT, df.phi < -PHI_CUT], ['cohesive', 'opposing'], 'neutral')
    df['class_psi'] = pd.qcut(df.psi, 4, labels=['psi Q1 (most repulsive)', 'psi Q2', 'psi Q3', 'psi Q4 (most attractive)'])
    rows = []
    for cl in ('class_phi', 'class_psi'):
        for k, x in df.groupby(cl, observed=True):
            row = {'model': tag, 'scheme': cl, 'class': k, 'n': len(x), 'nu_median': x.nu.median(),
                   'leftover_median': x[f'res_{tag}'].median(), 'leftover_mean': x[f'res_{tag}'].mean(),
                   'leftover_se': x[f'res_{tag}'].std(ddof=1) / np.sqrt(len(x))}
            for p in PARTS:
                row[f'res_{p}_mean'] = x[f'res_{tag}_{p}'].mean()
            for c in ('f_aro', 'FCR', 'NCPR', 'SCD', 'c_S', 'c_Q', 'phi', 'psi'):
                row[f'{c}_median'] = x[c].median()
            rows.append(row)
    return pd.DataFrame(rows)


def efficiency(df: pd.DataFrame, tag: str) -> pd.DataFrame:
    """Response per unit cohesion: S / c_S vs Q / c_Q, chains with net attractive charges, by nu bin."""
    x = df[(df.c_Q > 0.005) & (df.c_S > 0)].assign(bin=pd.cut(df.nu, BINS))
    x = x.assign(eff_S=x.S / x.c_S, eff_Q=x.Q / x.c_Q)
    g = x.groupby('bin', observed=True)
    out = pd.DataFrame({'n': g.size(), 'eff_S_median': g.eff_S.median(), 'eff_Q_median': g.eff_Q.median()})
    out['ratio_Q_over_S'] = out.eff_Q_median / out.eff_S_median
    out.index = out.index.astype(str)
    out.loc['all'] = [len(x), x.eff_S.median(), x.eff_Q.median(), x.eff_Q.median() / x.eff_S.median()]
    out.insert(0, 'model', tag)
    return out


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    pd.set_option('display.width', 250)
    pd.set_option('display.max_columns', 40)
    f3 = lambda v: f'{v:.3f}'  # noqa: E731
    bio = pd.read_csv(bc.BIO_CSV, low_memory=False).set_index('UniProt_ID')

    cal = average(parts_cal())
    mp_raw = parts_mpipi()
    mp_len = mp_raw.groupby('UniProt_ID')['n_residues'].first()
    mp = average(mp_raw)
    nu_m = pd.read_csv(os.path.join(MPIPI_DIR, 'mpipi_300K.csv'), usecols=['UniProt_ID', 'nu', 'nu_r2']).set_index('UniProt_ID')
    cal = cal.join(bio[['nu', 'nu_fit_r2', 'length']].rename(columns={'nu_fit_r2': 'nu_r2'}))
    mp = mp.join(nu_m)
    mp['length'] = mp_len.reindex(mp.index)   # from the runs: 35 Mpipi-GG proteins are not in BENDER_BIO
    for name, d in (('CALVADOS-2', cal), ('Mpipi-GG', mp)):
        print(f'{name}: {len(d)} proteins; max |rest| / |chi| = {(d.rest.abs() / d.chi.abs().clip(lower=1e-9)).median():.2e} (median)')

    print('reading 300 K npz for static energies ...', flush=True)
    cal = cal.join(static_energy(cal.index, 'C'))
    mp = mp.join(static_energy(mp.index, 'M'))
    print(f'  npz missing: CALVADOS-2 {cal.c_S.isna().sum()}, Mpipi-GG {mp.c_S.isna().sum()}')

    # A: classes by nu bin
    A = []
    for tag, d in (('CALVADOS-2', cal), ('Mpipi-GG', mp)):
        c = classes(d.dropna(subset=['nu', 'c_S']))
        print(f'\n== A {tag}: mechanism classes by nu bin (full set)')
        print(c.to_string(float_format=f3))
        A.append(c.assign(model=tag))
    A = pd.concat(A)
    A.index.name = 'nu_bin'
    pathA = os.path.join(OUT_DIR, 'A_classes_by_nu.csv')
    A.to_csv(pathA)

    # B: leftover decomposition
    seq_all, _ = t2.sequence_layers(bio)
    V, D, R2, C, E = [], [], [], [], []
    pairs, _ = t2b.load_pairs()
    for measure in ('nu', 'rg', 'both'):
        for set_name, sel in (('fit_r2', (pairs.nu_r2_C >= t2b.R2_MIN) & (pairs.nu_r2_M >= t2b.R2_MIN)),
                              ('all_paired', pairs.index == pairs.index)):
            j = pairs[sel]
            seq = seq_all.reindex(j.index)
            rn = lambda cols: np.column_stack([t1.rank_normal(seq[c].to_numpy()) for c in cols])  # noqa: E731
            jc = cal.reindex(j.index)[PARTS + ['chi']]
            jm = mp.reindex(j.index)[PARTS + ['chi']]
            if jc.isna().any().any() or jm.isna().any().any():
                raise ValueError('paired proteins missing parts')
            lc = leftover_parts(jc, t2b.compaction_design(j, 'C', measure), 'C')
            lm = leftover_parts(jm, t2b.compaction_design(j, 'M', measure), 'M')
            sdC, sdM = lc.res_C.std(), lm.res_M.std()
            shared = (lc.res_C / sdC + lm.res_M / sdM) / 2
            sh_parts = pd.DataFrame({p: (lc[f'res_C_{p}'] / sdC + lm[f'res_M_{p}'] / sdM) / 2 for p in PARTS})
            assert np.allclose(shared, sh_parts.sum(axis=1))
            for label, res, parts in (('shared', shared, sh_parts),
                                      ('CALVADOS-2', lc.res_C, lc[[f'res_C_{p}' for p in PARTS]].set_axis(PARTS, axis=1)),
                                      ('Mpipi-GG', lm.res_M, lm[[f'res_M_{p}' for p in PARTS]].set_axis(PARTS, axis=1))):
                v, dd, r2 = decompose(res, parts, seq, label, rn)
                tags = {'set': set_name, 'compaction': measure}
                V.append({**tags, **v})
                D += [{**tags, **r} for r in dd]
                R2 += [{**tags, **r} for r in r2]
            if set_name == 'fit_r2' and measure == 'nu':
                # check against the saved Test 2b leftovers
                t2b_res = pd.read_csv(os.path.join(t2b.OUT_DIR, 'residuals_fit_r2.csv'), index_col=0).reindex(j.index)
                print(f"\nTest 2b leftover reproduced: max |diff| C {np.abs(t2b_res.res_C - lc.res_C).max():.1e}, "
                      f"M {np.abs(t2b_res.res_M - lm.res_M).max():.1e}, shared {np.abs(t2b_res.shared - shared).max():.1e}")

    # B on each model's full set (nu measure), and C, E there
    for tag, d in (('CALVADOS-2', cal), ('Mpipi-GG', mp)):
        d = d.dropna(subset=['nu', 'c_S'])
        d = d[d.index.isin(seq_all.dropna().index)]   # descriptor regressions need sequence features
        print(f'{tag}: full-set decomposition on {len(d)} proteins with sequence features')
        d = d.assign(**{f'rg300_{tag}': d.rg300, f'nu_{tag}': d.nu})
        seq = seq_all.reindex(d.index)
        rn = lambda cols: np.column_stack([t1.rank_normal(seq[c].to_numpy()) for c in cols])  # noqa: E731
        lf = leftover_parts(d[PARTS + ['chi']], t2b.compaction_design(d, tag, 'nu'), tag)
        v, dd, r2 = decompose(lf[f'res_{tag}'], lf[[f'res_{tag}_{p}' for p in PARTS]].set_axis(PARTS, axis=1), seq, tag, rn)
        tags = {'set': 'full', 'compaction': 'nu'}
        V.append({**tags, **v})
        D += [{**tags, **r} for r in dd]
        R2 += [{**tags, **r} for r in r2]
        dd_ = d.join(lf).join(seq[['f_aro', 'FCR', 'NCPR', 'SCD']])
        C.append(class_compare(dd_, tag))
        E.append(efficiency(d, tag))

    V, D, R2, C, E = pd.DataFrame(V), pd.DataFrame(D), pd.DataFrame(R2), pd.concat(C), pd.concat(E)
    V.to_csv(os.path.join(OUT_DIR, 'B_variance_shares.csv'), index=False)
    D.to_csv(os.path.join(OUT_DIR, 'B_descriptor_split.csv'), index=False)
    R2.to_csv(os.path.join(OUT_DIR, 'B_group_r2.csv'), index=False)
    C.to_csv(os.path.join(OUT_DIR, 'C_class_compare.csv'), index=False)
    E.index.name = 'nu_bin'
    E.to_csv(os.path.join(OUT_DIR, 'C_efficiency.csv'))
    per = pd.concat([cal.add_suffix('_C'), mp.add_suffix('_M')], axis=1)
    per.to_csv(os.path.join(OUT_DIR, 'per_protein.csv'))

    print('\n== B variance shares of the leftover')
    print(V.to_string(index=False, float_format=f3))
    print('\n== B descriptor correlations with the leftover, split by part (fit_r2, nu; and full sets)')
    sel = D[((D.set == 'fit_r2') & (D.compaction == 'nu')) | (D.set == 'full')]
    print(sel.to_string(index=False, float_format=f3))
    print('\n== B group CV R^2 on each part (fit_r2, nu; full)')
    print(R2[((R2.set == 'fit_r2') & (R2.compaction == 'nu')) | (R2.set == 'full')].to_string(index=False, float_format=f3))
    print('\n== C leftover by mechanism class (full sets, nu)')
    print(C.to_string(index=False, float_format=f3))
    print('\n== C response per unit cohesion, chains with net attractive charges')
    print(E.to_string(float_format=f3))

    bc.write_provenance(pathA, __file__,
                        {'parts': 'S sticker, Q charge Boltzmann, dQ charge dH/dT, B bond, rest rounding; % <Rg>/K, '
                                  'mean over 278-320 K', 'phi': 'Q/(|S|+|Q|)', 'phi_cut': PHI_CUT,
                         'static': '300 K, -<E>/N in kT', 'residualisation': 'as Test 2b (own compaction spline + log N)',
                         'bins': BINS},
                        [bc.BIO_CSV] + [os.path.join(L.OUT_ROOT, f'fdt_results_{T}K.csv') for T in TEMPS]
                        + [os.path.join(MPIPI_DIR, f'mpipi_{T}K{s}.csv') for T in TEMPS for s in ('', '_fdt')])
    print(f'saved {OUT_DIR}')


if __name__ == '__main__':
    main()
