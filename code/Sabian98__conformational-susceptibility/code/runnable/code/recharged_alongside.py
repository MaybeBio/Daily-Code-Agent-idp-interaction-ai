"""
recharged_alongside.py — Mpipi-Recharged next to CALVADOS-2 and Mpipi-GG for the deck's susceptibility results, all at
300 K from stored data (Mpipi-Recharged has 300 K runs only; the deck's CALVADOS-2 / Mpipi-GG numbers average 278-320 K,
so the 300 K values of all three are given here for a like-for-like comparison).

  A  susceptibility chi(Rg) by each model's own nu bin, and its parts (slides 12-13)
  B  leftover of chi after own nu (spline) + log N (slide 15): composition correlations and rank partials (pi | D+E+K,
     D+E+K | pi), the pi x (D+E+K) tercile grid of the leftover (slide 16), thirds by nu (slides 17-18), and the
     sticker and charge parts of the leftover by decile of D+E+K and pi (slide 22)
  C  charge-held chains (slides 24-25): medians of the energy parts; the electrostatic dH/dT part relative to the charge
     part (CALVADOS-2 dielectric, Mpipi-Recharged screening)
  D  factor split of the sticker part (slide 26): chi_s / <Rg> = r(Rg, E_s) x CV(Rg) x sd(E_s)/RT / T, per own nu bin,
     from the per-frame energies in the npz
  E  sequence predictability of chi at 300 K (slide 28, Test 2-like): cross-validated R^2 of descriptor groups against a
     noise ceiling from the block-jackknife SE of chi
    python recharged_alongside.py [--out DIR] [--workers N]
Output: Analysis_FDT_RECHARGED/alongside/*.csv, summary.txt
"""
import sys

sys.dont_write_bytecode = True

import argparse  # noqa: E402
import os  # noqa: E402
from concurrent.futures import ThreadPoolExecutor  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

HERE = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, HERE)
import recharged_compare as RC  # noqa: E402
import total_response_evidence as E  # noqa: E402

TLA = RC.TLA
OUT = os.path.join(TLA, 'Analysis_FDT_RECHARGED', 'alongside')
NPZ = {'C': (f'{TLA}/Analysis_FDT/npz/300K', 'V_AH'), 'M': (f'{TLA}/Analysis_FDT_MPIPI/npz_300K', 'E_WF'),
       'R': (f'{TLA}/Analysis_FDT_RECHARGED/npz_300K', 'V_AH')}
RT_UNIT = {'C': 8.314462618e-3 * 300, 'R': 8.314462618e-3 * 300, 'M': 1.98720425864083e-3 * 300}
SE = {'C': (f'{TLA}/Analysis_FDT/fdt_results_300K.csv', None), 'M': (f'{TLA}/Analysis_FDT_MPIPI/mpipi_300K_fdt.csv', None),
      'R': (f'{TLA}/Analysis_FDT_RECHARGED/recharged_fdt_300K.csv', None)}
GROUPS = {'stickers': ['f_F', 'f_W', 'f_Y', 'f_R', 'lambda_avg', 'f_aro', 'f_pi', 'f_pi_extended', 'log1p_aromatic_pairs',
                       'log1p_cation_pi_pairs'],
          'charge': ['FCR', 'NCPR', 'uversky_R', 'kappa', 'SCD', 'f_D', 'f_E', 'f_K'],
          'hydropathy': ['hydropathy_avg', 'uversky_H', 'SHD', 'f_hydrophobic', 'f_polar']}


def rn(v):
    return stats.norm.ppf((stats.rankdata(v) - 0.5) / len(v))


def partial(y, x, z):
    Z = np.column_stack([np.ones(len(y)), rn(z)])
    ry, rx = rn(y) - Z @ np.linalg.lstsq(Z, rn(y), rcond=None)[0], rn(x) - Z @ np.linalg.lstsq(Z, rn(x), rcond=None)[0]
    return np.corrcoef(ry, rx)[0, 1]


def factors_one(args):
    k, uid = args
    path, key = NPZ[k]
    p = os.path.join(path, f'{uid}.npz')
    if k == 'M' and not os.path.exists(p):
        p = os.path.join(TLA, 'Analysis_FDT_MPIPI', 'compact', 'npz_300K', f'{uid}.npz')
    try:
        z = np.load(p)
        rg, es = z['rg'].astype(float), z[key].astype(float)
        ok = np.isfinite(rg) & np.isfinite(es)
        rg, es = rg[ok], es[ok]
        return {'model': E.MODELS[k], 'UniProt_ID': uid, 'coupling': -np.corrcoef(rg, es)[0, 1], 'breathing_CV': rg.std() / rg.mean(),
                'fluct_per_rootN': es.std() / RT_UNIT[k] / np.sqrt(len(str(z['sequence']))),
                'chi_s_rel': np.cov(rg, es)[0, 1] / (RT_UNIT[k] * 300.0) / rg.mean() * 100}
    except Exception as ex:  # noqa: BLE001
        return {'model': E.MODELS[k], 'UniProt_ID': uid, 'error': repr(ex)}


def main(out, workers):
    os.makedirs(out, exist_ok=True)
    pd.set_option('display.width', 220)
    pd.set_option('display.max_columns', 30)
    D = RC.load()
    bio = E.bio_table()
    lines = [f'{E.MODELS[k]}: {len(d)} chains at 300 K' for k, d in D.items()]
    # A
    A = RC.by_nu(D)
    A.to_csv(f'{out}/A_chi_by_nu.csv', index=False)
    lines += ['', '== A chi(Rg) by own nu bin, 300 K (% per K)', A[['model', 'own_nu_bin', 'n', 'chi_median', 'S_median', 'Q_median', 'dQ_median']].round(4).to_string(index=False)]
    # B
    comp, grids, thirds, dec = [], [], [], []
    for k, d in D.items():
        d = d[d.nu_r2 >= RC.R2_MIN]
        d = d.join(bio[[c for c in ('f_pi', 'f_DEK', 'FCR', 'hydropathy_avg') if c not in d]], how='inner').dropna(subset=['f_pi', 'f_DEK'])
        X = np.column_stack([RC.spline(d.nu.to_numpy()), np.log(d.N.to_numpy())])
        L = pd.Series(RC.residualise(d.chi.to_numpy(), X), index=d.index)
        Lz = L / L.std()
        parts = {p: pd.Series(RC.residualise(d[p].to_numpy(), X), index=d.index) / L.std() for p in ('S', 'Q', 'dQ', 'B')}
        comp.append({'model': E.MODELS[k], 'n': len(d), 'r2_nu_logN': 1 - L.var() / d.chi.var(),
                     'rho pi': stats.spearmanr(d.f_pi, L)[0], 'rho D+E+K': stats.spearmanr(d.f_DEK, L)[0],
                     'rho FCR': stats.spearmanr(d.FCR, L)[0], 'partial pi | D+E+K': partial(L.to_numpy(), d.f_pi.to_numpy(), d.f_DEK.to_numpy()),
                     'partial D+E+K | pi': partial(L.to_numpy(), d.f_DEK.to_numpy(), d.f_pi.to_numpy())})
        tp = pd.cut(d.f_pi, [-1, 0.068, 0.156, 1], labels=['low', 'mid', 'high'])
        tc = pd.cut(d.f_DEK, [-1, 0.18, 0.34, 1], labels=['low', 'mid', 'high'])
        g = Lz.groupby([tp, tc], observed=False).median().rename('median').reset_index()
        grids.append(g.assign(model=E.MODELS[k]))
        for col, t in (('pi', tp), ('D+E+K', tc)):
            for (b, th), x in d.groupby([pd.cut(d.nu, RC.BINS), t], observed=True):
                thirds.append({'model': E.MODELS[k], 'by': col, 'own_nu_bin': str(b), 'third': th, 'n': len(x),
                               'chi_median': x.chi.median() if len(x) >= 10 else np.nan})
        for col in ('f_DEK', 'f_pi'):
            q = pd.qcut(d[col].rank(method='first'), 10, labels=False)
            for de, idx in Lz.groupby(q).groups.items():
                dec.append({'model': E.MODELS[k], 'by': col, 'decile': int(de) + 1, 'L': Lz[idx].mean(),
                            **{f'{p}_part': parts[p][idx].mean() for p in parts}})
    Bc, Bg, Bt, Bd = pd.DataFrame(comp), pd.concat(grids), pd.DataFrame(thirds), pd.DataFrame(dec)
    for t, n in ((Bc, 'B_composition'), (Bg, 'B_grid'), (Bt, 'B_thirds'), (Bd, 'B_deciles')):
        t.to_csv(f'{out}/{n}.csv', index=False)
    tt = Bt.pivot_table(index=['by', 'own_nu_bin', 'third'], columns='model', values='chi_median').reset_index()
    ratio = []
    for (by, b), x in tt.groupby(['by', 'own_nu_bin']):
        x = x.set_index('third')
        if {'low', 'high'} <= set(x.index):
            hi, lo = ('high', 'low') if by == 'pi' else ('low', 'high')
            ratio.append({'by': by, 'own_nu_bin': b, **{m: x.at[hi, m] / x.at[lo, m] for m in E.MODELS.values() if m in x}})
    lines += ['', '== B leftover composition (300 K)', Bc.round(3).to_string(index=False), '',
              'grid (median leftover, SD): ' + Bg.pivot_table(index=['f_pi', 'f_DEK'], columns='model', values='median').round(2).to_string(),
              '', 'thirds: susceptible third / other third (high pi / low pi; low D+E+K / high D+E+K)', pd.DataFrame(ratio).round(2).to_string(index=False),
              '', 'parts by decile (lowest and highest tenth):',
              Bd[Bd.decile.isin([1, 10])].round(3).to_string(index=False)]
    # C
    rows = []
    for k, d in D.items():
        h = (d.Q > 0) & (d.Q > d.S)
        x = d[h]
        rows.append({'model': E.MODELS[k], 'n_held': int(h.sum()), 'frac': h.mean(), 'S_median': x.S.median(), 'Q_median': x.Q.median(),
                     'dQ_median': x.dQ.median(), 'chi_median': x.chi.median(),
                     'dQ/Q median (all chains, |Q| > 0.01)': (d.dQ / d.Q)[d.Q.abs() > 0.01].median(),
                     'chi held / chi others (same nu bin), median': (x.chi / d.groupby(pd.cut(d.nu, RC.BINS), observed=True).chi.transform('median')[h]).median()})
    Cc = pd.DataFrame(rows)
    Cc.to_csv(f'{out}/C_charge_held.csv', index=False)
    lines += ['', '== C charge-held chains (300 K)', Cc.round(3).to_string(index=False)]
    # D
    tasks = [(k, u) for k, d in D.items() for u in d.index]
    with ThreadPoolExecutor(workers) as ex:
        F = pd.DataFrame(list(ex.map(factors_one, tasks)))
    F = F.dropna(subset=['coupling'])
    nus = pd.concat([d.nu.rename('nu').to_frame().assign(model=E.MODELS[k]) for k, d in D.items()]).reset_index()
    F = F.merge(nus, on=['UniProt_ID', 'model'])
    Fd = F.groupby(['model', pd.cut(F.nu, RC.BINS)], observed=True)[['chi_s_rel', 'coupling', 'breathing_CV', 'fluct_per_rootN']].median().reset_index()
    Fd['nu'] = Fd.nu.astype(str)
    F.to_csv(f'{out}/D_factors_per_chain.csv', index=False)
    Fd.to_csv(f'{out}/D_factors_by_nu.csv', index=False)
    lines += ['', '== D sticker susceptibility factors by own nu bin (300 K)', Fd.round(3).to_string(index=False)]
    # E
    import test2_slope_vs_sequence as t2
    rows = []
    for k, d in D.items():
        se = pd.read_csv(SE[k][0])
        se = se[se.metric == 'rg'].set_index('UniProt_ID')
        d = d.join((100 * se.chi_total_se / se.A_mean).rename('se'), how='inner').dropna(subset=['se'])
        d = d.loc[d.index.intersection(bio.dropna(subset=['sequence']).index)]
        b = bio.reindex(d.index)
        seq, layers = t2.sequence_layers(b)
        y = d.chi.to_numpy()
        ceil = 1 - (d.se ** 2).mean() / y.var(ddof=1)
        r = {'model': E.MODELS[k], 'n': len(y), 'noise_ceiling': ceil}
        Xall = np.column_stack([rn(seq[c].to_numpy()) for _, cols in layers for c in cols])
        for g, gc in GROUPS.items():
            r[f'R2 {g} only'] = t2.cv_r2(np.column_stack([rn(seq[c].to_numpy()) for c in gc]), y, 'linear')
        r['R2 all (quadratic)'] = t2.cv_r2(Xall, y, 'quadratic')
        r['share of ceiling'] = r['R2 all (quadratic)'] / ceil
        rows.append(r)
        print(f'E {E.MODELS[k]} done', flush=True)
    Ee = pd.DataFrame(rows)
    Ee.to_csv(f'{out}/E_sequence.csv', index=False)
    lines += ['', '== E sequence predictability of chi at 300 K', Ee.round(3).to_string(index=False)]
    open(f'{out}/summary.txt', 'w').write('\n'.join(lines) + '\n')
    print('\n'.join(lines))


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', default=OUT)
    ap.add_argument('--workers', type=int, default=8)
    a = ap.parse_args()
    main(a.out, a.workers)
