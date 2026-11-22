"""
Test 2b: what sequence information do CALVADOS-2 and Mpipi-GG share in their susceptibility beyond compaction?

Per protein and model, chi(Rg) as % of <Rg> per K, averaged over the five temperatures (FDT chi_total).
Each model's chi is residualised on its own compaction at 300 K and log N; the residual is how much more or less
susceptible the protein is than its compaction and length predict in that model. Compaction measure (--compaction):
  nu     primary nu (cubic spline, 5 knots at quantiles)
  rg     swelling relative to an ideal chain, alpha = log(<Rg>_300 / N^0.5) (same spline); independent of the
         power-law fit, so it also holds where nu is ill-defined
  both   splines of nu and alpha together

Targets:
  shared      mean of the two z-scored residuals (what both models agree on beyond compaction)
  difference  z(CALVADOS residual) - z(Mpipi residual) (what separates the models beyond compaction)
  CALVADOS residual, Mpipi residual, and the raw chi of each model (reference)
Noise ceilings from the per-temperature block-jackknife SEs of chi (sampling noise only).

Predictors: Test 2's sequence layers (test2_slope_vs_sequence.sequence_layers), rank-normalised, cumulative
linear R^2 per layer, then quadratic and kNN on all; also each of three physical groups on its own
(stickers, charge, hydropathy), and univariate Spearman correlations with every descriptor.
Same CV code as Test 2 (5-fold, ridge penalty by inner 3-fold CV).

Sets: 'fit_r2' = proteins with status ok at all five Mpipi-GG temperatures and a power-law nu fit R^2 >= 0.99
in both models (nu is only a compaction measure where the fit holds); 'all_paired' = every protein ok in both.

Output: Traj Level Analysis/Analysis_FDT/test2b_cross_model_residual/ (r2_<set>[_<measure>].csv,
spearman_<set>[_<measure>].csv, residuals_<set>[_<measure>].csv; no suffix for nu)
Usage: python3 test2b_cross_model_residual.py [--compaction nu|rg|both]
"""
import sys

sys.dont_write_bytecode = True  # keep the NAS code folder free of __pycache__

import argparse  # noqa: E402
import os  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

HERE = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, HERE)
import test2_slope_vs_sequence as t2  # noqa: E402

t1, L, bc = t2.t1, t2.L, t2.bc
MPIPI_DIR = '/Volumes/User Homes/Nu Project/BENDER-MT/Traj Level Analysis/Analysis_FDT_MPIPI'
OUT_DIR = os.path.join(L.OUT_ROOT, 'test2b_cross_model_residual')
TEMPS = [278, 288, 300, 310, 320]
R2_MIN = 0.99
GROUPS = {
    'stickers': ['f_F', 'f_W', 'f_Y', 'f_R', 'lambda_avg', 'f_aro', 'f_pi', 'f_pi_extended', 'log1p_aromatic_pairs',
                 'log1p_cation_pi_pairs'],
    'charge': ['FCR', 'NCPR', 'uversky_R', 'kappa', 'SCD', 'f_D', 'f_E', 'f_K'],
    'hydropathy': ['hydropathy_avg', 'uversky_H', 'SHD', 'f_hydrophobic', 'f_polar'],
}


def chi_rel(frames: list) -> pd.DataFrame:
    """chi(Rg) and its sampling SE as % of <Rg> per K, averaged over the temperatures."""
    d = pd.concat(frames)
    g = d.groupby('UniProt_ID')
    A = g['A_mean'].mean()
    return pd.DataFrame({'chi': 100 * g['chi_total'].mean() / A,
                         'se': 100 * np.sqrt(g['chi_total_se'].apply(lambda s: (s ** 2).sum())) / g.size() / A,
                         'n_T': g.size()})


def spline(x: np.ndarray, knots: int = 5) -> np.ndarray:
    q = np.quantile(x, np.linspace(0, 1, knots + 2)[1:-1])
    return np.column_stack([x, x ** 2, x ** 3] + [np.clip(x - t, 0, None) ** 3 for t in q])


def residualise(y: np.ndarray, X: np.ndarray) -> np.ndarray:
    X1 = np.column_stack([np.ones(len(y)), X])
    return y - X1 @ np.linalg.lstsq(X1, y, rcond=None)[0]


def load_pairs() -> pd.DataFrame:
    rg = lambda df: df[df['metric'] == 'rg']  # noqa: E731
    cal_T = {T: rg(pd.read_csv(os.path.join(L.OUT_ROOT, f'fdt_results_{T}K.csv'),
                               usecols=['UniProt_ID', 'metric', 'chi_total', 'chi_total_se', 'A_mean'])) for T in TEMPS}
    cal = chi_rel(list(cal_T.values()))
    cal['rg300'] = cal_T[300].set_index('UniProt_ID')['A_mean']
    ok, mp = None, []
    for T in TEMPS:
        s = pd.read_csv(os.path.join(MPIPI_DIR, f'mpipi_{T}K.csv'), usecols=['UniProt_ID', 'status', 'nu', 'nu_r2'])
        ids = set(s.loc[s['status'] == 'ok', 'UniProt_ID'])
        ok = ids if ok is None else ok & ids
        if T == 300:
            nu_m = s.set_index('UniProt_ID')[['nu', 'nu_r2']].rename(columns={'nu': 'nu_M', 'nu_r2': 'nu_r2_M'})
        mp.append(rg(pd.read_csv(os.path.join(MPIPI_DIR, f'mpipi_{T}K_fdt.csv'),
                                 usecols=['UniProt_ID', 'metric', 'chi_total', 'chi_total_se', 'A_mean'])))
    mpi = chi_rel(mp)
    mpi['rg300'] = mp[TEMPS.index(300)].set_index('UniProt_ID')['A_mean']
    mpi = mpi[(mpi['n_T'] == len(TEMPS)) & mpi.index.isin(ok)]
    cal = cal[cal['n_T'] == len(TEMPS)]
    bio = pd.read_csv(bc.BIO_CSV, low_memory=False).set_index('UniProt_ID')
    j = cal.add_suffix('_C').join(mpi.add_suffix('_M'), how='inner').join(nu_m)
    j = j.join(bio[['nu', 'nu_fit_r2', 'length']].rename(columns={'nu': 'nu_C', 'nu_fit_r2': 'nu_r2_C'}))
    return j, bio


def compaction_design(j: pd.DataFrame, m: str, measure: str) -> np.ndarray:
    logN = np.log(j['length'].to_numpy())
    alpha = np.log(j[f'rg300_{m}'].to_numpy()) - 0.5 * logN
    parts = {'nu': [spline(j[f'nu_{m}'].to_numpy())], 'rg': [spline(alpha)],
             'both': [spline(j[f'nu_{m}'].to_numpy()), spline(alpha)]}[measure]
    return np.column_stack(parts + [logN])


def main(measure: str = 'nu'):
    os.makedirs(OUT_DIR, exist_ok=True)
    suffix = '' if measure == 'nu' else f'_{measure}'
    pairs, bio_all = load_pairs()
    sets = {'fit_r2': pairs[(pairs['nu_r2_C'] >= R2_MIN) & (pairs['nu_r2_M'] >= R2_MIN)], 'all_paired': pairs}
    for set_name, j in sets.items():
        j = j.copy()
        bio = bio_all.reindex(j.index)
        seq, layers = t2.sequence_layers(bio)
        if seq.isna().any().any():
            raise ValueError(f'missing sequence features: {seq.columns[seq.isna().any()].tolist()}')
        rn = lambda cols: np.column_stack([t1.rank_normal(seq[c].to_numpy()) for c in cols])  # noqa: E731
        for m in ('C', 'M'):
            j[f'res_{m}'] = residualise(j[f'chi_{m}'].to_numpy(), compaction_design(j, m, measure))
        sdC, sdM = j['res_C'].std(), j['res_M'].std()
        zC, zM = j['res_C'] / sdC, j['res_M'] / sdM
        nC, nM = (j['se_C'] / sdC) ** 2, (j['se_M'] / sdM) ** 2   # sampling-noise variance in z units
        j['shared'] = (zC + zM) / 2
        j['difference'] = zC - zM
        targets = {'shared': (j['shared'], (nC + nM) / 4), 'difference': (j['difference'], nC + nM),
                   'residual CALVADOS': (j['res_C'], j['se_C'] ** 2), 'residual Mpipi': (j['res_M'], j['se_M'] ** 2),
                   'chi CALVADOS (raw)': (j['chi_C'], j['se_C'] ** 2), 'chi Mpipi (raw)': (j['chi_M'], j['se_M'] ** 2)}
        rows = []
        all_cols = [c for _, cols in layers for c in cols]
        Xall = rn(all_cols)
        for name, (y, noise) in targets.items():
            y = y.to_numpy()
            row = {'target': name, 'n': len(y), 'noise_ceiling': 1 - noise.mean() / y.var(ddof=1)}
            cols = []
            for lname, lc in layers:
                cols += lc
                row[f'R2 {lname} (cum.)'] = t2.cv_r2(rn(cols), y, 'linear')
            for g, gc in GROUPS.items():
                row[f'R2 {g} only'] = t2.cv_r2(rn(gc), y, 'linear')
            row['R2 quadratic all'] = t2.cv_r2(Xall, y, 'quadratic')
            row['R2 knn all'] = t2.cv_r2(Xall, y, 'knn')
            row['R2 best'] = max(row[f'R2 {layers[-1][0]} (cum.)'], row['R2 quadratic all'], row['R2 knn all'])
            row['best share of explainable'] = row['R2 best'] / row['noise_ceiling']
            rows.append(row)
            print(f"[{set_name}/{measure}] {name:20s} n {len(y)}  ceiling {row['noise_ceiling']:.3f}  best {row['R2 best']:.3f}", flush=True)
        out = pd.DataFrame(rows)
        path = os.path.join(OUT_DIR, f'r2_{set_name}{suffix}.csv')
        out.to_csv(path, index=False)
        sp = pd.DataFrame({t: {c: stats.spearmanr(seq[c], j[t])[0] for c in all_cols} for t in ('shared', 'difference', 'res_C', 'res_M')})
        sp.index.name = 'descriptor'
        sp.to_csv(os.path.join(OUT_DIR, f'spearman_{set_name}{suffix}.csv'))
        j[['nu_C', 'nu_M', 'rg300_C', 'rg300_M', 'length', 'chi_C', 'chi_M', 'res_C', 'res_M', 'shared', 'difference']].to_csv(
            os.path.join(OUT_DIR, f'residuals_{set_name}{suffix}.csv'))
        bc.write_provenance(path, __file__, {'set': set_name, 'nu_fit_r2_min': R2_MIN if set_name == 'fit_r2' else None,
                                             'compaction': measure,
                                             'residualisation': 'own compaction (cubic spline(s), 5 knots) + log N, OLS',
                                             'groups': GROUPS, 'layers': {n: c for n, c in layers}},
                            [bc.BIO_CSV] + [os.path.join(L.OUT_ROOT, f'fdt_results_{T}K.csv') for T in TEMPS]
                            + [os.path.join(MPIPI_DIR, f'mpipi_{T}K{s}.csv') for T in TEMPS for s in ('', '_fdt')])
        with pd.option_context('display.width', 300, 'display.max_columns', 40):
            print(out.to_string(index=False, float_format=lambda v: f'{v:.3f}'))
            top = sp.reindex(sp['shared'].abs().sort_values(ascending=False).index).head(12)
            print('\nstrongest univariate Spearman with the shared residual:')
            print(top.round(3).to_string())
            top = sp.reindex(sp['difference'].abs().sort_values(ascending=False).index).head(8)
            print('\nstrongest with the difference (CALVADOS minus Mpipi):')
            print(top.round(3).to_string())
            print(f"\ncorr(res_C, res_M) = {np.corrcoef(j['res_C'], j['res_M'])[0, 1]:.3f}\n")
    print(f'saved {OUT_DIR}')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--compaction', choices=['nu', 'rg', 'both'], default='nu')
    main(ap.parse_args().compaction)
