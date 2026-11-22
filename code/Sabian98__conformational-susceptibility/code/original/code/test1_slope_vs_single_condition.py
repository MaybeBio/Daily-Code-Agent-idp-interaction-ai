"""
Test 1 (Paper 1, E1/E3): do the temperature-slope features correlate with single-condition (300 K) features?

Set: proteins whose trajectories are complete (not truncated) at all five temperatures, 278-320 K, taken from
the FDT tables (Traj Level Analysis/Analysis_FDT/fdt_results_<T>K.csv; 1,000 production frames per temperature).

Slope feature   slope_A = least-squares slope of <A>(T) over the five temperatures (per K), for Rg, Ree,
                asphericity, the five contact-graph metrics (Calpha < 0.8 nm, |i-j| >= 3) and nu (Analysis_<T>K_v2).
                Its sampling SE follows from the block SEs of <A> at each temperature.
State features  single condition, 300 K: nu, log N, Rg, Ree, asphericity and the five graph metrics.

Per slope feature:
  Spearman rho with every state feature
  5-fold CV R^2 from the state features: linear, quadratic + interactions (ridge), k-nearest neighbours (k=50)
  on rank-normalised features; also log N alone
  reliability = 1 - mean(SE^2) / var(slope): the most variance any predictor can explain, since the rest is
  sampling noise; R2 / reliability is the share of the explainable variance that the state features explain
  single-condition fluctuations: Spearman rho(chi_FDT at 300 K, slope) and CV R^2 of state + chi_FDT(300 K),
  chi_FDT = -Cov(A, G_T) from the 300 K trajectory alone

Output: Traj Level Analysis/Analysis_FDT/test1_slope_vs_single_condition/
  features.csv           one row per protein: state features, slopes and slope SEs, chi_FDT(300 K)
  spearman.csv           slope x state Spearman matrix
  r2_summary.csv         the per-slope summary printed at the end

Usage: python3 test1_slope_vs_single_condition.py
"""

import sys

sys.dont_write_bytecode = True  # keep the NAS script folder free of __pycache__

import os  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

HERE = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, HERE)
import TODO_T2_trajectory_loader as L  # noqa: E402

bc = L.bc
TEMPS = L.TEMPS
T_REF = 300
FDT_METRICS = L.GEOMETRY + L.GRAPH_METRICS
OUT_DIR = os.path.join(L.OUT_ROOT, 'test1_slope_vs_single_condition')
N_FOLDS, K_NN, RIDGE, SEED = 5, 50, 1.0, 0


def rank_normal(x: np.ndarray) -> np.ndarray:
    """Inverse-normal transform of ranks (robust to heavy tails and scale)."""
    return stats.norm.ppf((stats.rankdata(x) - 0.5) / len(x))


def slope_and_se(means: np.ndarray, ses: np.ndarray) -> tuple:
    """Row-wise least-squares slope over TEMPS and its SE from independent per-temperature SEs."""
    t = np.array(TEMPS, dtype=float)
    w = (t - t.mean()) / ((t - t.mean()) ** 2).sum()
    return means @ w, np.sqrt((ses ** 2) @ (w ** 2))


def build_features() -> pd.DataFrame:
    df = L.load_results(columns=['A_mean', 'A_mean_se', 'chi_total'])
    complete = df.groupby('UniProt_ID').agg(n_T=('T', 'nunique'), truncated=('truncated', 'any'))
    keep = complete.index[(complete['n_T'] == len(TEMPS)) & ~complete['truncated']]
    df = df[df['UniProt_ID'].isin(keep)]
    wide = df.pivot_table(index='UniProt_ID', columns=['metric', 'T'], values=['A_mean', 'A_mean_se', 'chi_total'])

    nu = {}
    for T in TEMPS:
        v = pd.read_csv(os.path.join(L.V2_DIR_FMT.format(T=T), 'scalar_features_v2.csv'),
                        usecols=['UniProt_ID', 'nu', 'nu_block_se', 'n_frames_dcd'])
        nu[T] = v.drop_duplicates('UniProt_ID', keep='last').set_index('UniProt_ID')

    f = pd.DataFrame(index=wide.index)
    f['n_residues'] = df.drop_duplicates('UniProt_ID').set_index('UniProt_ID')['n_residues']
    f['logN'] = np.log(f['n_residues'])
    f['nu_300'] = nu[T_REF]['nu'].reindex(f.index)
    for m in FDT_METRICS:
        f[f'{m}_300'] = wide[('A_mean', m, T_REF)]
        s, se = slope_and_se(np.column_stack([wide[('A_mean', m, T)] for T in TEMPS]),
                             np.column_stack([wide[('A_mean_se', m, T)] for T in TEMPS]))
        f[f'slope_{m}'], f[f'slope_{m}_se'] = s, se
        f[f'chiFDT300_{m}'] = wide[('chi_total', m, T_REF)]
    s, se = slope_and_se(np.column_stack([nu[T]['nu'].reindex(f.index) for T in TEMPS]),
                         np.column_stack([nu[T]['nu_block_se'].reindex(f.index) for T in TEMPS]))
    f['slope_nu'], f['slope_nu_se'] = s, se

    # the v2 nu values must come from the same DCDs the FDT tables used
    dcd = df[df['metric'] == 'rg'].pivot_table(index='UniProt_ID', columns='T', values='n_frames_dcd')
    for T in TEMPS:
        bad = (nu[T]['n_frames_dcd'].reindex(f.index) != dcd[T].reindex(f.index)).sum()
        if bad:
            raise ValueError(f'{bad} proteins at {T} K: v2 nu from a different DCD than the FDT tables')
    return f.dropna()


# ---------------------------------------------------------------------------
# cross-validated regressions (numpy only)
# ---------------------------------------------------------------------------

def _folds(n: int):
    idx = np.random.default_rng(SEED).permutation(n)
    return np.array_split(idx, N_FOLDS)


def _lstsq_predict(Xtr, ytr, Xte, ridge=0.0):
    mu, sd = Xtr.mean(0), Xtr.std(0) + 1e-12
    A, B = (Xtr - mu) / sd, (Xte - mu) / sd
    A1, B1 = np.column_stack([np.ones(len(A)), A]), np.column_stack([np.ones(len(B)), B])
    reg = ridge * np.eye(A1.shape[1])
    reg[0, 0] = 0.0
    coef = np.linalg.solve(A1.T @ A1 + reg, A1.T @ ytr)
    return B1 @ coef


def _quad(X):
    i, j = np.triu_indices(X.shape[1], k=1)
    return np.column_stack([X, X ** 2, X[:, i] * X[:, j]])


def _knn_predict(Xtr, ytr, Xte, k=K_NN, chunk=500):
    out = np.empty(len(Xte))
    sq = (Xtr ** 2).sum(1)
    for a in range(0, len(Xte), chunk):
        b = Xte[a:a + chunk]
        d = (b ** 2).sum(1)[:, None] + sq[None, :] - 2.0 * b @ Xtr.T
        nn = np.argpartition(d, k, axis=1)[:, :k]
        out[a:a + chunk] = ytr[nn].mean(1)
    return out


def cv_r2(X: np.ndarray, y: np.ndarray, model: str) -> float:
    pred = np.empty(len(y))
    for te in _folds(len(y)):
        tr = np.setdiff1d(np.arange(len(y)), te)
        if model == 'linear':
            pred[te] = _lstsq_predict(X[tr], y[tr], X[te])
        elif model == 'quadratic':
            pred[te] = _lstsq_predict(_quad(X[tr]), y[tr], _quad(X[te]), ridge=RIDGE)
        else:
            pred[te] = _knn_predict(X[tr], y[tr], X[te])
    return 1.0 - ((y - pred) ** 2).sum() / ((y - y.mean()) ** 2).sum()


# ---------------------------------------------------------------------------

def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    f = build_features()
    f.to_csv(os.path.join(OUT_DIR, 'features.csv'))
    state = ['nu_300', 'logN'] + [f'{m}_300' for m in FDT_METRICS]
    targets = [f'slope_{m}' for m in FDT_METRICS] + ['slope_nu']
    Z = np.column_stack([rank_normal(f[c].to_numpy()) for c in state])

    rho = pd.DataFrame({c: [stats.spearmanr(f[t], f[c])[0] for t in targets] for c in state}, index=targets)
    rho.to_csv(os.path.join(OUT_DIR, 'spearman.csv'))

    rows = []
    for t in targets:
        y = f[t].to_numpy()
        rel = 1.0 - (f[f'{t}_se'] ** 2).mean() / y.var(ddof=1)
        r2 = {m: cv_r2(Z, y, m) for m in ('linear', 'quadratic', 'knn')}
        best = max(r2.values())
        yr = rank_normal(y)
        row = {'slope': t, 'n': len(y), 'reliability': rel,
               'strongest_state_feature': rho.loc[t].abs().idxmax(), 'rho_strongest': rho.loc[t, rho.loc[t].abs().idxmax()],
               'R2_logN': cv_r2(Z[:, [1]], y, 'quadratic'),
               'R2_linear': r2['linear'], 'R2_quadratic': r2['quadratic'], 'R2_knn': r2['knn'], 'R2_best': best,
               'explained_share_of_explainable': best / rel if rel > 0 else np.nan,
               'R2_best_rank_target': max(cv_r2(Z, yr, m) for m in ('linear', 'quadratic', 'knn'))}
        m = t.replace('slope_', '')
        if m in FDT_METRICS:
            chi = f[f'chiFDT300_{m}'].to_numpy()
            row['rho_chiFDT300'] = stats.spearmanr(chi, y)[0]
            row['R2_state_plus_chiFDT300'] = cv_r2(np.column_stack([Z, rank_normal(chi)]), y, 'quadratic')
        rows.append(row)
    summ = pd.DataFrame(rows)
    summ.to_csv(os.path.join(OUT_DIR, 'r2_summary.csv'), index=False)
    bc.write_provenance(os.path.join(OUT_DIR, 'r2_summary.csv'), __file__,
                        {'set': 'complete at all five temperatures', 'T_ref': T_REF, 'folds': N_FOLDS,
                         'k_nn': K_NN, 'ridge': RIDGE, 'seed': SEED, 'state_features': state},
                        [os.path.join(L.OUT_ROOT, f'fdt_results_{T}K.csv') for T in TEMPS]
                        + [os.path.join(L.V2_DIR_FMT.format(T=T), 'scalar_features_v2.csv') for T in TEMPS])
    with pd.option_context('display.width', 250, 'display.max_columns', 30):
        print(f'{len(f)} proteins complete at all five temperatures\n')
        print(summ.to_string(index=False, float_format=lambda v: f'{v:.3f}'))
        print('\nSpearman rho, slope x 300 K state:')
        print(rho.to_string(float_format=lambda v: f'{v:+.2f}'))
    print(f'\nsaved {OUT_DIR}')


if __name__ == '__main__':
    main()
