"""
Test 2 (Paper 1, E2): how much of the temperature-slope features can sequence descriptors explain?

Same set and slopes as Test 1 (test1_slope_vs_single_condition.py): proteins complete at all five
temperatures, least-squares slope of <A>(T) over 278-320 K, noise ceiling from the per-temperature block SEs.

Sequence layers, added cumulatively (BENDER_BIO.csv, recomputed descriptors):
  L1 length        log N
  L2 composition   the 20 residue fractions, mean and SD of CALVADOS-2 lambda, mean Kyte-Doolittle
                   hydropathy, FCR, NCPR, aromatic / hydrophobic / polar fractions
  L3 patterning    kappa (Das-Pappu), SCD, SHD
  L4 Uversky       <H> (Kyte-Doolittle rescaled to 0-1), <R> = |NCPR|, dCH = <R> - (2.785 <H> - 1.151)
                   (distance above the charge-hydropathy boundary), Dunker disorder / order-promoting fractions
  L5 pi            f_pi, f_pi_extended, aromatic and cation-pi pair counts (log1p); BENDER_BIO pi_score is empty
Contact frequencies from the simulations are not sequence features and are not used.

Models (5-fold CV, ridge penalty chosen by inner 3-fold CV on each training split): linear on each cumulative
layer set; quadratic + interactions and k-nearest neighbours (k=50) on all layers; predictors rank-normalised.
Reported per target: cumulative linear R^2 per layer, best R^2 with all layers, noise ceiling, best / ceiling,
and for the slopes also R^2 from the 300 K state (Test 1 features) and from state + sequence (both quadratic).
The same layers are also fitted to the 300 K state features themselves, as a reference.

Output: Traj Level Analysis/Analysis_FDT/test2_slope_vs_sequence/  (layers_r2.csv, sequence_features.csv)
Usage: python3 test2_slope_vs_sequence.py
"""

import sys

sys.dont_write_bytecode = True  # keep the NAS code folder free of __pycache__

import os  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

HERE = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, HERE)
import test1_slope_vs_single_condition as t1  # noqa: E402

L, bc = t1.L, t1.bc
OUT_DIR = os.path.join(L.OUT_ROOT, 'test2_slope_vs_sequence')
ALPHAS = [0.1, 1.0, 10.0, 100.0, 1000.0]
AA = 'ACDEFGHIKLMNPQRSTVWY'
KD = dict(A=1.8, R=-4.5, N=-3.5, D=-3.5, C=2.5, Q=-3.5, E=-3.5, G=-0.4, H=-3.2, I=4.5, L=3.8, K=-3.9,
          M=1.9, F=2.8, P=-1.6, S=-0.8, T=-0.7, W=-0.9, Y=-1.3, V=4.2)


def sequence_layers(bio: pd.DataFrame) -> tuple:
    """Per-protein sequence descriptors and the cumulative layer definitions."""
    s = pd.DataFrame(index=bio.index)
    s['logN'] = np.log(bio['length'])
    for a in AA:
        s[f'f_{a}'] = bio[f'f_{a}']
    for c in ['lambda_avg', 'lambda_std', 'hydropathy_avg', 'FCR', 'NCPR', 'f_aro', 'f_hydrophobic', 'f_polar']:
        s[c] = bio[c]
    for c in ['kappa', 'SCD', 'SHD']:
        s[c] = bio[c].fillna(bio[c].median())  # kappa is undefined without both charge signs
    # mean over standard residues: 3 protist sequences carry X, which has no hydropathy
    s['uversky_H'] = bio['sequence'].map(lambda q: np.mean([(KD[c] + 4.5) / 9.0 for c in q if c in KD]))
    s['uversky_R'] = bio['NCPR'].abs()
    s['dCH'] = s['uversky_R'] - (2.785 * s['uversky_H'] - 1.151)
    s['f_disorder_promoting'] = bio['f_disorder_promoting']
    s['f_order_promoting'] = bio['f_order_promoting']
    for c in ['f_pi', 'f_pi_extended']:
        s[c] = bio[c]
    s['log1p_aromatic_pairs'] = np.log1p(bio['n_aromatic_pairs'])
    s['log1p_cation_pi_pairs'] = np.log1p(bio['n_cation_pi_pairs'])
    layers = [
        ('L1 length', ['logN']),
        ('L2 composition', [f'f_{a}' for a in AA] + ['lambda_avg', 'lambda_std', 'hydropathy_avg', 'FCR', 'NCPR',
                                                     'f_aro', 'f_hydrophobic', 'f_polar']),
        ('L3 patterning', ['kappa', 'SCD', 'SHD']),
        ('L4 Uversky', ['uversky_H', 'uversky_R', 'dCH', 'f_disorder_promoting', 'f_order_promoting']),
        ('L5 pi', ['f_pi', 'f_pi_extended', 'log1p_aromatic_pairs', 'log1p_cation_pi_pairs']),
    ]
    return s, layers


def _ridge_all_alphas(Xtr, ytr, Xte):
    """Ridge predictions for every penalty in ALPHAS from one Gram matrix (standardised features, free intercept)."""
    mu, sd = Xtr.mean(0), Xtr.std(0) + 1e-12
    A = np.column_stack([np.ones(len(Xtr)), (Xtr - mu) / sd])
    B = np.column_stack([np.ones(len(Xte)), (Xte - mu) / sd])
    G, b = A.T @ A, A.T @ ytr
    out = {}
    for a in ALPHAS:
        reg = a * np.eye(G.shape[0])
        reg[0, 0] = 0.0
        out[a] = B @ np.linalg.solve(G + reg, b)
    return out


def _ridge_nested(X, y, model, tr, te):
    """Fit on tr with the ridge penalty chosen by inner 3-fold CV; predict te."""
    feat = t1._quad if model == 'quadratic' else (lambda Z: Z)
    inner = np.array_split(np.random.default_rng(1).permutation(tr), 3)
    err = {a: 0.0 for a in ALPHAS}
    for v in inner:
        fit = np.setdiff1d(tr, v)
        for a, p in _ridge_all_alphas(feat(X[fit]), y[fit], feat(X[v])).items():
            err[a] += ((y[v] - p) ** 2).sum()
    best = min(err, key=err.get)
    return _ridge_all_alphas(feat(X[tr]), y[tr], feat(X[te]))[best]


def cv_r2(X, y, model):
    pred = np.empty(len(y))
    for te in t1._folds(len(y)):
        tr = np.setdiff1d(np.arange(len(y)), te)
        pred[te] = t1._knn_predict(X[tr], y[tr], X[te]) if model == 'knn' else _ridge_nested(X, y, model, tr, te)
    return 1.0 - ((y - pred) ** 2).sum() / ((y - y.mean()) ** 2).sum()


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    f = t1.build_features()
    bio = pd.read_csv(bc.BIO_CSV, low_memory=False).set_index('UniProt_ID').reindex(f.index)
    seq, layers = sequence_layers(bio)
    if seq.isna().any().any():
        raise ValueError(f'missing sequence features: {seq.columns[seq.isna().any()].tolist()}')
    seq.to_csv(os.path.join(OUT_DIR, 'sequence_features.csv'))
    rn = lambda cols: np.column_stack([t1.rank_normal(seq[c].to_numpy()) for c in cols])  # noqa: E731

    # noise ceilings: slopes from Test 1's SEs; 300 K state from the block SEs at 300 K
    res = L.load_results(columns=['A_mean_se'])
    se300 = res[res['T'] == t1.T_REF].pivot_table(index='UniProt_ID', columns='metric', values='A_mean_se').reindex(f.index)
    nu_se = pd.read_csv(os.path.join(L.V2_DIR_FMT.format(T=t1.T_REF), 'scalar_features_v2.csv'),
                        usecols=['UniProt_ID', 'nu_block_se']).drop_duplicates('UniProt_ID', keep='last') \
        .set_index('UniProt_ID')['nu_block_se'].reindex(f.index)
    state = ['nu_300', 'logN'] + [f'{m}_300' for m in t1.FDT_METRICS]
    Zstate = np.column_stack([t1.rank_normal(f[c].to_numpy()) for c in state])
    targets = [(f'slope_{m}', f[f'slope_{m}_se']) for m in t1.FDT_METRICS] + [('slope_nu', f['slope_nu_se'])]
    targets += [(f'{m}_300', se300[m]) for m in t1.FDT_METRICS] + [('nu_300', nu_se)]

    all_cols = [c for _, cols in layers for c in cols]
    Xall = rn(all_cols)
    rows = []
    for t, se in targets:
        y = f[t].to_numpy()
        row = {'target': t, 'kind': 'slope' if t.startswith('slope_') else 'state 300 K', 'n': len(y),
               'noise_ceiling': 1.0 - (se ** 2).mean() / y.var(ddof=1)}
        cols = []
        for name, lc in layers:
            cols += lc
            row[f'R2 {name} (cum.)'] = cv_r2(rn(cols), y, 'linear')
        row['R2 quadratic all'] = cv_r2(Xall, y, 'quadratic')
        row['R2 knn all'] = cv_r2(Xall, y, 'knn')
        row['R2 sequence best'] = max(row[f'R2 {layers[-1][0]} (cum.)'], row['R2 quadratic all'], row['R2 knn all'])
        row['sequence share of explainable'] = row['R2 sequence best'] / row['noise_ceiling']
        if row['kind'] == 'slope':
            row['R2 state (Test 1)'] = cv_r2(Zstate, y, 'quadratic')
            row['R2 state + sequence'] = cv_r2(np.column_stack([Zstate, Xall]), y, 'quadratic')
        rows.append(row)
        print(f"{t:28s} ceiling {row['noise_ceiling']:.3f}  sequence best {row['R2 sequence best']:.3f}", flush=True)
    out = pd.DataFrame(rows)
    path = os.path.join(OUT_DIR, 'layers_r2.csv')
    out.to_csv(path, index=False)
    bc.write_provenance(path, __file__, {'layers': {n: c for n, c in layers}, 'alphas': ALPHAS,
                                         'folds': t1.N_FOLDS, 'k_nn': t1.K_NN, 'seed': t1.SEED},
                        [bc.BIO_CSV, os.path.join(t1.OUT_DIR, 'features.csv')])
    with pd.option_context('display.width', 300, 'display.max_columns', 40):
        print(out.to_string(index=False, float_format=lambda v: f'{v:.3f}'))
    print(f'\nsaved {OUT_DIR}')


if __name__ == '__main__':
    main()
