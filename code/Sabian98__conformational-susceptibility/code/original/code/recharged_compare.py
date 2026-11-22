"""
recharged_compare.py — Mpipi-Recharged as a third model at 300 K, next to CALVADOS-2 and Mpipi-GG at 300 K.

All three at one temperature (Mpipi-Recharged has 300 K only), chi(Rg) as % of <Rg> per K. The deck's numbers average
278-320 K; at 300 K alone sampling noise is larger (cross-model leftover r 0.28 instead of 0.44 for CALVADOS-2 vs
Mpipi-GG), so compare these tables with each other, not with the deck.

  A  agreement: Spearman of chi between models, sign agreement; overall, per CALVADOS-2 nu bin and per kingdom
  B  chi by each model's own nu bin (deck slide 12)
  C  expanded chains (CALVADOS-2 nu > 0.63): Boltzmann part, dH/dT part (CALVADOS-2 dielectric, Mpipi-Recharged
     screening, Mpipi-GG none) and net, medians per model (deck slide 14)
  D  energy parts: sticker S, Boltzmann charge Q, dH/dT dQ, bond B; charge-held (Q > 0 and Q > S) share per model and
     nu bin, and the overlap of charge-held chains (deck slides 24-25)
  E  leftover (own nu spline + log N, as Test 2b): cross-model correlations, variance share of each part, for each
     model and for the three-model shared leftover (deck slides 19, 21)

Inputs: Analysis_FDT_RECHARGED/recharged_fdt_300K.csv and recharged_300K_nu.csv (mpipi_recharged_fdt.py,
recharged_nu.py), Analysis_FDT/fdt_results_300K.csv + BENDER_BIO nu, Analysis_FDT_MPIPI/mpipi_300K{,_fdt}.csv.
Output: Analysis_FDT_RECHARGED/comparison_300K/*.csv and summary.txt.
"""
import sys

sys.dont_write_bytecode = True

import os  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402

TLA = '/Volumes/User Homes/Nu Project/BENDER-MT/Traj Level Analysis'
BIO_CSV = '/Volumes/User Homes/Nu Project/BENDER-BIO/BENDER_BIO.csv'
RCH = os.path.join(TLA, 'Analysis_FDT_RECHARGED')
OUT = os.path.join(RCH, 'comparison_300K')
BINS = [0.30, 0.45, 0.50, 0.54, 0.58, 0.63, 1.0]
R2_MIN = 0.99
MODELS = {'C': 'CALVADOS-2', 'M': 'Mpipi-GG', 'R': 'Mpipi-Recharged'}
PARTS = ['S', 'Q', 'dQ', 'B', 'rest']


def spline(x: np.ndarray, knots: int = 5) -> np.ndarray:
    """Same basis as test2b_cross_model_residual.spline."""
    q = np.quantile(x, np.linspace(0, 1, knots + 2)[1:-1])
    return np.column_stack([x, x ** 2, x ** 3] + [np.clip(x - t, 0, None) ** 3 for t in q])


def residualise(y: np.ndarray, X: np.ndarray) -> np.ndarray:
    X1 = np.column_stack([np.ones(len(y)), X])
    return y - X1 @ np.linalg.lstsq(X1, y, rcond=None)[0]


def rel(d: pd.DataFrame, cols: dict) -> pd.DataFrame:
    """chi and parts as % of <Rg> per K."""
    d = d[d.metric == 'rg'].set_index('UniProt_ID')
    out = pd.DataFrame({k: 100 * v(d) / d.A_mean for k, v in cols.items()})
    out['rg'] = d.A_mean
    out['N'] = d.n_residues
    return out


def load() -> dict:
    c = pd.read_csv(os.path.join(TLA, 'Analysis_FDT', 'fdt_results_300K.csv'))
    C = rel(c, {'chi': lambda d: d.chi_total, 'boltz': lambda d: d.chi_Boltzmann, 'S': lambda d: d.chi_AH,
                'Q': lambda d: d.chi_DH - d.chi_dielectric, 'dQ': lambda d: d.chi_dielectric, 'B': lambda d: d.chi_bond})
    bio = pd.read_csv(BIO_CSV, usecols=['UniProt_ID', 'nu', 'nu_fit_r2', 'FCR', 'kingdom', 'organism'], low_memory=False).set_index('UniProt_ID')
    C = C.join(bio.rename(columns={'nu_fit_r2': 'nu_r2'}))

    m = pd.read_csv(os.path.join(TLA, 'Analysis_FDT_MPIPI', 'mpipi_300K_fdt.csv'))
    M = rel(m, {'chi': lambda d: d.chi_total, 'boltz': lambda d: d.chi_total, 'S': lambda d: d.chi_WF,
                'Q': lambda d: d.chi_DH, 'dQ': lambda d: 0 * d.chi_DH, 'B': lambda d: d.chi_bond})
    s = pd.read_csv(os.path.join(TLA, 'Analysis_FDT_MPIPI', 'mpipi_300K.csv'), usecols=['UniProt_ID', 'status', 'nu', 'nu_r2'])
    s = s[s.status == 'ok'].set_index('UniProt_ID')
    M = M.join(s[['nu', 'nu_r2']], how='inner')

    r = pd.read_csv(os.path.join(RCH, 'recharged_fdt_300K.csv'))
    r = r[r.max_bond_nm <= 1.0]   # a bond > 1 nm after unwrapping is a broken run (Q983T3: frames at 1e13 kJ/mol)
    R = rel(r, {'chi': lambda d: d.chi_total, 'boltz': lambda d: d.chi_Boltzmann, 'S': lambda d: d.chi_WF,
                'Q': lambda d: d.chi_DH - d.chi_screening, 'dQ': lambda d: d.chi_screening, 'B': lambda d: d.chi_bond})
    nu = pd.read_csv(os.path.join(RCH, 'recharged_300K_nu.csv')).set_index('UniProt_ID')
    R = R.join(nu[['nu', 'nu_r2']], how='inner')
    for d in (C, M, R):   # Mpipi-GG: per-term energies vs the logged total leave a small rest (as in energy_mechanism_split)
        d['rest'] = d.chi - d[PARTS[:-1]].sum(axis=1)
    return {'C': C.dropna(subset=['chi', 'nu']), 'M': M.dropna(subset=['chi', 'nu']), 'R': R.dropna(subset=['chi', 'nu'])}


def agreement(D: dict) -> pd.DataFrame:
    rows = []
    trio = D['C'].index.intersection(D['M'].index).intersection(D['R'].index)
    for a, b in (('R', 'C'), ('R', 'M'), ('C', 'M')):
        for set_name, ids in (('pair', D[a].index.intersection(D[b].index)), ('three-model set', trio)):
            x, y = D[a].chi.reindex(ids), D[b].chi.reindex(ids)
            rows.append({'pair': f'{MODELS[a]} vs {MODELS[b]}', 'set': set_name, 'bin_nu_C': 'all', 'n': len(ids),
                         'spearman': spearmanr(x, y)[0], 'same_sign': (np.sign(x) == np.sign(y)).mean()})
    ids = D['R'].index.intersection(D['C'].index)
    b = pd.cut(D['C'].nu.reindex(ids), BINS)
    for lab, g in pd.Series(ids, index=ids).groupby(b, observed=True):
        x, y = D['R'].chi.reindex(g.index), D['C'].chi.reindex(g.index)
        rows.append({'pair': 'Mpipi-Recharged vs CALVADOS-2', 'set': 'pair', 'bin_nu_C': str(lab), 'n': len(g),
                     'spearman': spearmanr(x, y)[0], 'same_sign': (np.sign(x) == np.sign(y)).mean()})
    return pd.DataFrame(rows)


def by_kingdom(D: dict) -> pd.DataFrame:
    """Agreement within each kingdom (BENDER_BIO kingdom of the CALVADOS-2 entry)."""
    rows = []
    ids = D['R'].index.intersection(D['C'].index)
    king, org = D['C'].kingdom.reindex(ids), D['C'].organism.reindex(ids)
    for k, g in pd.Series(ids, index=ids).groupby(king):
        i = g.index
        i3 = i.intersection(D['M'].index)
        row = {'kingdom': k, 'n': len(i), 'organisms': org[i].nunique(),
               'spearman_R_C': spearmanr(D['R'].chi[i], D['C'].chi[i])[0],
               'same_sign_R_C': (np.sign(D['R'].chi[i]) == np.sign(D['C'].chi[i])).mean(), 'n_three_model': len(i3)}
        for a, b in (('R', 'M'), ('C', 'M')):
            row[f'spearman_{a}_{b}'] = spearmanr(D[a].chi[i3], D[b].chi[i3])[0] if len(i3) > 20 else np.nan
        rows.append(row)
    return pd.DataFrame(rows)


def by_nu(D: dict) -> pd.DataFrame:
    rows = []
    for k, d in D.items():
        for lab, g in d.groupby(pd.cut(d.nu, BINS), observed=True):
            rows.append({'model': MODELS[k], 'own_nu_bin': str(lab), 'n': len(g), 'chi_median': g.chi.median(),
                         **{f'{p}_median': g[p].median() for p in PARTS}})
        rows.append({'model': MODELS[k], 'own_nu_bin': 'all', 'n': len(d), 'chi_median': d.chi.median(),
                     **{f'{p}_median': d[p].median() for p in PARTS}})
    return pd.DataFrame(rows)


def expanded(D: dict) -> pd.DataFrame:
    rows = []
    trio = D['C'].index.intersection(D['M'].index).intersection(D['R'].index)
    sel_all = D['C'].index[D['C'].nu > 0.63]
    for set_name, ids in (('CALVADOS-2 nu > 0.63, three-model set', trio.intersection(sel_all)),
                          ('CALVADOS-2 nu > 0.63, CALVADOS-2 and Mpipi-Recharged', sel_all.intersection(D['R'].index))):
        for k in ('C', 'M', 'R'):
            d = D[k].reindex(ids).dropna(subset=['chi'])
            if not len(d):
                continue
            rows.append({'set': set_name, 'model': MODELS[k], 'n': len(d), 'boltzmann_median': d.boltz.median(),
                         'dHdT_part_median': d.dQ.median(), 'net_median': d.chi.median(),
                         'frac_contract_on_heating': (d.chi < 0).mean()})
    return pd.DataFrame(rows)


def charge_held(D: dict) -> tuple:
    rows, held = [], {}
    for k, d in D.items():
        h = (d.Q > 0) & (d.Q > d.S)
        held[k] = set(d.index[h])
        for lab, g in h.groupby(pd.cut(d.nu, BINS), observed=True):
            rows.append({'model': MODELS[k], 'own_nu_bin': str(lab), 'n': len(g), 'frac_charge_held': g.mean()})
        rows.append({'model': MODELS[k], 'own_nu_bin': 'all', 'n': len(d), 'frac_charge_held': h.mean()})
    trio = D['C'].index.intersection(D['M'].index).intersection(D['R'].index)
    ov = []
    for set_name, ids in (('CALVADOS-2 and Mpipi-Recharged', D['C'].index.intersection(D['R'].index)),
                          ('three-model set', trio)):
        ids = set(ids)
        hc, hr, hm = held['C'] & ids, held['R'] & ids, held['M'] & ids
        ov.append({'set': set_name, 'n': len(ids), 'held_C': len(hc), 'held_R': len(hr),
                   'held_M': len(hm) if set_name == 'three-model set' else np.nan,
                   'held_C_and_R': len(hc & hr), 'held_all_three': len(hc & hr & hm) if set_name == 'three-model set' else np.nan})
    return pd.DataFrame(rows), pd.DataFrame(ov)


def leftovers(D: dict) -> tuple:
    trio = D['C'].index.intersection(D['M'].index).intersection(D['R'].index)
    trio = [u for u in trio if all(D[k].nu_r2.get(u, 0) >= R2_MIN for k in D)]
    full = {k: d[d.nu_r2 >= R2_MIN] for k, d in D.items()}
    corr, shares = [], []

    def resid(d):
        X = np.column_stack([spline(d.nu.to_numpy()), np.log(d.N.to_numpy())])
        out = pd.DataFrame({'L': residualise(d.chi.to_numpy(), X)}, index=d.index)
        for p in PARTS:
            out[p] = residualise(d[p].to_numpy(), X)
        return out

    def share(res, label, n):
        v = res.L.var(ddof=1)
        shares.append({'target': label, 'n': n, **{f'share_{p}': np.cov(res[p], res.L)[0, 1] / v for p in PARTS}})

    for k, d in full.items():
        share(resid(d), f'{MODELS[k]}, own full set', len(d))
    Rt = {k: resid(D[k].loc[trio]) for k in D}
    for k in D:
        share(Rt[k], f'{MODELS[k]}, three-model set', len(trio))
    z = {k: Rt[k] / Rt[k].L.std(ddof=1) for k in D}
    shared = sum(z.values()) / 3
    share(shared, 'shared leftover (three models)', len(trio))
    for a, b in (('R', 'C'), ('R', 'M'), ('C', 'M')):
        corr.append({'pair': f'{MODELS[a]} vs {MODELS[b]}', 'n': len(trio), 'r_leftover': np.corrcoef(Rt[a].L, Rt[b].L)[0, 1]})
    pairRC = D['C'].index.intersection(D['R'].index)
    pairRC = [u for u in pairRC if D['C'].nu_r2[u] >= R2_MIN and D['R'].nu_r2[u] >= R2_MIN]
    lc, lr = resid(D['C'].loc[pairRC]), resid(D['R'].loc[pairRC])
    corr.append({'pair': 'Mpipi-Recharged vs CALVADOS-2 (all paired)', 'n': len(pairRC), 'r_leftover': np.corrcoef(lc.L, lr.L)[0, 1]})
    zc, zr = lc / lc.L.std(ddof=1), lr / lr.L.std(ddof=1)
    share((zc + zr) / 2, 'shared leftover (CALVADOS-2 + Mpipi-Recharged, all paired)', len(pairRC))
    return pd.DataFrame(corr), pd.DataFrame(shares)


def main():
    os.makedirs(OUT, exist_ok=True)
    pd.set_option('display.width', 220)
    pd.set_option('display.max_columns', 30)
    D = load()
    tables = {'A_agreement': agreement(D), 'A_agreement_by_kingdom': by_kingdom(D), 'B_chi_by_nu': by_nu(D), 'C_expanded_dHdT': expanded(D)}
    tables['D_charge_held'], tables['D_charge_held_overlap'] = charge_held(D)
    tables['E_leftover_correlation'], tables['E_leftover_shares'] = leftovers(D)
    lines = [f'{MODELS[k]}: {len(d)} proteins at 300 K' for k, d in D.items()]
    for name, t in tables.items():
        t.to_csv(os.path.join(OUT, f'{name}.csv'), index=False)
        lines += ['', f'== {name}', t.to_string(index=False, float_format=lambda v: f'{v:.3f}')]
    open(os.path.join(OUT, 'summary.txt'), 'w').write('\n'.join(lines) + '\n')
    print('\n'.join(lines))


if __name__ == '__main__':
    main()
