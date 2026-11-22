"""
heldout_extrapolation.py — RETROSPECTIVE, post hoc (2026-10-20): are the misses of the frozen model in the held-out test
extrapolation errors? Reads the prepared held-out and calibration npz (heldout_analysis.py prepare / wuttke-types);
changes no frozen file.

Kernel rows. The calibrated model adds sum_j gamma_j E_j(x; T) to the energy, with ten forcing directions j: the
hydrophilic groups DE, K, R, pol (N Q S T H) and G, and the TEA classes aliph, pro, arom, arg and rest. For every series
(the Wuttke calibration constructs, lambda-repressor, and the eight held-out series) the kernel row is
    k_j = d s / d gamma_j  at gamma = 0,
with s = 100 d ln<A>/dT the same reweighted least-squares slope over the simulated temperatures as in the test
(heldout_analysis.slope without the jackknife), by central differences (+-delta_j, delta_j = 0.1 tau_j). At gamma = 0
the reweighting stays close to the plain ensemble, so the rows are reliable even where the frozen model's reweighting was
not (PNt). Rows are scaled by tau_j, the frozen strength of the term the direction belongs to (gamma_hyd for the
hydrophilic groups, gamma_tea for the TEA classes), so k~_j is the change of slope per prior standard deviation.

Extrapolation. With the four calibration rows K (4 x 10):
  outside span   d = |(I - P) k~|, P the projector onto the row space of K (in % per K), and d / |k~|
  predictive SD  sqrt(k~' S k~), S = (I + K'K / sigma^2)^-1 the posterior covariance of the scaled strengths under a
                 unit Gaussian prior after the four calibration slopes with noise sigma (0.02, 0.05, 0.1 % per K)
Miss: |s_A0 - s_exp| of the frozen model from primary.csv (seven measured series per force field; lambda-repressor from
the calibration fit, reported separately because it was excluded from the calibration as non-monotonic).
Statistic: Spearman between miss and each extrapolation measure, per force field and pooled (14 series), with an exact
permutation P within force fields.
    BENDER_REANALYSIS_DIR=... python heldout_extrapolation.py         (openmpipi env)
Output: heldout_temperature/results/exploratory/extrapolation_{rows,series}.csv, extrapolation_summary.txt
"""
import sys

sys.dont_write_bytecode = True

import itertools  # noqa: E402
import os  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

HERE = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, HERE)
import heldout_analysis as HA  # noqa: E402

OUT = os.path.join(HA.OUT, 'exploratory')
HYDG = ['DE', 'K', 'R', 'pol', 'G']
TEAC = ['aliph', 'pro', 'arom', 'arg', 'rest']
DIRS = [f'hyd_{g}' for g in HYDG] + [f'tea_{c}' for c in TEAC]
KEY = {**{f'hyd_{g}': f'E_hyd_{g}' for g in HYDG}, **{f'tea_{c}': f'E_{c}' for c in TEAC}}
SIGMAS = (0.02, 0.05, 0.1)
LAMBDA = 'lambda_P6C_I84C'


def tau(model, d):
    return HA.FROZEN[model]['hyd'] if d.startswith('hyd') else HA.FROZEN[model]['tea']


def value(z, obs, d, g, T):
    """Ensemble average of the observable with only direction d switched on at strength g."""
    u = g * z[KEY[d]] / (HA.R_GAS * T) if g else np.zeros(len(z['E_aliph']))
    w = np.exp(-(u - u.min()))
    w /= w.sum()
    if obs == 'Rg':
        return float((w * z['obs_Rg']).sum())
    if obs == 'Rh':
        return float(1.0 / (w * z['obs_inv']).sum())
    if obs == 'mean_loop':
        return float((w * z['obs_loop']).sum())
    raise ValueError(obs)


def slope(ser, obs, d, g):
    Ts = sorted(ser)
    return HA.lnslope(Ts, [value(ser[T], obs, d, g, float(T)) for T in Ts])


def kernel(model, ser, obs):
    row = {}
    for d in DIRS:
        h = 0.1 * tau(model, d)
        row[d] = (slope(ser, obs, d, h) - slope(ser, obs, d, -h)) / (2 * h) * tau(model, d)   # per prior SD
    row['s_plain'] = slope(ser, obs, DIRS[0], 0.0)
    return row


def wuttke(model, name):
    d = os.path.join(HA.OUT, f'wuttke_types_{model}')
    ser = {}
    for T in (278, 288, 300, 310, 320):
        z = np.load(os.path.join(d, f'{name}_{T}K.npz'))
        zz = {k: z[k] for k in z.files}
        zz['obs_Rg'] = zz['A']
        ser[T] = zz
    return ser


def perm_spearman(groups):
    """Spearman over pooled (x, y) pairs; exact permutation of y within each group."""
    xs = np.concatenate([g[0] for g in groups])
    ys = np.concatenate([g[1] for g in groups])
    obs = stats.spearmanr(xs, ys)[0]
    perms = [list(itertools.permutations(range(len(g[1])))) for g in groups]
    rng = np.random.default_rng(3)
    n_tot = np.prod([len(p) for p in perms])
    if n_tot <= 200000:
        combos = itertools.product(*perms)
    else:
        combos = ([p[rng.integers(len(p))] for p in perms] for _ in range(100000))
    null = []
    for c in combos:
        yy = np.concatenate([np.asarray(g[1])[list(pi)] for g, pi in zip(groups, c)])
        null.append(stats.spearmanr(xs, yy)[0])
    null = np.array(null)
    return obs, float((null >= obs - 1e-12).mean()), len(null)


def main():
    os.makedirs(OUT, exist_ok=True)
    rows, ser_rows, lines = [], [], ['RETROSPECTIVE, post hoc analysis of the held-out test (2026-10-20); no frozen file changed',
                                     'kernel rows d s / d gamma_j at gamma = 0, per prior SD (tau_j = frozen strength), % per K', '']
    for model in ('calvados', 'recharged'):
        prim = pd.read_csv(os.path.join(HA.OUT, f'primary_{model}', 'primary.csv'))
        a0 = prim[prim.variant == 'A0'].set_index('set')
        cal = {}
        for n in HA.FOUR + [LAMBDA]:
            r = kernel(model, wuttke(model, n), 'Rg')
            r.update(model=HA.LABEL[model], series=n, role='calibration' if n in HA.FOUR else 'lambda (not calibrated)')
            rows.append(r)
            cal[n] = r
        Kc = np.array([[cal[n][d] for d in DIRS] for n in HA.FOUR])
        U, s, Vt = np.linalg.svd(Kc, full_matrices=False)
        rank = int((s > 1e-9 * s[0]).sum())
        P = Vt[:rank].T @ Vt[:rank]
        lines.append(f'{HA.LABEL[model]}: calibration rows rank {rank} of {len(DIRS)}; singular values ' + ', '.join(f'{v:.3g}' for v in s))
        fit = pd.read_csv(f'{HA.WU}/calibration/{model}/fit_hyd_final_constructs.csv').set_index('construct')
        targets = [(st, obs) for st, obs in HA.OBS_OF.items()] + [(LAMBDA, 'Rg')]
        for st, obs in targets:
            if st == LAMBDA:
                k = np.array([cal[LAMBDA][d] for d in DIRS])
                miss = abs(fit.loc[LAMBDA, 'slope_fit'] - fit.loc[LAMBDA, 'slope_exp'])
                s_a0, s_exp, ess_min = fit.loc[LAMBDA, 'slope_fit'], fit.loc[LAMBDA, 'slope_exp'], np.nan
            else:
                ser = HA.series(model, st)
                r = kernel(model, ser, obs)
                r.update(model=HA.LABEL[model], series=st, role='held out')
                rows.append(r)
                k = np.array([r[d] for d in DIRS])
                s_a0, s_exp, ess_min = a0.loc[st, 's_model'], a0.loc[st, 's_exp'], a0.loc[st, 'min_ess']
                miss = abs(s_a0 - s_exp) if np.isfinite(s_exp) else np.nan
            out = k - P @ k
            s_pl = cal[LAMBDA]['s_plain'] if st == LAMBDA else r['s_plain']
            row = {'model': HA.LABEL[model], 'series': st, 'observable': obs, 's_exp': s_exp, 's_A0': s_a0, 'miss': miss,
                   's_plain': s_pl, 's_linear': s_pl + k.sum(),       # check: linear response at the frozen strengths

                   'A0_min_ess': ess_min, 'norm_k': float(np.linalg.norm(k)), 'outside_span': float(np.linalg.norm(out)),
                   'outside_frac': float(np.linalg.norm(out) / np.linalg.norm(k))}
            for sg in SIGMAS:
                S = np.linalg.inv(np.eye(len(DIRS)) + Kc.T @ Kc / sg ** 2)
                row[f'pred_sd_{sg:g}'] = float(np.sqrt(k @ S @ k))
            # direction carrying most of the part outside the span, and the part itself by direction
            row['main_outside_direction'] = DIRS[int(np.argmax(np.abs(out)))]
            row.update({f'out_{d}': float(v) for d, v in zip(DIRS, out)})
            ser_rows.append(row)
    R, Sx = pd.DataFrame(rows), pd.DataFrame(ser_rows)
    R.to_csv(os.path.join(OUT, 'extrapolation_rows.csv'), index=False)
    Sx.to_csv(os.path.join(OUT, 'extrapolation_series.csv'), index=False)
    pd.set_option('display.width', 250)
    lines += ['', 'kernel rows (% per K per prior SD):', R[['model', 'series', 'role'] + DIRS].round(3).to_string(index=False), '',
              'check: s_linear = s_plain + sum of the row (all scaled strengths 1) against the frozen model s_A0',
              Sx[['model', 'series', 's_plain', 's_linear', 's_A0', 'A0_min_ess']].round(3).to_string(index=False), '',
              Sx[['model', 'series', 's_exp', 's_A0', 'miss', 'A0_min_ess', 'norm_k', 'outside_span', 'outside_frac',
                  'pred_sd_0.05', 'main_outside_direction']].round(3).to_string(index=False), '']
    H = Sx[(Sx.series != LAMBDA) & Sx.miss.notna()]
    for meas in ['outside_span', 'outside_frac'] + [f'pred_sd_{sg:g}' for sg in SIGMAS]:
        per = []
        for m, g in H.groupby('model'):
            rho, p, n = perm_spearman([(g[meas].to_numpy(), g.miss.to_numpy())])
            per.append(f'{m} rho {rho:.2f} (P {p:.3f}, n {len(g)})')
        rho, p, n = perm_spearman([(g[meas].to_numpy(), g.miss.to_numpy()) for _, g in H.groupby('model')])
        lines.append(f'miss vs {meas}: pooled rho {rho:.2f} (one-sided permutation P {p:.3f} within force fields, {n} permutations); ' + '; '.join(per))
    rel = Sx[Sx.series == LAMBDA][['model', 'miss', 'outside_span', 'pred_sd_0.05']]
    lines += ['', 'lambda-repressor (excluded from the calibration):', rel.round(3).to_string(index=False)]
    open(os.path.join(OUT, 'extrapolation_summary.txt'), 'w').write('\n'.join(lines) + '\n')
    print('\n'.join(lines))


if __name__ == '__main__':
    main()
