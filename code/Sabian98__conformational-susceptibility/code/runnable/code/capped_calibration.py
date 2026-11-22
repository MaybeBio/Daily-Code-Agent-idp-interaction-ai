"""
capped_calibration.py — the Wuttke calibration of the capped hydration forcing (the main model from 2026-10-24) in the
form the manuscript reports for the linear calibration (tea_calibrate.py, Fig. 3): leave-one-construct-out slopes and
the model curves of all five constructs.

  prepare  per-residue contact numbers of the lambda-repressor runs (not a calibration construct, so not cached by
           saturating_forcing.py), cached as saturating_forcing does; check: the linear energies from the contact numbers
           equal the stored class and hydrophilic energies of the calibration npz to 1e-4 relative
  loo      for each of the four fitted constructs, the capped hydrophilic strength refitted on the other three
           (capped_model.fit's residuals and prior), and the left-out construct's slope at 300 K from its model curve,
           100 x the linear coefficient of a quadratic in T - 300 over the measured range up to 320 K, exactly as
           tea_calibrate.report; plain-model slopes alike (both strengths zero). Check: the same path with the linear A0
           form reproduces the stored leave-one-out slopes of tea_calibrate (cv_hyd_T320_four_constructs.csv) to 1e-4
  curves   relative <Rg>(T)/<Rg>(300) for the plain and capped models of all five constructs, each fitted construct with
           its leave-one-out strength and lambda-repressor with the full fit (capped_model.params), with the Kish
           effective sample fraction of the capped reweighting at each temperature (the plain curve needs none: 1)
    BENDER_REANALYSIS_DIR=... python capped_calibration.py      (openmpipi env; prepare runs only where the cache is missing)
Output: heldout_temperature/capped_forcing/calibration/{loo.csv, curves.csv, summary.txt}
"""
import sys

sys.dont_write_bytecode = True

import os  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

HERE = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, os.environ.get('BENDER_REANALYSIS_DIR', HERE))
sys.path.insert(0, HERE)
import capped_model as CM  # noqa: E402
import saturating_forcing as SF  # noqa: E402

HA, TC = SF.HA, SF.TCal
H = SF.HX.H
H.OUT = SF.HA.OUT
OUT = os.path.join(CM.BM, 'heldout_temperature', 'capped_forcing', 'calibration')
LAMBDA = 'lambda_P6C_I84C'


def prepare_lambda(model):
    msgs = []
    for T in SF.WT_TEMPS:
        out = os.path.join(SF.cache_dir(model), f'wuttke_{LAMBDA}_{T}K.npz')
        if os.path.exists(out):
            continue
        d = os.path.join(SF.WUTTKE_RUNS[model], f'{T}K', LAMBDA)
        if model == 'calvados':
            xyz, seq, _ = SF.bc.load_calvados(os.path.join(d, f'{LAMBDA}.dcd'), os.path.join(d, 'top.pdb'), k=SF.L.K_FRAMES)
        else:
            import mpipi_recharged_fdt as R
            xyz, seq, _ = R.load_run(d, LAMBDA)
        ref = np.load(os.path.join(HA.WU, 'calibration', model, f'{LAMBDA}_{T}K.npz'))
        c = SF.contacts(model, xyz, seq)
        wt, wh = SF.weights(model, seq, float(T), TC.HYD_RES)
        et, eh = SF.energy(model, c, wt, np.inf), SF.energy(model, c, wh, np.inf)
        st_t = sum(ref[f'E_{k}'] for k in HA.CLS)
        scale = max(np.abs(st_t).max(), np.abs(ref['E_hyd']).max(), 1e-9)
        dev = max(np.abs(et - st_t).max(), np.abs(eh - ref['E_hyd']).max()) / scale
        if dev > 1e-4:
            raise ValueError(f'{model} {LAMBDA} {T} K: linear energies from contacts differ from the stored ones ({dev:.1e})')
        np.savez_compressed(out, c=c, seq=seq, T=T)
        msgs.append(f'{model} {LAMBDA} {T} K: cached, energy check {dev:.1e}')
    return msgs


def lambda_series(model):
    ser = {}
    for T in SF.WT_TEMPS:
        z = dict(np.load(os.path.join(HA.WU, 'calibration', model, f'{LAMBDA}_{T}K.npz')))
        z['obs_Rg'] = z['A']
        ser[T] = SF.with_energies(model, z, *SF.cached(model, 'wuttke', LAMBDA, T), T, 'S1', CM.CSAT)
    return ser


def curve(data, gam):
    Ts = np.array(sorted(data))
    m = {T: H.ensemble_value(data[T], 'Rg', gam, float(T)) for T in Ts}
    return Ts, np.array([m[T] / m[300] for T in Ts])


def ess(data, gam):
    return np.array([H.ess(data[T], gam, float(T)) for T in sorted(data)])


def slope(unit, Ts, y):
    re_, _, inside = H.wuttke_target(unit, tuple(Ts))
    f = lambda v: 100 * np.polyfit(Ts[inside] - 300, v[inside], 2)[1]  # noqa: E731
    return f(re_), f(y)


def fit_on(model, units):
    from scipy.optimize import least_squares

    def resid(p):
        gam = SF.gam_of(model, p[0])
        r = [H.unit_resid(model, u, kind, d, gam, None) for u, (kind, d) in units.items()]
        return np.concatenate(r + [np.array([p[0] / 10.0])])
    return float(least_squares(resid, x0=[0.3 if model == 'recharged' else 1.0], method='trf', bounds=(0, np.inf)).x[0])


def main():
    os.makedirs(OUT, exist_ok=True)
    L = ['Capped hydration forcing (c_sat 2, D E K R), Wuttke calibration: leave-one-construct-out slopes at 300 K', '']
    loo, curves = [], []
    for k, model in CM.KEY.items():
        L += prepare_lambda(model)
        units = SF.wuttke_units(model, 'S1', CM.CSAT)
        g_full = CM.params()[0][k][0]
        if abs(fit_on(model, units) - g_full) > 1e-3 * g_full:
            raise ValueError(f'{model}: refit on all four constructs differs from capped_model_params')
        plain = {'g_tsat': 0.0, 'g_hsat': 0.0}
        a0 = SF.wuttke_units(model, 'A0', np.inf)
        st = pd.read_csv(f'{HA.WU}/calibration/{model}/cv_hyd_T320_four_constructs.csv').set_index('construct')
        for left in a0:
            g = fit_on(model, {u: v for u, v in a0.items() if u != left})
            Ts, y = curve(a0[left][1], SF.gam_of(model, g))
            if abs(slope(left, Ts, y)[1] - st.loc[left, 'slope_fit']) > 1e-4:
                raise ValueError(f'{model} {left}: the A0 path does not reproduce the stored leave-one-out slope')
        L.append(f'{model}: check passed, the A0 form through this path reproduces the stored leave-one-out slopes')
        for left in units:
            g = fit_on(model, {u: v for u, v in units.items() if u != left})
            data = units[left][1]
            Ts, yc = curve(data, SF.gam_of(model, g))
            _, yp = curve(data, plain)
            se, sc = slope(left, Ts, yc)
            _, sp = slope(left, Ts, yp)
            loo.append({'model': model, 'construct': left, 'gamma_loo': g, 'slope_exp': se, 'slope_fit': sc, 'slope_plain': sp})
            curves += [{'model': model, 'construct': left, 'curve': lab, 'T': T, 'rel': v, 'ess': e}
                       for lab, y, ee in (('plain', yp, np.ones(len(Ts))), ('capped', yc, ess(data, SF.gam_of(model, g))))
                       for T, v, e in zip(Ts, y, ee)]
        lam = lambda_series(model)
        Ts, yc = curve(lam, SF.gam_of(model, g_full))
        _, yp = curve(lam, plain)
        curves += [{'model': model, 'construct': LAMBDA, 'curve': lab, 'T': T, 'rel': v, 'ess': e}
                   for lab, y, ee in (('plain', yp, np.ones(len(Ts))), ('capped', yc, ess(lam, SF.gam_of(model, g_full))))
                   for T, v, e in zip(Ts, y, ee)]
    T_ = pd.DataFrame(loo)
    T_.to_csv(os.path.join(OUT, 'loo.csv'), index=False)
    pd.DataFrame(curves).to_csv(os.path.join(OUT, 'curves.csv'), index=False)
    for m, g in T_.groupby('model', sort=False):
        L.append(f'== {m}: strengths left one out {", ".join(f"{x:.3f}" for x in g.gamma_loo)} (full fit {CM.params()[0]["C" if m == "calvados" else "R"][0]:.3f})')
        for r in g.itertuples():
            L.append(f'  {r.construct:16s} measured {r.slope_exp:+.3f}, capped (left out) {r.slope_fit:+.3f}, plain {r.slope_plain:+.3f}')
        L.append(f'  median |error|: capped {np.median(np.abs(g.slope_fit - g.slope_exp)):.3f}, plain '
                 f'{np.median(np.abs(g.slope_plain - g.slope_exp)):.3f} % per K')
    open(os.path.join(OUT, 'summary.txt'), 'w').write('\n'.join(L) + '\n')
    print('\n'.join(L))


if __name__ == '__main__':
    main()
