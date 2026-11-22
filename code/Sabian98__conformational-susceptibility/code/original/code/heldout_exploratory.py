"""
heldout_exploratory.py — POST HOC, exploratory analyses of the held-out temperature test (not part of the pre-registered
analysis; heldout_analysis.py stays frozen and is only imported). Mpipi-Recharged only, while CALVADOS-2 is still blind.

gp   Does the over-contraction of the Gly/Pro-rich chains come from the forcing on Gly and Pro? Variants of the frozen
     model, each with gamma_hyd refitted on the four Wuttke constructs only (tea_calibrate-style residuals, as
     heldout_analysis.unit_resid), then applied to every held-out set:
       A0       frozen (hydrophilic set with Gly, TEA on all classes incl. Pro)
       noG      hydrophilic set without Gly, TEA on all classes
       noGP     hydrophilic set without Gly, TEA without the Pro class
     The TEA strength stays at the calibration prior (0.4). Slopes, jackknife SEs and effective samples use
     heldout_analysis's own estimators; only the bias energy is replaced for the noGP variant.
chph Literature-motivated form (Wuttke 2014 for charged groups, Baxa 2024 for hydrophobic side chains and an inert
     backbone): hydrophilic forcing on D, E, K, R only (one strength, groups g_DE, g_K, g_R), none on Gly or N, Q, S, T, H,
     plus the TEA (hydrophobic) term.
       scan  gamma_TEA = 0, 0.4, 1, 2, 4; for each, gamma_chg refitted on the four Wuttke constructs only
       fit   gamma_chg and gamma_TEA both fitted on the registered training split (Wuttke four, PNt-I025, PNt-I175,
             sfAFP-I175 end-point ratios; heldout_analysis.unit_resid, weak priors as A3), then the test sets
    BENDER_REANALYSIS_DIR=... python heldout_exploratory.py gp|chph       (openmpipi env)
Output: heldout_temperature/results/exploratory/{gp,chph}.csv, {gp,chph}_summary.txt
"""
import sys

sys.dont_write_bytecode = True

import os  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.optimize import least_squares  # noqa: E402

HERE = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, HERE)
import heldout_analysis as H  # noqa: E402

STAGE = os.path.join(H.OUT, 'primary_recharged')          # holds npz_recharged and wuttke_types_recharged only
OUT = os.path.join(H.OUT, 'exploratory')
MODEL = 'recharged'
TEA = H.FROZEN[MODEL]['tea']
_frozen_bias = H.bias


def bias(z, gam, T):
    """As heldout_analysis.bias, plus the key 'tea_noP' (TEA on every class except Pro)."""
    u = _frozen_bias(z, {k: v for k, v in gam.items() if k != 'tea_noP'}, T)
    if gam.get('tea_noP'):
        u = u + gam['tea_noP'] * sum(z[f'E_{c}'] for c in H.CLS if c != 'pro') / (H.R_GAS * T)
    return u


def gam_of(variant, g):
    return {'A0': {'tea': TEA, 'hyd_all': g}, 'noG': {'tea': TEA, 'hyd_noG': g},
            'noGP': {'tea_noP': TEA, 'hyd_noG': g}}[variant]


def fit_wuttke(variant, units):
    def resid(p):
        gam = gam_of(variant, p[0])
        r = [H.unit_resid(MODEL, u, 'wuttke', d, gam, None) for u, d in units.items()]
        return np.concatenate(r + [np.array([p[0] / 10.0])])           # the calibration's weak prior on gamma_hyd
    sol = least_squares(resid, x0=[H.FROZEN[MODEL]['hyd']], method='trf')
    rms = float(np.sqrt(np.mean(np.concatenate([H.unit_resid(MODEL, u, 'wuttke', d, gam_of(variant, sol.x[0]), None)
                                                for u, d in units.items()]) ** 2)))
    return float(sol.x[0]), rms


def gp():
    H.OUT = STAGE
    H.bias = bias
    os.makedirs(OUT, exist_ok=True)
    s_exp, _ = H.experimental()
    units = H.wuttke_units(MODEL)
    rows, lines = [], ['EXPLORATORY, post hoc (Mpipi-Recharged): forcing on Gly and Pro switched off',
                       'gamma_hyd refitted on the four Wuttke constructs only; held-out slopes in % per K', '']
    fits = {}
    for v in ('A0', 'noG', 'noGP'):
        g, rms = fit_wuttke(v, units)
        fits[v] = g
        lines.append(f'{v:5s} gamma_hyd {g:.3f} (Wuttke fit, rms standardized residual {rms:.2f})')
    lines.append('')
    for st, obs in H.OBS_OF.items():
        ser = H.series(MODEL, st)
        row = {'set': st, 'observable': obs, 's_exp': s_exp.get(st, np.nan)}
        for v, g in fits.items():
            gam = gam_of(v, g)
            s, se, Ts = H.slope(ser, obs, gam)
            row[v], row[f'{v}_se'] = s, se
            row[f'{v}_min_ess'] = min(H.ess(ser[T], gam, float(T)) for T in Ts)
        rows.append(row)
    t = pd.DataFrame(rows).set_index('set')
    t.to_csv(os.path.join(OUT, 'gp.csv'))
    show = t[['s_exp', 'A0', 'noG', 'noGP', 'noGP_se', 'noGP_min_ess']].round(3)
    lines.append(show.to_string())
    ctl = [s for s in H.CONTROLS]
    con = [s for s in H.CONTRACT]
    for v in fits:
        miss = {s: abs(t.loc[s, v] - t.loc[s, 's_exp']) for s in ctl}
        rat = {s: t.loc[s, v] / t.loc[s, 's_exp'] for s in con}
        lines.append(f'{v:5s} controls |model - exp| ' + ', '.join(f'{s} {m:.2f}' for s, m in miss.items())
                     + ' | ratios ' + ', '.join(f'{s} {r:.2f}' for s, r in rat.items()))
    open(os.path.join(OUT, 'gp_summary.txt'), 'w').write('\n'.join(lines) + '\n')
    print('\n'.join(lines))


def chg_gam(g_chg, g_tea):
    return {'tea': g_tea, 'g_DE': g_chg, 'g_K': g_chg, 'g_R': g_chg}


def fit_free(make, x0, units, ratio_exp, prior, model=MODEL):
    """Least squares over the given units for the parameters of make(*p); prior(p) gives extra residuals."""
    def resid(p):
        gam = make(*p)
        r = [H.unit_resid(model, u, kind, d, gam, ratio_exp) for u, (kind, d) in units.items()]
        return np.concatenate(r + [prior(p)])
    sol = least_squares(resid, x0=x0, method='trf', bounds=(0, np.inf))
    gam = make(*sol.x)
    rms = float(np.sqrt(np.mean(np.concatenate([H.unit_resid(model, u, kind, d, gam, ratio_exp)
                                                for u, (kind, d) in units.items()]) ** 2)))
    return sol.x, rms


def heldout_row(gam):
    s_exp, _ = H.experimental()
    out = {}
    for st, obs in H.OBS_OF.items():
        ser = H.series(MODEL, st)
        s, se, Ts = H.slope(ser, obs, gam)
        out[st] = (s, se, min(H.ess(ser[T], gam, float(T)) for T in Ts), s_exp.get(st, np.nan))
    return out


def chph():
    H.OUT = STAGE
    os.makedirs(OUT, exist_ok=True)
    s_exp, ratio_exp = H.experimental()
    wu = {n: ('wuttke', d) for n, d in H.wuttke_units(MODEL).items()}
    lines = ['EXPLORATORY, post hoc (Mpipi-Recharged): hydrophilic forcing on D, E, K, R only (no Gly, no N Q S T H)',
             'plus the TEA hydrophobic term; held-out slopes in % per K (exp in brackets), [min effective sample]', '']
    rows = []
    order = list(H.OBS_OF)

    def report(label, gam, extra):
        hr = heldout_row(gam)
        rows.append({'variant': label, **extra, **{f'{st}': hr[st][0] for st in order},
                     **{f'{st}_se': hr[st][1] for st in order}, **{f'{st}_ess': hr[st][2] for st in order}})
        lines.append(f'{label}: ' + ', '.join(f'{k} {v:.3f}' if isinstance(v, float) else f'{k} {v}' for k, v in extra.items()))
        for st in order:
            s, se, e, x = hr[st]
            lines.append(f'    {st:16s} {s:+.3f} +- {se:.3f} [{e:.3f}]  ({x:+.3f})')
        con = [hr[s][0] / hr[s][3] for s in H.CONTRACT]
        ctl = [abs(hr[s][0] - hr[s][3]) for s in H.CONTROLS]
        test = ['asyn4W26C-I030', 'asyn39W66C-I030', 'p53N93-I118']
        lines.append('    ratios ' + ', '.join(f'{s} {r:.2f}' for s, r in zip(H.CONTRACT, con))
                     + ' | controls |model-exp| ' + ', '.join(f'{s} {c:.2f}' for s, c in zip(H.CONTROLS, ctl))
                     + ' | test-set geo-mean ratio %.2f' % np.exp(np.mean(np.log([max(1e-6, hr[s][0] / hr[s][3]) for s in test]))))
        lines.append('')

    lines.append('--- scan: gamma_TEA fixed, gamma_chg refitted on the Wuttke four only')
    for gt in (0.0, 0.4, 1.0, 2.0, 4.0):
        (gc,), rms = fit_free(lambda c, gt=gt: chg_gam(c, gt), [0.3], wu, ratio_exp, lambda p: np.array([p[0] / 10.0]))
        report(f'scan gamma_TEA {gt:.1f}', chg_gam(gc, gt), {'gamma_chg': float(gc), 'gamma_TEA': gt, 'wuttke_rms': rms})

    lines.append('--- fit: gamma_chg and gamma_TEA on the registered training split (Wuttke four, PNt-I025, PNt-I175, sfAFP)')
    units = dict(wu)
    for st in H.TRAIN_BAXA:
        units[st] = ('baxa', H.series(MODEL, st))
    (gc, gt), rms = fit_free(chg_gam, [0.3, TEA], units, ratio_exp,
                             lambda p: np.array([p[0] / 10.0, (p[1] - TEA) / 1.0]))
    report('fit on training split', chg_gam(gc, gt), {'gamma_chg': float(gc), 'gamma_TEA': float(gt), 'train_rms': rms})

    pd.DataFrame(rows).to_csv(os.path.join(OUT, 'chph.csv'), index=False)
    open(os.path.join(OUT, 'chph_summary.txt'), 'w').write('\n'.join(lines) + '\n')
    print('\n'.join(lines))


# ----------------------------------------------------------------------------- Amendment 3, secondary analysis S1
S1_TEA = {'calvados': 2.0, 'recharged': 0.4}                 # calibration priors of the TEA term (frozen model)
S1_GAMMA = {'calvados': 0.469789, 'recharged': 0.294234}     # gamma_chg fitted on the Wuttke four (s1-fit), 2026-10-11


def s1_fit(model):
    """gamma_chg of S1, fitted on the four Wuttke constructs only (no held-out data are read)."""
    H.OUT = os.path.dirname(STAGE)                            # results/: wuttke_types_<model> live there
    wu = {n: ('wuttke', d) for n, d in H.wuttke_units(model).items()}
    (gc,), rms = fit_free(lambda c: chg_gam(c, S1_TEA[model]), [0.3], wu, None,
                          lambda p: np.array([p[0] / 10.0]), model=model)
    print(f'S1 {LABEL_OF[model]}: gamma_chg {gc:.6f} (gamma_TEA {S1_TEA[model]}), Wuttke rms standardized residual {rms:.3f}')
    return float(gc)


LABEL_OF = {'calvados': 'CALVADOS-2', 'recharged': 'Mpipi-Recharged'}


def s1(model, stage):
    """Pre-registered secondary analysis S1 (Amendment 3) for one force field, from a staging folder holding only
    npz_<model> (and wuttke_types_<model>)."""
    gc = s1_fit(model)
    if S1_GAMMA[model] is None or abs(gc - S1_GAMMA[model]) > 1e-6:
        raise ValueError(f'gamma_chg {gc} differs from the registered value {S1_GAMMA[model]}')
    H.OUT = stage
    s_exp, _ = H.experimental()
    gam = chg_gam(gc, S1_TEA[model])
    rows = []
    for st, obs in H.OBS_OF.items():
        ser = H.series(model, st)
        if len(ser) < 2:
            continue
        for lab, g in (('A0', H.frozen_gam(model)), ('S1', gam)):
            s, se, Ts = H.slope(ser, obs, g)
            rows.append({'model': LABEL_OF[model], 'variant': lab, 'set': st, 'observable': obs, 's_model': s, 'se': se,
                         's_exp': s_exp.get(st, np.nan), 'min_ess': min(H.ess(ser[T], g, float(T)) for T in Ts)})
    t = pd.DataFrame(rows)
    lines = [f'Amendment 3, secondary analysis S1, {LABEL_OF[model]}: hydrophilic forcing on D E K R only '
             f'(gamma_chg {gc:.4f}), TEA {S1_TEA[model]}; s in % per K', '']
    for lab in ('A0', 'S1'):
        tt = t[t.variant == lab].set_index('set')
        c = H.criteria(tt)
        m = tt.dropna(subset=['s_exp'])
        mae = float((m.s_model - m.s_exp).abs().mean())
        r = float(np.corrcoef(m.s_model, m.s_exp)[0, 1])
        sf = abs(tt.loc['sfAFP-I175', 's_model'] - tt.loc['sfAFP-I175', 's_exp']) if 'sfAFP-I175' in tt.index else np.nan
        lines += [f'{lab}: ' + ', '.join(f'{k} {v}' for k, v in c.items()) + f', MAE {mae:.3f}, r {r:+.2f}, sfAFP |d| {sf:.3f}',
                  tt[['observable', 's_model', 'se', 's_exp', 'min_ess']].round(3).to_string(), '']
    os.makedirs(stage, exist_ok=True)
    t.to_csv(os.path.join(stage, f's1_{model}.csv'), index=False)
    open(os.path.join(stage, f's1_{model}_summary.txt'), 'w').write('\n'.join(lines) + '\n')
    print('\n'.join(lines))


if __name__ == '__main__':
    if sys.argv[1] == 's1-fit':
        s1_fit(sys.argv[2])
    elif sys.argv[1] == 's1':
        s1(sys.argv[2], sys.argv[3])
    else:
        {'gp': gp, 'chph': chph}[sys.argv[1]]()
