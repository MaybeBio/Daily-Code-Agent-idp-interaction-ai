"""
noise_floor.py — how close the capped model is to the error that measurement and sampling noise alone would leave.
Exploratory (2026-10-24). For each held-out temperature series with a measured slope:
  measured SE   PNt and sfAFP: from the end-point Rg and their ranges in Baxa et al. 2024 (heldout_analysis.experimental);
                alpha-synuclein: least-squares slope of ln(peak distance) on T with each point's 66% interval taken as
                +-1 SD (Woodard et al. 2018 Table 1, make_heldout_inputs.WOODARD); Aap and p53: not reported (DLS fits
                shown only as plotted slopes), so 0 (a lower bound on their noise)
  simulated SE  10-block jackknife of the capped model's slope (held-out runs reweighted, heldout_analysis.slope; PNt from
                the direct calibrated runs reweighted to the capped model, saturating_forcing.pnt_direct)
A perfect model would miss each series by |N(0, sigma)| with sigma^2 = SE_meas^2 + SE_sim^2, on average sqrt(2/pi) sigma.
Reported per series: miss, sigma, z = miss / sigma, and the smallest Kish effective sample fraction over temperatures; per force field: observed mean |miss| against the expected noise floor,
with and without the two controls and with and without PNt (whose cap value was fitted).
Caveat stated, not corrected: for alpha-synuclein the model observable is the mean loop distance and the measured one the
peak of the distance distribution (Amendment 1), a possible systematic difference that is not noise.
    BENDER_REANALYSIS_DIR=... python noise_floor.py      (openmpipi env)
Output: heldout_temperature/capped_forcing/noise_floor.{csv,txt}
"""
import sys

sys.dont_write_bytecode = True

import math  # noqa: E402
import os  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

HERE = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, os.environ.get('BENDER_REANALYSIS_DIR', HERE))
sys.path.insert(0, HERE)
import capped_model as CM  # noqa: E402
import heldout_analysis as HA  # noqa: E402
import make_heldout_inputs as MH  # noqa: E402
import saturating_forcing as SF  # noqa: E402

OUT = f'{CM.BM}/heldout_temperature/capped_forcing'
SERIES = ['PNt-I025', 'PNt-I175', 'asyn4W26C-I030', 'asyn39W66C-I030', 'p53N93-I118', 'sfAFP-I175', 'AapPGR-I195']
CONTROLS = ['sfAFP-I175', 'AapPGR-I195']
FITTED = ['PNt-I025', 'PNt-I175']


def measured_se():
    _, ratio = HA.experimental()
    se = {k: 100 * v[1] / 48 for k, v in ratio.items()}           # ln end-point ratio over 48 K
    rows = [l.split() for l in MH.WOODARD.splitlines()]
    for loop, key in (('4-26', 'asyn4W26C-I030'), ('39-66', 'asyn39W66C-I030')):
        T = np.array([float(r[1]) for r in rows if r[0] == loop])
        lo = np.array([float(r[3]) for r in rows if r[0] == loop])
        hi = np.array([float(r[4]) for r in rows if r[0] == loop])
        sd = (np.log(hi) - np.log(lo)) / 2                         # 66% interval as +-1 SD of ln(distance)
        dT = T - T.mean()
        se[key] = 100 * math.sqrt((dT ** 2 * sd ** 2).sum()) / (dT ** 2).sum()
    se['p53N93-I118'] = 0.0
    se['AapPGR-I195'] = 0.0
    return se


def main():
    s_exp, _ = HA.experimental()
    mse = measured_se()
    rows = []
    for model in ('calvados', 'recharged'):
        g_h = CM.params()[0]['C' if model == 'calvados' else 'R'][0]
        gam = SF.gam_of(model, g_h)
        SF.load_direct_rg(model)
        for st in SERIES:
            if st.startswith('PNt'):
                s, se, me = SF.pnt_direct(model, st, 'S1', CM.CSAT, g_h)
            else:
                ser = SF.heldout_series(model, st, 'S1', CM.CSAT)
                s, se, _ = HA.slope(ser, HA.OBS_OF[st], gam)
                me = min(HA.ess(ser[T], gam, float(T)) for T in ser)
            sig = math.hypot(mse[st], se)
            rows.append({'model': HA.LABEL[model], 'series': st, 's_model': s, 'se_model': se, 's_exp': s_exp[st],
                         'se_exp': mse[st], 'miss': abs(s - s_exp[st]), 'sigma': sig,
                         'z': abs(s - s_exp[st]) / sig if sig > 0 else np.nan,
                         'expected_noise_miss': math.sqrt(2 / math.pi) * sig, 'min_ess': me})
    df = pd.DataFrame(rows)
    os.makedirs(OUT, exist_ok=True)
    df.to_csv(os.path.join(OUT, 'noise_floor.csv'), index=False)
    L = ['Capped model against the noise floor (slopes in % per K); expected miss of a perfect model = sqrt(2/pi) sigma,',
         'sigma^2 = SE_meas^2 + SE_sim^2 (Aap, p53 measured SE not reported: 0, a lower bound)', '']
    for m, d in df.groupby('model', sort=False):
        L.append(f'== {m}')
        for r in d.itertuples():
            L.append(f'  {r.series:16s} model {r.s_model:+.3f} +- {r.se_model:.3f} | measured {r.s_exp:+.3f} +- {r.se_exp:.3f} | '
                     f'miss {r.miss:.3f}, sigma {r.sigma:.3f}, z {r.z:.1f}')
        for lab, dd in (('all 7', d), ('without controls', d[~d.series.isin(CONTROLS)]),
                        ('not fitted (without PNt)', d[~d.series.isin(FITTED)]),
                        ('not fitted, without controls', d[~d.series.isin(FITTED + CONTROLS)])):
            L.append(f'  {lab:30s}: observed mean miss {dd.miss.mean():.3f}, noise floor {dd.expected_noise_miss.mean():.3f}, '
                     f'ratio {dd.miss.mean() / dd.expected_noise_miss.mean():.1f}')
        L.append('')
    L.append('caveat: alpha-synuclein model observable is the mean loop distance, the measured one the peak (Amendment 1)')
    open(os.path.join(OUT, 'noise_floor.txt'), 'w').write('\n'.join(L) + '\n')
    print('\n'.join(L))


if __name__ == '__main__':
    main()
