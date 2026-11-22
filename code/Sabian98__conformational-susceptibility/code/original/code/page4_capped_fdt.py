"""
page4_capped_fdt.py — the direct PAGE4 runs (page4_direct.py, 150 mM) under the capped hydration forcing (the main model
from 2026-10-24), in place of the linear A0 forcing of page4_direct.summary: the 300 K heating response of every form
and the heating curves of every form with a temperature series.

Per run, the stored frames and contact numbers of page4_capped_predictions.run_frames (Rg and the linear energies are
checked there against the stored run).
  300 K    the forcing per frame is the temperature derivative of the hydration energy, D = dE/dT = (E(301) - E(299)) / 2
           (exact: the per-residue weights are linear in T):
             capped  D = g_tea dE_tea^cap/dT + g_h dE_chg^cap/dT   (c_sat 2, strengths page4_capped_predictions.capped_gamma)
             check   the linear A0 forcing through the same path must equal the stored D300 to 1e-5 relative (storage
                     precision; 2e-6 seen)
           total = -Cov(Rg, G_T + D / RT) / <Rg> (% per K), 10-block jackknife SE (page4_direct.fdt_total)
  curves   each temperature reweighted to the capped model, w = exp(-(g_tea E_tea^cap(T) + g_h E_chg^cap(T)) / RT);
           <Rg>(T) / <Rg>(300) and the slope at 300 K from a quadratic in T - 300 (page4_direct.curve_slope), SEs by the
           10-block jackknife with the same block left out at every temperature; the smallest Kish effective sample
           fraction over temperatures. Check: the linear A0 curve through the same path reproduces the stored calibrated
           slope (page4_direct curves.csv) to 1e-6 % per K
Forms: every form with a stored 300 K run (page4_direct npz) of each model; curves for those with four or more
temperatures including 300 K.
    BENDER_REANALYSIS_DIR=... python page4_capped_fdt.py      (openmpipi env)
Output: heldout_temperature/capped_forcing/page4/{fdt_300K.csv, curves_capped.csv, fdt_300K.txt}
"""
import sys

sys.dont_write_bytecode = True

import os  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

HERE = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, os.environ.get('BENDER_REANALYSIS_DIR', HERE))
sys.path.insert(0, HERE)
import page4_capped_predictions as PC  # noqa: E402

PD, SF = PC.PD, PC.SF
TEMPS = (278, 288, 300, 310, 320)


def dEdT(model, z, form_hyd, csat):
    up, dn = PC.energies(model, z, 301.0, form_hyd, csat), PC.energies(model, z, 299.0, form_hyd, csat)
    return (up[0] - dn[0]) / 2.0, (up[1] - dn[1]) / 2.0


def stored(model, tag):
    out = {}
    for f in os.listdir(os.path.join(PD.out_dir(model), 'npz')):
        if f.startswith('PAGE4-') and f.endswith('K.npz') and f'{tag}_' in f:
            form, T = f[len('PAGE4-'):-len('K.npz')].split(f'{tag}_')
            out.setdefault(form, []).append(int(T))
    return {k: sorted(v) for k, v in out.items()}


def curve(data, us, temps, keep=None):
    m = {T: PC.value(data[T], 'Rg', us[T], keep)[0] for T in temps}
    rel = np.array([m[T] / m[300] for T in temps])
    return rel, 100 * np.polyfit(np.array(temps) - 300.0, rel, 2)[1]


def curve_se(data, us, temps):
    rel, s = curve(data, us, temps)
    n = len(data[300]['obs_Rg'])
    jk = [curve(data, us, temps, PC.FS.keep_mask(n, b)) for b in range(PC.NB)]
    f = (PC.NB - 1) / PC.NB
    r = np.array([x[0] for x in jk])
    v = np.array([x[1] for x in jk])
    ess = min(PC.value(data[T], 'Rg', us[T])[1] for T in temps)
    return rel, s, np.sqrt(f * ((r - r.mean(0)) ** 2).sum(0)), float(np.sqrt(f * ((v - v.mean()) ** 2).sum())), ess


def main():
    rows, crows = [], []
    for model in ('calvados', 'recharged'):
        M = PD.MODELS[model]
        g_tea, g_h = PC.capped_gamma(model)
        a0_tea, a0_hyd = M['gam']['aliph'], M['gam']['hyd']
        rt = 8.314462618e-3 * 300.0
        Rg = PD.fdt.R_GAS
        st_curves = pd.read_csv(os.path.join(PD.out_dir(model), 'curves.csv'))
        for form, temps in sorted(stored(model, PC.FS.TAG[model]).items()):
            if 300 not in temps:
                continue
            heat = [T for T in temps if T in TEMPS]
            temps = heat if len(heat) >= 4 else [300]
            data = {T: PC.run_frames(model, form, T)[0] for T in temps}
            z = data[300]
            lt, lh = dEdT(model, z, SF.HYDSET['A0'], np.inf)
            dev = np.abs(a0_tea * lt + a0_hyd * lh - z['D300']).max() / max(np.abs(z['D300']).max(), 1e-12)
            if dev > 1e-5:                       # storage precision (2e-6 seen); a convention error would be O(1)
                raise ValueError(f'{model} {form}: linear forcing differs from the stored D300 ({dev:.1e})')
            ct, ch = dEdT(model, z, SF.HYDSET['S1'], PC.CSAT)
            A = z['obs_Rg']
            tot, se = PD.fdt_total(A, z['G_T'] + (g_tea * ct + g_h * ch) / rt)
            tot_a0, se_a0 = PD.fdt_total(A, z['G_T'] + z['D300'] / rt)
            pl, pl_se = PD.fdt_total(A, z['G_T'])
            rows.append({'model': M['label'], 'form': form, 'total_capped': tot, 'total_capped_se': se, 'total_A0': tot_a0,
                         'total_A0_se': se_a0, 'plain': pl, 'plain_se': pl_se, 'linear_check': dev})
            if len(temps) < 4:
                continue
            u_a0, u_cap, u_0 = {}, {}, {T: np.zeros(len(data[T]['obs_Rg'])) for T in temps}
            for T in temps:
                et, eh = PC.energies(model, data[T], T, SF.HYDSET['A0'], np.inf)
                u_a0[T] = (a0_tea * et + a0_hyd * eh) / (Rg * T)
                et, eh = PC.energies(model, data[T], T, SF.HYDSET['S1'], PC.CSAT)
                u_cap[T] = (g_tea * et + g_h * eh) / (Rg * T)
            ref = st_curves[(st_curves.form == form) & (st_curves.curve == 'calibrated')]
            if len(ref):
                s_a0 = curve(data, u_a0, temps)[1]
                if abs(s_a0 - ref.slope300.iloc[0]) > 1e-6:
                    raise ValueError(f'{model} {form}: A0 curve slope {s_a0} differs from the stored {ref.slope300.iloc[0]}')
            for lab, us in (('plain', u_0), ('capped', u_cap)):
                rel, s, rse, sse, ess = curve_se(data, us, temps)
                crows.append({'model': M['label'], 'form': form, 'curve': lab, 'temps': ','.join(map(str, temps)),
                              'slope300': s, 'slope300_se': sse, 'min_ess': ess,
                              **{f'rel_{T}K': v for T, v in zip(temps, rel)}, **{f'rel_se_{T}K': v for T, v in zip(temps, rse)}})
            print(model, form, 'done', flush=True)
    df = pd.DataFrame(rows)
    os.makedirs(PC.OUT, exist_ok=True)
    df.to_csv(os.path.join(PC.OUT, 'fdt_300K.csv'), index=False)
    cv = pd.DataFrame(crows)
    cv.to_csv(os.path.join(PC.OUT, 'curves_capped.csv'), index=False)
    L = ['PAGE4 direct runs at 300 K, 150 mM: heating response of Rg (% per K), capped vs A0 forcing, 10-block SEs;',
         'curve slope = slope at 300 K of the directly simulated heating curve reweighted to the capped model', '']
    for m, g in df.groupby('model', sort=False):
        wt = g.set_index('form').loc['WT']
        c = cv[(cv.model == m) & (cv.curve == 'capped')].set_index('form')
        L.append(f'== {m}')
        for r in g.itertuples():
            d = r.total_capped - wt.total_capped
            se = np.hypot(r.total_capped_se, wt.total_capped_se)
            cs = f' | curve slope {c.loc[r.form, "slope300"]:+.3f} +- {c.loc[r.form, "slope300_se"]:.3f} [ESS {c.loc[r.form, "min_ess"]:.2f}]' \
                if r.form in c.index else ''
            L.append(f'  {r.form:10s} capped {r.total_capped:+.3f} +- {r.total_capped_se:.3f} (change from WT {d:+.3f} +- {se:.3f})'
                     f'{cs} | A0 {r.total_A0:+.3f} | plain {r.plain:+.3f} | linear check {r.linear_check:.1e}')
        L.append('')
    L.append('check passed: the A0 curve through this path reproduces every stored calibrated curve slope (1e-6)')
    open(os.path.join(PC.OUT, 'fdt_300K.txt'), 'w').write('\n'.join(L) + '\n')
    print('\n'.join(L))


if __name__ == '__main__':
    main()
