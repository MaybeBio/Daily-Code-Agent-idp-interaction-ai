"""
page4_registration.py — the predictions registered in Amendment 4 (heldout_temperature/amendment4/): the smFRET
temperature series of WT, HIPK1- and CLK2-PAGE4 (Kulkarni et al. PNAS 2017 constructs, Alexa 555/647, R0 5.24 nm)
under three forms of the solvation forcing, from the per-frame data page4_direct.py stored (no trajectories read):
  plain  no hydration term (reference)
  A0     the frozen calibration (hydrophilic forcing on D E K R N Q S T H G, gamma_hyd 0.358 / 0.209)
  S1     the Amendment 3 form (hydrophilic forcing on D E K R only, gamma_chg 0.469789 / 0.294234)
TEA term at the calibration prior in A0 and S1 (CALVADOS-2 2.0, Mpipi-Recharged 0.4). Each T is reweighted with
exp(-sum_c gamma_c E_c(T) / RT), as the Wuttke calibration and page4_fret_series.py (A0 here must reproduce its numbers).
Phosphosites carry the hydration energies of the unmodified Ser/Thr letters, as stored (so S1 puts no forcing on them).
Per point: mean, 10-block jackknife SE, Kish effective sample fraction. Slopes per 10 K from a linear fit:
primary over 288-310 K (the simulated temperatures inside the requested 10-40 C range, ESS >= 0.38), secondary over
278-320 K; jackknife leaving the same block out at every T. Per force field and forcing: the WT slope, the ratio
s(form) / s(WT) and the dampening 1 - s(form) / s(WT) for HIPK1 and CLK2 (SE by the jackknife of the ratio, WT and form
blocks paired by index; the runs are independent, so this is a resampling convenience, checked against propagation).
    BENDER_REANALYSIS_DIR=... python page4_registration.py       (openmpipi env)
Output: heldout_temperature/amendment4/predictions/{series.csv, slopes.csv, ratios.csv, summary.txt}
"""
import sys

sys.dont_write_bytecode = True

import os  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

HERE = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, HERE)
import page4_direct as PD  # noqa: E402
import page4_fret_series as FS  # noqa: E402

OUT = os.path.join(PD.BM, 'heldout_temperature/amendment4/predictions')
FORMS = ['WT', 'HIPK1occ', 'CLK2occ']
OBS = ['E18_63', 'E63_102', 'Rg']
TEA = ('aliph', 'pro', 'arom', 'arg', 'rest')
G_CHG = {'calvados': 0.469789, 'recharged': 0.294234}           # Amendment 3, frozen
WINDOWS = {'primary 288-310 K': (288, 300, 310), 'secondary 278-320 K': (278, 288, 300, 310, 320)}
NB = FS.NB


def forcings(model):
    g = PD.MODELS[model]['gam']
    g_tea = g['aliph']
    return {'plain': {c: 0.0 for c in g},
            'A0': dict(g),
            'S1': dict({c: g_tea for c in TEA}, **{f'hyd_{k}': G_CHG[model] for k in ('DE', 'K', 'R')})}


def slope(temps, y):
    return 10 * np.polyfit(np.array(temps, float), y, 1)[0]


def main():
    os.makedirs(OUT, exist_ok=True)
    rows, srows, rrows, store = [], [], [], {}
    for model in ('calvados', 'recharged'):
        lab_m = PD.MODELS[model]['label']
        data = {form: FS.load(model, form) for form in FORMS}
        for form, ser in data.items():                       # stored hydrophilic groups add up to the class
            for T, z in ser.items():
                d = np.abs(sum(z[f'E_hyd_{k}'] for k in ('DE', 'K', 'R', 'pol', 'G')) - z['E_hyd']).max()
                if d > 1e-6 * max(1.0, np.abs(z['E_hyd']).max()):
                    raise ValueError(f'{model} {form} {T} K: hydrophilic groups miss the total by {d}')
        for fname, gam in forcings(model).items():
            for form in FORMS:
                for obs in OBS:
                    full, ess, jk = FS.series(data[form], obs, gam)
                    sc = 10.0 if obs == 'Rg' else 1.0                      # Rg in Angstrom
                    full, jk = sc * full, sc * jk
                    se = FS.jk_se(jk)
                    for T, m, s, e in zip(FS.TEMPS, full, se, ess):
                        rows.append({'model': lab_m, 'forcing': fname, 'form': form, 'obs': obs, 'T': T,
                                     'mean': m, 'se': s, 'ess': e})
                    for wname, temps in WINDOWS.items():
                        ix = [FS.TEMPS.index(t) for t in temps]
                        s0 = slope(temps, full[ix])
                        rep = np.array([slope(temps, r[ix]) for r in jk])
                        store[(model, fname, form, obs, wname)] = (s0, rep)
                        srows.append({'model': lab_m, 'forcing': fname, 'form': form, 'obs': obs, 'window': wname,
                                      'per10K': s0, 'per10K_se': float(FS.jk_se(rep[:, None])[0]),
                                      'min_ess': float(ess[ix].min())})
        for fname in forcings(model):
            for obs in OBS:
                for wname in WINDOWS:
                    w0, wrep = store[(model, fname, 'WT', obs, wname)]
                    for form in ('HIPK1occ', 'CLK2occ'):
                        c0, crep = store[(model, fname, form, obs, wname)]
                        r0 = c0 / w0
                        se_jk = float(FS.jk_se((crep / wrep)[:, None])[0])
                        se_pr = abs(r0) * np.hypot(FS.jk_se(crep[:, None])[0] / c0, FS.jk_se(wrep[:, None])[0] / w0)
                        rrows.append({'model': lab_m, 'forcing': fname, 'form': form, 'obs': obs, 'window': wname,
                                      'ratio': r0, 'ratio_se': max(se_jk, float(se_pr)), 'damp_pct': 100 * (1 - r0),
                                      'damp_se': 100 * max(se_jk, float(se_pr)),
                                      'diff_per10K': c0 - w0,
                                      'diff_se': float(np.hypot(FS.jk_se(crep[:, None])[0], FS.jk_se(wrep[:, None])[0]))})
    df, sl, ra = pd.DataFrame(rows), pd.DataFrame(srows), pd.DataFrame(rrows)
    # A0 must reproduce page4_fret_series (same reweighting, same blocks)
    ref = pd.read_csv(os.path.join(FS.OUT, 'series.csv'))
    ref = ref[ref.hydration == 'calibrated'].set_index(['model', 'form', 'obs', 'T'])['mean']
    mine = df[df.forcing == 'A0'].set_index(['model', 'form', 'obs', 'T'])['mean']
    dev = float((mine - ref.reindex(mine.index)).abs().max())
    if not dev < 1e-9:
        raise ValueError(f'A0 differs from page4_fret_series by {dev}')
    df.to_csv(os.path.join(OUT, 'series.csv'), index=False)
    sl.to_csv(os.path.join(OUT, 'slopes.csv'), index=False)
    ra.to_csv(os.path.join(OUT, 'ratios.csv'), index=False)

    L = ['PAGE4 smFRET temperature series, predictions registered in Amendment 4 (150 mM, pH 7.5 runs)',
         f'check: A0 reproduces page4_fret_series calibrated means to {dev:.1e}', '',
         'slopes per 10 K (E absolute; Rg in Angstrom) +- jackknife SE [min Kish fraction]']
    for wname in WINDOWS:
        L.append(f'== {wname}')
        for (m, f), g in sl[sl.window == wname].groupby(['model', 'forcing'], sort=False):
            parts = []
            for obs in OBS:
                h = g[g.obs == obs].set_index('form')
                parts.append(f'{obs}: ' + ', '.join(f"{fm} {h.loc[fm, 'per10K']:+.4f}+-{h.loc[fm, 'per10K_se']:.4f}"
                                                    f"[{h.loc[fm, 'min_ess']:.2f}]" for fm in FORMS))
            L.append(f'  {m:16s} {f:5s} ' + ' | '.join(parts))
        L.append('  dampening 1 - s(form)/s(WT), % (+- SE)')
        for (m, f), g in ra[ra.window == wname].groupby(['model', 'forcing'], sort=False):
            if f == 'plain':
                continue
            parts = [f"{r.obs} {r.form[:5]} {r.damp_pct:+.0f}+-{r.damp_se:.0f}" for _, r in g.iterrows()]
            L.append(f'  {m:16s} {f:5s} ' + ', '.join(parts))
        L.append('')
    L.append('absolute means at 300 K (E; Rg in Angstrom)')
    h = df[df['T'] == 300]
    for (m, f), g in h.groupby(['model', 'forcing'], sort=False):
        L.append(f'  {m:16s} {f:5s} ' + ', '.join(f"{r.obs} {r.form[:5]} {r['mean']:.3f}" for _, r in g.iterrows()))
    open(os.path.join(OUT, 'summary.txt'), 'w').write('\n'.join(L) + '\n')
    print('\n'.join(L))


if __name__ == '__main__':
    main()
