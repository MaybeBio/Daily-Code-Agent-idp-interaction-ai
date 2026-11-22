"""
page4_nmr_predictions.py — what the PAGE4 direct simulations predict for two NMR temperature series (the proposal to the
Daughdrill lab): paramagnetic relaxation enhancement from a spin label at the native Cys63, and the hydrodynamic radius
from diffusion NMR, for WT, HIPK1- and CLK2-PAGE4 at 278-320 K (150 mM, pH 7.5 runs of page4_direct.py).

Per frame (the 3,000 frames page4_direct.py stored, read again from the DCDs by their frame_index; Rg checked against the
stored values):
  r6_i     r^-6 between the bead of residue i and the bead of Cys63 (nm^-6); the PRE rate is proportional to <r^-6> at a
           fixed correlation time, so its temperature dependence at fixed tau_c is the susceptibility of <r^-6>
  r6N      the mean of r6_i over the N-motif (residues 4-12)
  cNA      contacts between the N-motif (4-12) and the central acidic region (43-62), bead centres closer than 1.2 nm
  invR     (1/N^2) sum_{i != j} 1/r_ij; the Kirkwood hydrodynamic radius is Rh = 1/<invR>, as in the held-out test
Each temperature is reweighted with exp(-sum_c gamma_c E_c(T) / RT) under the three forcings of Amendment 4 (plain, A0
the frozen calibration, S1 the charged-group form; page4_registration.forcings), with 10-block jackknife SEs and the Kish
effective sample fraction. Slopes per 10 K from a linear fit over 288-310 K (primary, as Amendment 4) and 278-320 K.
Susceptibility here = relative slope, % per K, of each quantity; the phosphorylation effect = form minus WT.
Exploratory, after the fact; the values are to be frozen and hashed before any NMR measurement.
    BENDER_REANALYSIS_DIR=... python page4_nmr_predictions.py      (openmpipi env)
Output: Traj Level Analysis/Analysis_FDT/total_response/page4/nmr_predictions/{series.csv, profiles.csv, slopes.csv,
        summary.txt}
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
import page4_registration as PR  # noqa: E402

OUT = os.path.join(PD.BM, 'Traj Level Analysis/Analysis_FDT/total_response/page4/nmr_predictions')
FORMS = ['WT', 'HIPK1occ', 'CLK2occ']
LABEL = PD.MODELS
C63 = 62                          # 0-based index of Cys63
NMOTIF = np.arange(3, 12)         # residues 4-12
ACIDIC = np.arange(42, 62)        # residues 43-62
CUT = 1.2                         # nm, as the Fig. 7h contact maps
WINDOWS = PR.WINDOWS


def frame_obs(model, form, T):
    name = f'PAGE4-{form}{FS.TAG[model]}'
    z = dict(np.load(os.path.join(PD.out_dir(model), 'npz', f'{name}_{T}K.npz')))
    d = os.path.join(PD.MODELS[model]['root'], f'{T}K', name)
    xyz_all, _, _ = PD.HC.bc.load_calvados(os.path.join(d, f'{name}.dcd'), os.path.join(d, 'top.pdb'), k=None,
                                           production=False)
    xyz = xyz_all[z['frame_index']]
    if PD.HC.bc.max_bond(xyz) > 1.0:
        xyz = PD.HC.bc.unwrap_chain(xyz, 50.0)
    rg = np.sqrt(((xyz - xyz.mean(1, keepdims=True)) ** 2).sum(2).mean(1))
    dev = np.abs(rg - z['obs_Rg']).max()
    if dev > 1e-3 * z['obs_Rg'].mean():
        raise ValueError(f'{name} {T} K: Rg differs from the stored frames by {dev:.4f} nm')
    n = xyz.shape[1]
    r63 = np.linalg.norm(xyz - xyz[:, [C63]], axis=2)
    with np.errstate(divide='ignore'):
        r6 = np.where(np.arange(n) == C63, np.nan, r63 ** -6.0)
    for i in range(n):
        if i != C63:
            z[f'obs_r6_{i + 1}'] = r6[:, i]
    z['obs_r6N'] = r6[:, NMOTIF].mean(1)
    dNA = np.linalg.norm(xyz[:, NMOTIF, None, :] - xyz[:, None, ACIDIC, :], axis=3)
    z['obs_cNA'] = (dNA < CUT).sum((1, 2)).astype(float)
    iu = np.triu_indices(n, 1)
    inv = np.zeros(len(xyz))
    for a in range(0, len(xyz), 200):
        dd = np.linalg.norm(xyz[a:a + 200, iu[0]] - xyz[a:a + 200, iu[1]], axis=2)
        inv[a:a + 200] = 2.0 * (1.0 / dd).sum(1) / n ** 2
    z['obs_invR'] = inv
    return z, n


def slope_per10(temps, y):
    return 10 * np.polyfit(np.array(temps, float), y, 1)[0]


def main():
    os.makedirs(OUT, exist_ok=True)
    rows, prof, srows, store = [], [], [], {}
    for model in ('calvados', 'recharged'):
        lab = PD.MODELS[model]['label']
        data, n = {}, None
        for form in FORMS:
            data[form] = {}
            for T in FS.TEMPS:
                data[form][T], n = frame_obs(model, form, T)
            print(f'{lab} {form}: frames read and checked', flush=True)
        for fname, gam in PR.forcings(model).items():
            for form in FORMS:
                ser = data[form]
                for obs in ('r6N', 'cNA', 'invR'):
                    full, ess, jk = FS.series(ser, obs, gam)
                    if obs == 'invR':                        # Rh = 1 / <invR>, in Angstrom
                        full, jk, obs = 10.0 / full, 10.0 / jk, 'Rh'
                    se = FS.jk_se(jk)
                    for T, m, s, e in zip(FS.TEMPS, full, se, ess):
                        rows.append({'model': lab, 'forcing': fname, 'form': form, 'obs': obs, 'T': T, 'mean': m,
                                     'se': s, 'ess': e})
                    i300 = FS.TEMPS.index(300)
                    for wname, temps in WINDOWS.items():
                        ix = [FS.TEMPS.index(t) for t in temps]
                        s0 = slope_per10(temps, full[ix]) / full[i300] * 10       # % per K
                        rep = np.array([slope_per10(temps, r[ix]) / r[i300] * 10 for r in jk])
                        store[(model, fname, form, obs, wname)] = (s0, rep)
                        srows.append({'model': lab, 'forcing': fname, 'form': form, 'obs': obs, 'window': wname,
                                      'pct_per_K': s0, 'se': float(FS.jk_se(rep[:, None])[0]),
                                      'min_ess': float(ess[ix].min())})
                if fname != 'plain' and fname != 'A0':
                    continue
                for i in range(1, n + 1):                    # PRE profiles, plain and A0
                    if i == C63 + 1:
                        continue
                    full, ess, jk = FS.series(ser, f'r6_{i}', gam)
                    for T, m, s in zip(FS.TEMPS, full, FS.jk_se(jk)):
                        prof.append({'model': lab, 'forcing': fname, 'form': form, 'residue': i, 'T': T,
                                     'r6_mean_nm-6': m, 'r6_se': s, 'r_eff_nm': m ** (-1 / 6)})
    df, sl = pd.DataFrame(rows), pd.DataFrame(srows)
    for (model, fname, form, obs, wname), (s0, rep) in list(store.items()):
        if form == 'WT':
            continue
        w0, wrep = store[(model, fname, 'WT', obs, wname)]
        m = sl[(sl.model == PD.MODELS[model]['label']) & (sl.forcing == fname) & (sl.form == form) & (sl.obs == obs)
               & (sl.window == wname)].index
        sl.loc[m, 'minus_WT'] = s0 - w0
        sl.loc[m, 'minus_WT_se'] = float(np.hypot(FS.jk_se(rep[:, None])[0], FS.jk_se(wrep[:, None])[0]))
    df.to_csv(os.path.join(OUT, 'series.csv'), index=False)
    pd.DataFrame(prof).to_csv(os.path.join(OUT, 'profiles.csv'), index=False)
    sl.to_csv(os.path.join(OUT, 'slopes.csv'), index=False)
    L = ['PAGE4 NMR temperature series predicted from the direct runs (150 mM, pH 7.5), relative slopes in % per K',
         'r6N = <r^-6> from Cys63 over the N-motif (PRE at fixed tau_c); cNA = N-motif to acidic-region contacts; Rh = Kirkwood',
         'radius. Each cell: form slope +- SE; for HIPK1 and CLK2 also form minus WT +- SE; [min Kish fraction]', '']
    for wname in WINDOWS:
        L.append(f'== {wname}')
        g0 = sl[sl.window == wname]
        for (m, f), g in g0.groupby(['model', 'forcing'], sort=False):
            parts = []
            for obs in ('r6N', 'cNA', 'Rh'):
                h = g[g.obs == obs].set_index('form')
                txt = [f"WT {h.loc['WT', 'pct_per_K']:+.3f}+-{h.loc['WT', 'se']:.3f}"]
                for fm in ('HIPK1occ', 'CLK2occ'):
                    txt.append(f"{fm[:5]} {h.loc[fm, 'pct_per_K']:+.3f} (diff {h.loc[fm, 'minus_WT']:+.3f}+-{h.loc[fm, 'minus_WT_se']:.3f})")
                parts.append(f'{obs}: ' + ', '.join(txt) + f" [{h['min_ess'].min():.2f}]")
            L.append(f'  {m:16s} {f:5s}')
            L += [f'      {p}' for p in parts]
        L.append('')
    L.append('values at 300 K (r6N nm^-6; cNA contacts; Rh Angstrom)')
    h = df[df['T'] == 300]
    for (m, f), g in h.groupby(['model', 'forcing'], sort=False):
        L.append(f'  {m:16s} {f:5s} ' + ', '.join(f"{r.obs} {r.form[:5]} {r['mean']:.4g}" for _, r in g.iterrows()))
    open(os.path.join(OUT, 'summary.txt'), 'w').write('\n'.join(L) + '\n')
    print('\n'.join(L))


if __name__ == '__main__':
    main()
