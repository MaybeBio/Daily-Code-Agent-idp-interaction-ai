"""
page4_capped_predictions.py — predictions of the capped (saturating) hydration forcing for the PAGE4 temperature series,
recorded next to the Amendment 4 predictions of the current forms before any PAGE4 temperature data exist (2026-10-24).
Exploratory model; the Amendment 4 predictions remain the test of the current calibration.

Capped model: the charged-group form (hydrophilic forcing on D E K R) with both terms saturating at c_sat = 2 full
contacts per residue (saturating_forcing.py), TEA strength at its prior and the hydrophilic strength refitted on the
four Wuttke constructs only (read from saturating_<model>.csv: 2.893 CALVADOS-2, 0.927 Mpipi-Recharged). Chosen
because it brought the direct PNt slopes to the measured ones in both force fields while keeping the Wuttke fit.

Frames: the 3,000 stored frames of each direct PAGE4 run (page4_direct.py, 150 mM, pH 7.5, 278-320 K), read again
from the DCDs by frame_index (Rg checked against the stored values); per-residue contact numbers as saturating_forcing.
contacts (CALVADOS-2 on the unmodified sequence, phosphosites keep Ser/Thr hydration as stored; Mpipi-Recharged on the
run sequence, phosphate-like Glu beads), checked against the stored class and hydrophilic energies.
Observables as Amendment 4 (smFRET efficiencies E 18-63 and 63-102, Rg) and as the NMR proposal
(page4_nmr_predictions.py: Kirkwood Rh, N-motif to acidic-region contacts cNA). Each temperature reweighted to the
capped model; 10-block jackknife; Kish effective sample fraction. Slopes per 10 K over 288-310 K (primary) and 278-320 K;
dampening 1 - s(form)/s(WT) as page4_registration.py. Check: the linear S1 form computed through the same path
reproduces the registered S1 slopes.
Not for the experimental collaborators (predicted magnitudes stay with us until the data are in).
    BENDER_REANALYSIS_DIR=... python page4_capped_predictions.py      (openmpipi env)
Output: heldout_temperature/capped_forcing/page4/{series.csv, slopes.csv, ratios.csv, summary.txt}
"""
import sys

sys.dont_write_bytecode = True

import os  # noqa: E402

import mdtraj as md  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

HERE = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, os.environ.get('BENDER_REANALYSIS_DIR', HERE))
sys.path.insert(0, HERE)
import page4_direct as PD  # noqa: E402
import page4_fret_series as FS  # noqa: E402
import page4_nmr_predictions as NMR  # noqa: E402
import page4_registration as PR  # noqa: E402
import saturating_forcing as SF  # noqa: E402

OUT = os.path.join(PD.BM, 'heldout_temperature/capped_forcing/page4')
FORMS = ['WT', 'HIPK1occ', 'CLK2occ']
OBS = ['E18_63', 'E63_102', 'Rg', 'Rh', 'cNA']
CSAT = 2.0
NB = FS.NB
TEMPS = FS.TEMPS


def capped_gamma(model):
    t = pd.read_csv(os.path.join(SF.OUT, f'saturating_{model}.csv'))
    row = t[(t.variant == 'S1') & (t.c_sat == CSAT)]
    assert len(row) == 1
    return SF.HX.S1_TEA[model], float(row.gamma_h.iloc[0])


def run_frames(model, form, T):
    """Stored per-frame data plus contact numbers, Rh and cNA of one direct run."""
    name = f'PAGE4-{form}{FS.TAG[model]}'
    z = dict(np.load(os.path.join(PD.out_dir(model), 'npz', f'{name}_{T}K.npz')))
    d = os.path.join(PD.MODELS[model]['root'], f'{T}K', name)
    xyz_all, _, _ = PD.HC.bc.load_calvados(os.path.join(d, f'{name}.dcd'), os.path.join(d, 'top.pdb'), k=None,
                                           production=False)
    xyz = xyz_all[z['frame_index']].astype(float)
    if PD.HC.bc.max_bond(xyz) > 1.0:
        xyz = PD.HC.bc.unwrap_chain(xyz, 50.0)
    rg = np.sqrt(((xyz - xyz.mean(1, keepdims=True)) ** 2).sum(2).mean(1))
    if np.abs(rg - z['obs_Rg']).max() > 1e-3 * z['obs_Rg'].mean():
        raise ValueError(f'{name} {T} K: frames differ from the stored ones')
    # the sequence the stored energies used: CALVADOS-2 the unmodified chain (phosphosites as Ser/Thr), read from the WT
    # run; Mpipi-Recharged the run's own chain (phosphate-like Glu beads); the energy check below confirms either choice
    ref = d if model == 'recharged' else os.path.join(PD.MODELS[model]['root'], f'{T}K', f'PAGE4-WT{FS.TAG[model]}')
    seq = ''.join(r.code for r in md.load_topology(os.path.join(ref, 'top.pdb')).residues)
    c = SF.contacts(model, xyz, seq)
    wt, wh = SF.weights(model, seq, float(T), SF.TCal.HYD_RES)
    et, eh = SF.energy(model, c, wt, np.inf), SF.energy(model, c, wh, np.inf)
    st = sum(z[f'E_{k}'] for k in PD.CLS)
    # relative to the energy size; at 300 K the forcing is zero by construction, so the absolute difference is used
    dev = max(np.abs(et - st).max(), np.abs(eh - z['E_hyd']).max()) / max(np.abs(st).max(), np.abs(z['E_hyd']).max(), 1.0)
    if dev > 1e-4:
        raise ValueError(f'{name} {T} K: E_lin from contacts differs from the stored energies ({dev:.1e})')
    n = xyz.shape[1]
    iu = np.triu_indices(n, 1)
    inv = np.zeros(len(xyz))
    for a in range(0, len(xyz), 200):
        inv[a:a + 200] = 2.0 * (1.0 / np.linalg.norm(xyz[a:a + 200, iu[0]] - xyz[a:a + 200, iu[1]], axis=2)).sum(1) / n ** 2
    dNA = np.linalg.norm(xyz[:, NMR.NMOTIF, None, :] - xyz[:, None, NMR.ACIDIC, :], axis=3)
    z.update(obs_invR=inv, obs_cNA=(dNA < NMR.CUT).sum((1, 2)).astype(float), c=c, seq=seq)
    return z, dev


def energies(model, z, T, form_hyd, csat):
    wt, wh = SF.weights(model, z['seq'], float(T), form_hyd)
    return SF.energy(model, z['c'], wt, csat), SF.energy(model, z['c'], wh, csat)


def value(z, obs, u, keep=None):
    k = np.ones(len(u), bool) if keep is None else keep
    w = np.exp(-(u[k] - u[k].min()))
    w /= w.sum()
    if obs == 'Rh':
        return 10.0 / float((w * z['obs_invR'][k]).sum()), float(1.0 / (w ** 2).sum() / len(w))
    a = z[f'obs_{obs}'][k]
    return float((w * a).sum()) * (10.0 if obs == 'Rg' else 1.0), float(1.0 / (w ** 2).sum() / len(w))


def series(ser, obs, us):
    n = len(us[TEMPS[0]])
    full = np.array([value(ser[T], obs, us[T])[0] for T in TEMPS])
    ess = np.array([value(ser[T], obs, us[T])[1] for T in TEMPS])
    jk = np.array([[value(ser[T], obs, us[T], FS.keep_mask(n, b))[0] for T in TEMPS] for b in range(NB)])
    return full, ess, jk


def main():
    os.makedirs(OUT, exist_ok=True)
    reg = pd.read_csv(os.path.join(PR.OUT, 'slopes.csv'))
    rows, srows, rrows, store, checks = [], [], [], {}, []
    for model in ('calvados', 'recharged'):
        lab = PD.MODELS[model]['label']
        gt, gh = capped_gamma(model)
        g_s1 = PR.G_CHG[model]
        data, devs = {}, []
        for form in FORMS:
            data[form] = {}
            for T in TEMPS:
                data[form][T], dv = run_frames(model, form, T)
                devs.append(dv)
        checks.append(f'{lab}: contact numbers reproduce the stored energies to {max(devs):.1e}')
        forms = {'S1 (linear, check)': ('DEKR', np.inf, gt, g_s1), f'capped S1, c_sat {CSAT:g}': ('DEKR', CSAT, gt, gh)}
        for fname, (hyd, csat, g_t, g_h) in forms.items():
            for form in FORMS:
                us = {}
                for T in TEMPS:
                    et, eh = energies(model, data[form][T], T, hyd, csat)
                    us[T] = (g_t * et + g_h * eh) / (PD.fdt.R_GAS * T)
                for obs in OBS:
                    full, ess, jk = series(data[form], obs, us)
                    se = FS.jk_se(jk)
                    for T, m, s, e in zip(TEMPS, full, se, ess):
                        rows.append({'model': lab, 'forcing': fname, 'form': form, 'obs': obs, 'T': T, 'mean': m,
                                     'se': s, 'ess': e})
                    for wname, temps in PR.WINDOWS.items():
                        ix = [TEMPS.index(t) for t in temps]
                        s0 = PR.slope(temps, full[ix])
                        rep = np.array([PR.slope(temps, r[ix]) for r in jk])
                        store[(model, fname, form, obs, wname)] = (s0, rep)
                        srows.append({'model': lab, 'forcing': fname, 'form': form, 'obs': obs, 'window': wname,
                                      'per10K': s0, 'per10K_se': float(FS.jk_se(rep[:, None])[0]),
                                      'min_ess': float(ess[ix].min())})
            for obs in OBS:
                for wname in PR.WINDOWS:
                    w0, wrep = store[(model, fname, 'WT', obs, wname)]
                    for form in ('HIPK1occ', 'CLK2occ'):
                        c0, crep = store[(model, fname, form, obs, wname)]
                        r0 = c0 / w0
                        se_jk = float(FS.jk_se((crep / wrep)[:, None])[0])
                        se_pr = abs(r0) * np.hypot(FS.jk_se(crep[:, None])[0] / c0, FS.jk_se(wrep[:, None])[0] / w0)
                        rrows.append({'model': lab, 'forcing': fname, 'form': form, 'obs': obs, 'window': wname,
                                      'ratio': r0, 'damp_pct': 100 * (1 - r0), 'damp_se': 100 * max(se_jk, float(se_pr)),
                                      'diff_per10K': c0 - w0,
                                      'diff_se': float(np.hypot(FS.jk_se(crep[:, None])[0], FS.jk_se(wrep[:, None])[0]))})
        # check: linear S1 through this path = registered S1
        for form in FORMS:
            for obs in ('E18_63', 'E63_102', 'Rg'):
                for wname in PR.WINDOWS:
                    mine = store[(model, 'S1 (linear, check)', form, obs, wname)][0]
                    ref = reg[(reg.model == lab) & (reg.forcing == 'S1') & (reg.form == form) & (reg.obs == obs)
                              & (reg.window == wname)].per10K.iloc[0]
                    if abs(mine - ref) > 1e-6 * max(1.0, abs(ref)):
                        raise ValueError(f'{lab} {form} {obs} {wname}: S1 {mine} differs from the registered {ref}')
        checks.append(f'{lab}: linear S1 through this path reproduces every registered S1 slope (E18_63, E63_102, Rg)')
    df, sl, ra = pd.DataFrame(rows), pd.DataFrame(srows), pd.DataFrame(rrows)
    df.to_csv(os.path.join(OUT, 'series.csv'), index=False)
    sl.to_csv(os.path.join(OUT, 'slopes.csv'), index=False)
    ra.to_csv(os.path.join(OUT, 'ratios.csv'), index=False)
    L = ['PAGE4 temperature series: predictions of the capped hydration forcing (charged-group form, c_sat 2), recorded '
         '2026-10-24 before any PAGE4 temperature data', ''] + checks + ['']
    for m in ('CALVADOS-2', 'Mpipi-Recharged'):
        t = pd.read_csv(os.path.join(SF.OUT, f"saturating_{'calvados' if m == 'CALVADOS-2' else 'recharged'}.csv"))
        r = t[(t.variant == 'S1') & (t.c_sat == CSAT)].iloc[0]
        L.append(f'{m}: gamma_TEA {SF.HX.S1_TEA["calvados" if m == "CALVADOS-2" else "recharged"]}, gamma_chg {r.gamma_h:.3f} '
                 f'(refitted on Wuttke, rms {r.wuttke_rms:.2f})')
    for wname in PR.WINDOWS:
        L += ['', f'== slopes per 10 K, {wname} (E absolute; Rg, Rh in Angstrom; cNA contacts) +- SE [min Kish fraction]']
        for (m, f), g in sl[sl.window == wname].groupby(['model', 'forcing'], sort=False):
            for obs in OBS:
                h = g[g.obs == obs].set_index('form')
                L.append(f'  {m:16s} {f:22s} {obs:7s} ' + ', '.join(
                    f"{fm} {h.loc[fm, 'per10K']:+.4f}+-{h.loc[fm, 'per10K_se']:.4f}[{h.loc[fm, 'min_ess']:.2f}]" for fm in FORMS))
        L.append('  dampening 1 - s(form)/s(WT), % +- SE')
        for (m, f), g in ra[ra.window == wname].groupby(['model', 'forcing'], sort=False):
            L.append(f'  {m:16s} {f:22s} ' + ', '.join(f"{x.obs} {x.form.replace('occ', '')} {x.damp_pct:+.0f}+-{x.damp_se:.0f}"
                                                      for x in g.itertuples()))
    L += ['', 'registered Amendment 4 values for comparison: heldout_temperature/amendment4/predictions/summary.txt']
    open(os.path.join(OUT, 'summary.txt'), 'w').write('\n'.join(L) + '\n')
    print('\n'.join(L))


if __name__ == '__main__':
    main()
