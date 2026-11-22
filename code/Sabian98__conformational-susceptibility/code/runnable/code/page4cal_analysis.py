"""
page4cal_analysis.py — PAGE4 heating curves simulated directly with the saturating hydration forcing
(GAIVI_PAGE4CAL, Mpipi-Recharged 2026-10-25, CALVADOS-2 2026-11-02). The forcing vanishes at 300 K, so the plain 300 K runs
of page4_direct.py are also the 300 K point of the calibrated model; the new runs supply 278, 288, 310 and 320 K.
For each form:
  direct calibrated curve of <Rg> (and the smFRET efficiencies) relative to 300 K, 3,000 evenly spaced frames of the
  second half of each run (page4_direct.load_frames), slope at 300 K from a quadratic fit, 10-block jackknife SEs that
  leave the same block out at every temperature
  compared with (i) the fluctuation-response total response at 300 K of the calibrated model (page4_capped_fdt.py,
  fdt_300K.csv) and (ii) the reweighted calibrated curve of the plain runs (curves_capped.csv), whose effective sample
  fraction fell below 5% for every Mpipi-Recharged form
  change from WT of the direct slope (the CLK2 drop in heat sensitivity)
Checks: every run complete (DCD header = frames), chain intact after unwrapping, top.pdb sequence = package forms.fasta
(Mpipi-Recharged) or = residue names of the plain 300 K run (CALVADOS-2, whose custom phosphoresidues the loader does not
map back to the fasta letters).
Four CALVADOS-2 runs were interrupted and resumed from a checkpoint (278 K CLK2full, Nmotif; 288 K HIPK1occ, WT): 30,000
frames follow the 4,500 to 5,000 frames of the first attempt, so the second half analysed lies wholly after the restart.
    BENDER_REANALYSIS_DIR=... python page4cal_analysis.py [recharged|calvados]      (openmpipi env, default recharged)
Output: heldout_temperature/capped_forcing/page4/direct_cal_{recharged,calvados}/{npz/, curves.csv, summary.txt}
"""
import sys

sys.dont_write_bytecode = True

import os  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

HERE = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, os.environ.get('BENDER_REANALYSIS_DIR', HERE))
sys.path.insert(0, HERE)
import page4_direct as PD  # noqa: E402
import page4_phosphoforms as P4  # noqa: E402

BM = PD.BM
RUNS = os.path.join(BM, 'page4cal_MPIPI_Recharged')
PKG = os.path.expanduser((__import__('os').environ['CS_DATA'] + '/Scripts/GAIVI_PAGE4CAL/mpipi'))
CAPD = os.path.join(BM, 'heldout_temperature', 'capped_forcing', 'page4')
OUT = os.path.join(CAPD, 'direct_cal_recharged')
TEMPS = [278, 288, 300, 310, 320]
OBS = ['Rg', 'E18_63', 'E63_102']
NB = 10
LABEL = 'Mpipi-Recharged'
MODEL, TAG = 'recharged', '-I150'


def use(model):
    """Point the module at one force field (recharged is the import default)."""
    global RUNS, PKG, OUT, LABEL, MODEL, TAG
    MODEL = model
    if model == 'calvados':
        RUNS = os.path.join(BM, 'page4cal_CALVADOS2')
        PKG = os.path.expanduser((__import__('os').environ['CS_DATA'] + '/Scripts/GAIVI_PAGE4CAL/calvados'))
        OUT = os.path.join(CAPD, 'direct_cal_calvados')
        LABEL, TAG = 'CALVADOS-2', '-I150-pH75'


def top_residues(path):
    return [x[17:20].strip() for x in open(path) if x.startswith('ATOM')]


def forms_fasta():
    fa = {}
    for block in open(os.path.join(PKG, 'forms.fasta')).read().split('>')[1:]:
        h, *s = block.split('\n')
        fa[h.strip()] = ''.join(s).strip()
    return fa


def run_obs(form, T, seq):
    """Per-frame observables of one temperature, cached; 300 K from the plain run of page4_direct."""
    if T == 300:
        z = np.load(os.path.join(PD.out_dir(MODEL), 'npz', f'PAGE4-{form}{TAG}_300K.npz'))
        return {k: z[f'obs_{k}'] for k in OBS}, 'plain 300 K run (forcing zero at 300 K)'
    cache = os.path.join(OUT, 'npz', f'{form}_{T}K.npz')
    if os.path.exists(cache):
        z = np.load(cache)
        return {k: z[k] for k in OBS}, str(z['note'])
    name = f'PAGE4cal-{form}{TAG}'
    row = {'dir': os.path.join(RUNS, f'{T}K', name), 'name': name}
    xyz_all, xyz, idx, mb, top_seq = PD.load_frames(row)
    if MODEL == 'calvados':
        plain = os.path.join(PD.MODELS[MODEL]['root'], '300K', f'PAGE4-{form}{TAG}', 'top.pdb')
        if top_residues(os.path.join(row['dir'], 'top.pdb')) != top_residues(plain):
            raise ValueError(f'{name} {T} K: top.pdb residues differ from the plain 300 K run')
    elif top_seq != seq:
        raise ValueError(f'{name} {T} K: top.pdb sequence differs from forms.fasta')
    bond = float(PD.HC.bc.max_bond(xyz))
    if bond > 1.0:
        raise ValueError(f'{name} {T} K: chain broken after unwrapping (max bond {bond:.2f} nm)')
    obs = P4.observables(xyz)
    note = f'{len(xyz_all)} frames, max bond {bond:.3f} nm'
    os.makedirs(os.path.dirname(cache), exist_ok=True)
    np.savez_compressed(cache, note=note, frame_index=idx, **{k: obs[k] for k in OBS})
    return {k: obs[k] for k in OBS}, note


def curve(series, k, drop=None):
    m = {}
    for T, o in series.items():
        a = o[k]
        if drop is not None:
            a = np.delete(a, np.array_split(np.arange(len(a)), NB)[drop])
        m[T] = a.mean()
    Ts = np.array(sorted(m))
    rel = np.array([m[T] / m[300] for T in Ts])
    absl = np.array([m[T] for T in Ts])
    return Ts, rel, absl


def slope(Ts, y, relative):
    c = np.polyfit(Ts - 300.0, y, 2)[1]
    return 100 * c if relative else 10 * c          # Rg: % per K of the 300 K value; E: change per 10 K


def with_se(series, k):
    Ts, rel, absl = curve(series, k)
    y = rel if k == 'Rg' else absl
    s = slope(Ts, y, k == 'Rg')
    jk = [curve(series, k, b) for b in range(NB)]
    sj = np.array([slope(t, r if k == 'Rg' else a, k == 'Rg') for t, r, a in jk])
    yj = np.array([r if k == 'Rg' else a for _, r, a in jk])
    f = (NB - 1) / NB
    return Ts, y, s, float(np.sqrt(f * ((sj - sj.mean()) ** 2).sum())), np.sqrt(f * ((yj - yj.mean(0)) ** 2).sum(0)), sj


def main():
    os.makedirs(OUT, exist_ok=True)
    fa = forms_fasta()
    forms = [f for f in fa if os.path.exists(os.path.join(RUNS, '278K', f'PAGE4cal-{f}{TAG}'))]
    fdt = pd.read_csv(os.path.join(CAPD, 'fdt_300K.csv'))
    fdt = fdt[fdt.model == LABEL].set_index('form')
    rw = pd.read_csv(os.path.join(CAPD, 'curves_capped.csv'))
    rw = rw[(rw.model == LABEL)].set_index(['form', 'curve'])
    rows, notes, jk = [], [], {}
    for form in forms:
        series = {}
        for T in TEMPS:
            series[T], note = run_obs(form, T, fa[form])
            notes.append(f'  {form:9s} {T} K: {note}')
        for k in OBS:
            Ts, y, s, se, yse, sj = with_se(series, k)
            jk[(form, k)] = sj
            r = {'form': form, 'obs': k, 'slope300': s, 'slope300_se': se,
                 **{f'{"rel" if k == "Rg" else "mean"}_{T}K': v for T, v in zip(Ts, y)},
                 **{f'se_{T}K': v for T, v in zip(Ts, yse)}}
            if k == 'Rg':
                r['fdt_total'], r['fdt_total_se'] = fdt.loc[form, 'total_capped'], fdt.loc[form, 'total_capped_se']
                r['reweighted'], r['reweighted_se'] = rw.loc[(form, 'capped'), 'slope300'], rw.loc[(form, 'capped'), 'slope300_se']
                r['reweighted_min_ess'] = rw.loc[(form, 'capped'), 'min_ess']
                r['plain_curve'] = rw.loc[(form, 'plain'), 'slope300']
            rows.append(r)
    cv = pd.DataFrame(rows)
    cv.to_csv(os.path.join(OUT, 'curves.csv'), index=False)

    g = cv[cv.obs == 'Rg'].set_index('form')
    L = [f'PAGE4 heating curves simulated with the saturating hydration forcing, {LABEL} (300 ns per run, 3,000 frames of '
         'the second half; 300 K = plain run), 10-block jackknife SEs', '', 'runs'] + notes + ['']
    L += ['(1) slope of <Rg> at 300 K (% per K): direct calibrated curve, fluctuation-response at 300 K, reweighted curve',
          '    (min ESS = smallest Kish fraction of the reweighting over temperatures), plain direct curve']
    t1 = pd.DataFrame({'direct': g.slope300, 'se': g.slope300_se, 'FDT 300 K': g.fdt_total, 'se_f': g.fdt_total_se,
                       'direct-FDT (SE)': (g.slope300 - g.fdt_total) / np.hypot(g.slope300_se, g.fdt_total_se),
                       'reweighted': g.reweighted, 'se_r': g.reweighted_se, 'min ESS': g.reweighted_min_ess,
                       'plain curve': g.plain_curve})
    L += [t1.round(3).to_string(), '']
    n2 = int((t1['direct-FDT (SE)'].abs() <= 2).sum())
    L += [f'    direct and FDT agree within two combined SEs for {n2} of {len(t1)} forms; '
          f'mean |direct - FDT| {np.mean(np.abs(g.slope300 - g.fdt_total)):.3f} % per K', '']
    L += ['(2) change of the direct calibrated slope from WT (% per K; positive = less heat-sensitive), jackknife over blocks '
          'paired across forms is not possible (independent runs), so SEs are combined in quadrature']
    t2 = []
    for form in [f for f in forms if f != 'WT']:
        d = g.loc[form, 'slope300'] - g.loc['WT', 'slope300']
        se = np.hypot(g.loc[form, 'slope300_se'], g.loc['WT', 'slope300_se'])
        df = g.loc[form, 'fdt_total'] - g.loc['WT', 'fdt_total']
        dfse = np.hypot(g.loc[form, 'fdt_total_se'], g.loc['WT', 'fdt_total_se'])
        t2.append({'form': form, 'direct change': d, 'se': se, 'z': d / se, 'FDT change': df, 'se_f': dfse})
    L += [pd.DataFrame(t2).set_index('form').round(3).to_string(), '']
    L += ['(3) relative <Rg> at each temperature, direct calibrated']
    L += [g[[f'rel_{T}K' for T in TEMPS]].round(4).to_string(), '']
    e = cv[cv.obs == 'E18_63'].set_index('form')
    L += ['(4) smFRET efficiency 18-63, direct calibrated: slope per 10 K at 300 K']
    L += [e[['slope300', 'slope300_se']].round(4).to_string()]
    open(os.path.join(OUT, 'summary.txt'), 'w').write('\n'.join(L) + '\n')
    print('\n'.join(L))


if __name__ == '__main__':
    use(sys.argv[1] if len(sys.argv) > 1 else 'recharged')
    main()
