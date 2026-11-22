"""
pnt_cal_analysis.py — the direct calibrated PNt runs (GAIVI_PNTCAL, make_pnt_calibrated_inputs.py): relative slope of
<Rg> with temperature under A0 and S1, simulated with the forcing in the model, against the measured slope (Baxa et al.
2024) and the reweighted values of the held-out test. Exploratory, after the test; nothing frozen is changed.

Frames as heldout_analysis.load: 3,000 evenly spaced frames of the second half of the run, chain unwrapped in the box.
Each DCD must hold the 60,000 frames of one run. Two CALVADOS-2 files (A0, 303 K, both salts) also hold the first
frames of a job that stopped and was restarted from scratch onto the same file: there the frames before the single
discontinuity are dropped when exactly 60,000 continuous frames follow it (checked: the only jump larger than five
times the median frame-to-frame change is at that join).
Slope s = 100 d ln<Rg>/dT by least squares over the simulated temperatures, 10-block jackknife with the same block
removed at every temperature (heldout_analysis.slope with no reweighting). Also, per temperature, <Rg> of the direct
run against the plain held-out run and the reweighted value from it (with its Kish effective sample fraction).
    BENDER_REANALYSIS_DIR=... python pnt_cal_analysis.py [--model calvados|recharged|both]     (openmpipi env)
Output: heldout_temperature/results/exploratory/pnt_cal_<model>.csv, pnt_cal_<model>_temps.csv, pnt_cal_summary.txt
"""
import sys

sys.dont_write_bytecode = True

import argparse  # noqa: E402
import os  # noqa: E402

import mdtraj as md  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

HERE = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, os.environ.get('BENDER_REANALYSIS_DIR', HERE))
sys.path.insert(0, HERE)
import heldout_analysis as HA  # noqa: E402
import heldout_exploratory as HX  # noqa: E402
import make_heldout_inputs as MH  # noqa: E402
import make_pnt_calibrated_inputs as MP  # noqa: E402

ROOT = {'calvados': f'{HA.BM}/pntcal_CALVADOS2', 'recharged': f'{HA.BM}/pntcal_MPIPI_Recharged'}
NAME = {'calvados': lambda st, v: f'HC-{st}-{v}-pH74', 'recharged': lambda st, v: f'HC-{st}-{v}'}
NFR = 60000
OUT = os.path.join(HA.OUT, 'exploratory')


def frames(model, st, v, T):
    xyz, note = load_xyz(model, st, v, T)
    o = HA.observables(xyz, 'PNt')
    return {'obs_Rg': o['Rg'], 'obs_inv': o['inv'], 'E_aliph': np.zeros(len(xyz))}, note


def load_xyz(model, st, v, T):
    """The 3,000 analysed frames (nm, chain unwrapped) of one direct run, and a note when frames were dropped."""
    name = NAME[model](st, v)
    d = os.path.join(ROOT[model], f'{T}K', name)
    t = md.load(os.path.join(d, f'{name}.dcd'), top=os.path.join(d, 'top.pdb'))
    seq = ''.join(r.code for r in t.topology.residues)
    if seq != MH.PNT:
        raise ValueError(f'{name} {T} K: top.pdb is not PNt')
    x = t.xyz.astype(float)
    note = ''
    if len(x) != NFR:
        c = x - x.mean(1, keepdims=True)
        step = np.sqrt(((c[1:] - c[:-1]) ** 2).sum(2).mean(1))
        jumps = np.where(step > 5 * np.median(step))[0]
        if len(x) > NFR and list(jumps) == [len(x) - NFR - 1]:
            note = f'{len(x)} frames; first {len(x) - NFR} (stopped job) dropped at the only discontinuity'
            x = x[len(x) - NFR:]
        else:
            raise ValueError(f'{name} {T} K: {len(x)} frames, discontinuities at {jumps[:5]}')
    n = len(x)
    idx = n // 2 + HA.bc.even_indices(n - n // 2, HA.K)
    xyz = x[idx]
    if HA.bc.max_bond(xyz) > 1.0:
        xyz = HA.bc.unwrap_chain(xyz, HA.box_nm(len(seq)))
    return xyz, note


def reweighted_gam(model, v):
    return HA.frozen_gam(model) if v == 'A0' else HX.chg_gam(HX.S1_GAMMA[model], HX.S1_TEA[model])


def analyse(model):
    s_exp, _ = HA.experimental()
    none = {'tea': 0.0}
    rows, trows, notes = [], [], []
    for st, (ionic, temps) in MP.SETS.items():
        plain = HA.series(model, st)
        for v in ('A0', 'S1'):
            ser = {}
            for T in temps:
                ser[T], note = frames(model, st, v, T)
                if note:
                    notes.append(f'{NAME[model](st, v)} {T} K: {note}')
                g = reweighted_gam(model, v)
                trows.append({'model': HA.LABEL[model], 'set': st, 'variant': v, 'T': T,
                              'Rg_direct_nm': ser[T]['obs_Rg'].mean(),
                              'Rg_plain_nm': plain[T]['obs_Rg'].mean(),
                              'Rg_reweighted_nm': HA.ensemble_value(plain[T], 'Rg', g, float(T)),
                              'ess_reweighted': HA.ess(plain[T], g, float(T))})
            s, se, _ = HA.slope(ser, 'Rg', none)
            sp, sep, _ = HA.slope(plain, 'Rg', none)
            sr, ser_, _ = HA.slope(plain, 'Rg', reweighted_gam(model, v))
            rows.append({'model': HA.LABEL[model], 'set': st, 'variant': v, 's_exp': s_exp[st],
                         's_direct': s, 'se_direct': se, 's_reweighted': sr, 'se_reweighted': ser_,
                         's_plain': sp, 'se_plain': sep,
                         'min_ess_reweighted': min(HA.ess(plain[T], reweighted_gam(model, v), float(T)) for T in temps),
                         'direct_minus_exp_in_se': (s - s_exp[st]) / se})
    return pd.DataFrame(rows), pd.DataFrame(trows), notes


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--model', default='calvados', choices=['calvados', 'recharged', 'both'])
    a = ap.parse_args()
    models = ['calvados', 'recharged'] if a.model == 'both' else [a.model]
    os.makedirs(OUT, exist_ok=True)
    L = ['Direct calibrated PNt runs: relative slope of <Rg> (% per K), direct against measured and reweighted', '']
    for m in models:
        R, Tt, notes = analyse(m)
        tag = 'calvados' if m == 'calvados' else 'recharged'
        R.to_csv(os.path.join(OUT, f'pnt_cal_{tag}.csv'), index=False)
        Tt.to_csv(os.path.join(OUT, f'pnt_cal_{tag}_temps.csv'), index=False)
        L.append(f'== {HA.LABEL[m]}')
        L += [f'  note: {x}' for x in notes]
        for r in R.itertuples():
            L.append(f'  {r.set} {r.variant}: direct {r.s_direct:+.3f} +- {r.se_direct:.3f} | measured {r.s_exp:+.3f} '
                     f'({r.direct_minus_exp_in_se:+.1f} SE) | reweighted {r.s_reweighted:+.3f} +- {r.se_reweighted:.3f} '
                     f'(min ESS {r.min_ess_reweighted:.3f}) | plain {r.s_plain:+.3f} +- {r.se_plain:.3f}')
        L.append('  <Rg> per temperature (nm): direct / reweighted [ESS] / plain')
        for (st, v), g in Tt.groupby(['set', 'variant']):
            L.append(f'    {st} {v}: ' + '; '.join(f'{int(x.T)} K {x.Rg_direct_nm:.2f} / {x.Rg_reweighted_nm:.2f} '
                                                 f'[{x.ess_reweighted:.3f}] / {x.Rg_plain_nm:.2f}' for x in g.itertuples()))
        L.append('')
    path = os.path.join(OUT, 'pnt_cal_summary.txt')
    open(path, 'w').write('\n'.join(L) + '\n')
    print('\n'.join(L))


if __name__ == '__main__':
    main()
