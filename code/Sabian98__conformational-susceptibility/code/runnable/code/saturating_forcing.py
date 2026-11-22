"""
saturating_forcing.py — does a saturating (burial-limited) conversion of hydration free energy into contact stickiness
fix the PNt over-contraction while keeping the Wuttke calibration and the held-out series? Exploratory (2026-10-24),
after pnt_cal_analysis.py (PNt miss real) and residual_compactness_check.py (misses track contacts per residue,
suggestively). Nothing frozen is changed; PNt is used only as a test here, the strengths are refitted on Wuttke alone.

The calibrated forcing is pairwise additive. With B_ij the attractive pair shape that the hydration change scales
(CALVADOS-2: U_AH(r; 1) - U_AH(r; 0), between -eps and 0; Mpipi-Recharged: the attractive Wang-Frenkel shape, between
-1 and 0) and w_i(T) the per-residue weight per unit strength (TEA: -dG_E,i(T); hydrophilic: du_h,i(T) - du_h,i(T0) on
the hydrophilic set; x 4.184 for Mpipi-Recharged), the energy is
    E_lin = sum_{i<j} (w_i + w_j)/2 B_ij = 1/2 sum_i w_i beta_i,     beta_i = sum_j B_ij = -eps c_i,
so each residue's hydration change is counted in proportion to c_i, its contact number in units of a full contact.
The saturating form counts it at most c_sat times:
    E_sat = -eps/2 sum_i w_i h(c_i),     h(c) = c_sat (1 - exp(-c / c_sat)),
which equals E_lin for sparse contacts and as c_sat -> infinity. c_sat is not identifiable from the four Wuttke
constructs, whose contact numbers are similar, so it is scanned (infinity to 0.5 full contacts, finer below 2).
For each c_sat and each variant (A0-type: hydrophilic set D E K R N Q S T H G; S1-type: D E K R) the TEA strength stays at
its prior (2.0 CALVADOS-2, 0.4 Mpipi-Recharged), both terms take the saturating form, and the hydrophilic strength is
refitted on the four Wuttke constructs with the routine that fitted S1 (heldout_exploratory.fit_free, curves to 320 K).
Predictions: the seven held-out series with a measured slope by reweighting the plain held-out ensembles
(heldout_analysis.slope, observables of the test), and PNt by reweighting the direct calibrated runs of the same variant
(A0 or S1 direct ensembles, bias = new energy - that run's energy), with the Kish effective sample fraction.
Checks: E_lin from the per-residue contact numbers equals the stored class and hydrophilic energies on every cached
ensemble; with c_sat = infinity and the frozen strengths, the held-out slopes equal primary.csv.

    BENDER_REANALYSIS_DIR=... python saturating_forcing.py prepare --model calvados|recharged [--workers 4]
    BENDER_REANALYSIS_DIR=... python saturating_forcing.py fit --model calvados|recharged
Output: heldout_temperature/results/exploratory/saturating/{cache_<model>/, saturating_<model>.csv, summary_<model>.txt}
"""
import sys

sys.dont_write_bytecode = True

import argparse  # noqa: E402
import glob  # noqa: E402
import multiprocessing as mp  # noqa: E402
import os  # noqa: E402
import traceback  # noqa: E402

import mdtraj as md  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import sparse  # noqa: E402

HERE = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, os.environ.get('BENDER_REANALYSIS_DIR', HERE))
sys.path.insert(0, HERE)
import heldout_analysis as HA  # noqa: E402
import heldout_exploratory as HX  # noqa: E402
import make_pnt_calibrated_inputs as MP  # noqa: E402
import pnt_cal_analysis as PA  # noqa: E402

TCal, W, fdt, L, bc = HA.TCal, HA.W, HA.fdt, HA.L, HA.bc
OUT = os.path.join(HA.OUT, 'exploratory', 'saturating')
WUTTKE_RUNS = {'calvados': f'{HA.BM}/Wuttke Runs/wuttke_CALVADOS', 'recharged': f'{HA.BM}/Wuttke Runs/wuttke_MPIPI'}
HO_RUNS = {'calvados': f'{HA.BM}/heldout_CALVADOS2', 'recharged': f'{HA.BM}/heldout_MPIPI_Recharged'}
HO_NAME = {'calvados': lambda st, ph: f'HO-{st}-pH{ph}', 'recharged': lambda st, ph: f'HO-{st}'}
PH = {'PNt-I025': 74, 'PNt-I175': 74, 'sfAFP-I175': 75, 'asyn4W26C-I030': 75, 'asyn39W66C-I030': 75,
      'AapPGR-I195': 74, 'p53N93-I118': 70}
WT_TEMPS = [278, 288, 300, 310, 320]
CSAT = [np.inf, 8.0, 4.0, 3.0, 2.0, 1.5, 1.25, 1.0, 0.75, 0.5]   # finer below 2 (2026-10-24): typical c is ~3 in CALVADOS-2, ~1.7 in Mpipi-Recharged
HYDSET = {'A0': TCal.HYD_RES, 'S1': 'DEKR'}
EPS = {'calvados': fdt.EPS_AH, 'recharged': 1.0}
UNIT = {'calvados': 1.0, 'recharged': 4.184}


def cache_dir(model):
    return os.path.join(OUT, f'cache_{model}')


# ----------------------------------------------------------------------------- per-residue contact numbers
def incidence(n, i, j):
    p = len(i)
    return sparse.csr_matrix((np.ones(2 * p), (np.r_[np.arange(p), np.arange(p)], np.r_[i, j])), shape=(p, n))


def contacts(model, xyz, seq):
    """c_i per frame (frames x N): minus the summed attractive pair shape of residue i over its non-bonded partners, in
    units of a full contact, with the pair list and shape of the stored energies."""
    xyz = np.asarray(xyz, float)
    n = len(seq)
    if model == 'calvados':
        q, lam, sig, _ = fdt.calvados_residue_params(seq, L.load_residue_params())
        i, j, s, l, shift, _ = fdt.pair_setup(q, lam, sig)
        shape = lambda r: fdt.ah_pair_energy(r, s, np.ones_like(l), shift) - fdt.ah_pair_energy(r, s, np.zeros_like(l), shift)  # noqa: E731
        box = None
    else:
        import tea_crossff as X
        tb = X.rc_tables()
        ix = np.array([tb['index'][a] for a in seq])
        i, j = np.triu_indices(n, k=2)
        sg, mu = tb['sig'][ix[j], ix[i]], tb['mu'][ix[j], ix[i]]
        rc = 3 * sg
        alpha = 2 * 3 ** (2 * mu) * (3 / (2 * (3 ** (2 * mu) - 1))) ** 3
        shape = lambda r: TCal.wf_attractive_shape(r, sg, mu, rc, alpha)  # noqa: E731
        box = HA.box_nm(n)
    M = incidence(n, i, j)
    c = np.zeros((len(xyz), n), np.float32)
    step = max(1, int(4e6 // len(i)))
    for a in range(0, len(xyz), step):
        dv = xyz[a:a + step, i] - xyz[a:a + step, j]
        if box is not None:
            dv -= box * np.round(dv / box)
        c[a:a + step] = -(M.T @ shape(np.linalg.norm(dv, axis=2)).T).T / EPS[model]
    return c


def weights(model, seq, T, hydset):
    """Per-residue weights per unit strength: (TEA, hydrophilic), kJ/mol per unit pair shape."""
    u = UNIT[model]
    wt = np.array([-W.dGE(a, T) for a in seq]) * u
    wh = TCal.hyd_weights(seq, T, hydset) * u
    return wt, wh


def energy(model, c, w, csat):
    h = c if np.isinf(csat) else csat * (1.0 - np.exp(-c / csat))
    return -0.5 * EPS[model] * (h @ w)


# ----------------------------------------------------------------------------- prepare: cache c_i for every ensemble
def tasks(model):
    out = []
    for name in HA.FOUR:
        for T in WT_TEMPS:
            out.append(('wuttke', name, T))
    for st, ph in PH.items():
        for f in sorted(glob.glob(os.path.join(HA.OUT, f'npz_{model}', f'{HO_NAME[model](st, ph)}_*K.npz'))):
            out.append(('heldout', st, int(np.load(f)['T'])))
    if have_direct(model):
        for st, (_, temps) in MP.SETS.items():
            for v in ('A0', 'S1'):
                for T in temps:
                    out.append((f'direct_{v}', st, T))
    return out


def have_direct(model):
    """True when every direct calibrated PNt run of this model has its trajectory (GAIVI_PNTCAL)."""
    return all(os.path.exists(os.path.join(PA.ROOT[model], f'{T}K', PA.NAME[model](st, v), f'{PA.NAME[model](st, v)}.dcd'))
               for st, (_, temps) in MP.SETS.items() for v in ('A0', 'S1') for T in temps)


def prepare_one(args):
    model, kind, name, T = args
    out = os.path.join(cache_dir(model), f'{kind}_{name}_{T}K.npz')
    tag = f'{kind} {name} {T} K'
    if os.path.exists(out):
        return tag, 'exists'
    try:
        check = ''
        if kind == 'wuttke':
            d = os.path.join(WUTTKE_RUNS[model], f'{T}K', name)
            if model == 'calvados':
                xyz, seq, _ = bc.load_calvados(os.path.join(d, f'{name}.dcd'), os.path.join(d, 'top.pdb'), k=L.K_FRAMES)
            else:
                import mpipi_recharged_fdt as R
                xyz, seq, _ = R.load_run(d, name)
            ref = dict(np.load(os.path.join(HA.OUT, f'wuttke_types_{model}', f'{name}_{T}K.npz')))
        elif kind == 'heldout':
            nm = HO_NAME[model](name, PH[name])
            d = os.path.join(HO_RUNS[model], f'{T}K', nm)
            seq = ''.join(r.code for r in md.load_topology(os.path.join(d, 'top.pdb')).residues)
            xyz, _, _ = HA.load(d, nm, seq)
            ref = dict(np.load(os.path.join(HA.OUT, f'npz_{model}', f'{nm}_{T}K.npz')))
        else:
            v = kind.split('_')[1]
            xyz, check = PA.load_xyz(model, name, v, T)
            seq = MP.MH.PNT
            ref = None
        c = contacts(model, xyz, seq)
        if ref is not None:                                # E_lin from c_i must equal the stored energies
            wt, wh = weights(model, seq, float(T), TCal.HYD_RES)
            et, eh = energy(model, c, wt, np.inf), energy(model, c, wh, np.inf)
            st_t = sum(ref[f'E_{k}'] for k in HA.CLS)
            dt, dh = np.abs(et - st_t).max(), np.abs(eh - ref['E_hyd']).max()
            scale = max(np.abs(st_t).max(), np.abs(ref['E_hyd']).max(), 1e-9)
            if max(dt, dh) / scale > 1e-4:
                raise ValueError(f'E_lin differs from the stored energies (TEA {dt:.2e}, hyd {dh:.2e} kJ/mol)')
            check = f'E_lin matches stored energies to {max(dt, dh) / scale:.1e}'
        os.makedirs(cache_dir(model), exist_ok=True)
        np.savez_compressed(out, c=c, seq=seq, T=T)
        return tag, f'ok, {len(c)} frames, median c {np.median(c):.2f}, 95th pct {np.percentile(c, 95):.2f}; {check}'
    except Exception as ex:
        return tag, f'FAILED {ex!r} | {traceback.format_exc()[-400:]}'


def prepare(model, workers):
    with mp.get_context('spawn').Pool(workers) as pool:
        for tag, msg in pool.imap_unordered(prepare_one, [(model, *t) for t in tasks(model)]):
            print(f'{tag}: {msg}', flush=True)


# ----------------------------------------------------------------------------- fit and predict
def cached(model, kind, name, T):
    z = np.load(os.path.join(cache_dir(model), f'{kind}_{name}_{T}K.npz'))
    return z['c'], str(z['seq'])


def with_energies(model, z, c, seq, T, variant, csat, tag='sat'):
    wt, wh = weights(model, seq, float(T), HYDSET[variant])
    z = dict(z)
    z[f'E_hyd_t{tag}'] = energy(model, c, wt, csat)
    z[f'E_hyd_h{tag}'] = energy(model, c, wh, csat)
    return z


def gam_of(model, g_h):
    return {'g_tsat': HX.S1_TEA[model], 'g_hsat': g_h}


def wuttke_units(model, variant, csat):
    base = HA.wuttke_units(model)
    units = {}
    for n, ser in base.items():
        units[n] = ('wuttke', {T: with_energies(model, z, *cached(model, 'wuttke', n, T), T, variant, csat)
                               for T, z in ser.items()})
    return units


def heldout_series(model, st, variant, csat):
    ser = HA.series(model, st)
    return {T: with_energies(model, z, *cached(model, 'heldout', st, T), T, variant, csat) for T, z in ser.items()}


def frozen_gam_sat(model, variant):
    """Strengths of the existing linear models, in this script's keys (for the c_sat = infinity check)."""
    g = HA.FROZEN[model]['hyd'] if variant == 'A0' else HX.S1_GAMMA[model]
    return gam_of(model, g)


def pnt_direct(model, st, variant, csat, g_h):
    """PNt slope by reweighting the direct runs of the same variant to the saturating model."""
    g0t, g0h = HX.S1_TEA[model], (HA.FROZEN[model]['hyd'] if variant == 'A0' else HX.S1_GAMMA[model])
    ser, ess = {}, []
    for T in MP.SETS[st][1]:
        c, seq = cached(model, f'direct_{variant}', st, T)
        wt, wh = weights(model, seq, float(T), HYDSET[variant])
        u_run = g0t * energy(model, c, wt, np.inf) + g0h * energy(model, c, wh, np.inf)
        u_new = g0t * energy(model, c, wt, csat) + g_h * energy(model, c, wh, csat)
        ser[T] = {'obs_Rg': PA_RG[(variant, st, T)], 'E_aliph': np.zeros(len(c)), 'E_hyd_d': u_new - u_run}
        ess.append(HA.ess(ser[T], {'g_d': 1.0}, float(T)))
    s, se, _ = HA.slope(ser, 'Rg', {'g_d': 1.0})
    return s, se, min(ess)


PA_RG = {}


def load_direct_rg(model):
    for st, (_, temps) in MP.SETS.items():
        for v in ('A0', 'S1'):
            for T in temps:
                xyz, _ = PA.load_xyz(model, st, v, T)
                PA_RG[(v, st, T)] = HA.observables(xyz, 'PNt')['Rg']


def fit(model):
    s_exp, _ = HA.experimental()
    HX.H.OUT = HA.OUT
    prim = pd.read_csv(os.path.join(HA.OUT, f'primary_{model}', 'primary.csv'))
    prim = prim[prim.variant == 'A0'].set_index('set')
    L_ = [f'Saturating hydration forcing, {HA.LABEL[model]}: strengths refitted on the four Wuttke constructs, '
          'predictions in % per K', '']
    # check: c_sat = infinity with the frozen A0 strengths reproduces the primary held-out slopes
    for st in PH:
        s, _, _ = HA.slope(heldout_series(model, st, 'A0', np.inf), HA.OBS_OF[st], frozen_gam_sat(model, 'A0'))
        if abs(s - prim.loc[st, 's_model']) > 1e-6:
            raise ValueError(f'{st}: linear check {s:.6f} differs from primary.csv {prim.loc[st, "s_model"]:.6f}')
    L_.append('check passed: with c_sat = infinity and the frozen A0 strengths, all seven held-out slopes equal primary.csv')
    direct = have_direct(model)
    if direct:
        load_direct_rg(model)
    rows = []
    for variant in ('A0', 'S1'):
        for csat in CSAT:
            wu = wuttke_units(model, variant, csat)
            (g,), rms = HX.fit_free(lambda x: gam_of(model, x), [0.3], wu, None, lambda p: np.array([p[0] / 10.0]),
                                    model=model)
            gam = gam_of(model, g)
            row = {'model': HA.LABEL[model], 'variant': variant, 'c_sat': csat, 'gamma_h': g, 'wuttke_rms': rms}
            for st in PH:
                ser = heldout_series(model, st, variant, csat)
                s, se, Ts = HA.slope(ser, HA.OBS_OF[st], gam)
                row[f'{st}'] = s
                row[f'{st}_ess'] = min(HA.ess(ser[T], gam, float(T)) for T in Ts)
            if direct:
                for st in MP.SETS:
                    s, se, e = pnt_direct(model, st, variant, csat, g)
                    row[f'{st}_direct'] = s
                    row[f'{st}_direct_se'] = se
                    row[f'{st}_direct_ess'] = e
            contract = ['asyn4W26C-I030', 'asyn39W66C-I030', 'p53N93-I118']
            pnt = [f'{st}_direct' if direct else st for st in MP.SETS]
            row['MAE_aSyn_p53'] = np.mean([abs(row[k] - s_exp[k]) for k in contract])
            row['MAE_PNt'] = np.mean([abs(row[k] - s_exp[k.replace('_direct', '')]) for k in pnt])
            row['MAE_controls'] = np.mean([abs(row[k] - s_exp[k]) for k in HA.CONTROLS])
            rows.append(row)
            print(f'{variant} c_sat {csat}: gamma_h {g:.3f}, Wuttke rms {rms:.2f}', flush=True)
    R = pd.DataFrame(rows)
    os.makedirs(OUT, exist_ok=True)
    R.to_csv(os.path.join(OUT, f'saturating_{model}.csv'), index=False)
    L_ += ['', f'measured: ' + ', '.join(f'{k} {s_exp[k]:+.2f}' for k in PH), '']
    sfx = '_direct' if direct else ''
    for r in R.to_dict('records'):
        pn = ' '.join(f'{st} {r[st + sfx]:+.2f}' for st in MP.SETS)
        L_.append(f"{r['variant']} c_sat {r['c_sat']:>4}: gamma_h {r['gamma_h']:.3f} (Wuttke rms {r['wuttke_rms']:.2f}) | "
                  f"PNt {pn} | aSyn {r['asyn4W26C-I030']:+.2f} {r['asyn39W66C-I030']:+.2f} p53 {r['p53N93-I118']:+.2f} | "
                  f"sfAFP {r['sfAFP-I175']:+.2f} Aap {r['AapPGR-I195']:+.2f} | MAE PNt {r['MAE_PNt']:.2f}, "
                  f"aSyn+p53 {r['MAE_aSyn_p53']:.2f}, controls {r['MAE_controls']:.2f}")
    L_ += ['', 'effective sample fraction (min over temperatures): held-out reweighting from plain'
           + ('; PNt from the direct runs' if direct else '; PNt from plain (direct runs not yet available)')]
    for r in R.to_dict('records'):
        e = ', '.join(f"{st} {r[st + '_ess']:.3f}" for st in PH)
        if direct:
            e += ' | ' + ', '.join(f"{st} direct {r[st + '_direct_ess']:.3f}" for st in MP.SETS)
        L_.append(f"  {r['variant']} c_sat {r['c_sat']:>4}: {e}")
    open(os.path.join(OUT, f'summary_{model}.txt'), 'w').write('\n'.join(L_) + '\n')
    print('\n'.join(L_))



if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('mode', choices=['prepare', 'fit'])
    ap.add_argument('--model', choices=['calvados', 'recharged'], required=True)
    ap.add_argument('--workers', type=int, default=4)
    a = ap.parse_args()
    prepare(a.model, a.workers) if a.mode == 'prepare' else fit(a.model)
