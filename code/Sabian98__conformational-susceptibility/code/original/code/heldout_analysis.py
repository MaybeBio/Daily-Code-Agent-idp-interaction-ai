"""
heldout_analysis.py — the analysis fixed by heldout_temperature/PREREGISTRATION.md (registered 2026-10-11): does the
calibrated hydration term predict independent temperature series of disordered proteins (GAIVI_HELDOUT runs)?
Written and checked on the PAGE4 direct runs before any held-out trajectory existed.

Modes
  prepare --model calvados|recharged [--root DIR] [--workers N]
      every run of GAIVI_HELDOUT/<model>/runs.tsv found under DIR/<T>K/<name>/: 3,000 evenly spaced frames of the second
      half (DCD header checked), per-frame Rg, Kirkwood 1/r mean (Rh = 1/<mean 1/r_ij>), loop distances (alpha-synuclein
      4-26 and 39-66), the TEA class energies and the hydrophilic energies by residue group and by residue type
      (tea_calibrate's own functions, HYD_GROUPS extended by one entry per type) -> results/npz_<model>/<name>_<T>K.npz
  wuttke-types --model M --runs ROOT --fasta FASTA --out DIR     (on GAIVI, where the Wuttke trajectories are)
      the Wuttke calibration npz again, with per-type hydrophilic energies (needed for alternative A2)
  guinier      Rg of PNt and sfAFP from the raw SAXS profiles (Guinier, q Rg <= 1.0, iterated) -> results/guinier.csv
  primary      frozen model A0 and plain model, both force fields: relative slopes s = 100 d ln<A>/dT (least squares
               over the simulated temperatures, 10-block jackknife with the same block removed at every T; alpha-synuclein
               uses the mean loop distance, Amendment 1), criteria
               C1 (sign), C2 (magnitude), C3 (controls) -> results/primary.csv, primary_summary.txt
  alternatives A1-A4 fitted on the training units (Wuttke four constructs, PNt-I025, PNt-I175, sfAFP-I175; PNt-I525
               dropped by Amendment 2) with the tea_calibrate loss, units weighted equally; leave-one-unit-out
               CV picks one (ties within 10% -> fewer parameters); it is then tested once on alpha-synuclein, p53 and
               Aap -> results/alternatives_*.csv, alternatives_summary.txt
  selftest     on the PAGE4 direct runs only: per-type energies sum to the group and total energies; the reweighted
               <Rg> curves equal page4_direct's curves.csv; estimators and criteria run end to end on a mock table
    BENDER_REANALYSIS_DIR=... python heldout_analysis.py <mode> ...      (openmpipi env)
"""
import sys

sys.dont_write_bytecode = True

import argparse  # noqa: E402
import glob  # noqa: E402
import math  # noqa: E402
import multiprocessing as mp  # noqa: E402
import os  # noqa: E402
import struct  # noqa: E402
import traceback  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.optimize import least_squares  # noqa: E402

HERE = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, os.environ.get('BENDER_REANALYSIS_DIR', HERE))
sys.path.insert(0, HERE)
import tea_calibrate as TCal  # noqa: E402
import tea_reweight as W  # noqa: E402

fdt, bc, L = TCal.fdt, TCal.bc, TCal.L
BM = '/Volumes/User Homes/Nu Project/BENDER-MT'
HO = f'{BM}/heldout_temperature'
OUT = f'{HO}/results'
PKG = os.path.expanduser('~/Downloads/Scripts/GAIVI_HELDOUT')
WU = f'{BM}/Traj Level Analysis/Wuttke2014'
ROOTS = {'calvados': f'{BM}/heldout_CALVADOS', 'recharged': f'{BM}/heldout_MPIPI'}
PKGDIR = {'calvados': 'calvados', 'recharged': 'mpipi'}
LABEL = {'calvados': 'CALVADOS-2', 'recharged': 'Mpipi-Recharged'}
K = 3000
T0 = 300.0
R_GAS = fdt.R_GAS
CLS = TCal.CLS
FROZEN = {'calvados': {'tea': 2.0, 'hyd': 0.358}, 'recharged': {'tea': 0.4, 'hyd': 0.209}}
A4_PRIOR = {'calvados': 0.5, 'recharged': 0.3}
TYPES = list(TCal.HYD_RES)                                   # D E K R N Q S T H G
GROUPS = dict(TCal.HYD_GROUPS)                               # DE, K, R, pol (NQSTH), G
LOOPS = {'asyn4W26C': ('loop', 3, 25), 'asyn39W66C': ('loop', 38, 65)}
FOUR = TCal.FOUR
TRAIN_BAXA = ['PNt-I025', 'PNt-I175', 'sfAFP-I175']
TEST = ['asyn4W26C-I030', 'asyn39W66C-I030', 'p53N93-I118', 'AapPGR-I195']
CONTRACT = ['PNt-I175', 'PNt-I025', 'asyn4W26C-I030', 'asyn39W66C-I030', 'p53N93-I118']
CONTROLS = ['sfAFP-I175', 'AapPGR-I195']
OBS_OF = {'PNt-I025': 'Rg', 'PNt-I175': 'Rg', 'PNt-I525': 'Rg', 'sfAFP-I175': 'Rg', 'asyn4W26C-I030': 'mean_loop',
          'asyn39W66C-I030': 'mean_loop', 'AapPGR-I195': 'Rh', 'p53N93-I118': 'Rh'}   # Amendment 1: mean loop distance


# ----------------------------------------------------------------------------- experimental reference (pre-registration §4)
def lnslope(T, y):
    return 100 * np.polyfit(np.asarray(T, float), np.log(np.asarray(y, float)), 1)[0]


def experimental():
    """s_exp (% per K) and, for the training units, the ln end-point ratio with its error."""
    asyn = {'4-26': [33.8, 33.8, 32.3, 30.9, 29.3], '39-66': [38.9, 37.6, 34.8, 33.2, 32.9]}
    s = {'PNt-I175': 100 * math.log(46.98 / 52.72) / 48, 'PNt-I025': 100 * math.log(41.83 / 47.07) / 48,
         'sfAFP-I175': 100 * math.log(24.29 / 23.69) / 48,
         'asyn4W26C-I030': lnslope([5, 10, 20, 30, 40], asyn['4-26']),
         'asyn39W66C-I030': lnslope([5, 10, 20, 30, 40], asyn['39-66']),
         'AapPGR-I195': 100 * -0.015 / 38.39, 'p53N93-I118': 100 * -0.076 / 32.8}
    ends = {'PNt-I175': (52.72, 0.14, 46.98, 0.14), 'PNt-I025': (47.07, 0.12, 41.83, 0.21),
            'sfAFP-I175': (23.69, 0.12, 24.29, 0.08)}
    ratio = {k: (math.log(b / a), math.hypot(ea / a, eb / b)) for k, (a, ea, b, eb) in ends.items()}
    return s, ratio


# ----------------------------------------------------------------------------- per-run preparation
def runs(model):
    r = pd.read_csv(os.path.join(PKG, PKGDIR[model], 'runs.tsv'), sep='\t')
    fa = {}
    for block in open(os.path.join(PKG, PKGDIR[model], 'forms.fasta')).read().split('>')[1:]:
        h, *q = block.split('\n')
        fa[h.strip()] = ''.join(q).strip()
    r['seq'] = r.form.map(fa)
    return r


def box_nm(n):
    return float(int(n / 5) + 20 if n > 300 else 50)


def load(run_dir, name, seq):
    dcd = os.path.join(run_dir, f'{name}.dcd')
    xyz_all, top_seq, info = bc.load_calvados(dcd, os.path.join(run_dir, 'top.pdb'), k=None, production=False)
    n = len(xyz_all)
    nset = struct.unpack('<i', open(dcd, 'rb').read(12)[8:12])[0]
    if nset != n:
        raise ValueError(f'DCD header lists {nset} frames but the file holds {n}; rerun')
    if top_seq != seq:
        raise ValueError('top.pdb sequence differs from forms.fasta')
    idx = n // 2 + bc.even_indices(n - n // 2, K)
    xyz = xyz_all[idx]
    if bc.max_bond(xyz) > 1.0:
        xyz = bc.unwrap_chain(xyz, box_nm(len(seq)))
    return xyz, idx, n


def energies(model, xyz, seq, T):
    """TEA classes, hydrophilic total, groups and per-type parts, from tea_calibrate's own functions."""
    saved = TCal.HYD_GROUPS
    TCal.HYD_GROUPS = dict(saved, **{f'aa{a}': a for a in TYPES})
    try:
        if model == 'calvados':
            E = W.class_energies(xyz, seq, T)
            E.update(TCal.calvados_hyd_energies(xyz, seq, T))
        else:
            E = TCal.recharged_class_energies(xyz, seq, T, box_nm(len(seq)))
    finally:
        TCal.HYD_GROUPS = saved
    return E


def observables(xyz, form):
    x = np.asarray(xyz, float)
    rg = np.sqrt(((x - x.mean(1, keepdims=True)) ** 2).sum(2).mean(1))
    i, j = np.triu_indices(x.shape[1], k=1)
    inv = np.zeros(len(x))
    step = max(1, int(4e6 // len(i)))
    for a in range(0, len(x), step):
        inv[a:a + step] = (1.0 / np.linalg.norm(x[a:a + step, i] - x[a:a + step, j], axis=2)).mean(1)
    out = {'Rg': rg, 'inv': inv}
    if form in LOOPS:
        _, p, q = LOOPS[form]
        out['loop'] = np.linalg.norm(x[:, p] - x[:, q], axis=1)
    return out


def _prepare_one(args):
    model, row, root, od = args
    name, T = row['name'], int(row['temp'])
    out = os.path.join(od, f'{name}_{T}K.npz')
    if os.path.exists(out):
        return name, T, 'exists'
    try:
        xyz, idx, n = load(os.path.join(root, f'{T}K', name), name, row['seq'])
        E = energies(model, xyz, row['seq'], float(T))
        obs = observables(xyz, row['form'])
        np.savez_compressed(out, T=T, form=row['form'], set=row['set'], frame_index=idx, n_frames_total=n,
                            **{f'obs_{k}': v for k, v in obs.items()}, **{f'E_{k}': v for k, v in E.items()})
        return name, T, 'ok'
    except Exception as ex:  # noqa: BLE001
        return name, T, f'FAILED {ex!r} | {traceback.format_exc()[-300:]}'.replace('\n', ' ')


def prepare(model, root, workers):
    od = os.path.join(OUT, f'npz_{model}')
    os.makedirs(od, exist_ok=True)
    r = runs(model)
    r = r[[os.path.exists(os.path.join(root, f'{int(t)}K', n, f'{n}.dcd')) for n, t in zip(r['name'], r.temp)]]
    print(f'{model}: {len(r)} runs found under {root}')
    with mp.get_context('spawn').Pool(workers) as pool:
        for name, T, msg in pool.imap_unordered(_prepare_one, [(model, x, root, od) for x in r.to_dict('records')]):
            print(f'{name} {T} K: {msg}', flush=True)


def wuttke_types(model, root, fasta, out):
    """tea_calibrate.prepare for the Wuttke constructs, sequentially, with per-type hydrophilic energies."""
    os.makedirs(out, exist_ok=True)
    saved = TCal.HYD_GROUPS
    TCal.HYD_GROUPS = dict(saved, **{f'aa{a}': a for a in TYPES})
    try:
        dyes = TCal.read_fasta_dyes(fasta)
        for name, v in dyes.items():
            for T in TCal.WUTTKE_TEMPS:
                print(TCal._prepare_one(('calvados' if model == 'calvados' else 'recharged', root, name, v, T, out)), flush=True)
    finally:
        TCal.HYD_GROUPS = saved


# ----------------------------------------------------------------------------- estimators
def series(model, st):
    """{T: dict of arrays} for one set."""
    out = {}
    for f in glob.glob(os.path.join(OUT, f'npz_{model}', '*.npz')):
        z = np.load(f)
        if str(z['set']) == st:
            out[int(z['T'])] = {k: z[k] for k in z.files}
    return dict(sorted(out.items()))


def hyd_energy(z, variant):
    """Hydrophilic energy (per unit strength) of a variant: frozen set, without G, without G/S/T."""
    if variant == 'all':
        return z['E_hyd']
    if variant == 'noG':
        return z['E_hyd'] - z['E_hyd_G']
    if variant == 'noGST':
        return z['E_hyd'] - z['E_hyd_G'] - z['E_hyd_aaS'] - z['E_hyd_aaT']
    raise ValueError(variant)


def bias(z, gam, T):
    """Reduced bias energy sum_c gamma_c E_c / RT for a parameter dict: tea (all TEA classes), hyd_<variant>, or the
    five groups g_DE, g_K, g_R, g_pol, g_G."""
    u = np.zeros(len(z['E_aliph']))
    if gam.get('tea'):
        u += gam['tea'] * sum(z[f'E_{c}'] for c in CLS)
    for k, g in gam.items():
        if g == 0 or k == 'tea':
            continue
        if k.startswith('hyd_'):
            u += g * hyd_energy(z, k[4:])
        elif k.startswith('g_'):
            u += g * z[f'E_hyd_{k[2:]}']
    return u / (R_GAS * T)


def weights(z, gam, T, keep=None):
    u = bias(z, gam, T)
    if keep is not None:
        u = u[keep]
    w = np.exp(-(u - u.min()))
    return w / w.sum()


def kde_peak(r, w, bw=0.1):
    grid = np.arange(max(0.2, r.min() - 0.5), r.max() + 0.5, 0.005)
    dens = np.zeros_like(grid)
    for a in range(0, len(r), 500):
        dens += (w[a:a + 500, None] * np.exp(-0.5 * ((grid[None, :] - r[a:a + 500, None]) / bw) ** 2)).sum(0)
    return grid[np.argmax(dens)]


def ensemble_value(z, obs, gam, T, keep=None):
    w = weights(z, gam, T, keep)
    sel = (lambda a: a[keep]) if keep is not None else (lambda a: a)
    if obs == 'Rg':
        return float((w * sel(z['obs_Rg'])).sum())
    if obs == 'Rh':
        return float(1.0 / (w * sel(z['obs_inv'])).sum())
    if obs == 'peak':
        return float(kde_peak(sel(z['obs_loop']), w))
    if obs == 'mean_loop':
        return float((w * sel(z['obs_loop'])).sum())
    raise ValueError(obs)


def ess(z, gam, T):
    w = weights(z, gam, T)
    return float(1.0 / (w ** 2).sum() / len(w))


def slope(ser, obs, gam, nb=10):
    """s = 100 d ln<A>/dT by least squares over the temperatures, and its jackknife SE (same block removed at every T)."""
    Ts = np.array(sorted(ser))
    val = lambda keep: lnslope(Ts, [ensemble_value(ser[T], obs, gam, float(T), keep) for T in Ts])  # noqa: E731
    s = val(None)
    n = len(ser[Ts[0]]['obs_Rg'])
    blocks = np.array_split(np.arange(n), nb)
    v = []
    for b in blocks:
        keep = np.ones(n, bool)
        keep[b] = False
        v.append(val(keep))
    v = np.array(v)
    return s, float(np.sqrt((nb - 1) / nb * ((v - v.mean()) ** 2).sum())), Ts


def end_ratio(ser, obs, gam, nb=10):
    Ts = sorted(ser)
    lo, hi = ser[Ts[0]], ser[Ts[-1]]
    f = lambda keep: math.log(ensemble_value(hi, obs, gam, float(Ts[-1]), keep) / ensemble_value(lo, obs, gam, float(Ts[0]), keep))  # noqa: E731
    n = len(lo['obs_Rg'])
    v = []
    for b in np.array_split(np.arange(n), nb):
        keep = np.ones(n, bool)
        keep[b] = False
        v.append(f(keep))
    v = np.array(v)
    return f(None), float(np.sqrt((nb - 1) / nb * ((v - v.mean()) ** 2).sum()))


def criteria(tab):
    """tab: DataFrame indexed by set with s_model, se, s_exp. Pre-registration §5."""
    c1 = all((tab.loc[s, 's_model'] + 2 * tab.loc[s, 'se'] < 0) for s in CONTRACT if s in tab.index)
    ratios = np.array([tab.loc[s, 's_model'] / tab.loc[s, 's_exp'] for s in CONTRACT if s in tab.index])
    positive = len(ratios) > 0 and bool((ratios > 0).all())
    c2 = positive and 0.67 <= np.exp(np.log(ratios).mean()) <= 1.5 and bool(((ratios >= 0.5) & (ratios <= 2)).all())
    c3 = all(abs(tab.loc[s, 's_model'] - tab.loc[s, 's_exp']) <= 0.10 for s in CONTROLS if s in tab.index)
    complete = all(s in tab.index for s in CONTRACT + CONTROLS)
    gm = float(np.exp(np.log(ratios).mean())) if len(ratios) and (ratios > 0).all() else float('nan')
    return {'C1_sign': c1, 'C2_magnitude': c2, 'C3_controls': c3, 'geo_mean_ratio': gm,
            'PASS': bool(c1 and c2 and c3 and complete), 'complete': complete}


# ----------------------------------------------------------------------------- primary analysis (A0)
def frozen_gam(model, plain=False):
    return {'tea': 0.0, 'hyd_all': 0.0} if plain else {'tea': FROZEN[model]['tea'], 'hyd_all': FROZEN[model]['hyd']}


def primary():
    s_exp, _ = experimental()
    rows, summ = [], ['Held-out test, frozen model A0 (pre-registration §5); s in % per K']
    for model in ('calvados', 'recharged'):
        tab = []
        for st, obs in OBS_OF.items():
            ser = series(model, st)
            if len(ser) < 2:
                continue
            for lab, gam in (('plain', frozen_gam(model, True)), ('A0', frozen_gam(model))):
                s, se, Ts = slope(ser, obs, gam)
                e = min(ess(ser[T], gam, float(T)) for T in Ts)
                row = {'model': LABEL[model], 'variant': lab, 'set': st, 'observable': obs, 'temps': ','.join(map(str, Ts)),
                       's_model': s, 'se': se, 's_exp': s_exp.get(st, np.nan), 'min_ess': e}
                if obs == 'mean_loop':
                    row['s_peak_secondary'] = slope(ser, 'peak', gam)[0]
                if obs == 'Rh':
                    row['s_Rg'] = slope(ser, 'Rg', gam)[0]
                tab.append(row)
        if not tab:
            summ.append(f'{LABEL[model]}: no prepared runs')
            continue
        t = pd.DataFrame(tab)
        rows.append(t)
        for lab in ('plain', 'A0'):
            tt = t[t.variant == lab].set_index('set')
            c = criteria(tt)
            summ += ['', f'{LABEL[model]}, {lab}: ' + ', '.join(f'{k} {v}' for k, v in c.items()),
                     tt[['observable', 's_model', 'se', 's_exp', 'min_ess']].round(3).to_string()]
    if rows:
        os.makedirs(OUT, exist_ok=True)
        pd.concat(rows).to_csv(os.path.join(OUT, 'primary.csv'), index=False)
    open(os.path.join(OUT, 'primary_summary.txt'), 'w').write('\n'.join(summ) + '\n')
    print('\n'.join(summ))


# ----------------------------------------------------------------------------- alternatives (A1-A4)
ALTS = {'A0': ([], {}), 'A1': (['hyd_noG'], {}), 'A2': (['hyd_noGST'], {}), 'A3': (['hyd_noG', 'tea'], {}),
        'A4': (['g_DE', 'g_K', 'g_R', 'g_pol', 'g_G'], {})}


def alt_start(model, alt):
    f = FROZEN[model]
    if alt == 'A0':
        return {'tea': f['tea'], 'hyd_all': f['hyd']}
    base = {'tea': f['tea']}
    if alt in ('A1', 'A3'):
        base['hyd_noG'] = f['hyd']
    if alt == 'A2':
        base['hyd_noGST'] = f['hyd']
    if alt == 'A4':
        base.update({g: f['hyd'] for g in ALTS['A4'][0]})
    return base


def prior_resid(model, alt, gam):
    f = FROZEN[model]
    r = []
    for k in ALTS[alt][0]:
        if k == 'tea':
            r.append((gam[k] - f['tea']) / 1.0)                      # tea_calibrate tau for the TEA classes
        elif k.startswith('hyd_'):
            r.append((gam[k] - 0.0) / 10.0)                          # tea_calibrate tau_hyd, prior 0
        else:
            r.append((gam[k] - f['hyd']) / A4_PRIOR[model])          # pre-registration §6 (A4)
    return np.array(r)


def wuttke_units(model):
    """{construct: series {T: z}} with A = dye-segment Rg, from the calibration npz (per-type energies if present)."""
    d = os.path.join(WU, 'calibration', model)
    dt = os.path.join(OUT, f'wuttke_types_{model}')
    out = {}
    for n in FOUR:
        ser = {}
        for f in sorted(glob.glob(os.path.join(dt if os.path.isdir(dt) else d, f'{n}_*K.npz'))):
            z = np.load(f)
            T = int(z['T'])
            if T > 320:
                continue
            zz = {k: z[k] for k in z.files}
            zz['obs_Rg'] = zz['A']
            ser[T] = zz
        out[n] = ser
    return out


_WT = {}


def wuttke_target(unit, Ts):
    """Measured relative curve at Ts as in the calibration (cubic interpolation, temperatures up to 320 K)."""
    if (unit, Ts) not in _WT:
        ex = pd.read_csv(f'{WU}/wuttke2014_fig2b_Rg_digitized.csv')
        TCal.TMAX = 320
        _WT[(unit, Ts)] = TCal.exp_curve(ex[ex.construct == unit].sort_values('T_K'), np.array(Ts))
    return _WT[(unit, Ts)]


def unit_resid(model, unit, kind, data, gam, ratio_exp):
    """Standardized residuals of one training unit, scaled so that the unit's sum of squares is its mean."""
    if kind == 'wuttke':
        Ts = np.array(sorted(data))
        re_, ee, inside = wuttke_target(unit, tuple(Ts))
        m = {T: ensemble_value(data[T], 'Rg', gam, float(T)) for T in Ts}
        se = {T: TCal.block_se(data[T]['obs_Rg']) / m[300] / np.sqrt(max(ess(data[T], gam, float(T)), 1e-3)) for T in Ts}
        keep = (Ts != 300) & inside
        r = np.array([(m[T] / m[300] - re_[k]) / np.hypot(ee[k], se[T]) for k, T in enumerate(Ts)])[keep]
    else:
        val, sm = end_ratio(data, 'Rg', gam)
        ve, se_ = ratio_exp[unit]
        r = np.array([(val - ve) / math.hypot(se_, sm)])
    return r / np.sqrt(len(r))


def fit_alt(model, alt, units, ratio_exp):
    keys = ALTS[alt][0]
    start = alt_start(model, alt)

    def resid(p):
        gam = dict(start, **dict(zip(keys, p)))
        r = [unit_resid(model, u, kind, data, gam, ratio_exp) for u, (kind, data) in units.items()]
        r.append(prior_resid(model, alt, gam))
        return np.concatenate(r)

    if not keys:
        return start
    sol = least_squares(resid, x0=[start[k] for k in keys], method='trf')
    return dict(start, **dict(zip(keys, sol.x)))


def alternatives():
    s_exp, ratio_exp = experimental()
    summ, rows = ['Alternatives A1-A4 (pre-registration §6-7)'], []
    for model in ('calvados', 'recharged'):
        units = {n: ('wuttke', ser) for n, ser in wuttke_units(model).items()}
        for st in TRAIN_BAXA:
            ser = series(model, st)
            if len(ser) >= 2:
                units[st] = ('baxa', ser)
        # Amendment 2: PNt-I525 is secondary only and not a training unit (its 7 C Guinier point is inconsistent)
        has_types = all('E_hyd_aaS' in next(iter(d.values())) for k, (kind, d) in units.items() if d)
        alts = [a for a in ALTS if a != 'A0' and (a != 'A2' or has_types)]
        if 'A2' not in alts:
            summ.append(f'{LABEL[model]}: A2 needs per-type energies for the Wuttke runs (mode wuttke-types on GAIVI); skipped')
        cv = {}
        for alt in alts:
            errs = []
            for out in units:
                tr = {u: v for u, v in units.items() if u != out}
                gam = fit_alt(model, alt, tr, ratio_exp)
                kind, data = units[out]
                errs.append(float((unit_resid(model, out, kind, data, gam, ratio_exp) ** 2).sum()))
            cv[alt] = float(np.mean(errs))
        best = min(cv, key=cv.get)
        near = [a for a in cv if cv[a] <= 1.10 * cv[best]]
        chosen = min(near, key=lambda a: (len(ALTS[a][0]), cv[a]))
        gam = fit_alt(model, chosen, units, ratio_exp)
        tab = []
        for st in TEST:
            ser = series(model, st)
            if len(ser) >= 2:
                s, se, _ = slope(ser, OBS_OF[st], gam)
                tab.append({'set': st, 's_model': s, 'se': se, 's_exp': s_exp[st]})
        t = pd.DataFrame(tab).set_index('set') if tab else pd.DataFrame()
        c1 = all((t.loc[s, 's_model'] + 2 * t.loc[s, 'se'] < 0) for s in ['asyn4W26C-I030', 'asyn39W66C-I030', 'p53N93-I118'] if s in t.index)
        ratios = np.array([t.loc[s, 's_model'] / t.loc[s, 's_exp'] for s in ['asyn4W26C-I030', 'asyn39W66C-I030', 'p53N93-I118'] if s in t.index])
        c2 = bool(len(ratios) == 3 and (ratios > 0).all() and 0.67 <= np.exp(np.log(ratios).mean()) <= 1.5 and ((ratios >= 0.5) & (ratios <= 2)).all())
        c3 = 'AapPGR-I195' in t.index and abs(t.loc['AapPGR-I195', 's_model'] - t.loc['AapPGR-I195', 's_exp']) <= 0.10
        rows.append(pd.DataFrame([{'model': LABEL[model], 'alt': a, 'cv_error': v, 'chosen': a == chosen} for a, v in cv.items()]))
        summ += ['', f'{LABEL[model]}: units {list(units)}', f'  CV error {cv}; chosen {chosen}; fitted {gam}',
                 t.round(3).to_string() if len(t) else '  no test runs',
                 f'  test: C1 {c1}, C2 {c2}, C3 {c3}, PASS {bool(c1 and c2 and c3)}']
    if rows:
        pd.concat(rows).to_csv(os.path.join(OUT, 'alternatives_cv.csv'), index=False)
    open(os.path.join(OUT, 'alternatives_summary.txt'), 'w').write('\n'.join(summ) + '\n')
    print('\n'.join(summ))


# ----------------------------------------------------------------------------- Guinier (experimental side)
def guinier():
    raw = f'{HO}/raw/saxs_github'
    files = {'PNt-I025': 'PNt_0-KCl_T-dep', 'PNt-I175': 'PNt_150mM-KCl_T-dep', 'PNt-I525': 'PNt_500mM-KCl_T-dep',
             'sfAFP-I175': 'sfafpSNAC_150mM-KCl_T-dep'}
    rows = []
    for st, d in files.items():
        for f in sorted(glob.glob(os.path.join(raw, d, '*.dat'))):
            q, I, e = np.loadtxt(f, comments='#', unpack=True)
            tc = int(os.path.basename(f).split('_')[-1].split('C')[0])
            rg, se, ok, npts = 30.0, np.nan, False, 0
            for _ in range(30):
                m = (q * rg <= 1.0) & (I > 0)
                if m.sum() < 10:
                    break
                A = np.vstack([np.ones(m.sum()), q[m] ** 2]).T
                wgt = I[m] / e[m]
                coef, *_ = np.linalg.lstsq(A * wgt[:, None], np.log(I[m]) * wgt, rcond=None)
                if coef[1] >= 0:
                    break
                new = math.sqrt(-3 * coef[1])
                cov = np.linalg.pinv((A * wgt[:, None]).T @ (A * wgt[:, None]))
                se = 1.5 / new * math.sqrt(cov[1, 1])
                npts = int(m.sum())
                if abs(new - rg) < 1e-3:
                    rg, ok = new, True
                    break
                rg = new
            rows.append({'set': st, 'file': os.path.basename(f), 'T_C': tc, 'Rg_A': rg, 'se_A': se, 'n_points': npts,
                         'qmax_Rg': float(q[(q * rg <= 1.0)].max() * rg) if ok else np.nan, 'ok': ok})
    g = pd.DataFrame(rows).sort_values(['set', 'T_C', 'file'])
    g = g.groupby(['set', 'T_C'], as_index=False).agg({'Rg_A': 'mean', 'se_A': 'mean', 'n_points': 'min', 'ok': 'all', 'file': 'first'})
    os.makedirs(OUT, exist_ok=True)
    g.to_csv(os.path.join(OUT, 'guinier.csv'), index=False)
    print(g.round(3).to_string(index=False))
    for st, (a, b) in {'PNt-I175': (52.72, 46.98), 'PNt-I025': (47.07, 41.83), 'sfAFP-I175': (23.69, 24.29)}.items():
        x = g[g.set == st].sort_values('T_C')
        if len(x):
            print(f'{st}: Guinier {x.Rg_A.iloc[0]:.2f} -> {x.Rg_A.iloc[-1]:.2f} A; MFF (paper) {a} -> {b} A')


# ----------------------------------------------------------------------------- selftest on PAGE4 (no held-out data)
def selftest():
    import page4_direct as PD
    msgs = []
    for model in ('calvados', 'recharged'):
        r = PD.runs_table(model)
        r = r[(r.form == 'WT') & (r.set == 'heat') & r.present]
        cv = pd.read_csv(os.path.join(PD.out_dir(model), 'curves.csv'))
        ref = cv[(cv.form == 'WT') & (cv.curve == 'calibrated')].iloc[0]
        ser = {}
        for _, row in r.iterrows():
            T = int(row.temp)
            xyz, idx, n = load(row['dir'], row['name'], row['seq'] if model == 'recharged' else row['wt_seq'])
            pz = np.load(os.path.join(PD.out_dir(model), 'npz', f"{row['name']}_{T}K.npz"))
            assert np.array_equal(idx, pz['frame_index']), 'frames differ from page4_direct'
            E = energies(model, xyz, row['seq'] if model == 'recharged' else row['wt_seq'], float(T))
            s_types = sum(E[f'hyd_aa{a}'] for a in TYPES)
            s_groups = sum(E[f'hyd_{g}'] for g in GROUPS)
            d1 = np.abs(s_types - E['hyd']).max() / max(1e-12, np.abs(E['hyd']).max())
            d2 = np.abs(s_groups - E['hyd']).max() / max(1e-12, np.abs(E['hyd']).max())
            d3 = np.abs(E['hyd_pol'] - sum(E[f'hyd_aa{a}'] for a in 'NQSTH')).max() / max(1e-12, np.abs(E['hyd_pol']).max())
            obs = observables(xyz, 'asyn4W26C')
            z = {'obs_Rg': obs['Rg'], 'obs_inv': obs['inv'], 'obs_loop': obs['loop'], **{f'E_{k}': v for k, v in E.items()}}
            ser[T] = z
            msgs.append(f'{model} WT {T} K: per-type sum vs hyd {d1:.1e}, groups vs hyd {d2:.1e}, pol vs N+Q+S+T+H {d3:.1e}; '
                        f'Rh/Rg {1 / obs["inv"].mean() / obs["Rg"].mean():.3f}')
        gam = frozen_gam(model)
        Ts = sorted(ser)
        m = {T: ensemble_value(ser[T], 'Rg', gam, float(T)) for T in Ts}
        dev = max(abs(m[T] / m[300] - ref[f'rel_{T}K']) for T in Ts)
        msgs.append(f'{model}: reweighted <Rg>(T)/<Rg>(300) against page4_direct curves.csv, max |d| {dev:.2e}')
        for obs in ('Rg', 'Rh', 'peak', 'mean_loop'):
            s, se, _ = slope(ser, obs, gam)
            msgs.append(f'{model}: mock slope of {obs} (PAGE4 WT, 18 bead pairs stand-in) {s:+.3f} +- {se:.3f} % per K')
        hyd_variants = {v: hyd_energy(ser[300], v) for v in ('all', 'noG', 'noGST')}
        msgs.append(f'{model}: hydrophilic energy at 300 K is zero by construction: '
                    + ', '.join(f'{k} {np.abs(v).max():.1e}' for k, v in hyd_variants.items()))
        s320 = {v: float(np.abs(hyd_energy(ser[320], v)).mean()) for v in ('all', 'noG', 'noGST')}
        msgs.append(f'{model}: mean |E_hyd| at 320 K by variant {s320} (all >= noG >= noGST expected)')
    mock = pd.DataFrame({'s_model': [-0.25, -0.24, -0.40, -0.45, -0.20, 0.0, -0.05], 'se': [0.02] * 7,
                         's_exp': [-0.240, -0.246, -0.424, -0.503, -0.232, 0.052, -0.039]},
                        index=CONTRACT + CONTROLS)
    msgs.append(f'criteria on a passing mock table: {criteria(mock)}')
    mock2 = mock.copy()
    mock2.loc['sfAFP-I175', 's_model'] = -0.40
    msgs.append(f'criteria with sfAFP at -0.40: {criteria(mock2)}')
    s_exp, ratio = experimental()
    msgs.append('experimental slopes: ' + ', '.join(f'{k} {v:+.3f}' for k, v in s_exp.items()))
    os.makedirs(OUT, exist_ok=True)
    open(os.path.join(OUT, 'selftest.txt'), 'w').write('\n'.join(msgs) + '\n')
    print('\n'.join(msgs))


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('mode', choices=['prepare', 'wuttke-types', 'guinier', 'primary', 'alternatives', 'selftest'])
    ap.add_argument('--model', choices=['calvados', 'recharged'])
    ap.add_argument('--root')
    ap.add_argument('--workers', type=int, default=6)
    ap.add_argument('--runs')
    ap.add_argument('--fasta')
    ap.add_argument('--out')
    a = ap.parse_args()
    if a.mode == 'prepare':
        prepare(a.model, a.root or ROOTS[a.model], a.workers)
    elif a.mode == 'wuttke-types':
        wuttke_types(a.model, a.runs, a.fasta, a.out)
    else:
        globals()[a.mode.replace('-', '_')]()
