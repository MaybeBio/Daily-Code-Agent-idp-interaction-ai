"""
tea_calibrate.py — per-residue-class TEA hydration scales fitted to measured temperature series, one fit per force
field, by reweighting plain-model ensembles at each temperature (no TEA simulations needed).

Model (tea_reweight.py): at temperature T the TEA-augmented force field is H + sum_c gamma_c E_c(x; T), with classes
aliph (AVLIMC), pro (P), arom (FYW), arg (R), rest (GSTNQHKDE); E_c from TEA eqs 2-3 (Arg kappa CHARMM36m; a negative
fitted gamma_arg means the data prefer the opposite Arg sign, as with the CGenFF kappa). CALVADOS-2: gamma in mol/kcal
on lambda. Mpipi-Recharged: gamma' (dimensionless) on the Wang-Frenkel eps_ij (deps_ij = -gamma' dG_E,ij), as in
tea_crossff.py. Each plain ensemble at T is reweighted with exp(-sum_c gamma_c E_c / RT); at 300 K nothing changes.

Fit: relative curves r_k(T) = <A_k>(T) / <A_k>(300 K) (A = dye-to-dye segment Rg; relative, because the measured Rg
comes from FRET through Sanchez theory with its own absolute offset) against the measured curve interpolated to the
simulated temperatures (cubic in T), residuals divided by sqrt(sigma_exp^2 + sigma_model^2) (sigma_model: 10-block SE),
plus a prior (gamma_c - gamma0_c) / tau_c; scipy least_squares; uncertainties from the Jacobian. Calibration and
held-out constructs are reported separately.

Modes
  prepare --model calvados|recharged --runs ROOT --fasta FASTA --out DIR [--temps ...] [--workers N]
      ROOT/<T>K/<name>/<name>.dcd + top.pdb (GAIVI layout); dye positions from the FASTA headers ('dyes at Ca/Cb')
      -> DIR/<model>/<name>_<T>K.npz (A, E_<class>, T)
  fit --model M --npz DIR --exp CSV --calib a,b,c --test d,e [--gamma0 2] [--tau 1] [--tau-arg 4]
      -> DIR/<model>/fit_*.csv, fit_summary.txt
  fit ... --free hyd                      fit only the hydrophilic-solvation term (TEA classes held at gamma0)
  cv --model M --npz DIR --exp CSV --free LIST   leave-one-construct-out over all five constructs
  selftest --out DIR [--n 16] [--noise 0.003]

Hydrophilic-solvation term ('hyd', added 2026-10-09 after TEA alone failed on the Wuttke data): Wuttke et al. 2014
attribute the collapse of hydrophilic IDPs to the solvation free energy of charged/polar groups and the backbone becoming
less favourable with temperature, which TEA leaves out (kappa ~ 0 for charged residues). Here the contact free energy of
residues in HYD_RES (D E K R N Q S T H, and G as the backbone analogue N-methylacetamide) changes by
    dG_contact,i(T) - dG_contact,i(T0) = -[du_h,i(T) - du_h,i(T0)]        (TEA eq 2, Table S1; i.e. kappa = -1)
mapped like TEA: dlambda_ij = gamma_hyd (w_i + w_j)/2 (CALVADOS-2, gamma_hyd in mol/kcal) or deps_ij = gamma_hyd
4.184 (w_i + w_j)/2 kJ/mol (Mpipi-Recharged), w_i = du_h,i(T) - du_h,i(T0) > 0 above 300 K (stickier on heating).
One strength gamma_hyd per force field, prior 0 with a wide width.
      synthetic data from BENDER CALVADOS-2 ensembles at 278-320 K (whole-chain Rg) with known gammas; recovery for the
      designed protein set (tea_300K/design/selected.csv) vs random sets of the same size
"""
import sys

sys.dont_write_bytecode = True

import argparse  # noqa: E402
import glob  # noqa: E402
import multiprocessing as mp  # noqa: E402
import os  # noqa: E402
import re  # noqa: E402
import traceback  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.optimize import least_squares  # noqa: E402

HERE = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, os.environ.get('BENDER_REANALYSIS_DIR', HERE))
sys.path.insert(0, HERE)
import tea_reweight as W  # noqa: E402

L, fdt, bc = W.L, W.fdt, W.bc
CLS = list(W.CLASSES)
T0 = W.T0
WUTTKE_TEMPS = [278, 288, 300, 310, 320, 335, 350]
HYD_RES = 'DEKRNQSTHG'
HYD_GROUPS = {'DE': 'DE', 'K': 'K', 'R': 'R', 'pol': 'NQSTH', 'G': 'G'}   # 'hyd' = sum of the groups (for sensitivity)
ALL = CLS + ['hyd']
INTERP_DEG = 3   # polynomial degree for interpolating the measured curves


def hyd_weights(seq, T, res=HYD_RES):
    """w_i = du_h,i(T) - du_h,i(T0) (kcal/mol) for residues in res, 0 otherwise."""
    return np.array([(W.du_h(a, T) - W.du_h(a, T0)) if a in res else 0.0 for a in seq])
R_KJ = fdt.R_GAS


# ----------------------------------------------------------------------------- inputs
def read_fasta_dyes(path):
    """{name: (sequence, i, j)} with 1-based dye residue positions from 'dyes at Ca/Cb' in the header."""
    out, name, seq, dyes = {}, None, [], None
    for line in open(path):
        line = line.strip()
        if line.startswith('>'):
            if name:
                out[name] = (''.join(seq), *dyes)
            name = line[1:].split('|')[0].strip().split()[0]
            m = re.search(r'dyes at C(\d+)/C(\d+)', line)
            dyes, seq = ((int(m.group(1)), int(m.group(2))) if m else (None, None)), []
        elif line:
            seq.append(line)
    if name:
        out[name] = (''.join(seq), *dyes)
    for k, (s, i, j) in out.items():
        if i is None:
            cys = [p + 1 for p, a in enumerate(s) if a == 'C']
            out[k] = (s, cys[0], cys[-1])
    return out


def segment_rg(xyz, i, j):
    """Rg of residues i..j (1-based, inclusive) per frame, nm."""
    x = np.asarray(xyz[:, i - 1:j], dtype=float)
    c = x - x.mean(axis=1, keepdims=True)
    return np.sqrt((c ** 2).sum(axis=2).mean(axis=1))


def wf_attractive_shape(r, sig, mu, rc, alpha):
    """Attractive part of the Wang-Frenkel shape (nu = 1), split as CALVADOS splits Ashbaugh-Hatch: the well depth -1
    inside the minimum r_min, the tail beyond it, 0 past rc. A hydration change scales this part only; scaling the full
    potential (eps multiplies the repulsive wall too, mu up to 12) put the forcing on near-neighbour wall contacts and
    made the reweighting collapse (effective sample ~0)."""
    rmin = rc * (3 / (1 + 2 * (rc / sig) ** (2 * mu))) ** (1 / (2 * mu))
    full = alpha * ((sig / r) ** (2 * mu) - 1) * ((rc / r) ** (2 * mu) - 1) ** 2
    return np.where(r < rmin, -1.0, np.where(r < rc, full, 0.0))


def recharged_class_energies(xyz, seq, T, box):
    """Per-frame E_c (kJ/mol per unit gamma') on the attractive part of the Mpipi-Recharged Wang-Frenkel term."""
    import tea_crossff as X
    tb = X.rc_tables()
    ix = np.array([tb['index'][a] for a in seq])
    i, j = np.triu_indices(len(seq), k=2)
    sig, mu = tb['sig'][ix[j], ix[i]], tb['mu'][ix[j], ix[i]]
    rc = 3 * sig
    alpha = 2 * 3 ** (2 * mu) * (3 / (2 * (3 ** (2 * mu) - 1))) ** 3
    de_res = np.array([-W.dGE(a, T) for a in seq]) * 4.184               # deps_i per unit gamma', kJ/mol
    cls = np.array([next(c for c, r in W.CLASSES.items() if a in r) for a in seq])
    w = {c: 0.5 * (np.where(cls == c, de_res, 0.0)[i] + np.where(cls == c, de_res, 0.0)[j]) for c in CLS}
    hw = hyd_weights(seq, T) * 4.184
    w['hyd'] = 0.5 * (hw[i] + hw[j])
    for g, res in HYD_GROUPS.items():
        hg = hyd_weights(seq, T, res) * 4.184
        w[f'hyd_{g}'] = 0.5 * (hg[i] + hg[j])
    E = {c: np.zeros(len(xyz)) for c in w}
    step = max(1, int(4e6 // max(1, len(i))))
    for a in range(0, len(xyz), step):
        dv = xyz[a:a + step, i].astype(float) - xyz[a:a + step, j].astype(float)
        dv -= box * np.round(dv / box)
        r = np.linalg.norm(dv, axis=2)
        shape = wf_attractive_shape(r, sig, mu, rc, alpha)
        for c in w:
            E[c][a:a + step] = (shape * w[c]).sum(axis=1)
    return E


def calvados_hyd_energies(xyz, seq, T):
    """Per-frame E_hyd and its residue-group parts (kJ/mol per unit gamma_hyd, mol/kcal) on CALVADOS-2 Ashbaugh-Hatch."""
    q, lam, sig, _ = fdt.calvados_residue_params(seq, L.load_residue_params())
    i, j, s, l, shift, _ = fdt.pair_setup(q, lam, sig)
    sets = dict({'hyd': HYD_RES}, **{f'hyd_{g}': r for g, r in HYD_GROUPS.items()})
    wij = {k: 0.5 * (hyd_weights(seq, T, r)[i] + hyd_weights(seq, T, r)[j]) for k, r in sets.items()}
    E = {k: np.zeros(len(xyz)) for k in sets}
    step = max(1, int(fdt.MAX_ELEMS // max(1, len(i))))
    for a in range(0, len(xyz), step):
        r = np.linalg.norm(np.asarray(xyz[a:a + step], dtype=float)[:, i] - np.asarray(xyz[a:a + step], dtype=float)[:, j], axis=2)
        B = fdt.ah_pair_energy(r, s, np.ones_like(l), shift) - fdt.ah_pair_energy(r, s, np.zeros_like(l), shift)
        for k in sets:
            E[k][a:a + step] = (B * wij[k]).sum(axis=1)
    return E


def calvados_hyd_energy(xyz, seq, T):
    """Per-frame E_hyd (kJ/mol per unit gamma_hyd, mol/kcal) on the CALVADOS-2 Ashbaugh-Hatch term."""
    return calvados_hyd_energies(xyz, seq, T)['hyd']


def _prepare_one(args):
    model, root, name, seq_i_j, T, out = args
    seq0, i, j = seq_i_j
    try:
        d = os.path.join(root, f'{T}K', name)
        if model == 'calvados':
            xyz, seq, info = bc.load_calvados(os.path.join(d, f'{name}.dcd'), os.path.join(d, 'top.pdb'), k=L.K_FRAMES)
            E = W.class_energies(xyz, seq, float(T))
            E.update(calvados_hyd_energies(xyz, seq, float(T)))
        else:
            import mpipi_recharged_fdt as R
            xyz, seq, info = R.load_run(d, name)
            E = recharged_class_energies(xyz, seq, float(T), R.box_nm(len(seq)))
        if seq != seq0:
            raise ValueError('simulated sequence differs from the FASTA')
        A = segment_rg(xyz, i, j)
        np.savez_compressed(os.path.join(out, f'{name}_{T}K.npz'), T=T, A=A, i=i, j=j,
                            **{f'E_{c}': E[c] for c in E})
        return name, T, None
    except Exception as ex:
        return name, T, f'{ex!r} | {traceback.format_exc()[-300:]}'.replace('\n', ' ')


def prepare(model, root, fasta, out, temps=WUTTKE_TEMPS, workers=4):
    dyes = read_fasta_dyes(fasta)
    od = os.path.join(out, model)
    os.makedirs(od, exist_ok=True)
    tasks = [(model, root, n, v, T, od) for n, v in dyes.items() for T in temps]
    with mp.get_context('spawn').Pool(workers) as pool:
        for name, T, err in pool.imap_unordered(_prepare_one, tasks):
            print(f'{name} {T} K: ' + ('ok' if not err else f'FAILED {err[:200]}'), flush=True)


# ----------------------------------------------------------------------------- model curves and fit
def load_npz(d):
    """{name: {T: (A, {c: E_c})}} from DIR/<name>_<T>K.npz."""
    out = {}
    for f in glob.glob(os.path.join(d, '*K.npz')):
        m = re.match(r'(.+)_(\d+)K\.npz$', os.path.basename(f))
        z = np.load(f)
        out.setdefault(m.group(1), {})[int(m.group(2))] = (z['A'].astype(float),
                                                            {k[2:]: z[k] for k in z.files if k.startswith('E_')})
    return out


def block_se(a, n=10):
    b = np.array_split(a, n)
    return np.std([x.mean() for x in b], ddof=1) / np.sqrt(n)


TMAX = None   # optional upper temperature for fits (set by --tmax)


def model_curve(series, gam):
    """Relative curve <A>(T)/<A>(300) for per-class gammas; also the relative sampling SE per T, inflated by
    1/sqrt(Kish effective sample fraction) so that temperatures the reweighting barely reaches count less."""
    Ts = sorted(series)
    m, ess = {}, {}
    for T in Ts:
        A, E = series[T]
        m[T], ess[T] = W.reweight(A, E, gam, float(T))
    ref = m[300]
    se = {T: block_se(series[T][0]) / ref / np.sqrt(max(ess[T], 1e-3)) for T in Ts}
    return np.array(Ts), np.array([m[T] / ref for T in Ts]), np.array([se[T] for T in Ts])


def exp_curve(df, Ts):
    """Measured relative curve interpolated to Ts (cubic in T), and its relative error."""
    p = np.polyfit(df.T_K - T0, df.Rg_nm, INTERP_DEG)
    y = np.polyval(p, np.asarray(Ts) - T0)
    ref = np.polyval(p, 0.0)
    err = np.interp(Ts, df.T_K, df.err_nm) / ref
    inside = (np.asarray(Ts) >= df.T_K.min() - 0.5) & (np.asarray(Ts) <= df.T_K.max() + 3.0)   # no extrapolation
    if TMAX is not None:
        inside &= np.asarray(Ts) <= TMAX
    return y / ref, err, inside


def fit(data, targets, calib, gamma0=2.0, tau=1.0, tau_arg=4.0, free=CLS, tau_hyd=10.0):
    """data: {name: series}; targets: {name: (Ts, r_exp, err, inside)}. Returns gammas, table with SEs, solver."""
    g0 = dict({c: gamma0 for c in CLS}, hyd=0.0)
    taus = dict({c: (tau_arg if c == 'arg' else tau) for c in CLS}, hyd=tau_hyd)

    def resid(p):
        gam = dict(g0, **dict(zip(free, p)))
        r = []
        for n in calib:
            Ts, rm, sm = model_curve(data[n], gam)
            Te, re_, ee, inside = targets[n]
            keep = (Ts != 300) & inside
            r.append((rm[keep] - re_[keep]) / np.sqrt(ee[keep] ** 2 + sm[keep] ** 2))
        r.append(np.array([(gam[c] - g0[c]) / taus[c] for c in free]))
        return np.concatenate(r)

    sol = least_squares(resid, x0=[g0[c] for c in free], method='trf')
    J = sol.jac
    cov = np.linalg.pinv(J.T @ J)
    gam = dict(g0, **dict(zip(free, sol.x)))
    return gam, pd.DataFrame({'gamma': [gam[c] for c in free], 'se': np.sqrt(np.diag(cov))}, index=free), sol


def report(data, targets, gam, names, tag):
    rows = []
    for n in names:
        Ts, rm, sm = model_curve(data[n], gam)
        _, r0, _ = model_curve(data[n], {c: 0.0 for c in ALL})
        Te, re_, ee, inside = targets[n]
        k = (Ts != 300) & inside
        slope = lambda y: 100 * np.polyfit(Ts[inside] - T0, y[inside], 2)[1]  # noqa: E731  % per K at 300 K
        rows.append({'construct': n, 'set': tag, 'rms_rel_fit': float(np.sqrt(np.mean((rm[k] - re_[k]) ** 2))),
                     'rms_rel_plain': float(np.sqrt(np.mean((r0[k] - re_[k]) ** 2))),
                     'slope_exp': slope(re_), 'slope_fit': slope(rm), 'slope_plain': slope(r0)})
    return pd.DataFrame(rows)


def fit_mode(model, npz_dir, exp_csv, calib, test, gamma0, tau, tau_arg, fixed=('rest',), free=None, tag=''):
    d = os.path.join(npz_dir, model)
    data = load_npz(d)
    ex = pd.read_csv(exp_csv)
    targets = {}
    for n in calib + test:
        Ts = np.array(sorted(data[n]))
        targets[n] = (Ts, *exp_curve(ex[ex.construct == n].sort_values('T_K'), Ts))
    free = free or [c for c in CLS if c not in fixed]
    gam, tab, sol = fit(data, targets, calib, gamma0, tau, tau_arg, free=free)
    rep = pd.concat([report(data, targets, gam, calib, 'calibration'), report(data, targets, gam, test, 'held out')])
    tab.to_csv(os.path.join(d, f'fit{tag}_gammas.csv'))
    rep.to_csv(os.path.join(d, f'fit{tag}_constructs.csv'), index=False)
    txt = [f'{model}: free {free}; others held at prior (TEA classes gamma0 = {gamma0}, hyd 0); tau = {tau}, '
           f'tau_arg = {tau_arg}; calibration {calib}; cost {sol.cost:.3f}', tab.round(4).to_string(), '',
           rep.round(4).to_string(index=False)]
    open(os.path.join(d, f'fit{tag}_summary.txt'), 'w').write('\n'.join(txt) + '\n')
    print('\n'.join(txt))


def cv_mode(model, npz_dir, exp_csv, free, gamma0, tau, tau_arg, names=None, tag=''):
    """Leave one construct out: fit on the other four, predict the left-out one."""
    d = os.path.join(npz_dir, model)
    data = load_npz(d)
    ex = pd.read_csv(exp_csv)
    names = names or ['lambda_P6C_I84C', 'IN_NTD', 'CspM34', 'ProTa_C2_C56', 'ProTa_C56_C110']
    targets = {}
    for n in names:
        Ts = np.array(sorted(data[n]))
        targets[n] = (Ts, *exp_curve(ex[ex.construct == n].sort_values('T_K'), Ts))
    rows, gams = [], []
    for out in names:
        calib = [n for n in names if n != out]
        gam, tab, sol = fit(data, targets, calib, gamma0, tau, tau_arg, free=free)
        rows.append(report(data, targets, gam, [out], 'left out'))
        gams.append({'left_out': out, **{f'gamma_{c}': gam[c] for c in free}, **{f'se_{c}': tab.loc[c, 'se'] for c in free}})
    rep, gt = pd.concat(rows), pd.DataFrame(gams)
    rep.to_csv(os.path.join(d, f'cv{tag}_constructs.csv'), index=False)
    gt.to_csv(os.path.join(d, f'cv{tag}_gammas.csv'), index=False)
    txt = [f'{model}: leave-one-construct-out, free {free}', gt.round(4).to_string(index=False), '',
           rep.round(4).to_string(index=False),
           f'held-out slope: r = {np.corrcoef(rep.slope_exp, rep.slope_fit)[0, 1]:.3f}, '
           f'median |fit - exp| = {np.median(np.abs(rep.slope_fit - rep.slope_exp)):.3f} %/K '
           f'(plain: {np.median(np.abs(rep.slope_plain - rep.slope_exp)):.3f})']
    open(os.path.join(d, f'cv{tag}_summary.txt'), 'w').write('\n'.join(txt) + '\n')
    print('\n'.join(txt))


VARIANTS = {'all (DEKRNQSTH+G)': ['DE', 'K', 'R', 'pol', 'G'], 'without G (backbone)': ['DE', 'K', 'R', 'pol'],
            'without R': ['DE', 'K', 'pol', 'G'], 'charged only (DEKR)': ['DE', 'K', 'R'],
            'polar + G only': ['pol', 'G'], 'acidic only (DE)': ['DE']}
FIVE = ['lambda_P6C_I84C', 'IN_NTD', 'CspM34', 'ProTa_C2_C56', 'ProTa_C56_C110']
FOUR = ['IN_NTD', 'CspM34', 'ProTa_C2_C56', 'ProTa_C56_C110']


def sens_mode(model, npz_dir, exp_csv, gamma0):
    """Calibration sensitivity: hydrophilic residue set, interpolation degree, temperature range. For each setting: the
    final fit on the four constructs (gamma_hyd +/- SE) and leave-one-out over the four (held-out slope / curve errors)."""
    global INTERP_DEG, TMAX
    d = os.path.join(npz_dir, model)
    raw = load_npz(d)
    ex = pd.read_csv(exp_csv)
    settings = [(v, 3, 320.0) for v in VARIANTS] + [('all (DEKRNQSTH+G)', 2, 320.0), ('all (DEKRNQSTH+G)', 3, 335.0),
                                                    ('all (DEKRNQSTH+G)', 3, None)]
    rows = []
    for var, deg, tmax in settings:
        INTERP_DEG, TMAX = deg, tmax
        data = {n: {T: (A, dict(E, hyd=sum(E[f'hyd_{g}'] for g in VARIANTS[var]))) for T, (A, E) in s_.items()}
                for n, s_ in raw.items()}
        targets = {}
        for n in FIVE:
            Ts = np.array(sorted(data[n]))
            targets[n] = (Ts, *exp_curve(ex[ex.construct == n].sort_values('T_K'), Ts))
        gam, tab, _ = fit(data, targets, FOUR, gamma0, 1.0, 2.0, free=['hyd'])
        loo, gs = [], []
        for out in FOUR:
            g2, _, _ = fit(data, targets, [n for n in FOUR if n != out], gamma0, 1.0, 2.0, free=['hyd'])
            loo.append(report(data, targets, g2, [out], 'left out'))
            gs.append(g2['hyd'])
        loo = pd.concat(loo)
        lam = report(data, targets, gam, ['lambda_P6C_I84C'], 'held out')
        rows.append({'model': model, 'hyd residues': var, 'interp_deg': deg, 'tmax': tmax or 'all',
                     'gamma_hyd': gam['hyd'], 'se': tab.loc['hyd', 'se'], 'gamma_loo_min': min(gs), 'gamma_loo_max': max(gs),
                     'loo_slope_err': float(np.median(np.abs(loo.slope_fit - loo.slope_exp))),
                     'loo_slope_err_plain': float(np.median(np.abs(loo.slope_plain - loo.slope_exp))),
                     'loo_curve_rms': float(loo.rms_rel_fit.median()), 'loo_curve_rms_plain': float(loo.rms_rel_plain.median()),
                     'lambda_curve_rms': float(lam.rms_rel_fit.iloc[0]), 'lambda_curve_rms_plain': float(lam.rms_rel_plain.iloc[0])})
        print(rows[-1], flush=True)
    INTERP_DEG, TMAX = 3, None
    t = pd.DataFrame(rows)
    t.to_csv(os.path.join(d, 'sensitivity.csv'), index=False)
    open(os.path.join(d, 'sensitivity.txt'), 'w').write(t.round(4).to_string(index=False) + '\n')
    print(t.round(4).to_string(index=False))


# ----------------------------------------------------------------------------- self-test on BENDER ensembles
def _bender_series(args):
    uid, kingdom = args
    try:
        res = L.load_residue_params()
        dirs = W._dirs(uid, kingdom)
        s = {}
        for T, d in dirs.items():
            z = np.load(os.path.join(W.FDT_DIR, 'npz', f'{T}K', f'{uid}.npz'))
            xyz, seq, info = bc.load_calvados(os.path.join(d, f'{uid}.dcd'), os.path.join(d, 'top.pdb'), k=L.K_FRAMES)
            s[T] = (z['rg'].astype(float), W.class_energies(xyz, seq, float(T), residues=res))
        return uid, s, None
    except Exception as ex:
        return uid, None, repr(ex)


def selftest(out, n=16, noise=0.003, n_random=5, seed=1):
    os.makedirs(out, exist_ok=True)
    truth = {'aliph': 3.0, 'pro': 1.0, 'arom': 2.5, 'arg': -1.0, 'rest': 0.5}
    design = pd.read_csv(os.path.join(W.FDT_DIR, 'tea_300K', 'design', 'selected.csv')).UniProt_ID.tolist()[:n]
    bio = pd.read_csv(L.BIO_CSV, usecols=['UniProt_ID', 'kingdom', 'length']).set_index('UniProt_ID')
    pool_ids = pd.read_csv(os.path.join(W.FDT_DIR, 'tea_300K', 'design', 'candidates.csv')).UniProt_ID.tolist()
    rng = np.random.default_rng(seed)
    sets = {'designed': design}
    for r in range(n_random):
        sets[f'random{r + 1}'] = list(rng.choice([u for u in pool_ids if u not in design], n, replace=False))
    need = sorted({u for s in sets.values() for u in s})
    data = {}
    with mp.get_context('spawn').Pool(8) as pool:
        for uid, s, err in pool.imap_unordered(_bender_series, [(u, bio.kingdom.get(u)) for u in need]):
            if err:
                print(f'skip {uid}: {err[:150]}', flush=True)
            else:
                data[uid] = s
    rows = []
    for name, ids in sets.items():
        ids = [u for u in ids if u in data]
        targets = {}
        for u in ids:
            Ts, r, se = model_curve(data[u], truth)
            targets[u] = (Ts, r + rng.normal(0, noise, len(r)) * (Ts != 300), np.full(len(r), noise), np.ones(len(r), bool))
        gam, tab, _ = fit(data, targets, ids, gamma0=2.0, tau=10.0, tau_arg=10.0)   # weak prior: test the data
        for c in CLS:
            rows.append({'set': name, 'n': len(ids), 'class': c, 'truth': truth[c], 'fit': tab.loc[c, 'gamma'],
                         'se': tab.loc[c, 'se'], 'error': tab.loc[c, 'gamma'] - truth[c]})
    d = pd.DataFrame(rows)
    d.to_csv(os.path.join(out, 'selftest.csv'), index=False)
    piv = d.pivot_table(index='set', columns='class', values='se')
    err = d.groupby('set').apply(lambda x: np.sqrt((x.error ** 2).mean()), include_groups=False).rename('rms_error')
    txt = [f'synthetic data: BENDER CALVADOS-2 ensembles 278-320 K, whole-chain Rg, truth {truth}, noise {noise} '
           '(relative), weak prior (tau = 10)', '', 'fitted gamma per class:',
           d.pivot_table(index='set', columns='class', values='fit').round(2).to_string(), '', 'standard errors:',
           piv.round(2).to_string(), '', err.round(3).to_string()]
    open(os.path.join(out, 'selftest_summary.txt'), 'w').write('\n'.join(txt) + '\n')
    print('\n'.join(txt))


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest='mode', required=True)
    p = sub.add_parser('prepare')
    p.add_argument('--model', required=True, choices=['calvados', 'recharged'])
    p.add_argument('--runs', required=True)
    p.add_argument('--fasta', required=True)
    p.add_argument('--out', required=True)
    p.add_argument('--temps', type=int, nargs='+', default=WUTTKE_TEMPS)
    p.add_argument('--workers', type=int, default=4)
    f = sub.add_parser('fit')
    f.add_argument('--model', required=True, choices=['calvados', 'recharged'])
    f.add_argument('--npz', required=True)
    f.add_argument('--exp', required=True)
    f.add_argument('--calib', default='lambda_P6C_I84C,IN_NTD,CspM34')
    f.add_argument('--test', default='ProTa_C2_C56,ProTa_C56_C110')
    f.add_argument('--gamma0', type=float, default=2.0)
    f.add_argument('--tau', type=float, default=1.0)
    f.add_argument('--tau-arg', type=float, default=4.0)
    f.add_argument('--fixed', default='rest', help='classes held at gamma0 (comma-separated; the self-test shows '
                   '"rest" is unidentifiable because TEA gives polar/charged residues almost no forcing)')
    f.add_argument('--free', default=None, help='comma-separated scales to fit (default: TEA classes not --fixed); '
                   'e.g. hyd or aliph,pro,arom,arg,hyd')
    f.add_argument('--tag', default='')
    f.add_argument('--tmax', type=float, default=None, help='only temperatures up to this (K)')
    v = sub.add_parser('cv')
    v.add_argument('--model', required=True, choices=['calvados', 'recharged'])
    v.add_argument('--npz', required=True)
    v.add_argument('--exp', required=True)
    v.add_argument('--free', required=True)
    v.add_argument('--gamma0', type=float, default=2.0)
    v.add_argument('--tau', type=float, default=1.0)
    v.add_argument('--tau-arg', type=float, default=4.0)
    v.add_argument('--tag', default='')
    v.add_argument('--tmax', type=float, default=None, help='only temperatures up to this (K)')
    v.add_argument('--names', default=None, help='constructs in the leave-one-out (default: all five; the calibration '
                   'uses IN_NTD,CspM34,ProTa_C2_C56,ProTa_C56_C110)')
    z = sub.add_parser('sens')
    z.add_argument('--model', required=True, choices=['calvados', 'recharged'])
    z.add_argument('--npz', required=True)
    z.add_argument('--exp', required=True)
    z.add_argument('--gamma0', type=float, default=2.0)
    s = sub.add_parser('selftest')
    s.add_argument('--out', required=True)
    s.add_argument('--n', type=int, default=16)
    s.add_argument('--noise', type=float, default=0.003)
    a = ap.parse_args()
    if getattr(a, 'tmax', None) is not None:
        TMAX = a.tmax
    if a.mode == 'prepare':
        prepare(a.model, a.runs, a.fasta, a.out, a.temps, a.workers)
    elif a.mode == 'fit':
        fit_mode(a.model, a.npz, a.exp, a.calib.split(','), a.test.split(','), a.gamma0, a.tau, a.tau_arg,
                 tuple(x for x in a.fixed.split(',') if x), a.free.split(',') if a.free else None, a.tag)
    elif a.mode == 'sens':
        sens_mode(a.model, a.npz, a.exp, a.gamma0)
    elif a.mode == 'cv':
        cv_mode(a.model, a.npz, a.exp, a.free.split(','), a.gamma0, a.tau, a.tau_arg,
                names=a.names.split(',') if a.names else None, tag=a.tag)
    else:
        selftest(a.out, a.n, a.noise)
