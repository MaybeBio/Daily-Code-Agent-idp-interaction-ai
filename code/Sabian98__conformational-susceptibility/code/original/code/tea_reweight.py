"""
tea_reweight.py — TEA-augmented CALVADOS-2 at any temperature by reweighting plain CALVADOS-2 frames (no new
simulations), and a check of that against the FDT + TEA slope at 300 K.

TEA (Chen & Zeng, JACS Au 2026) changes only the Ashbaugh-Hatch stickiness with temperature:
    eq 2   du_h,i(T) = [du_h,i(T0) - dh_i] T/T0 + dh_i + dc_p,i [T (1 - ln T/T0) - T0]        (Table S1)
    eq 3   G_E,i(T) - G_E,i(T0) = kappa_i [du_h,i(T) - du_h,i(T0)]                            (Table 1)
    eq 6   lambda_i(T) - lambda_i(T0) = -gamma_c(i) [G_E,i(T) - G_E,i(T0)]      (sign: see tea_fdt.py)
with one scale per residue class c (aliph AVLIMC, pro P, arom FYW, arg R, rest GSTNQHKDE). lambda_ij is the mean of
lambda_i and lambda_j and U_AH is affine in lambda_ij, so at temperature T
    H_TEA(x; T) = H(x; T) + sum_c gamma_c E_c(x; T),   E_c = sum_ij dlambda_ij^c(T) [U_AH(r; 1) - U_AH(r; 0)]
and a plain CALVADOS-2 ensemble at T reweights to TEA at T with weights exp(-sum_c gamma_c E_c / RT). At T0 = 300 K
all E_c = 0: the 300 K ensemble is the TEA ensemble, as tea_fdt.py uses.

check mode: for a random sample of proteins with stored FDT frames at 278-320 K, reweight each temperature to TEA
(all classes at one gamma), fit <Rg>(T) with a quadratic, and compare its 300 K slope with chi + gamma * dchi from
fdt_results_300K.csv and tea_dchi_300K.csv; Kish effective sample size per temperature. The same comparison for
gamma = 0 checks the plain FDT slope against the plain five-temperature curve.
    python tea_reweight.py check --out DIR [--n 300] [--workers 8]
"""
import sys

sys.dont_write_bytecode = True

import argparse  # noqa: E402
import multiprocessing as mp  # noqa: E402
import os  # noqa: E402
import time  # noqa: E402
import traceback  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

HERE = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, os.environ.get('BENDER_REANALYSIS_DIR', HERE))
sys.path.insert(0, HERE)
import tea_fdt as TF  # noqa: E402

L, fdt, bc = TF.L, TF.fdt, TF.bc
T0 = TF.T0
TEMPS = [278, 288, 300, 310, 320]
CLASSES = {'aliph': 'AVLIMC', 'pro': 'P', 'arom': 'FYW', 'arg': 'R', 'rest': 'GSTNQHKDE'}
# Table S1 heat capacities of hydration, cal/mol/K (residue analogues; Val and Pro share propane)
DCP = {'A': 35.30, 'V': 40.50, 'P': 40.50, 'L': 50.49, 'I': 62.68, 'M': 43.66, 'F': 63.67, 'C': 39.10, 'Y': 41.91,
       'W': 50.57, 'S': 21.47, 'T': 37.74, 'N': 24.54, 'Q': 36.41, 'H': 41.86, 'G': 41.81, 'R': 58.02, 'K': 41.93,
       'D': 3.49, 'E': 4.68}
FDT_DIR = os.path.join(L.BENDER_MT, 'Traj Level Analysis', 'Analysis_FDT')


def du_h(aa: str, T: float) -> float:
    """Hydration free energy of residue aa's analogue at T (kcal/mol), eq 2."""
    _, du0, dh = TF.TEA[aa]
    cp = DCP[aa] / 1000.0
    return (du0 - dh) * T / T0 + dh + cp * (T * (1 - np.log(T / T0)) - T0)


def dGE(aa: str, T: float, arg: str = 'charmm36') -> float:
    """G_E,i(T) - G_E,i(T0), kcal/mol, eq 3."""
    kap = TF.TEA[aa][0] if not (aa == 'R' and arg == 'charmm36') else TF.KAPPA_R_CHARMM36
    return kap * (du_h(aa, T) - du_h(aa, T0))


def class_energies(xyz, seq, T, arg='charmm36', residues=None):
    """Per-frame E_c (kJ/mol per unit gamma_c, gamma in mol/kcal) for the CALVADOS-2 Ashbaugh-Hatch term."""
    residues = residues if residues is not None else L.load_residue_params()
    q, lam, sig, _ = fdt.calvados_residue_params(seq, residues)
    i, j, s, l, shift, _ = fdt.pair_setup(q, lam, sig)
    dl_res = np.array([-dGE(a, T, arg) for a in seq])            # dlambda_i per unit gamma at T
    cls = np.array([next(c for c, r in CLASSES.items() if a in r) for a in seq])
    w = {c: 0.5 * (np.where(cls == c, dl_res, 0.0)[i] + np.where(cls == c, dl_res, 0.0)[j]) for c in CLASSES}
    n = len(xyz)
    E = {c: np.zeros(n) for c in CLASSES}
    step = max(1, int(fdt.MAX_ELEMS // max(1, len(i))))
    for a in range(0, n, step):
        x = np.asarray(xyz[a:a + step], dtype=float)
        r = np.linalg.norm(x[:, i] - x[:, j], axis=2)
        B = fdt.ah_pair_energy(r, s, np.ones_like(l), shift) - fdt.ah_pair_energy(r, s, np.zeros_like(l), shift)
        for c in CLASSES:
            E[c][a:a + step] = (B * w[c]).sum(axis=1)
    return E


def reweight(A, E, gammas, T):
    """<A> under TEA at T for per-class gammas (dict), and the Kish effective sample size fraction."""
    u = sum(gammas.get(c, 0.0) * E[c] for c in E) / (fdt.R_GAS * T)
    w = np.exp(-(u - u.min()))
    return float((w * A).sum() / w.sum()), float(w.sum() ** 2 / (w ** 2).sum() / len(w))


def _dirs(uid, kingdom):
    out = {}
    for T in TEMPS:
        base = L.TRAJ_ROOT_FMT.format(T=T)
        d = L.find_traj_dir(uid, base, hint=kingdom, flat=(T != 300))
        if d is None:
            raise FileNotFoundError(f'{uid} {T} K')
        out[T] = d
    return out


GAMMAS = [0.0, 1.0, 2.0, 3.0]


def check_one(args):
    uid, kingdom = args
    t0 = time.time()
    try:
        res = L.load_residue_params()
        dirs = _dirs(uid, kingdom)
        rows = []
        for T, d in dirs.items():
            z = np.load(os.path.join(FDT_DIR, 'npz', f'{T}K', f'{uid}.npz'))
            xyz, seq, info = bc.load_calvados(os.path.join(d, f'{uid}.dcd'), os.path.join(d, 'top.pdb'), k=L.K_FRAMES)
            fi = info['production_start'] + info['stride'] * np.arange(len(xyz))
            if not np.array_equal(fi, z['frame_index']):
                raise ValueError(f'{T} K frames differ from the stored FDT frames')
            E = class_energies(xyz, seq, float(T), residues=res)
            rg = z['rg'].astype(float)
            for g in GAMMAS:
                m, ess = reweight(rg, E, {c: g for c in CLASSES}, float(T))
                rows.append({'T': T, 'gamma': g, 'rg': m, 'ess_frac': ess})
        d = pd.DataFrame(rows)
        out = {'UniProt_ID': uid, 'n_residues': len(seq), 'seconds': round(time.time() - t0, 1)}
        for g, gd in d.groupby('gamma'):
            gd = gd.sort_values('T')
            p = np.polyfit(gd['T'] - T0, gd.rg, 2)
            out[f'slope_rel_g{g:g}'] = 100 * p[1] / p[2]
            out[f'ess_min_g{g:g}'] = gd.ess_frac.min()
        return out, None
    except Exception as ex:
        return {'UniProt_ID': uid}, f'{ex!r} | {traceback.format_exc()[-300:]}'.replace('\n', ' ')


def check(out, n=300, workers=8, seed=0):
    os.makedirs(out, exist_ok=True)
    fdt300 = pd.read_csv(os.path.join(FDT_DIR, 'fdt_results_300K.csv'))
    fdt300 = fdt300[fdt300.metric == 'rg'].set_index('UniProt_ID')
    tea = pd.read_csv(os.path.join(FDT_DIR, 'tea_300K', 'tea_dchi_300K.csv')).set_index('UniProt_ID')
    bio = pd.read_csv(L.BIO_CSV, usecols=['UniProt_ID', 'kingdom']).set_index('UniProt_ID')
    have = set(fdt300.index)
    for T in TEMPS:
        have &= {f[:-4] for f in os.listdir(os.path.join(FDT_DIR, 'npz', f'{T}K')) if f.endswith('.npz') and not f.startswith('._')}
    ids = sorted(have & set(tea.index))
    rng = np.random.default_rng(seed)
    pick = list(rng.choice(ids, min(n, len(ids)), replace=False))
    path = os.path.join(out, 'tea_reweight_check.csv')
    print(f'{len(ids)} proteins with all five temperatures; checking {len(pick)}', flush=True)
    rows = []
    with mp.get_context('spawn').Pool(workers) as pool:
        for k, (row, err) in enumerate(pool.imap_unordered(check_one, [(u, bio.kingdom.get(u)) for u in pick]), 1):
            if err:
                print(f"failed {row['UniProt_ID']}: {err[:200]}", flush=True)
                continue
            rows.append(row)
            if k % 50 == 0:
                print(f'{k}/{len(pick)}', flush=True)
    d = pd.DataFrame(rows).set_index('UniProt_ID')
    d['chi_fdt_rel'] = 100 * fdt300.chi_total / fdt300.A_mean
    d['dchi_tea_rel'] = tea.dchi_rg_rel_charmm36
    d.to_csv(path)
    lines = [f'{len(d)} proteins; slope of <Rg>(T) (quadratic over 278-320 K, % of <Rg> per K at 300 K) from reweighted '
             'plain frames vs FDT chi + gamma * TEA dchi at 300 K']
    for g in GAMMAS:
        pred = d.chi_fdt_rel + g * d.dchi_tea_rel
        obs = d[f'slope_rel_g{g:g}']
        lines.append(f'gamma {g:g}: r = {np.corrcoef(pred, obs)[0, 1]:.3f}, median(obs - pred) = {(obs - pred).median():+.4f}, '
                     f'median |obs - pred| = {(obs - pred).abs().median():.4f}, median pred = {pred.median():+.4f}, '
                     f'min ESS fraction median {d[f"ess_min_g{g:g}"].median():.3f} (5th pct {d[f"ess_min_g{g:g}"].quantile(0.05):.3f})')
    open(os.path.join(out, 'summary.txt'), 'w').write('\n'.join(lines) + '\n')
    print('\n'.join(lines))


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest='mode', required=True)
    c = sub.add_parser('check')
    c.add_argument('--out', required=True)
    c.add_argument('--n', type=int, default=300)
    c.add_argument('--workers', type=int, default=8)
    a = ap.parse_args()
    check(a.out, a.n, a.workers)
