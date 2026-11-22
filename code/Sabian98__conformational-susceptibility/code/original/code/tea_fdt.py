"""
tea_fdt.py — hydration-driven dH/dT for CALVADOS-2 at 300 K from TEA (Chen & Zeng, JACS Au 2026, 6, 4478), added to
the FDT susceptibility without new simulations.

TEA makes the Ashbaugh-Hatch stickiness depend on temperature while keeping the model's own lambda at T0 = 300 K, so
the CALVADOS-2 300 K trajectories are also TEA-CALVADOS-2 trajectories at 300 K. From the paper:
    eq 2   du_h/dT at T0       = (du_h(T0) - dh) / T0                     (Table S1; the dc_p term vanishes at T0)
    eq 3   dG_E/dT             = kappa_i du_h/dT                           (Table 1)
    eq 6   lambda(T) - lambda(T0) = gamma [G_E(T) - G_E(T0)], gamma ~ 2 mol/kcal (best for TEA-HPS-Urry)
Sign: as printed, eq 6 with gamma > 0 would make Leu/Ile/Phe *less* sticky on heating, against the paper's own Fig. 1
("G_E decreases ... indicating stronger attractions"), its LCST results and its text. Here
    dlambda_i/dT = -gamma kappa_i (du_h(T0) - dh_i) / T0      (hydrophobic residues stickier on heating)
Pairs: lambda_ij = (lambda_i + lambda_j)/2 (CALVADOS and TEA's combination rule), so dlambda_ij/dT is the mean.
U_AH is affine in lambda_ij, so per frame
    D = sum_ij (dlambda_ij/dT per unit gamma) [U_AH(r; lambda=1) - U_AH(r; lambda=0)]      kJ/mol/K per unit gamma
and the extra susceptibility is linear in gamma:  chi_TEA(A) = chi(A) + gamma * dchi(A),  dchi(A) = -Cov(A, D)/(R T).
Arg: Table 1 lists kappa for the CGenFF analogue (-0.014) and the CHARMM36 one (0.059); both are computed.

Frames are the stored FDT frames (npz/300K frame_index, 1,000 per protein); the recomputed V_AH is checked against the
stored one for every protein.

    python tea_fdt.py --out DIR [--workers N] [--limit N]
Output: DIR/tea_dchi_300K.csv, one row per protein: dchi per unit gamma for Rg and Ree (absolute and % of <A> per K,
both Arg choices), mean/sd of D, the V_AH check.
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
import TODO_T2_trajectory_loader as L  # noqa: E402

fdt, bc = L.fdt, L.bc
T0 = 300.0
NPZ = os.path.join(L.BENDER_MT, 'Traj Level Analysis', 'Analysis_FDT', 'npz', '300K')
BASE = L.TRAJ_ROOT_FMT.format(T=300)

# Table 1 (kappa_i) and Table S1 (du_h(300 K), dh in kcal/mol), Chen & Zeng, JACS Au 2026
TEA = {  # one-letter: (kappa, du_h(T0), dh)
    'A': (-0.103, 2.45, 0.02), 'V': (-0.119, 2.54, -1.13), 'P': (-0.119, 2.54, -1.13), 'L': (-0.178, 2.90, -1.17),
    'I': (-0.170, 2.74, -1.93), 'M': (-0.120, 1.03, -3.75), 'F': (-0.100, 0.37, -5.71), 'C': (-0.054, -0.02, -3.69),
    'Y': (-0.039, -4.00, -10.60), 'W': (-0.039, -5.25, -12.94), 'S': (0.003, -4.38, -8.83), 'T': (-0.009, -4.24, -10.10),
    'N': (-0.002, -7.53, -13.01), 'Q': (-0.009, -7.27, -13.80), 'H': (0.001, -10.55, -18.72),
    'G': (-0.004, -7.15, -14.06), 'R': (-0.014, -45.74, -57.99), 'K': (-0.007, -53.54, -62.55),
    'D': (-0.002, -97.11, -111.95), 'E': (0.001, -97.05, -112.90),
}
KAPPA_R_CHARMM36 = 0.059


def dlambda_dT(arg: str = 'cgenff') -> dict:
    """dlambda/dT at T0 per unit gamma (gamma in mol/kcal), per residue."""
    out = {}
    for aa, (kap, du0, dh) in TEA.items():
        if aa == 'R' and arg == 'charmm36':
            kap = KAPPA_R_CHARMM36
        out[aa] = -kap * (du0 - dh) / T0
    return out


DL = {'cgenff': dlambda_dT('cgenff'), 'charmm36': dlambda_dT('charmm36')}
# --classes: the CHARMM36m-Arg forcing split by residue class (classes add up to 'charmm36'); Arg also with CGenFF kappa
CLASSES = {'aliph': 'AVLIMC', 'pro': 'P', 'arom': 'FYW', 'argC36': 'R', 'rest': 'GSTNQHKDE'}
DL_CLASSES = {c: {a: (v if a in res else 0.0) for a, v in DL['charmm36'].items()} for c, res in CLASSES.items()}
DL_CLASSES['argCGenFF'] = {a: (v if a == 'R' else 0.0) for a, v in DL['cgenff'].items()}
# hydrophilic-solvation forcing (Wuttke et al. 2014 mechanism; tea_calibrate.py 'hyd'): contact free energy of charged,
# polar and Gly (backbone) residues follows -du_h(T), i.e. kappa = -1: dlambda_i/dT = gamma_hyd (du_h(T0) - dh)/T0.
# Not part of the TEA classes (they still add up to 'charmm36'); a separate term with its own scale.
HYD_RES = 'DEKRNQSTHG'
DL_CLASSES['hyd'] = {a: ((du0 - dh) / T0 if a in HYD_RES else 0.0) for a, (_, du0, dh) in TEA.items()}
# its residue-group parts (they add up to 'hyd'), for sensitivity to which residues carry the hydrophilic forcing
HYD_GROUPS = {'DE': 'DE', 'K': 'K', 'R': 'R', 'pol': 'NQSTH', 'G': 'G'}
for _g, _res in HYD_GROUPS.items():
    DL_CLASSES[f'hyd_{_g}'] = {a: ((du0 - dh) / T0 if a in _res else 0.0) for a, (_, du0, dh) in TEA.items()}
DL_MODES = {'base': DL, 'classes': DL_CLASSES}
# --all-observables: every per-frame observable stored in the FDT npz, not only Rg and Ree (added 2026-10-10)
OBSERVABLES = ['rg', 'ree', 'delta', 'graph_fragmentation_index', 'graph_global_efficiency', 'graph_avg_clustering',
               'graph_transitivity', 'graph_degree_assortativity']
ALL_OBS = os.environ.get('TEA_ALL_OBS') == '1'   # set in the workers too (spawn re-imports this module)


def one(args):
    uid, d, mode = args
    dls = DL_MODES[mode]
    t0 = time.time()
    try:
        z = np.load(os.path.join(NPZ, f'{uid}.npz'))
        xyz, seq, info = bc.load_calvados(os.path.join(d, f'{uid}.dcd'), os.path.join(d, 'top.pdb'), k=L.K_FRAMES)
        fi = info['production_start'] + info['stride'] * np.arange(len(xyz))
        if not np.array_equal(fi, z['frame_index']):
            raise ValueError('frames differ from the stored FDT frames')
        q, lam, sig, _ = fdt.calvados_residue_params(seq, L.load_residue_params())
        i, j, s, l, shift, _ = fdt.pair_setup(q, lam, sig)
        dl = {k: 0.5 * (np.array([v[a] for a in seq])[i] + np.array([v[a] for a in seq])[j]) for k, v in dls.items()}
        n = len(xyz)
        D = {k: np.zeros(n) for k in dls}
        V_AH = np.zeros(n)
        step = max(1, int(fdt.MAX_ELEMS // max(1, len(i))))
        for a in range(0, n, step):
            x = xyz[a:a + step].astype(float)
            r = np.linalg.norm(x[:, i] - x[:, j], axis=2)
            u1, u0 = fdt.ah_pair_energy(r, s, np.ones_like(l), shift), fdt.ah_pair_energy(r, s, np.zeros_like(l), shift)
            V_AH[a:a + step] = (u0 + l * (u1 - u0)).sum(axis=1)
            for k in D:
                # dl: dlambda/dT per unit gamma (gamma mol/kcal x kcal/mol/K = 1/K); u1 - u0 is already kJ/mol, so no unit
                # factor. (Fixed 2026-10-09: an extra x4.184 here made every dchi 4.184x too large; outputs rescaled.)
                D[k][a:a + step] = ((u1 - u0) * dl[k]).sum(axis=1)
        dv = np.abs(V_AH - z['V_AH']).max() / max(1.0, np.abs(z['V_AH']).max())
        RT = fdt.R_GAS * T0
        row = {'UniProt_ID': uid, 'n_residues': len(seq), 'V_AH_check_rel': dv, 'seconds': round(time.time() - t0, 1)}
        for m in [x for x in OBSERVABLES if x in z.files] if ALL_OBS else ('rg', 'ree'):
            A = z[m].astype(float)
            ok = np.isfinite(A)           # degree assortativity is undefined for some frames
            for k in D:
                dchi = -fdt._cov(A[ok], D[k][ok]) / RT
                row[f'dchi_{m}_{k}'] = dchi
                row[f'dchi_{m}_rel_{k}'] = 100 * dchi / A[ok].mean()
            row[f'{m}_mean'] = A[ok].mean()
        for k in D:
            row[f'D_mean_{k}'], row[f'D_sd_{k}'] = D[k].mean(), D[k].std(ddof=1)
        return row, None
    except Exception as ex:
        return {'UniProt_ID': uid}, f'{ex!r} | {traceback.format_exc()[-300:]}'.replace('\n', ' ')


def index_dirs(ids: list) -> dict:
    """BENDER_300K/<group>/<ID>: list each group once instead of searching per protein."""
    out = {}
    for g in sorted(os.listdir(BASE)):
        gd = os.path.join(BASE, g)
        if not os.path.isdir(gd):
            continue
        if os.path.isfile(os.path.join(gd, f'{g}.dcd')):   # flat layout entry
            out.setdefault(g, gd)
            continue
        for u in os.listdir(gd):
            if u not in out and os.path.isfile(os.path.join(gd, u, f'{u}.dcd')):
                out[u] = os.path.join(gd, u)
    return {u: out[u] for u in ids if u in out}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', required=True)
    ap.add_argument('--workers', type=int, default=8)
    ap.add_argument('--limit', type=int, default=None)
    ap.add_argument('--classes', action='store_true', help='split the TEA term by residue class')
    ap.add_argument('--ids', default=None, help='file with one UniProt ID per line (default: all stored npz)')
    ap.add_argument('--all-observables', action='store_true', help='every stored per-frame observable (sets TEA_ALL_OBS=1)')
    a = ap.parse_args()
    if a.all_observables:
        os.environ['TEA_ALL_OBS'] = '1'
        global ALL_OBS
        ALL_OBS = True
    os.makedirs(a.out, exist_ok=True)
    mode = 'classes' if a.classes else 'base'
    path = os.path.join(a.out, 'tea_dchi_classes_300K.csv' if a.classes else 'tea_dchi_300K.csv')
    ids = sorted(f[:-4] for f in os.listdir(NPZ) if f.endswith('.npz') and not f.startswith('._'))
    if a.ids:
        want = {l.strip() for l in open(a.ids) if l.strip()}
        ids = [u for u in ids if u in want]
    done = set(pd.read_csv(path, usecols=['UniProt_ID']).UniProt_ID) if os.path.exists(path) else set()
    dirs = index_dirs([u for u in ids if u not in done])
    todo = [(u, d, mode) for u, d in dirs.items()][:a.limit]
    print(f'{len(ids)} npz, {len(done)} done, {len(dirs)} trajectories found, {len(todo)} to run', flush=True)
    t0, nf = time.time(), 0
    with mp.get_context('spawn').Pool(a.workers) as pool:
        for k, (row, err) in enumerate(pool.imap_unordered(one, todo), 1):
            if err:
                nf += 1
                print(f"failed {row['UniProt_ID']}: {err[:200]}", flush=True)
                continue
            pd.DataFrame([row]).to_csv(path, mode='a', header=not os.path.exists(path), index=False)
            if k % 250 == 0:
                el = (time.time() - t0) / 3600
                print(f'{k}/{len(todo)} {el:.2f} h, ~{el / k * (len(todo) - k):.1f} h left, {nf} failed', flush=True)
    bc.write_provenance(path, __file__, {'T0_K': T0, 'source': 'Chen & Zeng, JACS Au 2026, 6, 4478 (Table 1, Table S1, eqs 2, 3, 6)',
                                         'sign': 'dlambda/dT = -gamma kappa (du_h(T0) - dh)/T0',
                                         'frames': 'stored FDT frames (npz/300K)'}, [NPZ, BASE])
    print('finished', path, flush=True)


if __name__ == '__main__':
    main()
