"""
tea_crossff.py — TEA's hydration dH/dT (tea_fdt.py) applied to the Mpipi-GG and Mpipi-Recharged 300 K ensembles, to
test whether the three force fields agree once the same solvation forcing is added.

TEA gives each residue's change of pair interaction free energy with temperature (Chen & Zeng, JACS Au 2026, eqs 2-3,
Table 1, Table S1, combination rule = arithmetic mean):
    f_i = dG_E,i/dT sign-flipped = -kappa_i (du_h,i(T0) - dh_i) / T0      kcal/mol/K, > 0 = stickier on heating
For the Wang-Frenkel models the stickiness is the pair well depth eps_ij (an energy), so the forcing is mapped as
    deps_ij/dT = gamma' (f_i + f_j) / 2,       gamma' dimensionless (one scale per model; Mpipi-Recharged calibrated
                                               against Wuttke 2014, Mpipi-GG set at the others' relative size)
applied to the attractive part of U_WF only (-1 inside the minimum, the Wang-Frenkel tail beyond it, 0 past rc; fixed
2026-10-10, the full potential including the repulsive wall had been scaled), so per frame D = sum_ij deps_ij/dT *
att_ij(r) and, per unit gamma',
    dchi(A) = -Cov(A, D) / (R T).
TEA keeps each model's parameters at T0 = 300 K, so the stored 300 K ensembles are the TEA-augmented ensembles at
300 K and dchi is that model's exact extra slope there. Because dchi is linear in gamma', the cross-model test uses
scale-free quantities (rank correlations of dchi; switch points relative to each other).

Frames: the stored FDT frames (Mpipi-GG npz 'step'; Mpipi-Recharged npz 'frame_index'); the recomputed WF energy is
checked against the stored one for every protein.

    python tea_crossff.py --model gg|recharged --out DIR [--workers N] [--limit N]
Output: DIR/tea_dchi_<model>_300K.csv (dchi per unit gamma' for Rg and Ree, absolute and % of <A> per K, Arg kappa
CHARMM36m and CGenFF; the WF check).
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
import tea_fdt as TF  # noqa: E402  (TEA tables, T0)

bc = TF.bc
T0 = TF.T0
R_KJ, R_KCAL = 8.314462618e-3, 1.98720425864083e-3
TLA = '/Volumes/User Homes/Nu Project/BENDER-MT/Traj Level Analysis'
GG_NPZ = os.path.join(TLA, 'Analysis_FDT_MPIPI', 'npz_300K')
GG_RUNS = ['/Volumes/User Homes/Nu Project/BENDER-BIO/MPIPI/mpipi_sim/300k_runs/Mpipi 300k ',
           '/Volumes/User Homes/Nu Project/BENDER-BIO/MPIPI/mpipi_sim/300k_runs/Mpipi SAXS']
RC_NPZ = os.path.join(TLA, 'Analysis_FDT_RECHARGED', 'npz_300K')
RC_RUNS = '/Volumes/User Homes/Nu Project/MPIPI Recharged/300K'
MAX_ELEMS = 4e6


def forcing(arg: str) -> dict:
    """f_i in kcal/mol/K per residue (positive = stickier on heating); TF.DL already holds -kappa (du - dh)/T0."""
    return TF.DL[arg]


F = {'charmm36': forcing('charmm36'), 'cgenff': forcing('cgenff')}
F_MODES = {'base': F, 'classes': TF.DL_CLASSES}


def wf_attractive(r, sig, mu, nu, rc, alpha):
    """Attractive part of the Wang-Frenkel shape (well depth -1 inside the minimum, tail beyond, 0 past rc), the split
    CALVADOS uses for Ashbaugh-Hatch. The hydration forcing scales this part only: scaling the full potential (eps also
    multiplies the steep repulsive wall, mu up to 12) puts the forcing on near-neighbour wall contacts (fixed
    2026-10-10; earlier Mpipi TEA terms used the full potential)."""
    rmin = rc * ((1 + 2 * nu) / (1 + 2 * nu * (rc / sig) ** (2 * mu))) ** (1 / (2 * mu))
    full = alpha * ((sig / r) ** (2 * mu) - 1) * ((rc / r) ** (2 * mu) - 1) ** (2 * nu)
    return np.where(r < rmin, -1.0, np.where(r < rc, full, 0.0))


OBSERVABLES = ['rg', 'ree', 'delta', 'fragmentation_index', 'global_efficiency', 'avg_clustering', 'transitivity',
               'degree_assortativity']


def observables(z):
    """Rg and Ree; with TEA_ALL_OBS=1 (--all-observables) every per-frame observable stored in the npz."""
    keys = [m for m in OBSERVABLES if m in z.files] if os.environ.get('TEA_ALL_OBS') == '1' else ['rg', 'ree']
    return {m: z[m] for m in keys}


def finish(uid, seq, A_obs, D, check, t0, energy_unit_R):
    row = {'UniProt_ID': uid, 'n_residues': len(seq), 'WF_check_rel': check, 'seconds': round(time.time() - t0, 1)}
    RT = energy_unit_R * T0
    for m, A in A_obs.items():
        A = np.asarray(A, dtype=float)
        ok = np.isfinite(A)              # degree assortativity is undefined for some frames
        for k in D:
            c = np.cov(A[ok], D[k][ok], ddof=1)[0, 1]
            row[f'dchi_{m}_{k}'] = -c / RT
            row[f'dchi_{m}_rel_{k}'] = 100 * (-c / RT) / A[ok].mean()
        row[f'{m}_mean'] = A[ok].mean()
    for k in D:
        row[f'D_mean_{k}'], row[f'D_sd_{k}'] = D[k].mean(), D[k].std(ddof=1)
    return row


# ------------------------------- Mpipi-GG (LAMMPS) -------------------------------
def gg_one(args):
    uid, d, mode = args
    F = F_MODES[mode]
    t0 = time.time()
    try:
        import mdtraj as md
        import mpipi_temperature_hpc as MH
        z = np.load(os.path.join(GG_NPZ, f'{uid}.npz'))
        seq = bc.read_mpipi_sequence(f'{d}/seq.dat')
        with md.formats.XTCTrajectoryFile(f'{d}/traj.xtc') as fh:
            xyz, _, step, box = fh.read()
        pos = {int(s): k for k, s in enumerate(step)}
        idx = np.array([pos[int(s)] for s in z['step']])
        x = xyz[idx].astype(float)
        L_nm = float(np.median(box[:, 0, 0]))
        sysd = MH.read_mpipi_system(d)
        pp = MH.mpipi_pairs(sysd)
        if len(seq) != len(sysd['types']):
            raise ValueError('sequence and seq.dat atoms differ')
        nu, mu, sig, rc = pp['nu'], pp['mu'], pp['sig'], pp['rc']
        alpha = 2 * nu * (rc / sig) ** (2 * mu) * ((1 + 2 * nu) / (2 * nu * ((rc / sig) ** (2 * mu) - 1))) ** (2 * nu + 1)
        aa = np.array(list(seq))
        de = {k: 0.5 * (np.array([v[a] for a in aa])[pp['i']] + np.array([v[a] for a in aa])[pp['j']]) for k, v in F.items()}
        n = len(x)
        D = {k: np.zeros(n) for k in F}
        E_WF = np.zeros(n)
        xa, L = x * 10.0, L_nm * 10.0
        stepn = max(1, int(MAX_ELEMS // max(1, len(pp['i']))))
        for a in range(0, n, stepn):
            dv = xa[a:a + stepn, pp['i']] - xa[a:a + stepn, pp['j']]
            dv -= L * np.round(dv / L)
            r = np.linalg.norm(dv, axis=2)
            shape = np.where(r < rc, alpha * ((sig / r) ** (2 * mu) - 1) * ((rc / r) ** (2 * mu) - 1) ** (2 * nu), 0.0)
            E_WF[a:a + stepn] = (shape * (pp['eps_alpha'] / alpha)).sum(1)      # full potential: the stored-energy check
            att = wf_attractive(r, sig, mu, nu, rc, alpha)
            for k in D:
                D[k][a:a + stepn] = (att * de[k]).sum(1)            # kcal/mol/K per unit gamma'
        check = np.abs(E_WF - z['E_WF']).max() / max(1.0, np.abs(z['E_WF']).max())
        return finish(uid, seq, observables(z), D, check, t0, R_KCAL), None
    except Exception as ex:
        return {'UniProt_ID': uid}, f'{ex!r} | {traceback.format_exc()[-300:]}'.replace('\n', ' ')


def gg_dirs():
    ids = {f[:-4] for f in os.listdir(GG_NPZ) if f.endswith('.npz') and not f.startswith('._')}
    out = {}
    for root in GG_RUNS:
        if not os.path.isdir(root):
            continue
        for u in os.listdir(root):
            if u in ids and u not in out and os.path.isfile(os.path.join(root, u, 'traj.xtc')):
                out[u] = os.path.join(root, u)
    return out


# --------------------------- Mpipi-Recharged (OpenMpipi) ---------------------------
_RC = None


def rc_tables():
    global _RC
    if _RC is None:
        import OpenMpipi
        P = np.loadtxt(os.path.join(os.path.dirname(OpenMpipi.__file__), 'data_files', 'recharged_params.txt'))
        w = P[:-21 * 21].reshape(3, 21, 21)          # Discrete3DFunction: value[x + 21 y + 441 z] -> w[z, y, x]
        order = 'MGKTRADEYVLQWFSHNPCI'                  # mapping_dict indices 0..19 (20 = RNA U)
        _RC = {'eps': w[0], 'sig': w[1], 'mu': np.floor(w[2]), 'index': {a: k for k, a in enumerate(order)}}
    return _RC


def rc_one(args):
    uid, d, mode = args
    F = F_MODES[mode]
    t0 = time.time()
    try:
        import mpipi_recharged_fdt as R
        z = np.load(os.path.join(RC_NPZ, f'{uid}.npz'))
        xyz, seq, info = R.load_run(d, uid)
        fi = info['production_start'] + info['stride'] * np.arange(len(xyz))
        if not np.array_equal(fi, z['frame_index']):
            raise ValueError('frames differ from the stored FDT frames')
        tb = rc_tables()
        ix = np.array([tb['index'][a] for a in seq])
        i, j = np.triu_indices(len(seq), k=2)              # bonded (i, i+1) pairs excluded, as in OpenMpipi
        eps, sig, mu = tb['eps'][ix[j], ix[i]], tb['sig'][ix[j], ix[i]], tb['mu'][ix[j], ix[i]]
        rc = 3 * sig
        alpha = 2 * 3 ** (2 * mu) * (3 / (2 * (3 ** (2 * mu) - 1))) ** 3
        aa = np.array(list(seq))
        de = {k: 0.5 * (np.array([v[a] for a in aa])[i] + np.array([v[a] for a in aa])[j]) * 4.184 for k, v in F.items()}
        n, L = len(xyz), R.box_nm(len(seq))
        D = {k: np.zeros(n) for k in F}
        E_WF = np.zeros(n)
        stepn = max(1, int(MAX_ELEMS // max(1, len(i))))
        for a in range(0, n, stepn):
            dv = xyz[a:a + stepn, i].astype(float) - xyz[a:a + stepn, j].astype(float)
            dv -= L * np.round(dv / L)
            r = np.linalg.norm(dv, axis=2)
            shape = np.where(r < rc, alpha * ((sig / r) ** (2 * mu) - 1) * ((rc / r) ** (2 * mu) - 1) ** 2, 0.0)
            E_WF[a:a + stepn] = (shape * eps).sum(1)                             # full potential: the stored-energy check
            att = wf_attractive(r, sig, mu, 1, rc, alpha)
            for k in D:
                D[k][a:a + stepn] = (att * de[k]).sum(1)            # kJ/mol/K per unit gamma'
        check = np.abs(E_WF - z['V_AH']).max() / max(1.0, np.abs(z['V_AH']).max())
        return finish(uid, seq, observables(z), D, check, t0, R_KJ), None
    except Exception as ex:
        return {'UniProt_ID': uid}, f'{ex!r} | {traceback.format_exc()[-300:]}'.replace('\n', ' ')


def rc_dirs():
    ids = sorted(f[:-4] for f in os.listdir(RC_NPZ) if f.endswith('.npz') and not f.startswith('._'))
    return {u: os.path.join(RC_RUNS, u) for u in ids}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--model', required=True, choices=['gg', 'recharged'])
    ap.add_argument('--out', required=True)
    ap.add_argument('--workers', type=int, default=8)
    ap.add_argument('--limit', type=int, default=None)
    ap.add_argument('--classes', action='store_true', help='split the TEA term by residue class')
    ap.add_argument('--ids', default=None, help='file with one UniProt ID per line (default: all)')
    ap.add_argument('--all-observables', action='store_true', help='every stored per-frame observable (sets TEA_ALL_OBS=1)')
    a = ap.parse_args()
    if a.all_observables:
        os.environ['TEA_ALL_OBS'] = '1'   # inherited by the spawned workers
    os.makedirs(a.out, exist_ok=True)
    mode = 'classes' if a.classes else 'base'
    path = os.path.join(a.out, f'tea_dchi_{a.model}' + ('_classes' if a.classes else '') + '_300K.csv')
    dirs, fn = (gg_dirs(), gg_one) if a.model == 'gg' else (rc_dirs(), rc_one)
    if a.ids:
        want = {l.strip() for l in open(a.ids) if l.strip()}
        dirs = {u: d for u, d in dirs.items() if u in want}
    done = set(pd.read_csv(path, usecols=['UniProt_ID']).UniProt_ID) if os.path.exists(path) else set()
    todo = [(u, d, mode) for u, d in dirs.items() if u not in done][:a.limit]
    print(f'{a.model}: {len(dirs)} runs, {len(done)} done, {len(todo)} to run', flush=True)
    t0, nf = time.time(), 0
    with mp.get_context('spawn').Pool(a.workers) as pool:
        for k, (row, err) in enumerate(pool.imap_unordered(fn, todo), 1):
            if err:
                nf += 1
                print(f"failed {row['UniProt_ID']}: {err[:220]}", flush=True)
                continue
            pd.DataFrame([row]).to_csv(path, mode='a', header=not os.path.exists(path), index=False)
            if k % 250 == 0:
                el = (time.time() - t0) / 3600
                print(f'{k}/{len(todo)} {el:.2f} h, ~{el / k * (len(todo) - k):.1f} h left, {nf} failed', flush=True)
    bc.write_provenance(path, __file__, {'model': a.model, 'T0_K': T0, 'forcing': 'deps_ij/dT = gamma\' (f_i + f_j)/2, '
                                         'f_i = -kappa_i (du_h(T0) - dh_i)/T0 (TEA, JACS Au 2026)'}, [GG_NPZ if a.model == 'gg' else RC_NPZ])
    print('finished', path, flush=True)


if __name__ == '__main__':
    main()
