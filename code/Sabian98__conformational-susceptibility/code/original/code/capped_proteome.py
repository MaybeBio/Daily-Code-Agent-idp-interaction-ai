"""
capped_proteome.py — the capped hydration forcing (saturating_forcing.py, adopted as the main model 2026-10-24) for every
chain of the three force fields at 300 K, from the stored FDT frames, as tea_fdt.py (CALVADOS-2) and tea_crossff.py
(Mpipi-GG, Mpipi-Recharged) do for the linear forcing.

Per frame, residue i's contact number c_i = -sum_j S_ij / s0 in full contacts, with S the pair shape the linear forcing
scales (CALVADOS-2: U_AH(r; 1) - U_AH(r; 0), s0 = eps_AH; Mpipi: the attractive Wang-Frenkel shape, s0 = 1) over the
model's non-bonded pairs. With f_i the per-residue temperature derivative of the stickiness at 300 K per unit strength
(TEA: tea_fdt.DL['charmm36']; hydrophilic, charged-group form: DL_CLASSES hyd_DE + hyd_K + hyd_R; Mpipi-Recharged in kJ,
x 4.184; Mpipi-GG in kcal), the per-frame forcing is
    linear   D = 1/2 sum_i f_i beta_i = -s0/2 sum_i f_i c_i               (= tea_fdt / tea_crossff, checked per chain)
    capped   D = -s0/2 sum_i f_i h(c_i),   h(c) = c_sat (1 - exp(-c / c_sat)),   c_sat = 2
(the cap does not depend on temperature, so the 300 K derivative of the capped energy only differentiates f).
Per chain and observable (Rg, Ree): dchi = -Cov(A, D) / RT and its % of <A> per K, for the terms
    tea_lin, hyd_lin (full hydrophilic set, the stored A0 term), chg_lin (charged-group, S1),
    tea_cap, chg_cap (the capped model's two terms, per unit strength)
Check: tea_lin and hyd_lin must equal the stored per-chain values (tea_dchi*_classes_300K.csv) to 1e-6 relative.
    BENDER_REANALYSIS_DIR=... python capped_proteome.py --model calvados|gg|recharged [--workers 6] [--limit N]
Output: Traj Level Analysis/Analysis_FDT/capped_300K/capped_dchi_<model>_300K.csv (one row per chain, resumable)
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
from scipy import sparse  # noqa: E402

HERE = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, os.environ.get('BENDER_REANALYSIS_DIR', HERE))
sys.path.insert(0, HERE)
import tea_crossff as X  # noqa: E402
import tea_fdt as TF  # noqa: E402

fdt, bc, L = TF.fdt, TF.bc, TF.L
T0 = TF.T0
CSAT = 2.0
OUT = os.path.join(TF.L.BENDER_MT, 'Traj Level Analysis', 'Analysis_FDT', 'capped_300K')
STORED = {'calvados': 'tea_dchi_classes_300K.csv', 'gg': 'tea_dchi_gg_classes_300K.csv',
          'recharged': 'tea_dchi_recharged_classes_300K.csv'}
TEA_CLASSES = ['aliph', 'pro', 'arom', 'argC36', 'rest']
# per-residue temperature derivatives at T0 per unit strength (kcal/mol/K x mol/kcal = 1/K for CALVADOS-2 lambda units)
F_TEA = TF.DL['charmm36']
F_HYD = TF.DL_CLASSES['hyd']
F_CHG = {a: TF.DL_CLASSES['hyd_DE'][a] + TF.DL_CLASSES['hyd_K'][a] + TF.DL_CLASSES['hyd_R'][a] for a in F_TEA}
TERMS = {'tea': F_TEA, 'hyd': F_HYD, 'chg': F_CHG}


def incidence(n, i, j):
    p = len(i)
    return sparse.csr_matrix((np.ones(2 * p), (np.r_[np.arange(p), np.arange(p)], np.r_[i, j])), shape=(p, n))


def h(c):
    return CSAT * (1.0 - np.exp(-c / CSAT))


def forcings(seq, beta, s0, unit):
    """{term_lin, term_cap: per-frame D} from per-residue shape sums beta (frames x N)."""
    c = -beta / s0
    hc = h(c)
    out = {}
    for k, F in TERMS.items():
        f = np.array([F[a] for a in seq]) * unit
        out[f'{k}_lin'] = 0.5 * beta @ f
        if k != 'hyd':
            out[f'{k}_cap'] = -0.5 * s0 * (hc @ f)
    return out, c


def calvados_one(args):
    uid, d = args
    t0 = time.time()
    try:
        z = np.load(os.path.join(TF.NPZ, f'{uid}.npz'))
        xyz, seq, info = bc.load_calvados(os.path.join(d, f'{uid}.dcd'), os.path.join(d, 'top.pdb'), k=L.K_FRAMES)
        fi = info['production_start'] + info['stride'] * np.arange(len(xyz))
        if not np.array_equal(fi, z['frame_index']):
            raise ValueError('frames differ from the stored FDT frames')
        q, lam, sig, _ = fdt.calvados_residue_params(seq, L.load_residue_params())
        i, j, s, l, shift, _ = fdt.pair_setup(q, lam, sig)
        M = incidence(len(seq), i, j)
        n = len(xyz)
        beta = np.zeros((n, len(seq)))
        step = max(1, int(fdt.MAX_ELEMS // max(1, len(i))))
        for a in range(0, n, step):
            x = xyz[a:a + step].astype(float)
            r = np.linalg.norm(x[:, i] - x[:, j], axis=2)
            S = fdt.ah_pair_energy(r, s, np.ones_like(l), shift) - fdt.ah_pair_energy(r, s, np.zeros_like(l), shift)
            beta[a:a + step] = (M.T @ S.T).T
        D, c = forcings(seq, beta, fdt.EPS_AH, 1.0)
        return finish(uid, seq, {'rg': z['rg'], 'ree': z['ree']}, D, c, t0, fdt.R_GAS), None
    except Exception as ex:
        return {'UniProt_ID': uid}, f'{ex!r} | {traceback.format_exc()[-300:]}'.replace('\n', ' ')


def recharged_one(args):
    uid, d = args
    t0 = time.time()
    try:
        import mpipi_recharged_fdt as R
        z = np.load(os.path.join(X.RC_NPZ, f'{uid}.npz'))
        xyz, seq, info = R.load_run(d, uid)
        fi = info['production_start'] + info['stride'] * np.arange(len(xyz))
        if not np.array_equal(fi, z['frame_index']):
            raise ValueError('frames differ from the stored FDT frames')
        tb = X.rc_tables()
        ix = np.array([tb['index'][a] for a in seq])
        i, j = np.triu_indices(len(seq), k=2)
        sig, mu = tb['sig'][ix[j], ix[i]], tb['mu'][ix[j], ix[i]]
        rc = 3 * sig
        alpha = 2 * 3 ** (2 * mu) * (3 / (2 * (3 ** (2 * mu) - 1))) ** 3
        M = incidence(len(seq), i, j)
        n, Lb = len(xyz), R.box_nm(len(seq))
        beta = np.zeros((n, len(seq)))
        step = max(1, int(X.MAX_ELEMS // max(1, len(i))))
        for a in range(0, n, step):
            dv = xyz[a:a + step, i].astype(float) - xyz[a:a + step, j].astype(float)
            dv -= Lb * np.round(dv / Lb)
            S = X.wf_attractive(np.linalg.norm(dv, axis=2), sig, mu, 1, rc, alpha)
            beta[a:a + step] = (M.T @ S.T).T
        D, c = forcings(seq, beta, 1.0, 4.184)
        return finish(uid, seq, {'rg': z['rg'], 'ree': z['ree']}, D, c, t0, X.R_KJ), None
    except Exception as ex:
        return {'UniProt_ID': uid}, f'{ex!r} | {traceback.format_exc()[-300:]}'.replace('\n', ' ')


def gg_one(args):
    uid, d = args
    t0 = time.time()
    try:
        import mdtraj as md
        import mpipi_temperature_hpc as MH
        z = np.load(os.path.join(X.GG_NPZ, f'{uid}.npz'))
        seq = bc.read_mpipi_sequence(f'{d}/seq.dat')
        with md.formats.XTCTrajectoryFile(f'{d}/traj.xtc') as fh:
            xyz, _, stp, box = fh.read()
        pos = {int(s): k for k, s in enumerate(stp)}
        x = xyz[np.array([pos[int(s)] for s in z['step']])].astype(float)
        L_nm = float(np.median(box[:, 0, 0]))
        pp = MH.mpipi_pairs(MH.read_mpipi_system(d))
        nu, mu, sig, rc = pp['nu'], pp['mu'], pp['sig'], pp['rc']
        alpha = 2 * nu * (rc / sig) ** (2 * mu) * ((1 + 2 * nu) / (2 * nu * ((rc / sig) ** (2 * mu) - 1))) ** (2 * nu + 1)
        M = incidence(len(seq), pp['i'], pp['j'])
        n = len(x)
        beta = np.zeros((n, len(seq)))
        xa, La = x * 10.0, L_nm * 10.0
        step = max(1, int(X.MAX_ELEMS // max(1, len(pp['i']))))
        for a in range(0, n, step):
            dv = xa[a:a + step, pp['i']] - xa[a:a + step, pp['j']]
            dv -= La * np.round(dv / La)
            S = X.wf_attractive(np.linalg.norm(dv, axis=2), sig, mu, nu, rc, alpha)
            beta[a:a + step] = (M.T @ S.T).T
        D, c = forcings(seq, beta, 1.0, 1.0)
        return finish(uid, seq, {'rg': z['rg'], 'ree': z['ree']}, D, c, t0, X.R_KCAL), None
    except Exception as ex:
        return {'UniProt_ID': uid}, f'{ex!r} | {traceback.format_exc()[-300:]}'.replace('\n', ' ')


def finish(uid, seq, obs, D, c, t0, Rgas):
    row = {'UniProt_ID': uid, 'n_residues': len(seq), 'c_median': float(np.median(c)), 'c_p95': float(np.percentile(c, 95)),
           'seconds': round(time.time() - t0, 1)}
    RT = Rgas * T0
    for m, A in obs.items():
        A = np.asarray(A, float)
        ok = np.isfinite(A)
        for k, v in D.items():
            cv = np.cov(A[ok], v[ok], ddof=1)[0, 1]
            row[f'dchi_{m}_{k}'] = -cv / RT
            row[f'dchi_{m}_rel_{k}'] = 100 * (-cv / RT) / A[ok].mean()
        row[f'{m}_mean'] = A[ok].mean()
    for k, v in D.items():                     # spread of the per-frame forcing (factor split, total_response_evidence I)
        row[f'D_mean_{k}'], row[f'D_sd_{k}'] = float(v.mean()), float(v.std(ddof=1))
    return row


RUNNERS = {'calvados': calvados_one, 'recharged': recharged_one, 'gg': gg_one}


def dirs(model):
    if model == 'calvados':
        ids = sorted(f[:-4] for f in os.listdir(TF.NPZ) if f.endswith('.npz') and not f.startswith('._'))
        return TF.index_dirs(ids)
    return X.rc_dirs() if model == 'recharged' else X.gg_dirs()


def check(model, path):
    """tea_lin and hyd_lin against the stored per-chain values."""
    new = pd.read_csv(path).set_index('UniProt_ID')
    old = pd.read_csv(os.path.join(TF.L.BENDER_MT, 'Traj Level Analysis', 'Analysis_FDT', 'tea_300K', STORED[model])) \
        .set_index('UniProt_ID')
    ids = new.index.intersection(old.index)
    tea_old = sum(old.loc[ids, f'dchi_rg_rel_{k}'] for k in TEA_CLASSES)
    d_tea = (new.loc[ids, 'dchi_rg_rel_tea_lin'] - tea_old).abs().max() / max(tea_old.abs().max(), 1e-12)
    d_hyd = (new.loc[ids, 'dchi_rg_rel_hyd_lin'] - old.loc[ids, 'dchi_rg_rel_hyd']).abs().max() / \
        max(old.loc[ids, 'dchi_rg_rel_hyd'].abs().max(), 1e-12)
    return len(ids), float(d_tea), float(d_hyd)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--model', choices=list(RUNNERS), required=True)
    ap.add_argument('--workers', type=int, default=6)
    ap.add_argument('--limit', type=int, default=None)
    ap.add_argument('--check-only', action='store_true')
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    path = os.path.join(OUT, f'capped_dchi_{a.model}_300K.csv')
    if not a.check_only:
        done = set(pd.read_csv(path, usecols=['UniProt_ID']).UniProt_ID) if os.path.exists(path) else set()
        todo = [(u, d) for u, d in dirs(a.model).items() if u not in done][:a.limit]
        print(f'{a.model}: {len(done)} done, {len(todo)} to run', flush=True)
        t0, nf = time.time(), 0
        with mp.get_context('spawn').Pool(a.workers) as pool:
            for k, (row, err) in enumerate(pool.imap_unordered(RUNNERS[a.model], todo), 1):
                if err:
                    nf += 1
                    print(f"failed {row['UniProt_ID']}: {err[:200]}", flush=True)
                    continue
                pd.DataFrame([row]).to_csv(path, mode='a', header=not os.path.exists(path), index=False)
                if k % 250 == 0:
                    el = (time.time() - t0) / 3600
                    print(f'{k}/{len(todo)} {el:.2f} h, ~{el / k * (len(todo) - k):.1f} h left, {nf} failed', flush=True)
    n, dt, dh = check(a.model, path)
    print(f'check on {n} chains: linear TEA term matches the stored values to {dt:.1e}, hydrophilic to {dh:.1e} (relative)',
          flush=True)


if __name__ == '__main__':
    main()
