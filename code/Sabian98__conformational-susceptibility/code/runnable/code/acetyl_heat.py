"""
acetyl_heat.py — heat x acetylation: how neutralizing one Lys changes the calibrated total heating response, by exact
reweighting of the 300 K frames (part C of heat_crosses.py / heat_crosses_gg.py with Lys acetylation in place of
Ser/Thr phosphorylation).

Charge model, as for phosphorylation: acetylation removes the side-chain charge of Lys i (dq = -1 in CALVADOS-2, -0.75 in
Mpipi-GG, whose charges are scaled by 0.75) and keeps its size, stickiness and hydration. The Debye-Hueckel energy is
bilinear in the charges, so per frame
  CALVADOS-2   dH = dq phi_i,  dG = -dH / (R T^2) + dq (dphi_i/dT) / (R T)     (temperature-dependent screening)
  Mpipi-GG     dH = dq phi_i,  dG = -dH / (R T^2)
new total = -Cov_w(Rg, G_tot + dG) / <Rg>_w with w = exp(-dH/RT), G_tot = G + (gamma_tea D_tea + gamma_hyd D_hyd) / RT;
reported as the change in % per K with a 10-block jackknife SE and the Kish effective sample fraction of w.
Mpipi-Recharged has pair-specific Yukawa amplitudes A_ij instead of charges, and every neutral residue has A = 0 with all
partners, so acetyl-Lys is Lys with its Yukawa amplitudes set to zero and its Wang-Frenkel parameters and hydration kept
(energy terms and solvation forcing as mutation_scan_recharged.py):
  Mpipi-Recharged   dH = -Y_i,  dG = -dH / (R T^2) + (dkappa/dT) E_i / (R T)
with Y_i = sum_j A_ij exp(-kappa r) / r and E_i = sum_j A_ij exp(-kappa r) over the partners j of Lys i (r < 3.5 nm).
Setting MR_MUT to a residue letter instead substitutes that residue (validation against mutation_scan_recharged.py).
Chains: the 405 CALVADOS-2 chains of heat_crosses.py (its sample.csv; also used for Mpipi-Recharged) and the 211 shared
chains in Mpipi-GG.
    python acetyl_heat.py --model C|M|R --out DIR [--workers N]
Output: DIR/sites_{C,M,R}.csv, DIR/chains_{C,M,R}.csv (append; resumable); since 2026-10-19 also the change of the
        fixed-solvent susceptibility (acetyl_dchi, chains: chi_rel), first written to heat_acetyl/fixed_solvent/
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
import heat_crosses as HC  # noqa: E402
import capped_frames as CF  # noqa: E402

TR = os.path.join(HC.TLA, 'Analysis_FDT', 'total_response')
IDS = {'C': os.path.join(TR, 'heat_crosses', 'sample.csv'), 'M': os.path.join(TR, 'heat_crosses_shared', 'ids.txt'),
       'R': os.path.join(TR, 'heat_crosses', 'sample.csv')}
DQ = {'C': -1.0, 'M': -0.75}
MR_MUT = None   # Mpipi-Recharged: None = acetylation (Yukawa of the Lys off); a letter = that substitution instead
TARGET = 'K'   # residues perturbed in turn (set to 'ST' with DQ -2 / -1.5 to reproduce the phosphorylation results)


def sites_from(uid, seq, rg, Gtot, phi, dG_of, RT, Gfix=None):
    """Per-Lys change of the total response; the baseline uses ddof 0 as wcov does (heat_crosses.py fix of 2026-10-10).
    With Gfix (the thermal force without the hydration forcings) also the change of the fixed-solvent susceptibility,
    by the same reweighting (acetylation keeps the hydration of Lys, so dG has no forcing part)."""
    nf, m = len(rg), rg.mean()

    def perturbed(dH, dG, Gb):
        x = -dH / RT
        w = np.exp(x - x.max())
        f = lambda ix: (-HC.wcov(rg[ix], Gb[ix] + dG[ix], w[ix]) / ((w[ix] * rg[ix]).sum() / w[ix].sum())  # noqa: E731
                        + np.cov(rg[ix], Gb[ix], ddof=0)[0, 1] / rg[ix].mean()) * 100
        return f(np.arange(nf)), HC.jk(f, nf), float(w.sum() ** 2 / (w ** 2).sum() / nf), \
            float(((w * rg).sum() / w.sum() / m - 1) * 100)
    rows = []
    for p_, c in enumerate(seq):
        if c in TARGET:
            dH, dG = dG_of(p_)
            dt, se, ess, drg = perturbed(dH, dG, Gtot)
            r = {'UniProt_ID': uid, 'pos': p_ + 1, 'residue': c, 'acetyl_dtotal': dt, 'acetyl_dtotal_se': se,
                 'ess': ess, 'acetyl_dRg_pct': drg}
            if Gfix is not None:
                r['acetyl_dchi'], r['acetyl_dchi_se'] = perturbed(dH, dG, Gfix)[:2]
            rows.append(r)
    base = -np.cov(rg, Gtot)[0, 1] / m * 100
    row = {'UniProt_ID': uid, 'N': len(seq), 'n_K': seq.count('K'), 'total_rel': base}
    if Gfix is not None:
        row['chi_rel'] = -np.cov(rg, Gfix)[0, 1] / m * 100
    return row, pd.DataFrame(rows)


def one_calvados(args):
    uid, d = args
    t0 = time.time()
    try:
        TF, fdt, L, bc = HC.TF, HC.fdt, HC.L, HC.bc
        z = np.load(os.path.join(TF.NPZ, f'{uid}.npz'))
        xyz_all, seq, info = bc.load_calvados(os.path.join(d, f'{uid}.dcd'), os.path.join(d, 'top.pdb'), k=None)
        idx = bc.even_indices(len(xyz_all), L.K_FRAMES)
        if not np.array_equal(info['production_start'] + idx, z['frame_index']):
            raise ValueError('frames differ from the stored FDT frames')
        xyz, rg = xyz_all[idx], z['rg'].astype(float)
        rg_chk = np.sqrt(((xyz - xyz.mean(1, keepdims=True)) ** 2).sum(2).mean(1))
        if np.abs(rg_chk - rg).max() > 1e-3 * rg.mean():
            raise ValueError('Rg differs from the npz')
        N, nf = len(seq), len(xyz)
        q, lam, sig, _ = fdt.calvados_residue_params(seq, L.load_residue_params())
        i, j, s, l, shift, _ = fdt.pair_setup(q, lam, sig)
        npair = len(i)
        qmat = sparse.csr_matrix((np.concatenate([q[j], q[i]]), (np.concatenate([i, j]), np.concatenate([np.arange(npair)] * 2))), shape=(N, npair))
        dl = {k: np.array([v[a] for a in seq]) for k, v in (('tea', TF.DL['charmm36']), ('hyd', TF.DL_CLASSES['hyd']))}
        fpair = HC.G_TEA * 0.5 * (dl['tea'][i] + dl['tea'][j]) + HC.G_HYD * 0.5 * (dl['hyd'][i] + dl['hyd'][j])
        CAL, rc = HC.CAL, fdt.RC_YU
        k1, dk1, eps, deps = CAL['kappa'], CAL['dkappa_dT'], CAL['eps_yu'], CAL['deps_yu_dT']
        qq = q[i] * q[j]
        D, V1 = np.zeros(nf), np.zeros(nf)
        phi, dphi = np.zeros((nf, N)), np.zeros((nf, N))
        step = max(1, int(fdt.MAX_ELEMS // max(1, npair)))
        for a in range(0, nf, step):
            r = np.linalg.norm(xyz[a:a + step, i] - xyz[a:a + step, j], axis=2)
            u1, u0 = fdt.ah_pair_energy(r, s, np.ones_like(l), shift), fdt.ah_pair_energy(r, s, np.zeros_like(l), shift)
            inside = r < rc
            er, ec = np.exp(-k1 * r), np.exp(-k1 * rc)
            Y = np.where(inside, er / r - ec / rc, 0.0)
            dY = deps * Y + eps * dk1 * np.where(inside, ec - er, 0.0)
            V1[a:a + step] = (qq * eps * Y).sum(1)
            D[a:a + step] = ((u1 - u0) * fpair).sum(1)
            phi[a:a + step] = (qmat @ (eps * Y).T).T
            dphi[a:a + step] = (qmat @ dY.T).T
        vz = z['V_DH'].astype(float)
        if np.abs(V1 - vz).max() > 1e-3 * max(1.0, np.abs(vz).max()):
            raise ValueError('V_DH differs from the npz')
        Gtot = z['G_T'].astype(float) + (CF.term('C', uid, rg, d) if HC.FORCING == 'capped' else D / HC.RT)
        dq = DQ['C']

        def dG_of(p_):
            dH = dq * phi[:, p_]
            return dH, -dH / (fdt.R_GAS * HC.T0 ** 2) + dq * dphi[:, p_] / HC.RT
        row, S = sites_from(uid, seq, rg, Gtot, phi, dG_of, HC.RT, Gfix=z['G_T'].astype(float))
        row['seconds'] = round(time.time() - t0, 1)
        return row, S, None
    except Exception as ex:  # noqa: BLE001
        return {'UniProt_ID': uid}, None, f'{uid}: {ex!r} | {traceback.format_exc()[-300:]}'.replace('\n', ' ')


def one_gg(uid):
    t0 = time.time()
    try:
        import heat_crosses_gg as HG
        import total_response_noise as N
        SP, bs, mh = HG.SP, HG.bs, HG.mh
        xyz, G, rg = bs.frames_mpipi(uid, 300)
        rg2, G2, Dh, Dt, _ = N.mpipi(uid, N.TC.gg_dirs()[uid], 'M')
        if not (np.allclose(rg, rg2) and np.allclose(G, G2)):
            raise ValueError('frames differ between loaders')
        Gtot = G + (CF.term('M', uid, rg) if HC.FORCING == 'capped' else (N.GAMMA['M'] * Dh + N.GAMMA_TEA['M'] * Dt) / HG.RT)
        z = np.load(os.path.join(bs.MP_DIR, 'npz_300K', f'{uid}.npz'))
        seq = str(z['sequence'])
        sysd = mh.read_mpipi_system(os.path.join(bs.MP_RUNS[300], uid))
        pp = mh.mpipi_pairs(sysd)
        i, j = pp['i'], pp['j']
        q = np.array([sysd['qtype'].get(int(t), 0.0) for t in sysd['types']])
        qK = q[np.array(list(seq)) == 'K']
        if TARGET == 'K' and len(qK) and not np.allclose(qK, -DQ['M']):
            raise ValueError(f'unexpected Lys charge {sorted(set(np.round(qK, 3)))}')
        rc_all, rc_c = max(sysd['coul'].values()), pp['rc_c']
        M = SP.site_matrix(i, j, q, len(q))
        k, pref = sysd['kappa'], mh.QQR2E / sysd['diel']
        nf = len(xyz)
        E1, phi = np.zeros(nf), np.zeros((nf, len(q)))
        step = max(1, SP.MAX_ELEMS // max(1, len(i)))
        for a in range(0, nf, step):
            r = 10.0 * np.linalg.norm(xyz[a:a + step, i] - xyz[a:a + step, j], axis=2)
            E1[a:a + step] = np.where(r < rc_c, pref * pp['qq'] * np.exp(-k * r) / r, 0.0).sum(1)
            phi[a:a + step] = (M @ np.where(r < rc_all, pref * np.exp(-k * r) / r, 0.0).T).T
        ez = z['E_DH'].astype(float)
        if np.abs(ez).max() > 1e-6 and (np.corrcoef(E1, ez)[0, 1] < 0.999 or np.abs(E1 - ez).max() > 1e-4 * np.abs(ez).max()):
            raise ValueError('recomputed E_DH does not match the npz')
        dq = DQ['M']

        def dG_of(p_):
            dH = dq * phi[:, p_]
            return dH, -dH / (HG.R_KCAL * HG.T0 ** 2)
        row, S = sites_from(uid, seq, rg, Gtot, phi, dG_of, HG.RT, Gfix=G)
        row['seconds'] = round(time.time() - t0, 1)
        return row, S, None
    except Exception as ex:  # noqa: BLE001
        return {'UniProt_ID': uid}, None, f'{uid}: {ex!r} | {traceback.format_exc()[-300:]}'.replace('\n', ' ')


def one_recharged(args):
    uid, d = args
    t0 = time.time()
    try:
        import mpipi_recharged_fdt as R
        import mutation_scan_recharged as MS
        TC = MS.TC
        z = np.load(os.path.join(TC.RC_NPZ, f'{uid}.npz'))
        xyz, seq, info = R.load_run(d, uid)
        if not np.array_equal(info['production_start'] + info['stride'] * np.arange(len(xyz)), z['frame_index']):
            raise ValueError('frames differ from the stored FDT frames')
        tb = MS.tables()
        k = MS.kappa_of(MS.T0)
        dk = (MS.kappa_of(MS.T0 + 0.5) - MS.kappa_of(MS.T0 - 0.5)) / 1.0
        rg, G = z['rg'].astype(float), z['G_T'].astype(float)
        N, nf, L = len(seq), len(xyz), R.box_nm(len(seq))
        ix = np.array([tb['ix'][a] for a in seq])
        ft = {a: MS.G_TEA * MS.TF.DL['charmm36'][a] + MS.G_HYD * MS.TF.DL_CLASSES['hyd'][a] for a in MS.AA}
        fv = np.array([ft[a] for a in seq])
        i, j = np.triu_indices(N, k=2)
        D = np.zeros(nf)
        step = max(1, int(TC.MAX_ELEMS // max(1, len(i))))
        for a in range(0, nf, step):
            dv = xyz[a:a + step, i].astype(float) - xyz[a:a + step, j].astype(float)
            dv -= L * np.round(dv / L)
            D[a:a + step] = MS.pair_terms(np.linalg.norm(dv, axis=2), ix[i], ix[j], tb, k, fv[i], fv[j])[3]
        Gtot = G + (CF.term('R', uid, rg, d) if HC.FORCING == 'capped' else D / MS.RT)
        rgk = {}

        def dG_of(p_):
            other = np.array([q_ for q_ in range(N) if abs(q_ - p_) >= 2])
            dv = xyz[:, [p_]].astype(float) - xyz[:, other].astype(float)
            dv -= L * np.round(dv / L)
            r = np.linalg.norm(dv, axis=2)
            W0, Y0, E0, D0 = MS.pair_terms(r, np.full(len(other), ix[p_]), ix[other], tb, k, np.full(len(other), fv[p_]), fv[other])
            if MR_MUT is not None and HC.FORCING == 'capped':
                raise NotImplementedError('MR_MUT validation uses the linear forcing change; use mutation_scan_recharged.py')
            if MR_MUT is None:
                W1, Y1, E1, D1 = W0, 0.0 * Y0, 0.0 * E0, D0
            else:
                W1, Y1, E1, D1 = MS.pair_terms(r, np.full(len(other), tb['ix'][MR_MUT]), ix[other], tb, k,
                                               np.full(len(other), ft[MR_MUT]), fv[other])
            dH = (W1 - W0) + (Y1 - Y0)
            return dH, -dH / (MS.R_KJ * MS.T0 ** 2) + (-dk * (E1 - E0)) / MS.RT + (D1 - D0) / MS.RT
        row, S = sites_from(uid, seq, rg, Gtot, None, dG_of, MS.RT, Gfix=G)
        row['seconds'] = round(time.time() - t0, 1)
        return row, S, None
    except Exception as ex:  # noqa: BLE001
        return {'UniProt_ID': uid}, None, f'{uid}: {ex!r} | {traceback.format_exc()[-300:]}'.replace('\n', ' ')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--model', choices=['C', 'M', 'R'], required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--workers', type=int, default=6)
    ap.add_argument('--limit', type=int, default=None)
    ap.add_argument('--forcing', choices=['A0', 'capped'], default='A0')
    a = ap.parse_args()
    os.environ['BENDER_SCAN_FORCING'] = a.forcing           # capped baseline in the spawned workers (heat_crosses)
    os.makedirs(a.out, exist_ok=True)
    src = IDS[a.model]
    ids = sorted(pd.read_csv(src).UniProt_ID) if src.endswith('.csv') else sorted({l.strip() for l in open(src) if l.strip()})
    pc = os.path.join(a.out, f'chains_{a.model}.csv')
    done = set(pd.read_csv(pc, usecols=['UniProt_ID']).UniProt_ID) if os.path.exists(pc) else set()
    if a.model == 'C':
        todo, fn = list(HC.TF.index_dirs([u for u in ids if u not in done]).items()), one_calvados
    elif a.model == 'R':
        import tea_crossff as TC
        dirs = TC.rc_dirs()
        todo, fn = [(u, dirs[u]) for u in ids if u not in done and u in dirs], one_recharged
        print(f'{len([u for u in ids if u not in dirs])} chains without a Mpipi-Recharged run', flush=True)
    else:
        todo, fn = [u for u in ids if u not in done], one_gg
    todo = todo[:a.limit]
    print(f'{len(ids)} chains, {len(done)} done, {len(todo)} to run', flush=True)
    nfail = 0
    with mp.get_context('spawn').Pool(a.workers) as pool:
        for k, (row, S, err) in enumerate(pool.imap_unordered(fn, todo), 1):
            if err:
                nfail += 1
                print('failed', err[:300], flush=True)
                continue
            pd.DataFrame([row]).to_csv(pc, mode='a', header=not os.path.exists(pc), index=False)
            ps = os.path.join(a.out, f'sites_{a.model}.csv')
            if len(S):
                S.to_csv(ps, mode='a', header=not os.path.exists(ps), index=False)
            if k % 25 == 0:
                print(f'{k}/{len(todo)} done, {nfail} failed', flush=True)
    print('finished', flush=True)


if __name__ == '__main__':
    main()
