"""
forcing_role.py — what the hydration forcings change in the compaction rule, for phosphorylation and doubled salt
(acetylation: acetyl_heat.py, whose fixed-solvent output is in heat_acetyl/fixed_solvent/).

For every perturbation the same exact reweighting of the 300 K frames is applied twice, once with the calibrated thermal
force G_tot = G + (gamma_tea D_tea + gamma_hyd D_hyd)/RT (the total response, as in the paper) and once with the plain
thermal force G (the fixed-solvent susceptibility, the force fields without the forcings):
  change = -Cov_w(Rg, G_base + dG) / <Rg>_w + Cov(Rg, G_base) / <Rg>     (% per K), w = exp(-dH/RT)
The perturbation energies and dG are those of the paper (they contain no forcing term, since neither a phosphate charge
nor a change of screening changes the hydration of any residue):
  phosphorylation  CALVADOS-2 dq -2, Mpipi-GG dq -1.5 (per Ser/Thr; acetyl_heat.py with TARGET 'ST')
  doubled salt     kappa -> sqrt(2) kappa: CALVADOS-2 as heat_crosses.py, Mpipi-GG as heat_crosses_gg.py,
                   Mpipi-Recharged as recharged_traj_pass.py
The calibrated changes reproduce the stored ones (heat_crosses/, heat_crosses_shared/, traj_pass/), which is checked in
forcing_role_analysis.py. Chains: the 405 CALVADOS-2 chains of heat_crosses (also used for Mpipi-Recharged salt), the
211 shared Mpipi-GG chains.
    python forcing_role.py --what phospho|salt --model C|M|R --out DIR [--workers N]
Output: DIR/{phospho,salt}_{C,M,R}.csv (append; resumable)
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
import acetyl_heat as AH  # noqa: E402
import capped_frames as CF  # noqa: E402

HC = AH.HC
PHOS_DQ = {'C': -2.0, 'M': -1.5}


def change(rg, Gb, dH, dG, RT):
    """Change of -Cov(Rg, G_base)/<Rg> on reweighting (% per K), jackknife SE, Kish fraction, change of <Rg> (%)."""
    nf, m = len(rg), rg.mean()
    x = -dH / RT
    w = np.exp(x - x.max())
    f = lambda ix: (-HC.wcov(rg[ix], Gb[ix] + dG[ix], w[ix]) / ((w[ix] * rg[ix]).sum() / w[ix].sum())  # noqa: E731
                    + np.cov(rg[ix], Gb[ix], ddof=0)[0, 1] / rg[ix].mean()) * 100
    return f(np.arange(nf)), HC.jk(f, nf), float(w.sum() ** 2 / (w ** 2).sum() / nf), float(((w * rg).sum() / w.sum() / m - 1) * 100)


def chain_row(uid, rg, Gtot, G):
    m = rg.mean()
    return {'UniProt_ID': uid, 'total_rel': -np.cov(rg, Gtot)[0, 1] / m * 100, 'chi_rel': -np.cov(rg, G)[0, 1] / m * 100}


def salt_row(uid, rg, Gtot, G, dH, dG, RT, t0):
    row = chain_row(uid, rg, Gtot, G)
    row['salt_dtotal'], row['salt_dtotal_se'], row['salt_ess'], row['salt_dRg_pct'] = change(rg, Gtot, dH, dG, RT)
    row['salt_dchi'], row['salt_dchi_se'] = change(rg, G, dH, dG, RT)[:2]
    row['seconds'] = round(time.time() - t0, 1)
    return row


# ---------------------------------------------------------------- phosphorylation (reuses acetyl_heat.py per chain)
def _phospho(fn, model, args):
    AH.TARGET, AH.DQ[model] = 'ST', PHOS_DQ[model]
    row, S, err = fn(args)
    if S is not None and len(S):
        S = S.rename(columns=lambda c: c.replace('acetyl_', 'phospho_'))
    return row, S, err


def phospho_C(args):
    return _phospho(AH.one_calvados, 'C', args)


def phospho_M(uid):
    return _phospho(AH.one_gg, 'M', uid)


# ---------------------------------------------------------------- doubled salt
def salt_C(args):
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
        nf = len(xyz)
        q, lam, sig, _ = fdt.calvados_residue_params(seq, L.load_residue_params())
        i, j, s, l, shift, _ = fdt.pair_setup(q, lam, sig)
        npair = len(i)
        dl = {k: np.array([v[a] for a in seq]) for k, v in (('tea', TF.DL['charmm36']), ('hyd', TF.DL_CLASSES['hyd']))}
        fpair = HC.G_TEA * 0.5 * (dl['tea'][i] + dl['tea'][j]) + HC.G_HYD * 0.5 * (dl['hyd'][i] + dl['hyd'][j])
        CAL = HC.CAL
        eps, deps, rc = CAL['eps_yu'], CAL['deps_yu_dT'], fdt.RC_YU
        qq = q[i] * q[j]
        D, V1, V2, dV1, dV2 = (np.zeros(nf) for _ in range(5))
        step = max(1, int(fdt.MAX_ELEMS // max(1, npair)))
        for a in range(0, nf, step):
            r = np.linalg.norm(xyz[a:a + step, i] - xyz[a:a + step, j], axis=2)
            u1, u0 = fdt.ah_pair_energy(r, s, np.ones_like(l), shift), fdt.ah_pair_energy(r, s, np.zeros_like(l), shift)
            inside = r < rc
            D[a:a + step] = ((u1 - u0) * fpair).sum(1)
            for kk, dkk, V, dV in ((CAL['kappa'], CAL['dkappa_dT'], V1, dV1),
                                   (CAL['kappa'] * np.sqrt(2.0), CAL['dkappa_dT'] * np.sqrt(2.0), V2, dV2)):
                er, ec = np.exp(-kk * r), np.exp(-kk * rc)
                Y = np.where(inside, er / r - ec / rc, 0.0)
                V[a:a + step] = (qq * eps * Y).sum(1)
                dV[a:a + step] = (qq * (deps * Y + eps * dkk * np.where(inside, ec - er, 0.0))).sum(1)
        vz = z['V_DH'].astype(float)
        if np.abs(V1 - vz).max() > 1e-3 * max(1.0, np.abs(vz).max()):
            raise ValueError('V_DH differs from the npz')
        G = z['G_T'].astype(float)
        dH = V2 - V1
        dG = -dH / (fdt.R_GAS * HC.T0 ** 2) + (dV2 - dV1) / HC.RT
        return salt_row(uid, rg, G + (CF.term('C', uid, rg, d) if HC.FORCING == 'capped' else D / HC.RT), G, dH, dG, HC.RT, t0), None, None
    except Exception as ex:  # noqa: BLE001
        return {'UniProt_ID': uid}, None, f'{uid}: {ex!r} | {traceback.format_exc()[-300:]}'.replace('\n', ' ')


def salt_M(uid):
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
        sysd = mh.read_mpipi_system(os.path.join(bs.MP_RUNS[300], uid))
        pp = mh.mpipi_pairs(sysd)
        i, j = pp['i'], pp['j']
        k, pref = sysd['kappa'], mh.QQR2E / sysd['diel']
        k2, rc_c = k * np.sqrt(2.0), pp['rc_c']
        nf = len(xyz)
        E1, E2 = np.zeros(nf), np.zeros(nf)
        step = max(1, SP.MAX_ELEMS // max(1, len(i)))
        for a in range(0, nf, step):
            r = 10.0 * np.linalg.norm(xyz[a:a + step, i] - xyz[a:a + step, j], axis=2)
            E1[a:a + step] = np.where(r < rc_c, pref * pp['qq'] * np.exp(-k * r) / r, 0.0).sum(1)
            E2[a:a + step] = np.where(r < rc_c, pref * pp['qq'] * np.exp(-k2 * r) / r, 0.0).sum(1)
        ez = z['E_DH'].astype(float)
        if np.abs(ez).max() > 1e-6 and (np.corrcoef(E1, ez)[0, 1] < 0.999 or np.abs(E1 - ez).max() > 1e-4 * np.abs(ez).max()):
            raise ValueError('recomputed E_DH does not match the npz')
        dH = E2 - E1
        return salt_row(uid, rg, Gtot, G, dH, -dH / (HG.R_KCAL * HG.T0 ** 2), HG.RT, t0), None, None
    except Exception as ex:  # noqa: BLE001
        return {'UniProt_ID': uid}, None, f'{uid}: {ex!r} | {traceback.format_exc()[-300:]}'.replace('\n', ' ')


def salt_R(args):
    uid, d = args
    t0 = time.time()
    try:
        import mpipi_recharged_fdt as R
        import recharged_traj_pass as RTP
        MS, TC = RTP.MS, RTP.TC
        z = np.load(os.path.join(TC.RC_NPZ, f'{uid}.npz'))
        xyz, seq, info = R.load_run(d, uid)
        if not np.array_equal(info['production_start'] + info['stride'] * np.arange(len(xyz)), z['frame_index']):
            raise ValueError('frames differ from the stored FDT frames')
        tb = MS.tables()
        T0 = RTP.T0
        k = MS.kappa_of(T0)
        dk = MS.kappa_of(T0 + 0.5) - MS.kappa_of(T0 - 0.5)
        k2, dk2 = k * np.sqrt(2.0), dk * np.sqrt(2.0)
        rg, G = z['rg'].astype(float), z['G_T'].astype(float)
        N, nf, L = len(seq), len(xyz), R.box_nm(len(seq))
        ix = np.array([tb['ix'][a] for a in seq])
        ft = {a: MS.G_TEA * RTP.TF.DL['charmm36'][a] + MS.G_HYD * RTP.TF.DL_CLASSES['hyd'][a] for a in MS.AA}
        fv = np.array([ft[a] for a in seq])
        i, j = np.triu_indices(N, k=2)
        W, Y1, dY1, D, Y2, dY2 = (np.zeros(nf) for _ in range(6))
        step = max(1, int(TC.MAX_ELEMS // max(1, len(i))))
        for a in range(0, nf, step):
            dv = xyz[a:a + step, i].astype(float) - xyz[a:a + step, j].astype(float)
            dv -= L * np.round(dv / L)
            r = np.linalg.norm(dv, axis=2)
            w_, y_, dy_, d_ = RTP.per_pair(r, ix[i], ix[j], tb, k, dk, fv[i], fv[j])
            _, y2_, dy2_, _ = RTP.per_pair(r, ix[i], ix[j], tb, k2, dk2, fv[i], fv[j])
            W[a:a + step], Y1[a:a + step], dY1[a:a + step], D[a:a + step] = w_.sum(1), y_.sum(1), dy_.sum(1), d_.sum(1)
            Y2[a:a + step], dY2[a:a + step] = y2_.sum(1), dy2_.sum(1)
        if np.abs(W - z['V_AH']).max() > 1e-3 * max(1.0, np.abs(z['V_AH']).max()):
            raise ValueError('WF energy differs from the npz')
        dH = Y2 - Y1
        dG = -dH / (RTP.R_KJ * T0 ** 2) + (dY2 - dY1) / RTP.RT
        return salt_row(uid, rg, G + (CF.term('R', uid, rg, d) if HC.FORCING == 'capped' else D / RTP.RT), G, dH, dG, RTP.RT, t0), None, None
    except Exception as ex:  # noqa: BLE001
        return {'UniProt_ID': uid}, None, f'{uid}: {ex!r} | {traceback.format_exc()[-300:]}'.replace('\n', ' ')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--what', choices=['phospho', 'salt'], required=True)
    ap.add_argument('--model', choices=['C', 'M', 'R'], required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--workers', type=int, default=6)
    ap.add_argument('--limit', type=int, default=None)
    ap.add_argument('--forcing', choices=['A0', 'capped'], default='A0')
    a = ap.parse_args()
    os.environ['BENDER_SCAN_FORCING'] = a.forcing           # capped baseline in the spawned workers (heat_crosses)
    if a.what == 'phospho' and a.model == 'R':
        raise SystemExit('phosphorylation in Mpipi-Recharged is a Glu/Asp substitution; not part of this check')
    os.makedirs(a.out, exist_ok=True)
    src = AH.IDS[a.model]
    ids = sorted(pd.read_csv(src).UniProt_ID) if src.endswith('.csv') else sorted({l.strip() for l in open(src) if l.strip()})
    pc = os.path.join(a.out, f'{a.what}_{a.model}_chains.csv')
    done = set(pd.read_csv(pc, usecols=['UniProt_ID']).UniProt_ID) if os.path.exists(pc) else set()
    fn = {('phospho', 'C'): phospho_C, ('phospho', 'M'): phospho_M,
          ('salt', 'C'): salt_C, ('salt', 'M'): salt_M, ('salt', 'R'): salt_R}[(a.what, a.model)]
    if a.model == 'C':
        todo = list(HC.TF.index_dirs([u for u in ids if u not in done]).items())
    elif a.model == 'R':
        import tea_crossff as TC
        dirs = TC.rc_dirs()
        todo = [(u, dirs[u]) for u in ids if u not in done and u in dirs]
    else:
        todo = [u for u in ids if u not in done]
    todo = todo[:a.limit]
    print(f'{a.what} {a.model}: {len(ids)} chains, {len(done)} done, {len(todo)} to run', flush=True)
    nfail = 0
    with mp.get_context('spawn').Pool(a.workers) as pool:
        for k, (row, S, err) in enumerate(pool.imap_unordered(fn, todo), 1):
            if err:
                nfail += 1
                print('failed', err[:300], flush=True)
                continue
            pd.DataFrame([row]).to_csv(pc, mode='a', header=not os.path.exists(pc), index=False)
            if S is not None and len(S):
                ps = os.path.join(a.out, f'{a.what}_{a.model}_sites.csv')
                S.to_csv(ps, mode='a', header=not os.path.exists(ps), index=False)
            if k % 50 == 0:
                print(f'{k}/{len(todo)} done, {nfail} failed', flush=True)
    print('finished', flush=True)


if __name__ == '__main__':
    main()
