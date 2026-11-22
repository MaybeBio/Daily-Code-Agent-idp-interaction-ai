"""
heat_crosses.py — one trajectory pass per chain (CALVADOS-2, 300 K) for three questions about the calibrated total
heating response, total = -Cov(Rg, G_tot) / <Rg> (% per K), G_tot = G_T + (gamma_tea D_tea + gamma_hyd D_hyd) / RT
(G_T: the stored FDT thermal force; D: the TEA and hydrophilic forcings as tea_fdt.py; gamma_tea 2.0, gamma_hyd 0.358).

A  Dynamics at full resolution: every production frame (DCD spacing), integrated autocorrelation time of Rg (ns of
   simulation time; Langevin model time, compare chains only).
B  Thermal hotspots: the total split exactly by residue. Each pair's energy and forcing go half to each partner:
     solvation share  -Cov(Rg, D_r) / RT,  Boltzmann share  Cov(Rg, H_r) / (R T^2),  dielectric share  -Cov(Rg, dH_r/dT) / RT
   with H_r the residue's nonbonded (Ashbaugh-Hatch + Debye-Hueckel) energy; the shares add to the total minus the bond
   term (checked against the stored chain values).
C  Heat x phosphorylation and heat x salt, by exact reweighting of the same frames:
     phosphate at Ser/Thr i (charge -2, as salt_ph_phospho.py): dH = dq phi_i, dG = -dH/(R T^2) + dq (dphi_i/dT)/(R T)
     doubling the ionic strength: dH = V_DH(sqrt(2) kappa) - V_DH(kappa), dG likewise from dV/dT at both screenings
   new total = -Cov_w(Rg, G_tot + dG) / <Rg>_w with w = exp(-dH/RT); reported as the change in % per K, with a
   10-block jackknife SE and the Kish effective sample fraction of w.
Sample: proteins with annotated phosphosites (phosphosite_test/sites.csv), the Arg shortlist (arg_obs/A_arg_shortlist.csv)
and 300 random chains.
    python heat_crosses.py --out DIR [--workers N] [--limit N] [--forcing A0|capped] [--sample earlier/sample.csv]
Output: DIR/chains.csv, residues.csv, sites.csv (append; resumable)
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
import tea_fdt as TF  # noqa: E402

L, fdt, bc = TF.L, TF.fdt, TF.bc
T0 = 300.0
RT = fdt.R_GAS * T0
CAL = fdt.get_dielectric_params(T0)
G_TEA, G_HYD = 2.0, 0.358
# --forcing capped (2026-10-25): the main model's forcing (capped_proteome.py: charged-group and TEA terms, each residue's
# contact number capped by h(c) = c_sat (1 - exp(-c / c_sat)), strengths from capped_model.py) in place of the linear A0
# forcing. Passed to spawned workers through BENDER_SCAN_FORCING. The capped energy is a sum over residues, so the
# residue shares of part B take each residue's own term instead of half of each pair's.
FORCING = os.environ.get('BENDER_SCAN_FORCING', 'A0')
if FORCING == 'capped':
    import capped_model as CM  # noqa: E402
    import capped_proteome as CP  # noqa: E402
    _gh, _gt = CM.params()
    G_TEA, G_HYD = _gt['C'], _gh['C'][0]
DQ = -2.0
N_BLOCKS = 10
TLA = os.path.join(L.BENDER_MT, 'Traj Level Analysis')


def jk(f, n):
    blocks = np.array_split(np.arange(n), N_BLOCKS)
    v = np.array([f(np.setdiff1d(np.arange(n), b)) for b in blocks])
    return float(np.sqrt((N_BLOCKS - 1) / N_BLOCKS * ((v - v.mean()) ** 2).sum()))


def wcov(a, b, w):
    w = w / w.sum()
    return float((w * (a - (w * a).sum()) * (b - (w * b).sum())).sum())


def tau_int_series(x, c=5.0):
    x = np.asarray(x, float) - np.mean(x)
    n = len(x)
    f = np.fft.rfft(x, 2 * n)
    acf = np.fft.irfft(f * np.conj(f))[:n]
    acf /= acf[0]
    tau = 0.5
    for w in range(1, n):
        tau += acf[w]
        if w >= c * tau:
            break
    return tau


def one(args):
    uid, d = args
    t0 = time.time()
    try:
        z = np.load(os.path.join(TF.NPZ, f'{uid}.npz'))
        xyz_all, seq, info = bc.load_calvados(os.path.join(d, f'{uid}.dcd'), os.path.join(d, 'top.pdb'), k=None)
        idx = bc.even_indices(len(xyz_all), L.K_FRAMES)
        fi = info['production_start'] + idx
        if not np.array_equal(fi, z['frame_index']):
            raise ValueError('frames differ from the stored FDT frames')
        rg_all = np.sqrt(((xyz_all - xyz_all.mean(1, keepdims=True)) ** 2).sum(2).mean(1))
        xyz = xyz_all[idx]
        rg = z['rg'].astype(float)
        if np.abs(rg_all[idx] - rg).max() > 1e-3 * rg.mean():
            raise ValueError('Rg differs from the npz')
        N, nf = len(seq), len(xyz)
        run_ns = 100.0 * max(1.0, (N / 150.0) ** 1.5)
        tau_ns = tau_int_series(rg_all) * run_ns / info['n_frames_total']
        q, lam, sig, _ = fdt.calvados_residue_params(seq, L.load_residue_params())
        i, j, s, l, shift, _ = fdt.pair_setup(q, lam, sig)
        npair = len(i)
        half = sparse.csr_matrix((np.full(2 * npair, 0.5), (np.concatenate([i, j]), np.concatenate([np.arange(npair)] * 2))), shape=(N, npair))
        qmat = sparse.csr_matrix((np.concatenate([q[j], q[i]]), (np.concatenate([i, j]), np.concatenate([np.arange(npair)] * 2))), shape=(N, npair))
        dl = {k: np.array([v[a] for a in seq]) for k, v in (('tea', TF.DL['charmm36']), ('hyd', TF.DL_CLASSES['hyd']))}
        fpair = G_TEA * 0.5 * (dl['tea'][i] + dl['tea'][j]) + G_HYD * 0.5 * (dl['hyd'][i] + dl['hyd'][j])
        k1, k2 = CAL['kappa'], CAL['kappa'] * np.sqrt(2.0)
        dk1, dk2 = CAL['dkappa_dT'], CAL['dkappa_dT'] * np.sqrt(2.0)
        eps, deps, rc = CAL['eps_yu'], CAL['deps_yu_dT'], fdt.RC_YU
        qq = q[i] * q[j]
        Dres, Hres, dHres = np.zeros((nf, N)), np.zeros((nf, N)), np.zeros((nf, N))
        Sres = np.zeros((nf, N)) if FORCING == 'capped' else None
        V1, V2, dV1, dV2 = np.zeros(nf), np.zeros(nf), np.zeros(nf), np.zeros(nf)
        phi, dphi = np.zeros((nf, N)), np.zeros((nf, N))
        step = max(1, int(fdt.MAX_ELEMS // max(1, npair)))
        for a in range(0, nf, step):
            r = np.linalg.norm(xyz[a:a + step, i] - xyz[a:a + step, j], axis=2)
            u1, u0 = fdt.ah_pair_energy(r, s, np.ones_like(l), shift), fdt.ah_pair_energy(r, s, np.zeros_like(l), shift)
            uah = u0 + l * (u1 - u0)
            inside = r < rc
            out = {}
            for tag, kk, dkk in (('1', k1, dk1), ('2', k2, dk2)):
                er, ec = np.exp(-kk * r), np.exp(-kk * rc)
                Y = np.where(inside, er / r - ec / rc, 0.0)
                dYdk = np.where(inside, ec - er, 0.0)
                out[tag] = (Y, deps * Y + eps * dkk * dYdk)
            Y1, dY1 = out['1']
            Y2, dY2 = out['2']
            V1[a:a + step], V2[a:a + step] = (qq * eps * Y1).sum(1), (qq * eps * Y2).sum(1)
            dV1[a:a + step], dV2[a:a + step] = (qq * dY1).sum(1), (qq * dY2).sum(1)
            Dres[a:a + step] = (half @ ((u1 - u0) * fpair).T).T
            if Sres is not None:
                Sres[a:a + step] = 2.0 * (half @ (u1 - u0).T).T
            Hres[a:a + step] = (half @ (uah + qq * eps * Y1).T).T
            dHres[a:a + step] = (half @ (qq * dY1).T).T
            phi[a:a + step] = (qmat @ (eps * Y1).T).T
            dphi[a:a + step] = (qmat @ dY1.T).T
        vz = z['V_DH'].astype(float)
        if np.abs(V1 - vz).max() > 1e-3 * max(1.0, np.abs(vz).max()):
            raise ValueError('V_DH differs from the npz')
        if Sres is not None:                       # capped: D_i = -eps/2 f_i h(c_i), c_i = -sum_j S_ij / eps
            fres = G_TEA * np.array([CP.F_TEA[a] for a in seq]) + G_HYD * np.array([CP.F_CHG[a] for a in seq])
            Dres = -0.5 * fdt.EPS_AH * CP.h(-Sres / fdt.EPS_AH) * fres
        G = z['G_T'].astype(float)
        Dtot = Dres.sum(1)
        Gtot = G + Dtot / RT
        m = rg.mean()
        base = -fdt._cov(rg, Gtot) / m * 100
        chi = -fdt._cov(rg, G) / m * 100
        solv = -fdt._cov(rg, Dtot / RT) / m * 100
        row = {'UniProt_ID': uid, 'N': N, 'tau_int_ns': tau_ns, 'tau_int_frames_full': tau_int_series(rg_all),
               'frames_full': len(rg_all), 'chi_rel': chi, 'solvation_rel': solv, 'total_rel': base,
               'dV_dT_check': float(np.corrcoef(dV1, z['dVdT_DH'].astype(float))[0, 1])}
        # B residues
        res = []
        rgc = rg - m
        for r_ in range(N):
            res.append({'UniProt_ID': uid, 'pos': r_ + 1, 'residue': seq[r_],
                        'solv_share': -float((rgc * (Dres[:, r_] - Dres[:, r_].mean())).sum() / (nf - 1)) / RT / m * 100,
                        'boltz_share': float((rgc * (Hres[:, r_] - Hres[:, r_].mean())).sum() / (nf - 1)) / (RT * T0) / m * 100,
                        'diel_share': -float((rgc * (dHres[:, r_] - dHres[:, r_].mean())).sum() / (nf - 1)) / RT / m * 100})
        R_ = pd.DataFrame(res)
        R_['total_share'] = R_[['solv_share', 'boltz_share', 'diel_share']].sum(1)
        row['residue_sum_minus_total'] = float(R_.total_share.sum() - (base - (-(fdt._cov(rg, -z['V_bond'].astype(float) / (RT * T0))) / m * 100)))
        row['solv_residue_sum_check'] = float(R_.solv_share.sum() - solv)

        # C perturbations
        def perturbed(dH, dG):
            x = -dH / RT
            w = np.exp(x - x.max())
            f = lambda ix: (-wcov(rg[ix], Gtot[ix] + dG[ix], w[ix]) / ((w[ix] * rg[ix]).sum() / w[ix].sum())  # noqa: E731
                            + fdt._cov(rg[ix], Gtot[ix]) * (len(ix) - 1) / len(ix) / rg[ix].mean()) * 100
            # fixed 2026-10-10: the baseline used ddof 1 against wcov's ddof 0, shifting every change by +|total|/n;
            # outputs written before the fix were corrected in place by heat_crosses_ddof_fix.py
            return f(np.arange(nf)), jk(f, nf), float(w.sum() ** 2 / (w ** 2).sum() / nf), \
                float(((w * rg).sum() / w.sum() / m - 1) * 100)
        dHs = V2 - V1
        dGs = -dHs / (fdt.R_GAS * T0 ** 2) + (dV2 - dV1) / RT
        row['salt_dtotal'], row['salt_dtotal_se'], row['salt_ess'], row['salt_dRg_pct'] = perturbed(dHs, dGs)
        sites = []
        for p_, c in enumerate(seq):
            if c in 'ST':
                dH = DQ * phi[:, p_]
                dG = -dH / (fdt.R_GAS * T0 ** 2) + DQ * dphi[:, p_] / RT
                dt, se, ess, drg = perturbed(dH, dG)
                sites.append({'UniProt_ID': uid, 'pos': p_ + 1, 'residue': c, 'phospho_dtotal': dt, 'phospho_dtotal_se': se,
                              'ess': ess, 'phospho_dRg_pct': drg})
        row['seconds'] = round(time.time() - t0, 1)
        return row, R_, pd.DataFrame(sites), None
    except Exception as ex:  # noqa: BLE001
        return {'UniProt_ID': uid}, None, None, f'{uid}: {ex!r} | {traceback.format_exc()[-300:]}'.replace('\n', ' ')


def sample(seed=11):
    sites = pd.read_csv(os.path.join(TLA, 'Analysis_FDT', 'phosphosite_test', 'sites.csv'))
    ph = set(sites.loc[(sites.model == 'CALVADOS-2') & (sites.annotated), 'UniProt_ID'])
    arg = set(pd.read_csv(os.path.join(TLA, 'Analysis_FDT', 'total_response', 'arg_obs', 'A_arg_shortlist.csv')).UniProt_ID)
    allids = sorted(f[:-4] for f in os.listdir(TF.NPZ) if f.endswith('.npz') and not f.startswith('._'))
    rnd = set(np.random.default_rng(seed).choice(allids, 300, replace=False))
    return sorted(ph | arg | rnd), {'phosphosite': ph, 'arg_shortlist': arg, 'random': rnd}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', required=True)
    ap.add_argument('--workers', type=int, default=6)
    ap.add_argument('--limit', type=int, default=None)
    ap.add_argument('--ids', default=None, help='file with one UniProt ID per line instead of the default sample')
    ap.add_argument('--forcing', choices=['A0', 'capped'], default='A0')
    ap.add_argument('--sample', default=None, help='sample.csv of an earlier run: same chains and groups')
    a = ap.parse_args()
    os.environ['BENDER_SCAN_FORCING'] = a.forcing           # read by the spawned workers on import
    os.makedirs(a.out, exist_ok=True)
    if a.ids:
        ids = sorted({l.strip() for l in open(a.ids) if l.strip()})
        groups = {'given': set(ids)}
    else:
        ids, groups = sample()
    if a.sample:                                            # the chains and groups of an earlier run, unchanged
        t = pd.read_csv(a.sample)
        ids = sorted(t.UniProt_ID)
        groups = {g: set(t.loc[t[g].astype(bool), 'UniProt_ID']) for g in t.columns if g != 'UniProt_ID'}
    pd.DataFrame([{'UniProt_ID': u, **{g: u in s for g, s in groups.items()}} for u in ids]).to_csv(os.path.join(a.out, 'sample.csv'), index=False)
    pc = os.path.join(a.out, 'chains.csv')
    done = set(pd.read_csv(pc, usecols=['UniProt_ID']).UniProt_ID) if os.path.exists(pc) else set()
    dirs = TF.index_dirs([u for u in ids if u not in done])
    todo = [(u, d) for u, d in dirs.items()][:a.limit]
    print(f'{len(ids)} sampled, {len(done)} done, {len(todo)} to run', flush=True)
    nf = 0
    with mp.get_context('spawn').Pool(a.workers) as pool:
        for k, (row, R_, S_, err) in enumerate(pool.imap_unordered(one, todo), 1):
            if err:
                nf += 1
                print('failed', err[:300], flush=True)
                continue
            for df, name in ((pd.DataFrame([row]), 'chains.csv'), (R_, 'residues.csv'), (S_, 'sites.csv')):
                p = os.path.join(a.out, name)
                if len(df):
                    df.to_csv(p, mode='a', header=not os.path.exists(p), index=False)
            if k % 25 == 0:
                print(f'{k}/{len(todo)} done, {nf} failed', flush=True)
    print('finished', flush=True)


if __name__ == '__main__':
    main()
