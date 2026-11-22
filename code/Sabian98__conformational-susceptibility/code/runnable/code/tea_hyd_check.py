"""
tea_hyd_check.py — second-route check of the hydrophilic-solvation term (tea_fdt.py / tea_crossff.py --classes, 'hyd').

  calvados   for a random sample with stored FDT frames at 278-320 K: reweight each temperature with the hydrophilic
             energy (tea_calibrate.calvados_hyd_energy) at gamma_hyd, fit <Rg>(T) with a quadratic and compare the change
             of the 300 K slope with gamma_hyd * dchi_hyd from tea_dchi_classes_300K.csv (independent route: finite
             temperatures and the full eq 2 incl. dc_p, vs the 300 K linear-response term); Kish ESS per temperature
  recharged  Mpipi-Recharged has 300 K runs only: finite-difference derivative of the hydrophilic energy at
             300 +/- 0.5 K on the stored 300 K frames (tea_calibrate.recharged_class_energies), dchi = -Cov(Rg, dE/dT)/RT,
             vs the stored dchi_hyd (checks the implementation and the attractive-only mapping, not linear response)
    python tea_hyd_check.py calvados --out DIR [--n 300] [--gamma 0.1 0.358]
  gg         (added 2026-10-10) Mpipi-GG has stored 278-320 K runs: as 'calvados', reweighting each temperature with the
             hydrophilic energy E = sum_ij (w_i + w_j)/2 att_ij(r) (kcal/mol; w = du_h(T) - du_h(T0) for DEKRNQSTHG,
             tea_calibrate.hyd_weights; att = tea_crossff.wf_attractive with the run's LAMMPS pair parameters) at the
             Mpipi-GG strength, vs gamma * dchi_hyd from tea_dchi_gg_classes_300K.csv
    python tea_hyd_check.py recharged --out DIR [--n 100]
    python tea_hyd_check.py gg --out DIR [--n 100] [--gamma 0.1059]
"""
import sys

sys.dont_write_bytecode = True

import argparse  # noqa: E402
import multiprocessing as mp  # noqa: E402
import os  # noqa: E402
import traceback  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

HERE = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, os.environ.get('BENDER_REANALYSIS_DIR', HERE))
sys.path.insert(0, HERE)
import tea_calibrate as K  # noqa: E402

W, L, fdt, bc = K.W, K.L, K.fdt, K.bc
TEA_DIR = os.path.join(W.FDT_DIR, 'tea_300K')
RC_NPZ = os.path.join(L.BENDER_MT, 'Traj Level Analysis', 'Analysis_FDT_RECHARGED', 'npz_300K')
RC_RUNS = (__import__('os').environ['CS_DATA'] + '/MPIPI Recharged/300K')


def _cal_one(args):
    uid, kingdom, gammas = args
    try:
        dirs = W._dirs(uid, kingdom)
        rows = []
        for T, d in dirs.items():
            z = np.load(os.path.join(W.FDT_DIR, 'npz', f'{T}K', f'{uid}.npz'))
            xyz, seq, info = bc.load_calvados(os.path.join(d, f'{uid}.dcd'), os.path.join(d, 'top.pdb'), k=L.K_FRAMES)
            if not np.array_equal(info['production_start'] + info['stride'] * np.arange(len(xyz)), z['frame_index']):
                raise ValueError('frames differ')
            E = {'hyd': K.calvados_hyd_energy(xyz, seq, float(T))}
            rg = z['rg'].astype(float)
            for g in [0.0] + list(gammas):
                m, ess = W.reweight(rg, E, {'hyd': g}, float(T))
                rows.append({'T': T, 'gamma': g, 'rg': m, 'ess': ess})
        d = pd.DataFrame(rows)
        out = {'UniProt_ID': uid}
        for g, gd in d.groupby('gamma'):
            gd = gd.sort_values('T')
            p = np.polyfit(gd['T'] - W.T0, gd.rg, 2)
            out[f'slope_g{g:g}'] = 100 * p[1] / p[2]
            out[f'ess_min_g{g:g}'] = gd.ess.min()
        return out, None
    except Exception as ex:
        return {'UniProt_ID': uid}, f'{ex!r} | {traceback.format_exc()[-250:]}'.replace('\n', ' ')


GG_SIM = (__import__('os').environ['CS_DATA'] + '/BENDER-BIO/MPIPI/mpipi_sim')
GG_NPZ = os.path.join(L.BENDER_MT, 'Traj Level Analysis', 'Analysis_FDT_MPIPI')
R_KCAL = 1.98720425864083e-3


def _gg_dir(uid, T):
    if T == 300:
        for sub in ('Mpipi 300k ', 'Mpipi SAXS'):
            d = os.path.join(GG_SIM, '300k_runs', sub, uid)
            if os.path.isfile(os.path.join(d, 'traj.xtc')):
                return d
        return None
    d = os.path.join(GG_SIM, f'{T}k_runs', uid)
    return d if os.path.isfile(os.path.join(d, 'traj.xtc')) else None


def _gg_one(args):
    uid, gammas = args
    try:
        import mdtraj as md
        import mpipi_temperature_hpc as MH
        import tea_crossff as TC
        rows = []
        for T in (278, 288, 300, 310, 320):
            d = _gg_dir(uid, T)
            z = np.load(os.path.join(GG_NPZ, f'npz_{T}K', f'{uid}.npz'))
            seq = bc.read_mpipi_sequence(f'{d}/seq.dat')
            with md.formats.XTCTrajectoryFile(f'{d}/traj.xtc') as fh:
                xyz, _, step, box = fh.read()
            pos = {int(s_): k for k, s_ in enumerate(step)}
            x = xyz[np.array([pos[int(s_)] for s_ in z['step']])].astype(float) * 10.0      # A, as LAMMPS units
            Lb = float(np.median(box[:, 0, 0])) * 10.0
            pp = MH.mpipi_pairs(MH.read_mpipi_system(d))
            nu, mu, sig, rc = pp['nu'], pp['mu'], pp['sig'], pp['rc']
            alpha = 2 * nu * (rc / sig) ** (2 * mu) * ((1 + 2 * nu) / (2 * nu * ((rc / sig) ** (2 * mu) - 1))) ** (2 * nu + 1)
            w = K.hyd_weights(seq, float(T))
            wij = 0.5 * (w[pp['i']] + w[pp['j']])
            E = np.zeros(len(x))
            stepn = max(1, int(TC.MAX_ELEMS // max(1, len(pp['i']))))
            for a in range(0, len(x), stepn):
                dv = x[a:a + stepn, pp['i']] - x[a:a + stepn, pp['j']]
                dv -= Lb * np.round(dv / Lb)
                r = np.linalg.norm(dv, axis=2)
                E[a:a + stepn] = (TC.wf_attractive(r, sig, mu, nu, rc, alpha) * wij).sum(1)
            rg = z['rg'].astype(float)
            bv = np.diff(x, axis=1)
            bv -= Lb * np.round(bv / Lb)                              # LAMMPS wraps: unwrap along the chain for Rg
            xu = np.concatenate([x[:, :1], x[:, :1] + np.cumsum(bv, axis=1)], axis=1)
            rg_x = np.sqrt(((xu - xu.mean(1, keepdims=True)) ** 2).sum(2).mean(1)) / 10.0
            if np.abs(rg_x - rg).max() > 1e-3 * rg.mean():
                raise ValueError(f'{T} K: frames differ from the npz')
            for g in [0.0] + list(gammas):
                u = -g * E / (R_KCAL * T)
                wt = np.exp(u - u.max())
                rows.append({'T': T, 'gamma': g, 'rg': (wt * rg).sum() / wt.sum(), 'ess': wt.sum() ** 2 / (wt ** 2).sum() / len(wt)})
        dd = pd.DataFrame(rows)
        out = {'UniProt_ID': uid}
        for g, gd in dd.groupby('gamma'):
            gd = gd.sort_values('T')
            p_ = np.polyfit(gd['T'] - 300.0, gd.rg, 2)
            out[f'slope_g{g:g}'] = 100 * p_[1] / p_[2]
            out[f'ess_min_g{g:g}'] = gd.ess.min()
        return out, None
    except Exception as ex:
        return {'UniProt_ID': uid}, f'{ex!r} | {traceback.format_exc()[-250:]}'.replace('\n', ' ')


def check_gg(out, n=100, gammas=(0.1059,), seed=5, workers=8):
    os.makedirs(out, exist_ok=True)
    cls = pd.read_csv(os.path.join(TEA_DIR, 'tea_dchi_gg_classes_300K.csv')).set_index('UniProt_ID')
    have = set(cls.index)
    for T in (278, 288, 310, 320):
        have &= {f[:-4] for f in os.listdir(os.path.join(GG_NPZ, f'npz_{T}K')) if f.endswith('.npz') and not f.startswith('._')}
    have = sorted(u for u in have if all(_gg_dir(u, T) for T in (278, 288, 300, 310, 320)))
    pick = list(np.random.default_rng(seed).choice(have, min(n, len(have)), replace=False))
    rows = []
    with mp.get_context('spawn').Pool(workers) as pool:
        for row, err in pool.imap_unordered(_gg_one, [(u, gammas) for u in pick]):
            if err:
                print('failed', row['UniProt_ID'], err[:200], flush=True)
            else:
                rows.append(row)
    d = pd.DataFrame(rows).set_index('UniProt_ID')
    d['dchi_hyd'] = cls.dchi_rg_rel_hyd
    d.to_csv(os.path.join(out, 'hyd_check_gg.csv'))
    lines = [f'Mpipi-GG, {len(d)} proteins (of {len(have)} with all five temperatures): change of the 300 K slope (quadratic '
             'over 278-320 K, % per K) from reweighting with the hydrophilic energy vs gamma * dchi_hyd (FDT, 300 K)']
    for g in gammas:
        obs, pred = d[f'slope_g{g:g}'] - d['slope_g0'], g * d.dchi_hyd
        lines.append(f'gamma {g}: r = {np.corrcoef(pred, obs)[0, 1]:.3f}, regression slope {np.polyfit(pred, obs, 1)[0]:.3f}, '
                     f'median obs/pred {np.median(obs / pred):.3f}, median pred {pred.median():+.3f}, '
                     f'min ESS median {d[f"ess_min_g{g:g}"].median():.3f} (5th pct {d[f"ess_min_g{g:g}"].quantile(0.05):.3f})')
    open(os.path.join(out, 'hyd_check_gg.txt'), 'w').write('\n'.join(lines) + '\n')
    print('\n'.join(lines))


def check_calvados(out, n=300, gammas=(0.1, 0.358), seed=2, workers=8):
    os.makedirs(out, exist_ok=True)
    cls = pd.read_csv(os.path.join(TEA_DIR, 'tea_dchi_classes_300K.csv')).set_index('UniProt_ID')
    bio = pd.read_csv(L.BIO_CSV, usecols=['UniProt_ID', 'kingdom']).set_index('UniProt_ID')
    have = set(cls.index)
    for T in W.TEMPS:
        have &= {f[:-4] for f in os.listdir(os.path.join(W.FDT_DIR, 'npz', f'{T}K')) if f.endswith('.npz') and not f.startswith('._')}
    pick = list(np.random.default_rng(seed).choice(sorted(have), n, replace=False))
    rows = []
    with mp.get_context('spawn').Pool(workers) as pool:
        for row, err in pool.imap_unordered(_cal_one, [(u, bio.kingdom.get(u), gammas) for u in pick]):
            if err:
                print('failed', row['UniProt_ID'], err[:160], flush=True)
            else:
                rows.append(row)
    d = pd.DataFrame(rows).set_index('UniProt_ID')
    d['dchi_hyd'] = cls.dchi_rg_rel_hyd
    d.to_csv(os.path.join(out, 'hyd_check_calvados.csv'))
    lines = [f'CALVADOS-2, {len(d)} proteins: change of the 300 K slope (quadratic over 278-320 K, % per K) from reweighting '
             'with the hydrophilic energy vs gamma_hyd * dchi_hyd (FDT, 300 K)']
    for g in gammas:
        obs, pred = d[f'slope_g{g:g}'] - d['slope_g0'], g * d.dchi_hyd
        lines.append(f'gamma_hyd {g}: r = {np.corrcoef(pred, obs)[0, 1]:.3f}, regression slope {np.polyfit(pred, obs, 1)[0]:.3f}, '
                     f'median obs/pred {np.median(obs / pred):.3f}, median pred {pred.median():+.3f}, '
                     f'min ESS median {d[f"ess_min_g{g:g}"].median():.3f} (5th pct {d[f"ess_min_g{g:g}"].quantile(0.05):.3f})')
    open(os.path.join(out, 'hyd_check_calvados.txt'), 'w').write('\n'.join(lines) + '\n')
    print('\n'.join(lines))


def _rc_one(uid):
    try:
        import mpipi_recharged_fdt as R
        z = np.load(os.path.join(RC_NPZ, f'{uid}.npz'))
        xyz, seq, info = R.load_run(os.path.join(RC_RUNS, uid), uid)
        if not np.array_equal(info['production_start'] + info['stride'] * np.arange(len(xyz)), z['frame_index']):
            raise ValueError('frames differ')
        box = R.box_nm(len(seq))
        Ep, Em = (K.recharged_class_energies(xyz, seq, T, box)['hyd'] for T in (300.5, 299.5))
        rg = z['rg'].astype(float)
        dchi = -np.cov(rg, Ep - Em, ddof=1)[0, 1] / (fdt.R_GAS * 300.0)
        return {'UniProt_ID': uid, 'dchi_fd': 100 * dchi / rg.mean()}, None
    except Exception as ex:
        return {'UniProt_ID': uid}, f'{ex!r} | {traceback.format_exc()[-250:]}'.replace('\n', ' ')


def check_recharged(out, n=100, seed=3, workers=8):
    os.makedirs(out, exist_ok=True)
    cls = pd.read_csv(os.path.join(TEA_DIR, 'tea_dchi_recharged_classes_300K.csv')).set_index('UniProt_ID')
    pick = list(np.random.default_rng(seed).choice(sorted(cls.index), n, replace=False))
    rows = []
    with mp.get_context('spawn').Pool(workers) as pool:
        for row, err in pool.imap_unordered(_rc_one, pick):
            if err:
                print('failed', row['UniProt_ID'], err[:160], flush=True)
            else:
                rows.append(row)
    d = pd.DataFrame(rows).set_index('UniProt_ID')
    d['dchi_stored'] = cls.dchi_rg_rel_hyd
    d.to_csv(os.path.join(out, 'hyd_check_recharged.csv'))
    rel = (d.dchi_fd - d.dchi_stored).abs() / d.dchi_stored.abs()
    txt = (f'Mpipi-Recharged, {len(d)} proteins: finite-difference hydrophilic term (300 +/- 0.5 K, calibration code) vs '
           f'stored (tea_crossff.py): r = {np.corrcoef(d.dchi_fd, d.dchi_stored)[0, 1]:.6f}, median relative difference '
           f'{rel.median():.2e}, max {rel.max():.2e}')
    open(os.path.join(out, 'hyd_check_recharged.txt'), 'w').write(txt + '\n')
    print(txt)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('model', choices=['calvados', 'recharged', 'gg'])
    ap.add_argument('--out', required=True)
    ap.add_argument('--n', type=int, default=None)
    ap.add_argument('--gamma', type=float, nargs='+', default=[0.1, 0.358])
    a = ap.parse_args()
    if a.model == 'calvados':
        check_calvados(a.out, a.n or 300, a.gamma)
    elif a.model == 'gg':
        check_gg(a.out, a.n or 100, a.gamma if a.gamma != [0.1, 0.358] else [0.1059])
    else:
        check_recharged(a.out, a.n or 100)
