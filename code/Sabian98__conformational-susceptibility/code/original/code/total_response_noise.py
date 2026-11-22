"""
total_response_noise.py — sampling noise of the calibrated total heating response, for noise ceilings.

The total response at 300 K is chi_total + gamma_tea * dchi_tea + gamma_hyd * dchi_hyd = -Cov(Rg, G_tot) with the
per-frame force
    G_tot = G + (gamma_tea * D_tea + gamma_hyd * D_hyd) / (R T0)
(G = d(beta H)/dT of the plain model: CALVADOS-2 and Mpipi-Recharged npz 'G_T', Mpipi-GG -PotEng/(R T^2);
D_hyd, D_tea = the hydrophilic and TEA (CHARMM36m Arg) forcings per frame from tea_fdt.py / tea_crossff.py; strengths
as in total_response_evidence.py: the Wuttke calibration with TEA at its prior; Mpipi-GG at the others' relative size).
Standard error by the same delete-one-block jackknife (10 contiguous blocks) as chi_total_se in the FDT runs, for chi,
the hydrophilic term and the total. (No check is made in this script; an external check on 400 chains per model
reproduced the stored chi to 1.5e-5 relative and the hydrophilic and TEA terms exactly, 2026-10-10.)
A random subset per model (the stored proteome files hold no SE for the hydrophilic term).
    python total_response_noise.py --model C|M|R --n 400 --out DIR [--workers N]
Output: DIR/total_response_se_<model>.csv (relative values: % of <Rg> per K).
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
import tea_crossff as TC  # noqa: E402
import tea_fdt as TF  # noqa: E402

GAMMA = {'C': 0.358, 'R': 0.209, 'M': 0.1059}   # Wuttke calibration (C, R); Mpipi-GG at the others' relative size
GAMMA_TEA = {'C': 2.0, 'R': 0.4, 'M': 0.3000}   # the calibration's TEA prior; Mpipi-GG likewise by relative size
TEA_KEYS = ['aliph', 'pro', 'arom', 'argC36', 'rest']   # tea_fdt CLASSES: they add up to the CHARMM36m-Arg TEA forcing
N_BLOCKS = 10
T0 = TF.T0


def jk_se(f, n):
    blocks = np.array_split(np.arange(n), N_BLOCKS)
    jk = np.array([f(np.setdiff1d(np.arange(n), b)) for b in blocks])
    return float(np.sqrt((N_BLOCKS - 1) / N_BLOCKS * ((jk - jk.mean()) ** 2).sum()))


def summarise(uid, model, rg, G, D, Dtea, RT):
    g = GAMMA[model]
    Gt = G + (g * D + GAMMA_TEA[model] * Dtea) / RT

    def cv(x):
        return lambda idx: -float(np.cov(rg[idx], x[idx], ddof=1)[0, 1])

    n, m = len(rg), rg.mean()
    return {'UniProt_ID': uid, 'model': model, 'gamma_hyd': g, 'n_frames': n,
            'chi_rel': 100 * cv(G)(np.arange(n)) / m, 'chi_rel_se': 100 * jk_se(cv(G), n) / m,
            'hyd_rel': 100 * cv(D / RT)(np.arange(n)) / m, 'hyd_rel_se': 100 * jk_se(cv(D / RT), n) / m,
            'tea_rel': 100 * cv(Dtea / RT)(np.arange(n)) / m,
            'total_rel': 100 * cv(Gt)(np.arange(n)) / m, 'total_rel_se': 100 * jk_se(cv(Gt), n) / m}


# ---- per-frame D_hyd: CALVADOS-2 recomputed as in tea_fdt.one; Mpipi via tea_crossff with its finish() intercepted
def calvados(uid, d):
    z = np.load(os.path.join(TF.NPZ, f'{uid}.npz'))
    xyz, seq, info = TF.bc.load_calvados(os.path.join(d, f'{uid}.dcd'), os.path.join(d, 'top.pdb'), k=TF.L.K_FRAMES)
    fi = info['production_start'] + info['stride'] * np.arange(len(xyz))
    if not np.array_equal(fi, z['frame_index']):
        raise ValueError('frames differ from the stored FDT frames')
    q, lam, sig, _ = TF.fdt.calvados_residue_params(seq, TF.L.load_residue_params())
    i, j, s, l, shift, _ = TF.fdt.pair_setup(q, lam, sig)
    dls = {}
    for key, v in (('hyd', TF.DL_CLASSES['hyd']), ('tea', TF.DL['charmm36'])):
        dls[key] = 0.5 * (np.array([v[a] for a in seq])[i] + np.array([v[a] for a in seq])[j])
    D = {key: np.zeros(len(xyz)) for key in dls}
    step = max(1, int(TF.fdt.MAX_ELEMS // max(1, len(i))))
    for a in range(0, len(xyz), step):
        r = np.linalg.norm(xyz[a:a + step, i].astype(float) - xyz[a:a + step, j].astype(float), axis=2)
        u1, u0 = TF.fdt.ah_pair_energy(r, s, np.ones_like(l), shift), TF.fdt.ah_pair_energy(r, s, np.zeros_like(l), shift)
        for key in D:
            D[key][a:a + step] = ((u1 - u0) * dls[key]).sum(axis=1)
    return z['rg'].astype(float), z['G_T'].astype(float), D['hyd'], D['tea'], TF.fdt.R_GAS * T0


_CAUGHT = {}


def _catch(uid, seq, A_obs, D, check, t0, energy_unit_R):
    _CAUGHT[uid] = (np.asarray(A_obs['rg'], dtype=float), D['hyd'], sum(D[k] for k in TEA_KEYS), energy_unit_R * T0)
    return {'UniProt_ID': uid}


TC.finish = _catch   # applied in every worker too (module-level)


def mpipi(uid, d, model):
    row, err = (TC.gg_one if model == 'M' else TC.rc_one)((uid, d, 'classes'))
    if err:
        raise RuntimeError(err)
    rg, D, Dtea, RT = _CAUGHT.pop(uid)
    if model == 'M':
        z = np.load(os.path.join(TC.GG_NPZ, f'{uid}.npz'))
        G = -z['PotEng_kcal'].astype(float) / (TC.R_KCAL * T0 ** 2)
    else:
        z = np.load(os.path.join(TC.RC_NPZ, f'{uid}.npz'))
        G = z['G_T'].astype(float)
    return rg, G, D, Dtea, RT


def one(task):
    uid, d, model = task
    try:
        rg, G, D, Dtea, RT = calvados(uid, d) if model == 'C' else mpipi(uid, d, model)
        return summarise(uid, model, rg, G, D, Dtea, RT), None
    except Exception as ex:
        return None, f'{uid}: {ex!r} | {traceback.format_exc()[-300:]}'.replace('\n', ' ')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--model', choices=['C', 'M', 'R'], required=True)
    ap.add_argument('--n', type=int, default=400)
    ap.add_argument('--out', required=True)
    ap.add_argument('--workers', type=int, default=6)
    ap.add_argument('--seed', type=int, default=7)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    tea = os.path.join(TC.TLA, 'Analysis_FDT', 'tea_300K')
    f = {'C': 'tea_dchi_classes_300K.csv', 'M': 'tea_dchi_gg_classes_300K.csv', 'R': 'tea_dchi_recharged_classes_300K.csv'}[a.model]
    ids = sorted(pd.read_csv(os.path.join(tea, f), usecols=['UniProt_ID']).UniProt_ID)
    if a.model == 'R':
        ids = [u for u in ids if u != 'Q983T3']
    pick = list(np.random.default_rng(a.seed).choice(ids, size=min(a.n, len(ids)), replace=False))
    dirs = TF.index_dirs(pick) if a.model == 'C' else (TC.gg_dirs() if a.model == 'M' else TC.rc_dirs())
    tasks = [(u, dirs[u], a.model) for u in pick if u in dirs]
    print(f'{a.model}: {len(tasks)} proteins', flush=True)
    rows, fails = [], []
    with mp.get_context('spawn').Pool(a.workers) as pool:
        for k, (r, err) in enumerate(pool.imap_unordered(one, tasks), 1):
            (fails.append(err) if err else rows.append(r))
            if k % 100 == 0:
                print(f'  {k}/{len(tasks)}', flush=True)
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(a.out, f'total_response_se_{a.model}.csv'), index=False)
    print(f'finished {a.model}: {len(rows)} ok, {len(fails)} failed', flush=True)
    for e in fails[:3]:
        print('  ', e[:300])


if __name__ == '__main__':
    main()
