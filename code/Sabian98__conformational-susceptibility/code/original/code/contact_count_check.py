"""
contact_count_check.py — where does the sticker susceptibility sit in space? Close contacts vs the attractive tail
of the pair potential, for Mpipi-GG (Wang-Frenkel, rc = 3 sigma) and CALVADOS-2 (Ashbaugh-Hatch, rc = 2 nm),
at 300 K, by own-nu bin around and below the theta point.

Per protein, on 200 of the 1,000 FDT frames (every 5th; same frame selection as the FDT runners), every nonbonded
pair's sticker energy is recomputed from the model's own expression and split at r = 1.5 sigma_ij: first-shell
(close) pairs vs tail pairs (1.5 sigma_ij to the cutoff). Both wells sit at ~1.1 sigma, so the split treats the two
models alike (a fixed 0.8 nm split would put Mpipi-GG's larger beads' well on the boundary). At 1.5 sigma the pair
energy is 19% (WF) and ~32% (AH, lambda 0.5) of the well depth; beyond it lies 21% (WF) and ~39% (AH) of the
r^2-weighted attraction. Reported per protein:
  contacts_per_res     pairs with r < 1.5 sigma_ij and |i-j| >= 3, per residue (frame mean)
  E_per_res_kT         sticker energy per residue in kT; close_share = share from close pairs
  E_per_contact_kT     close-pair energy per close pair
  r_E_contacts         per-frame correlation of the sticker energy with the contact count (switch-like contacts
                       give ~ -1)
  chi_close_share      share of Cov(Rg, E_s) from close pairs (the rest from the tail); the covariance is what
                       sets the sticker susceptibility, chi_s = Cov(Rg, E_s) / (R T^2)
  var_close_share      share of Var(E_s) from close pairs, (Var close + Cov close,tail) / Var total
  r_rg_contacts        per-frame correlation of Rg with the contact count
Sample: up to 40 proteins per own-nu bin (random, seed 0) from chi_factor_split/per_protein.csv (fit R^2 >= 0.99).

Output: Traj Level Analysis/Analysis_FDT_MPIPI/compact/contact_check/ (per_protein.csv, bins.csv)
Usage: python3 contact_count_check.py [--workers 6]
"""
import sys

sys.dont_write_bytecode = True  # keep the NAS code folder free of __pycache__

import argparse  # noqa: E402
import multiprocessing as mp  # noqa: E402
import os  # noqa: E402
import traceback  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

HERE = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, HERE)
import TODO_T2_trajectory_loader as L  # noqa: E402
import mpipi_temperature_hpc as mh  # noqa: E402

fdt, bc = L.fdt, L.bc
MT = '/Volumes/User Homes/Nu Project/BENDER-MT/Traj Level Analysis'
MP = f'{MT}/Analysis_FDT_MPIPI'
SRC = f'{MP}/compact/chi_factor_split/per_protein.csv'
OUT_DIR = f'{MP}/compact/contact_check'
MPIPI_RUNS = {'main': '/Volumes/User Homes/Nu Project/BENDER-BIO/MPIPI/mpipi_sim/300k_runs/Mpipi 300k ',
              'compact': '/Volumes/User Homes/Nu Project/BENDER-MT/mpipi_compact'}
BINS = ['(0.3, 0.45]', '(0.45, 0.5]', '(0.5, 0.54]', '(0.54, 0.58]']
PER_BIN, SEED, T, STRIDE, CUT_SIGMA, MIN_SEP = 40, 0, 300.0, 5, 1.5, 3
MAX_ELEMS = 4e6   # pairs x frames per chunk
KT = {'Mpipi-GG': 1.98720425864083e-3 * T, 'CALVADOS-2': 8.314462618e-3 * T}
_RES = None


def mpipi_dir(uid):
    for root in MPIPI_RUNS.values():
        if os.path.isfile(os.path.join(root, uid, 'traj.xtc')):
            return os.path.join(root, uid)
    return None


def pair_energies_mpipi(uid):
    """Frames (nm), pair indices and u(r_nm) -> Wang-Frenkel pair energies (kcal/mol)."""
    d = mpipi_dir(uid)
    if d is None:
        raise FileNotFoundError(f'no Mpipi-GG run folder for {uid}')
    seq = bc.read_mpipi_sequence(f'{d}/seq.dat')
    xyz, _ = bc.load_mpipi(f'{d}/traj.xtc', n_res=len(seq))                     # nm, unwrapped
    xyz = np.asarray(xyz, dtype=float)[bc.even_indices(len(xyz), mh.K_FRAMES)][::STRIDE]
    pp = mh.mpipi_pairs(mh.read_mpipi_system(d))

    def u(r_nm):
        r = r_nm * 10.0
        v = pp['eps_alpha'] * ((pp['sig'] / r) ** (2 * pp['mu']) - 1) * ((pp['rc'] / r) ** (2 * pp['mu']) - 1) ** (2 * pp['nu'])
        return np.where(r < pp['rc'], v, 0.0)
    return xyz, pp['i'], pp['j'], u, pp['sig'] / 10.0


def pair_energies_calvados(uid, hint):
    global _RES
    if _RES is None:
        _RES = L.load_residue_params()
    d = L.find_traj_dir(uid, L.TRAJ_ROOT_FMT.format(T=300), hint, flat=False)
    if d is None:
        raise FileNotFoundError(f'no CALVADOS-2 300 K trajectory for {uid}')
    xyz, seq, _ = bc.load_calvados(os.path.join(d, f'{uid}.dcd'), os.path.join(d, 'top.pdb'), k=L.K_FRAMES)
    xyz = np.asarray(xyz, dtype=float)[::STRIDE]
    q, lam, sig, _ = fdt.calvados_residue_params(seq, _RES)
    i, j, s, l, shift, _ = fdt.pair_setup(q, lam, sig)
    return xyz, i, j, lambda r: fdt.ah_pair_energy(r, s, l, shift), s


def analyse(args):
    uid, model, hint = args
    try:
        xyz, i, j, ufun, sig = pair_energies_mpipi(uid) if model == 'Mpipi-GG' else pair_energies_calvados(uid, hint)
        rcut = CUT_SIGMA * sig
        kT, n, nf = KT[model], xyz.shape[1], xyz.shape[0]
        rg = np.sqrt(((xyz - xyz.mean(1, keepdims=True)) ** 2).sum(2).mean(1))
        Ec, Et, ncont, nclose = (np.zeros(nf) for _ in range(4))
        sep = (np.abs(j - i) >= MIN_SEP)[None, :]
        step = max(1, int(MAX_ELEMS // len(i)))
        for a in range(0, nf, step):
            r = np.linalg.norm(xyz[a:a + step, i] - xyz[a:a + step, j], axis=2)
            u, close = ufun(r), r < rcut
            Ec[a:a + step] = np.where(close, u, 0.0).sum(1)
            Et[a:a + step] = np.where(close, 0.0, u).sum(1)
            ncont[a:a + step] = (close & sep).sum(1)
            nclose[a:a + step] = close.sum(1)
        E = Ec + Et
        cov = lambda a, b: np.cov(a, b)[0, 1]  # noqa: E731
        cE = cov(rg, E)
        return {'UniProt_ID': uid, 'model': model, 'N': n, 'frames': len(rg),
                'contacts_per_res': ncont.mean() / n, 'E_per_res_kT': E.mean() / kT / n, 'close_share': Ec.mean() / E.mean(),
                'E_per_contact_kT': Ec.mean() / kT / max(nclose.mean(), 1e-9),
                'r_E_contacts': np.corrcoef(E, ncont)[0, 1], 'r_rg_contacts': np.corrcoef(rg, ncont)[0, 1],
                'r_rg_E': np.corrcoef(rg, E)[0, 1], 'chi_close_share': cov(rg, Ec) / cE,
                'var_close_share': (Ec.var(ddof=1) + cov(Ec, Et)) / E.var(ddof=1),
                'sdE_kT_per_sqrtN': E.std(ddof=1) / kT / np.sqrt(n), 'status': 'ok'}
    except Exception as e:
        return {'UniProt_ID': uid, 'model': model, 'status': f'failed: {e!r}'[:200], 'trace': traceback.format_exc()[-300:]}


def main(workers):
    os.makedirs(OUT_DIR, exist_ok=True)
    src = pd.read_csv(SRC, index_col=0)
    src = src[src['b'].isin(BINS)]
    pick = pd.concat([g.sample(min(len(g), PER_BIN), random_state=SEED) for _, g in src.groupby(['model', 'b'])])
    kingdom = pd.read_csv(bc.BIO_CSV, usecols=['UniProt_ID', 'kingdom']).set_index('UniProt_ID')['kingdom']
    tasks = [(u, m, kingdom.get(u)) for u, m in zip(pick.index, pick['model'])]
    print(f'{len(tasks)} proteins ({pick.groupby(["model", "b"]).size().to_dict()}), {workers} workers', flush=True)
    with mp.Pool(workers) as pool:
        res = []
        for k, r in enumerate(pool.imap_unordered(analyse, tasks), 1):
            res.append(r)
            if r['status'] != 'ok' or k % 40 == 0:
                print(f"[{k}/{len(tasks)}] {r['model']} {r['UniProt_ID']} {r['status']}", flush=True)
    per = pd.DataFrame(res).merge(pick[['model', 'nu', 'b']].reset_index(), on=['UniProt_ID', 'model'])
    per.to_csv(f'{OUT_DIR}/per_protein.csv', index=False)
    ok = per[per['status'] == 'ok']
    cols = ['contacts_per_res', 'E_per_res_kT', 'close_share', 'E_per_contact_kT', 'r_E_contacts', 'r_rg_contacts',
            'r_rg_E', 'chi_close_share', 'var_close_share', 'sdE_kT_per_sqrtN', 'N']
    g = ok.groupby(['model', 'b'])
    bins = g[cols].median()
    bins.insert(0, 'n', g.size())
    bins.to_csv(f'{OUT_DIR}/bins.csv')
    bc.write_provenance(f'{OUT_DIR}/bins.csv', __file__, {'T': T, 'frames': f'every {STRIDE}th of the 1,000 FDT frames',
                                                        'contact_cutoff': f'{CUT_SIGMA} sigma_ij', 'min_sep': MIN_SEP,
                                                        'per_bin': PER_BIN, 'seed': SEED}, [SRC])
    with pd.option_context('display.width', 250, 'display.max_columns', 20):
        print(bins.round(3).to_string())
    print(f"{(per['status'] != 'ok').sum()} failed; saved {OUT_DIR}")


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--workers', type=int, default=6)
    main(ap.parse_args().workers)
