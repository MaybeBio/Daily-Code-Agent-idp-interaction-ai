"""
mutation_scan_recharged.py — the in-silico mutational scan of mutation_scan.py in Mpipi-Recharged (300 K, 150 mM), on
the same chains, for a mutant-by-mutant comparison between the two calibrated models.

Energy terms as OpenMpipi.system_building.get_mpipi_system (recharged_params.txt; residue order MGKTRADEYVLQWFSHNPCI):
  Wang-Frenkel  eps_ij alpha ((s/r)^(2mu) - 1)((3s/r)^(2mu) - 1)^2 for r < 3 s, pair-specific eps, s, mu
  Yukawa        A_ij exp(-kappa r) / r for r < 3.5 nm, pair-specific A_ij; kappa from the salt and the temperature-dependent
                water dielectric (calculate_debye_length), so dH/dT = -(dkappa/dT) sum A_ij exp(-kappa r) for r < 3.5 nm
  bonded neighbours excluded; minimum image in the run's box (mpipi_recharged_fdt.box_nm).
Solvation forcing on the attractive Wang-Frenkel part (tea_crossff.wf_attractive), TEA at 0.4 and hydrophilic at 0.209
(the Wuttke calibration). For every substitution a -> b at position p: exact per-frame dH and change of the thermal force
dG = -dH/(R T^2) + d(dH/dT)/(R T) + dD/(R T); new total heating response by exact reweighting (as mutation_scan.py).
The unmutated total is not checked in this script; an external check (2026-10-10) found all 50 chains equal to the
stored proteome totals to 1e-8, and mutants equal to OpenMpipi full-sequence energies to 1e-6.
    python mutation_scan_recharged.py --out DIR [--workers N] [--forcing A0|capped] [--sample FILE]
  --forcing capped (2026-10-25): the capped main model, as mutation_scan.py (contact measure: the attractive WF shape).
Output: DIR/mutants_recharged.csv, DIR/chains_recharged.csv (append; resumable)
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
import heat_crosses as HC  # noqa: E402
import tea_crossff as TC  # noqa: E402

TF = TC.TF
T0 = 300.0
R_KJ = TC.R_KJ
RT = R_KJ * T0
G_TEA, G_HYD = 0.4, 0.209
if HC.FORCING == 'capped':                         # --forcing capped: the capped main model (heat_crosses.FORCING)
    _gh, _gt = HC.CM.params()
    G_TEA, G_HYD = _gt['R'], _gh['R'][0]
CSX = 150.0
AA = 'ACDEFGHIKLMNPQRSTVWY'
ORDER = 'MGKTRADEYVLQWFSHNPCI'


def tables():
    import OpenMpipi
    P = np.loadtxt(os.path.join(os.path.dirname(OpenMpipi.__file__), 'data_files', 'recharged_params.txt'))
    w = P[:-21 * 21].reshape(3, 21, 21)            # Discrete3DFunction: value[x + 21 y + 441 z] -> w[z, y, x]
    A = P[-21 * 21:].reshape(21, 21)               # Discrete2DFunction: value[x + 21 y] -> A[y, x]
    return {'eps': w[0], 'sig': w[1], 'mu': np.floor(w[2]), 'A': A, 'ix': {a: k for k, a in enumerate(ORDER)}}


def kappa_of(T):
    from OpenMpipi.system_building import calculate_debye_length
    return 1.0 / calculate_debye_length(T, CSX)


def pair_terms(r, ia, ib, tb, k, fa, fb):
    """Per-frame sums over pairs (a, b) with distances r [frames x pairs]: WF energy, Yukawa energy, Yukawa exp sum (for
    dH/dT), solvation forcing D (kJ/mol/K at the strengths above)."""
    eps, sig, mu = tb['eps'][ib, ia], tb['sig'][ib, ia], tb['mu'][ib, ia]
    rc = 3 * sig
    alpha = 2 * 3 ** (2 * mu) * (3 / (2 * (3 ** (2 * mu) - 1))) ** 3
    shape = np.where(r < rc, alpha * ((sig / r) ** (2 * mu) - 1) * ((rc / r) ** (2 * mu) - 1) ** 2, 0.0)
    Awf = (shape * eps).sum(1)
    Aij = tb['A'][ib, ia]
    ex = np.where(r < 3.5, np.exp(-k * r), 0.0)
    Y = (Aij * ex / r).sum(1)
    Yexp = (Aij * ex).sum(1)
    att = TC.wf_attractive(r, sig, mu, 1, rc, alpha)
    D = (att * 0.5 * (fa + fb) * 4.184).sum(1)
    return Awf, Y, Yexp, D


def wf_att(r, ia, ib, tb):
    """Attractive Wang-Frenkel shape per pair (s0 = 1), the contact measure of the capped forcing."""
    sig, mu = tb['sig'][ib, ia], tb['mu'][ib, ia]
    rc = 3 * sig
    alpha = 2 * 3 ** (2 * mu) * (3 / (2 * (3 ** (2 * mu) - 1))) ** 3
    return TC.wf_attractive(r, sig, mu, 1, rc, alpha)


def one(args):
    uid, d = args
    t0 = time.time()
    try:
        import mpipi_recharged_fdt as R
        z = np.load(os.path.join(TC.RC_NPZ, f'{uid}.npz'))
        xyz, seq, info = R.load_run(d, uid)
        fi = info['production_start'] + info['stride'] * np.arange(len(xyz))
        if not np.array_equal(fi, z['frame_index']):
            raise ValueError('frames differ from the stored FDT frames')
        tb = tables()
        k = kappa_of(T0)
        dk = (kappa_of(T0 + 0.5) - kappa_of(T0 - 0.5)) / 1.0
        rg, G = z['rg'].astype(float), z['G_T'].astype(float)
        N, nf, L = len(seq), len(xyz), R.box_nm(len(seq))
        ix = np.array([tb['ix'][a] for a in seq])
        ft = {a: G_TEA * TF.DL['charmm36'][a] + G_HYD * TF.DL_CLASSES['hyd'][a] for a in AA}
        fv = np.array([ft[a] for a in seq])
        i, j = np.triu_indices(N, k=2)
        D = np.zeros(nf)
        capped = HC.FORCING == 'capped'
        if capped:
            Mi = HC.CP.incidence(N, i, j)
            Sres = np.zeros((nf, N))
        step = max(1, int(TC.MAX_ELEMS // max(1, len(i))))
        for a in range(0, nf, step):
            dv = xyz[a:a + step, i].astype(float) - xyz[a:a + step, j].astype(float)
            dv -= L * np.round(dv / L)
            r = np.linalg.norm(dv, axis=2)
            D[a:a + step] = pair_terms(r, ix[i], ix[j], tb, k, fv[i], fv[j])[3]
            if capped:
                Sres[a:a + step] = (Mi.T @ wf_att(r, ix[i], ix[j], tb).T).T
        if capped:                                 # D = -1/2 sum_i f_i h(c_i), c_i = -S_i (kcal coefficients x 4.184)
            fcap = {a: 4.184 * (G_TEA * HC.CP.F_TEA[a] + G_HYD * HC.CP.F_CHG[a]) for a in AA}
            fr = np.array([fcap[a] for a in seq])
            cres = -Sres
            hres = HC.CP.h(cres)
            D = -0.5 * (hres @ fr)
        Gtot = G + D / RT
        m = rg.mean()
        base = -HC.fdt._cov(rg, Gtot) / m * 100
        rows = []
        for p in range(N):
            other = np.array([q for q in range(N) if abs(q - p) >= 2])
            dv = xyz[:, [p]].astype(float) - xyz[:, other].astype(float)
            dv -= L * np.round(dv / L)
            r = np.linalg.norm(dv, axis=2)
            o_ix = ix[other]
            W0, Y0, E0, D0 = pair_terms(r, np.full(len(other), ix[p]), o_ix, tb, k, np.full(len(other), fv[p]), fv[other])
            for b in AA:
                if b == seq[p]:
                    continue
                W1, Y1, E1, D1 = pair_terms(r, np.full(len(other), tb['ix'][b]), o_ix, tb, k, np.full(len(other), ft[b]), fv[other])
                dH = (W1 - W0) + (Y1 - Y0)
                if capped:                         # p's new coefficient; p and partners' contacts with p's new shapes
                    dS = wf_att(r, np.full(len(other), tb['ix'][b]), o_ix, tb) - wf_att(r, np.full(len(other), ix[p]), o_ix, tb)
                    dD = -0.5 * (fcap[b] * HC.CP.h(cres[:, p] - dS.sum(1)) - fr[p] * hres[:, p]
                                 + ((HC.CP.h(cres[:, other] - dS) - hres[:, other]) * fr[other]).sum(1))
                else:
                    dD = D1 - D0
                dG = -dH / (R_KJ * T0 ** 2) + (-dk * (E1 - E0)) / RT + dD / RT
                x = -dH / RT
                w = np.exp(x - x.max())
                f = lambda ixs: (-HC.wcov(rg[ixs], Gtot[ixs] + dG[ixs], w[ixs]) / ((w[ixs] * rg[ixs]).sum() / w[ixs].sum())  # noqa: E731
                                 + HC.fdt._cov(rg[ixs], Gtot[ixs]) * (len(ixs) - 1) / len(ixs) / rg[ixs].mean()) * 100
                rows.append({'UniProt_ID': uid, 'pos': p + 1, 'wt': seq[p], 'mut': b, 'dtotal': f(np.arange(nf)),
                             'dtotal_se': HC.jk(f, nf), 'ess': float(w.sum() ** 2 / (w ** 2).sum() / nf),
                             'dRg_pct': float(((w * rg).sum() / w.sum() / m - 1) * 100)})
        return {'UniProt_ID': uid, 'N': N, 'total_rel': base, 'seconds': round(time.time() - t0, 1)}, pd.DataFrame(rows), None
    except Exception as ex:  # noqa: BLE001
        return {'UniProt_ID': uid}, None, f'{uid}: {ex!r} | {traceback.format_exc()[-300:]}'.replace('\n', ' ')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', required=True)
    ap.add_argument('--workers', type=int, default=6)
    ap.add_argument('--limit', type=int, default=None)
    ap.add_argument('--forcing', choices=['A0', 'capped'], default='A0')
    ap.add_argument('--sample', default=None, help='sample.csv to use (default: the one in --out)')
    a = ap.parse_args()
    os.environ['BENDER_SCAN_FORCING'] = a.forcing           # read by heat_crosses in the spawned workers
    ids = list(pd.read_csv(a.sample or os.path.join(a.out, 'sample.csv')).UniProt_ID)
    os.makedirs(a.out, exist_ok=True)
    pc = os.path.join(a.out, 'chains_recharged.csv')
    done = set(pd.read_csv(pc, usecols=['UniProt_ID']).UniProt_ID) if os.path.exists(pc) else set()
    dirs = TC.rc_dirs()
    todo = [(u, dirs[u]) for u in ids if u not in done and u in dirs][:a.limit]
    print(f'{len(ids)} chains, {len(done)} done, {len(todo)} to run', flush=True)
    nf = 0
    with mp.get_context('spawn').Pool(a.workers) as pool:
        for kk, (row, M_, err) in enumerate(pool.imap_unordered(one, todo), 1):
            if err:
                nf += 1
                print('failed', err[:300], flush=True)
                continue
            pd.DataFrame([row]).to_csv(pc, mode='a', header=not os.path.exists(pc), index=False)
            pm = os.path.join(a.out, 'mutants_recharged.csv')
            M_.to_csv(pm, mode='a', header=not os.path.exists(pm), index=False)
            print(f'{kk}/{len(todo)} done ({row["seconds"]} s), {nf} failed', flush=True)
    print('finished', flush=True)


if __name__ == '__main__':
    main()
