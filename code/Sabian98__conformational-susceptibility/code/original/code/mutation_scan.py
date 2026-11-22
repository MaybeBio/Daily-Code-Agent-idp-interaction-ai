"""
mutation_scan.py — in-silico mutational scan of the heating response (CALVADOS-2, 300 K), from the stored frames.

For every position p and every substitution a -> b, the energy change per frame is exact for the frames' coordinates:
  Ashbaugh-Hatch pairs of p with the mutant's sigma and lambda (pair rule as fdt.pair_setup), Debye-Hueckel with the
  mutant's charge (His at pH 7.5 as calvados_residue_params), and the solvation forcing (TEA, gamma 2.0, and hydrophilic,
  gamma 0.358, per residue type as tea_fdt.py) with the mutant's pair energies.
  dG = -dH / (R T^2) + d(dV_DH/dT) / (R T) + dD / (R T)     (change of the thermal force)
  new total = -Cov_w(Rg, G_tot + dG) / <Rg>_w with w = exp(-dH / RT); change in % per K, 10-block jackknife SE, Kish
  effective sample fraction (reweighting is reliable only where it is large; mutations that reshape the chain are flagged).
Sample: the Arg shortlist (arg_obs/A_arg_shortlist.csv) and 25 random chains of <= 150 residues from heat_crosses' random set.
    python mutation_scan.py --out DIR [--workers N] [--forcing A0|capped] [--sample earlier/sample.csv]
  --forcing capped (2026-10-25): the capped main model (heat_crosses.FORCING); a substitution at p changes p's coefficient
  and, through p's new pair shapes, the contact numbers of p and of each partner, all inside the cap.
Output: DIR/mutants.csv (append; resumable), DIR/chains.csv
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

TF, L, fdt, bc = HC.TF, HC.L, HC.fdt, HC.bc
T0, RT, CAL = HC.T0, HC.RT, HC.CAL
AA = 'ACDEFGHIKLMNPQRSTVWY'


def one(args):
    uid, d = args
    t0 = time.time()
    try:
        z = np.load(os.path.join(TF.NPZ, f'{uid}.npz'))
        xyz, seq, info = bc.load_calvados(os.path.join(d, f'{uid}.dcd'), os.path.join(d, 'top.pdb'), k=L.K_FRAMES)
        fi = info['production_start'] + info['stride'] * np.arange(len(xyz))
        if not np.array_equal(fi, z['frame_index']):
            raise ValueError('frames differ from the stored FDT frames')
        res = L.load_residue_params()
        rg, G = z['rg'].astype(float), z['G_T'].astype(float)
        N, nf = len(seq), len(xyz)
        q, lam, sig, _ = fdt.calvados_residue_params(seq, res)
        i, j, s, l, shift, _ = fdt.pair_setup(q, lam, sig)
        dlt = {a: TF.DL['charmm36'][a] for a in AA}
        dlh = {a: TF.DL_CLASSES['hyd'][a] for a in AA}
        def fres(sq):
            return np.array([HC.G_TEA * dlt[a] for a in sq]), np.array([HC.G_HYD * dlh[a] for a in sq])
        ft, fh = fres(seq)
        fpair = 0.5 * (ft[i] + ft[j]) + 0.5 * (fh[i] + fh[j])
        eps, deps, k, dk, rc = CAL['eps_yu'], CAL['deps_yu_dT'], CAL['kappa'], CAL['dkappa_dT'], fdt.RC_YU
        ec = np.exp(-k * rc)
        # base solvation forcing per frame
        D = np.zeros(nf)
        capped = HC.FORCING == 'capped'
        Sres = np.zeros((nf, N)) if capped else None
        step = max(1, int(fdt.MAX_ELEMS // max(1, len(i))))
        for a in range(0, nf, step):
            r = np.linalg.norm(xyz[a:a + step, i] - xyz[a:a + step, j], axis=2)
            S_ = fdt.ah_pair_energy(r, s, np.ones_like(l), shift) - fdt.ah_pair_energy(r, s, np.zeros_like(l), shift)
            D[a:a + step] = (S_ * fpair).sum(1)
            if capped:
                Sres[a:a + step] = np.stack([np.bincount(i, S_[t], N) + np.bincount(j, S_[t], N) for t in range(len(S_))])
        if capped:                                 # capped main model: D = -eps/2 sum_i f_i h(c_i), c_i = -S_i / eps
            e0 = fdt.EPS_AH
            fcap = {a: HC.G_TEA * HC.CP.F_TEA[a] + HC.G_HYD * HC.CP.F_CHG[a] for a in AA}
            fr = np.array([fcap[a] for a in seq])
            cres = -Sres / e0
            hres = HC.CP.h(cres)
            D = -0.5 * e0 * (hres @ fr)
        Gtot = G + D / RT
        m = rg.mean()
        base = -fdt._cov(rg, Gtot) / m * 100
        rows = []
        for p in range(N):
            idx = np.flatnonzero((i == p) | (j == p))
            other = np.where(i[idx] == p, j[idx], i[idx])
            r = np.linalg.norm(xyz[:, [p]] - xyz[:, other], axis=2)           # frames x partners
            inside = r < rc
            er = np.exp(-k * r)
            Y = np.where(inside, er / r - ec / rc, 0.0)
            dY = np.where(inside, deps * Y + eps * dk * (ec - er), 0.0)
            u1o = fdt.ah_pair_energy(r, s[idx], np.ones(len(idx)), shift[idx])
            u0o = fdt.ah_pair_energy(r, s[idx], np.zeros(len(idx)), shift[idx])
            AHo = (u0o + l[idx] * (u1o - u0o)).sum(1)
            Do = ((u1o - u0o) * fpair[idx]).sum(1)
            qY, qdY = (q[other] * eps * Y).sum(1), (q[other] * dY).sum(1)
            for b in AA:
                if b == seq[p]:
                    continue
                mseq = seq[:p] + b + seq[p + 1:]
                qm, lm, sm, _ = fdt.calvados_residue_params(mseq, res)
                _, _, s2, l2, sh2, _ = fdt.pair_setup(qm, lm, sm)
                u1n = fdt.ah_pair_energy(r, s2[idx], np.ones(len(idx)), sh2[idx])
                u0n = fdt.ah_pair_energy(r, s2[idx], np.zeros(len(idx)), sh2[idx])
                ftm, fhm = fres(mseq)
                fpn = 0.5 * (ftm[p] + ftm[other]) + 0.5 * (fhm[p] + fhm[other])
                dq = qm[p] - q[p]
                dH = (u0n + l2[idx] * (u1n - u0n)).sum(1) - AHo + dq * qY
                if capped:                         # p's own term with its new coefficient; p and partners' new contacts
                    dS = (u1n - u0n) - (u1o - u0o)                                    # frames x partners
                    hp = HC.CP.h(cres[:, p] - dS.sum(1) / e0)
                    dD = -0.5 * e0 * (fcap[b] * hp - fr[p] * hres[:, p]
                                      + ((HC.CP.h(cres[:, other] - dS / e0) - hres[:, other]) * fr[other]).sum(1))
                else:
                    dD = ((u1n - u0n) * fpn).sum(1) - Do
                dG = -dH / (fdt.R_GAS * T0 ** 2) + dq * qdY / RT + dD / RT
                x = -dH / RT
                w = np.exp(x - x.max())
                f = lambda ix: (-HC.wcov(rg[ix], Gtot[ix] + dG[ix], w[ix]) / ((w[ix] * rg[ix]).sum() / w[ix].sum())  # noqa: E731
                                + fdt._cov(rg[ix], Gtot[ix]) * (len(ix) - 1) / len(ix) / rg[ix].mean()) * 100
                rows.append({'UniProt_ID': uid, 'pos': p + 1, 'wt': seq[p], 'mut': b, 'dtotal': f(np.arange(nf)),
                             'dtotal_se': HC.jk(f, nf), 'ess': float(w.sum() ** 2 / (w ** 2).sum() / nf),
                             'dRg_pct': float(((w * rg).sum() / w.sum() / m - 1) * 100)})
        return {'UniProt_ID': uid, 'N': N, 'total_rel': base, 'seconds': round(time.time() - t0, 1)}, pd.DataFrame(rows), None
    except Exception as ex:  # noqa: BLE001
        return {'UniProt_ID': uid}, None, f'{uid}: {ex!r} | {traceback.format_exc()[-300:]}'.replace('\n', ' ')


def sample(seed=13):
    TLA = HC.TLA
    arg = list(pd.read_csv(os.path.join(TLA, 'Analysis_FDT', 'total_response', 'arg_obs', 'A_arg_shortlist.csv')).UniProt_ID)
    hx = os.path.join(TLA, 'Analysis_FDT', 'total_response', 'heat_crosses')
    smp = pd.read_csv(os.path.join(hx, 'sample.csv'))
    ch = pd.read_csv(os.path.join(hx, 'chains.csv')).drop_duplicates('UniProt_ID').set_index('UniProt_ID')
    rnd = [u for u in smp.loc[smp.random, 'UniProt_ID'] if u in ch.index and ch.at[u, 'N'] <= 150 and u not in arg]
    rnd = list(np.random.default_rng(seed).choice(rnd, 25, replace=False))
    return arg + rnd, set(arg)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', required=True)
    ap.add_argument('--workers', type=int, default=6)
    ap.add_argument('--limit', type=int, default=None)
    ap.add_argument('--forcing', choices=['A0', 'capped'], default='A0')
    ap.add_argument('--sample', default=None, help='sample.csv of an earlier run: same chains')
    a = ap.parse_args()
    os.environ['BENDER_SCAN_FORCING'] = a.forcing           # read by heat_crosses in the spawned workers
    os.makedirs(a.out, exist_ok=True)
    ids, arg = sample()
    if a.sample:
        t = pd.read_csv(a.sample)
        ids, arg = list(t.UniProt_ID), set(t.loc[t.arg_shortlist.astype(bool), 'UniProt_ID'])
    pd.DataFrame({'UniProt_ID': ids, 'arg_shortlist': [u in arg for u in ids]}).to_csv(os.path.join(a.out, 'sample.csv'), index=False)
    pc = os.path.join(a.out, 'chains.csv')
    done = set(pd.read_csv(pc, usecols=['UniProt_ID']).UniProt_ID) if os.path.exists(pc) else set()
    dirs = TF.index_dirs([u for u in ids if u not in done])
    todo = [(u, d) for u, d in dirs.items()][:a.limit]
    print(f'{len(ids)} chains, {len(done)} done, {len(todo)} to run', flush=True)
    nf = 0
    with mp.get_context('spawn').Pool(a.workers) as pool:
        for k, (row, M_, err) in enumerate(pool.imap_unordered(one, todo), 1):
            if err:
                nf += 1
                print('failed', err[:300], flush=True)
                continue
            pd.DataFrame([row]).to_csv(pc, mode='a', header=not os.path.exists(pc), index=False)
            pm = os.path.join(a.out, 'mutants.csv')
            M_.to_csv(pm, mode='a', header=not os.path.exists(pm), index=False)
            print(f'{k}/{len(todo)} done ({row["seconds"]} s), {nf} failed', flush=True)
    print('finished', flush=True)


if __name__ == '__main__':
    main()
