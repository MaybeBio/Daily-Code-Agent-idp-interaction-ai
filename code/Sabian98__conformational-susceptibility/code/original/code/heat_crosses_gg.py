"""
heat_crosses_gg.py — heat x phosphorylation and heat x salt in Mpipi-GG (part C of heat_crosses.py in a second force
field), to test the 'compaction amplifies the heating response' rule across models.

Mpipi-GG at 300 K. Thermal force G = -H/(R T^2) (its energy function has no T dependence) plus the solvation terms at
the strengths used throughout (gamma_tea 0.30, gamma_hyd 0.106; per-frame forcings as tea_crossff.py, attractive part):
total = -Cov(Rg, G_tot) / <Rg> (% per K).
  phosphate at Ser/Thr i: charge -1.5 (Mpipi-GG charges are scaled by 0.75; as salt_ph_phospho.py), dH = dq phi_i,
  dG = -dH / (R T^2)
  doubling the ionic strength: kappa -> sqrt(2) kappa in the Debye-Hueckel term, dH = E_DH(sqrt(2) kappa) - E_DH(kappa)
new total by exact reweighting; 10-block jackknife SE; Kish effective sample fraction.
    python heat_crosses_gg.py --ids FILE --out DIR [--workers N]
Output: DIR/chains_gg.csv, DIR/sites_gg.csv (append; resumable)
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
import salt_ph_phospho as SP  # noqa: E402
import total_response_noise as N  # noqa: E402

bs, mh = SP.bs, SP.mh
T0 = 300.0
R_KCAL = bs.R_KCAL
RT = R_KCAL * T0
DQ = -1.5


def one(uid):
    t0 = time.time()
    try:
        xyz, G, rg = bs.frames_mpipi(uid, 300)
        rg2, G2, Dh, Dt, RTn = N.mpipi(uid, N.TC.gg_dirs()[uid], 'M')
        if not (np.allclose(rg, rg2) and np.allclose(G, G2)):
            raise ValueError('frames differ between loaders')
        Gtot = G + (N.GAMMA['M'] * Dh + N.GAMMA_TEA['M'] * Dt) / RT
        z = np.load(os.path.join(bs.MP_DIR, 'npz_300K', f'{uid}.npz'))
        seq = str(z['sequence'])
        sysd = mh.read_mpipi_system(os.path.join(bs.MP_RUNS[300], uid))
        pp = mh.mpipi_pairs(sysd)
        i, j = pp['i'], pp['j']
        q = np.array([sysd['qtype'].get(int(t), 0.0) for t in sysd['types']])
        rc_all, rc_c = max(sysd['coul'].values()), pp['rc_c']
        M = SP.site_matrix(i, j, q, len(q))
        k, pref = sysd['kappa'], mh.QQR2E / sysd['diel']
        k2 = k * np.sqrt(2.0)
        nf = len(xyz)
        E1, E2, phi = np.zeros(nf), np.zeros(nf), np.zeros((nf, len(q)))
        step = max(1, SP.MAX_ELEMS // max(1, len(i)))
        for a in range(0, nf, step):
            r = 10.0 * np.linalg.norm(xyz[a:a + step, i] - xyz[a:a + step, j], axis=2)
            E1[a:a + step] = np.where(r < rc_c, pref * pp['qq'] * np.exp(-k * r) / r, 0.0).sum(1)
            E2[a:a + step] = np.where(r < rc_c, pref * pp['qq'] * np.exp(-k2 * r) / r, 0.0).sum(1)
            phi[a:a + step] = (M @ np.where(r < rc_all, pref * np.exp(-k * r) / r, 0.0).T).T
        ez = z['E_DH'].astype(float)
        if np.abs(ez).max() > 1e-6 and (np.corrcoef(E1, ez)[0, 1] < 0.999 or np.abs(E1 - ez).max() > 1e-4 * np.abs(ez).max()):
            raise ValueError('recomputed E_DH does not match the npz')   # absolute tolerance added 2026-10-10
        m = rg.mean()
        base = -np.cov(rg, Gtot)[0, 1] / m * 100

        def perturbed(dH):
            dG = -dH / (R_KCAL * T0 ** 2)
            x = -dH / RT
            w = np.exp(x - x.max())
            f = lambda ix: (-HC.wcov(rg[ix], Gtot[ix] + dG[ix], w[ix]) / ((w[ix] * rg[ix]).sum() / w[ix].sum())  # noqa: E731
                            + np.cov(rg[ix], Gtot[ix], ddof=0)[0, 1] / rg[ix].mean()) * 100
            return f(np.arange(nf)), HC.jk(f, nf), float(w.sum() ** 2 / (w ** 2).sum() / nf), float(((w * rg).sum() / w.sum() / m - 1) * 100)
        row = {'UniProt_ID': uid, 'N': len(seq), 'total_rel': base, 'seconds': 0.0}
        row['salt_dtotal'], row['salt_dtotal_se'], row['salt_ess'], row['salt_dRg_pct'] = perturbed(E2 - E1)
        sites = []
        for p_, c in enumerate(seq):
            if c in 'ST':
                dt, se, ess, drg = perturbed(DQ * phi[:, p_])
                sites.append({'UniProt_ID': uid, 'pos': p_ + 1, 'residue': c, 'phospho_dtotal': dt, 'phospho_dtotal_se': se,
                              'ess': ess, 'phospho_dRg_pct': drg})
        row['seconds'] = round(time.time() - t0, 1)
        return row, pd.DataFrame(sites), None
    except Exception as ex:  # noqa: BLE001
        return {'UniProt_ID': uid}, None, f'{uid}: {ex!r} | {traceback.format_exc()[-300:]}'.replace('\n', ' ')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--ids', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--workers', type=int, default=6)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    ids = sorted({l.strip() for l in open(a.ids) if l.strip()})
    pc = os.path.join(a.out, 'chains_gg.csv')
    done = set(pd.read_csv(pc, usecols=['UniProt_ID']).UniProt_ID) if os.path.exists(pc) else set()
    todo = [u for u in ids if u not in done]
    print(f'{len(ids)} chains, {len(done)} done, {len(todo)} to run', flush=True)
    nf = 0
    with mp.get_context('spawn').Pool(a.workers) as pool:
        for k, (row, S_, err) in enumerate(pool.imap_unordered(one, todo), 1):
            if err:
                nf += 1
                print('failed', err[:300], flush=True)
                continue
            pd.DataFrame([row]).to_csv(pc, mode='a', header=not os.path.exists(pc), index=False)
            ps = os.path.join(a.out, 'sites_gg.csv')
            if len(S_):
                S_.to_csv(ps, mode='a', header=not os.path.exists(ps), index=False)
            if k % 25 == 0:
                print(f'{k}/{len(todo)} done, {nf} failed', flush=True)
    print('finished', flush=True)


if __name__ == '__main__':
    main()
