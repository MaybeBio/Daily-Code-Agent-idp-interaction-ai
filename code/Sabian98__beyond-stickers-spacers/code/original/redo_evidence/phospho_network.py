"""
phospho_network.py — step 1 of phospho_network_plan.md (hash 5b113f57). Per-frame contact networks on BENDER-MT's 1,000
CALVADOS-2 frames, charge-model reweighting for every Ser/Thr of the phosphosite-test proteins (dq = -2), and the PAGE4
validation against direct phosphoform runs. No annotation is used here except to choose the proteins.

  BENDER_REANALYSIS_DIR=/Users/taseef/Downloads/Scripts/bender_reanalysis \
  ~/.venvs/finches/bin/python phospho_network.py
      ->  phospho_network/sites_network.csv, phospho_network/profiles/<uid>.npz, phospho_network/page4_validation.txt,
          phospho_network/page4_profiles.npz
"""
import os
import sys

os.environ.setdefault('BENDER_REANALYSIS_DIR', '/Users/taseef/Downloads/Scripts/bender_reanalysis')
MT = '/Volumes/User Homes/Nu Project/BENDER-MT'
DYN = '/Volumes/User Homes/Nu Project/BENDER-DYN/redo_evidence'
sys.path.insert(0, os.path.join(MT, 'code'))
sys.path.insert(0, DYN)

import multiprocessing as mp  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402

import salt_ph_phospho as sp  # noqa: E402
import hub_frames as hf  # noqa: E402

OUT = os.path.join(DYN, 'phospho_network')
SITES = os.path.join(MT, 'Traj Level Analysis', 'Analysis_FDT', 'phosphosite_test', 'sites.csv')
DQ = -2.0
Q = {'K': 1, 'R': 1, 'D': -1, 'E': -1}


def networks(xyz, n, key):
    """Degree (frames x n) and hub indicator (frames x n) with the paper's rule; xyz in nm."""
    rng = hf.tie_rng(key)
    n_top = max(1, round(0.1 * n))
    D = np.zeros((len(xyz), n), np.int16)
    H = np.zeros((len(xyz), n), bool)
    for f in range(len(xyz)):
        X = (10.0 * np.asarray(xyz[f], float)).T
        adj = hf.contact_lists([list(X[0]), list(X[1]), list(X[2])], n)
        D[f] = [len(a) for a in adj]
        for i in hf.top_hubs(hf.betweenness(adj, n), n_top, rng):
            H[f, i] = True
    return D, H


def qloc(seq):
    return np.array([sum(Q.get(seq[j], 0) for j in range(max(0, i - 5), min(len(seq), i + 6)) if j != i)
                     for i in range(len(seq))])


def one(task):
    uid, hint, positions = task
    try:
        seq, rg, _, _, phi, beta = sp.calvados_terms(uid, hint)
        xyz = sp.bs.frames_calvados(uid, 300, hint)[0]
        n = len(seq)
        D, H = networks(xyz, n, uid)
        deg0, occ0 = D.mean(0), H.mean(0)
        Hf = H.astype(float)
        ql = qloc(seq)
        ncpr = sum(Q.get(c, 0) for c in seq) / n
        rows, dprof, oprof = [], [], []
        for pos in positions:
            i = pos - 1
            if seq[i] not in 'ST':
                raise ValueError(f'position {pos} is {seq[i]}, not S/T')
            x = -beta * DQ * phi[:, i]
            w = np.exp(x - x.max())
            w /= w.sum()
            ess = float(1.0 / (w ** 2).sum() / len(w))
            dw, ow = w @ D, w @ Hf
            dd, do = dw - deg0, ow - occ0
            far = np.abs(np.arange(n) - i) >= 3
            kr = far & np.isin(np.array(list(seq)), list('KR'))
            de = far & np.isin(np.array(list(seq)), list('DE'))
            rows.append({'UniProt_ID': uid, 'pos': pos, 'residue': seq[i], 'N': n, 'NCPR': ncpr, 'ess_frac_net': ess,
                         'dRg_pct_rw': (float(w @ rg) / rg.mean() - 1) * 100,
                         'R_hub': float(np.abs(do).sum()), 'R_mag': float(np.abs(dd[far]).sum()),
                         'R_shape': float(np.abs(dw / dw.sum() - deg0 / deg0.sum()).sum()),
                         'deg0_site': float(deg0[i]), 'occ0_site': float(occ0[i]),
                         's1_krde': float(np.median(dd[kr]) - np.median(dd[de])) if kr.any() and de.any() else np.nan,
                         's1_rho_charge': float(spearmanr(dd[far], ql[far]).correlation) if far.sum() > 5 else np.nan,
                         's1_total': float(dd.sum())})
            dprof.append(dd)
            oprof.append(do)
        np.savez_compressed(os.path.join(OUT, 'profiles', f'{uid}.npz'), seq=seq, positions=np.array(positions),
                            deg0=deg0, occ0=occ0, ddeg=np.array(dprof), docc=np.array(oprof))
        return rows, None
    except Exception as e:  # noqa: BLE001
        return None, f'{uid}: {e!r}'


# ------------------------------------------------------------------ PAGE4 validation (S4)
def page4_phi(xyz, seq):
    q, lam, sig, _ = sp.fdt.calvados_residue_params(seq, sp.L.load_residue_params())
    i, j, *_ = sp.fdt.pair_setup(q, lam, sig)
    M = sp.site_matrix(i, j, q, len(seq))
    k, eps, rc = sp.CAL['kappa'], sp.CAL['eps_yu'], sp.fdt.RC_YU
    ec = np.exp(-k * rc)
    xyz = np.asarray(xyz, float)
    phi = np.zeros((len(xyz), len(seq)))
    step = max(1, sp.MAX_ELEMS // max(1, len(i)))
    for a in range(0, len(xyz), step):
        r = np.linalg.norm(xyz[a:a + step, i] - xyz[a:a + step, j], axis=2)
        Y = np.where(r < rc, np.exp(-k * r) / r - ec / rc, 0.0)
        phi[a:a + step] = (M @ (eps * Y).T).T
    return phi, (k, eps, rc, ec)


def page4():
    base = os.path.join(MT, 'page4_CALVADOS', '300K')
    forms = {'HIPK1T51': [51], 'CLK2full': [7, 9, 51, 85, 73, 94, 71, 79]}

    def load(run):
        d = os.path.join(base, f'PAGE4-{run}-I150-pH75')
        xyz = sp.bc.load_calvados(os.path.join(d, f'PAGE4-{run}-I150-pH75.dcd'), os.path.join(d, 'top.pdb'),
                                  k=sp.L.K_FRAMES)[0]
        seq = ''.join(l.strip() for l in open(os.path.join(d, f'PAGE4-{run}-I150-pH75.fasta')) if not l.startswith('>'))
        return np.asarray(xyz, float), seq
    xw, sw = load('WT')
    n = len(sw)
    Dw, Hw = networks(xw, n, 'PAGE4-WT')
    phi, (k, eps, rc, ec) = page4_phi(xw, sw)
    rg = np.sqrt(((xw - xw.mean(1, keepdims=True)) ** 2).sum(2).mean(1))
    lines = ['S4  PAGE4: reweighted (unmodified direct run, charge model dq = -2, mutual terms) against direct phosphoform runs',
             f'    unmodified run: {len(xw)} frames, {n} residues', '']
    out = {'seq': sw, 'deg_wt': Dw.mean(0), 'occ_wt': Hw.mean(0)}
    for form, sites in forms.items():
        idx = [s - 1 for s in sites]
        assert all(sw[t] in 'ST' for t in idx), (form, [sw[t] for t in idx])
        dH = DQ * phi[:, idx].sum(1)
        for a in range(len(idx)):
            for b in range(a + 1, len(idx)):
                r = np.linalg.norm(xw[:, idx[a]] - xw[:, idx[b]], axis=1)
                dH += DQ * DQ * eps * np.where(r < rc, np.exp(-k * r) / r - ec / rc, 0.0)
        x = -sp.BETA_C * dH
        w = np.exp(x - x.max())
        w /= w.sum()
        ess = 1.0 / (w ** 2).sum() / len(w)
        xf, sf = load(form)
        Df, Hf = networks(xf, len(sf), f'PAGE4-{form}')
        rgf = np.sqrt(((xf - xf.mean(1, keepdims=True)) ** 2).sum(2).mean(1))
        rw_d, rw_o = w @ Dw - Dw.mean(0), w @ Hw.astype(float) - Hw.mean(0)
        di_d, di_o = Df.mean(0) - Dw.mean(0), Hf.mean(0) - Hw.mean(0)
        lines.append(f'  {form:9s} sites {sites}: effective sample fraction {ess:.3f}')
        lines.append(f'     Rg change: reweighted {(w @ rg / rg.mean() - 1) * 100:+.1f}%   direct {(rgf.mean() / rg.mean() - 1) * 100:+.1f}%')
        lines.append(f'     degree profile change, reweighted vs direct: Spearman {spearmanr(rw_d, di_d).correlation:+.3f}')
        lines.append(f'     hub occupancy change, reweighted vs direct:  Spearman {spearmanr(rw_o, di_o).correlation:+.3f}')
        out[f'{form}_rw_ddeg'], out[f'{form}_direct_ddeg'] = rw_d, di_d
        out[f'{form}_rw_docc'], out[f'{form}_direct_docc'] = rw_o, di_o
        out[f'{form}_sites'] = np.array(sites)
    np.savez_compressed(os.path.join(OUT, 'page4_profiles.npz'), **out)
    open(os.path.join(OUT, 'page4_validation.txt'), 'w').write('\n'.join(lines) + '\n')
    print('\n'.join(lines), flush=True)


def main():
    os.makedirs(os.path.join(OUT, 'profiles'), exist_ok=True)
    S = pd.read_csv(SITES)
    S = S[S.model == 'CALVADOS-2']
    kingdom = pd.read_csv(sp.bc.BIO_CSV, usecols=['UniProt_ID', 'kingdom']).set_index('UniProt_ID')['kingdom']
    tasks = [(u, kingdom.get(u), sorted(int(p) for p in g.pos)) for u, g in S.groupby('UniProt_ID')]
    print(f'{len(tasks)} proteins, {len(S)} Ser/Thr sites', flush=True)
    rows, fails = [], []
    with mp.Pool(8) as pool:
        for k, (r, err) in enumerate(pool.imap_unordered(one, tasks), 1):
            if err:
                fails.append(err)
            else:
                rows += r
            if k % 10 == 0:
                print(f'  {k}/{len(tasks)} proteins, {len(fails)} failed', flush=True)
    pd.DataFrame(rows).to_csv(os.path.join(OUT, 'sites_network.csv'), index=False)
    open(os.path.join(OUT, 'failed.txt'), 'w').write('\n'.join(fails) + '\n')
    print(f'{len(rows)} sites written, {len(fails)} proteins failed', flush=True)
    page4()


if __name__ == '__main__':
    main()
