"""
recharged_traj_pass.py — Mpipi-Recharged (300 K, 150 mM) counterparts of the CALVADOS-2 / Mpipi-GG trajectory analyses,
one read of the stored 1,000 frames per chain. Energies as mutation_scan_recharged.py (validated against OpenMpipi to
~1e-6); thermal force G_T from the npz; total = -Cov(Rg, G_T + D/RT) with the calibrated solvation forcing (TEA 0.4,
hydrophilic 0.209).
  P  pair types (slide 23, pair_type_split.py types): share of Cov(Rg, E_WF) per type, share of pairs, per-pair strength
  H  thermal hotspots: per residue (each pair half to each partner) the Boltzmann share Cov(Rg, H_r)/(R T^2) with H_r the
     Wang-Frenkel + Yukawa energy, the screening share -Cov(Rg, dH_r/dT)/RT, the solvation share -Cov(Rg, D_r)/RT
  S  heat x salt: ionic strength doubled (kappa and dkappa/dT x sqrt 2); exact reweighting; change of the total
  F  phosphomimetics: every Ser/Thr -> Glu and -> Asp (Mpipi-Recharged has no phosphoserine); change of the total and of
     <Rg> (the leverage counterpart of slides 34-35)
  B  basins (basin_susceptibility.Landscape): number, both-halves rule, between-basin share of the plain and total response
Sample: heat_crosses' chains, the slide-23 sample (pair_type_split) and the common set (basin_susceptibility).
    python recharged_traj_pass.py --out DIR [--workers N] [--limit N]
Output: DIR/{chains,pairs,residues,phosphomimetic}.csv (append; resumable)
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
import mutation_scan_recharged as MS  # noqa: E402

TC, TF = MS.TC, MS.TF
T0, R_KJ, RT = MS.T0, MS.R_KJ, MS.RT
TYPES = ['aro-aro', 'R-aro', 'K-aro', 'R-other', 'K-other', 'aro-other', 'other-other']
ARO = set('FWY')
TLA = HC.TLA


def pair_types(seq, i, j):
    a, b = np.array(list(seq))[i], np.array(list(seq))[j]
    ia, ib = np.isin(a, list(ARO)), np.isin(b, list(ARO))
    ra, rb, ka, kb = a == 'R', b == 'R', a == 'K', b == 'K'
    t = np.full(len(i), 'other-other', dtype=object)
    t[ia | ib] = 'aro-other'
    t[ka | kb] = 'K-other'
    t[ra | rb] = 'R-other'
    t[(ka & ib) | (kb & ia)] = 'K-aro'
    t[(ra & ib) | (rb & ia)] = 'R-aro'
    t[ia & ib] = 'aro-aro'
    return t


def per_pair(r, ia, ib, tb, k, dk, fa, fb):
    """Per-pair (frames x pairs) WF energy, Yukawa energy, Yukawa dH/dT, solvation forcing (kJ/mol, kJ/mol/K)."""
    eps, sig, mu = tb['eps'][ib, ia], tb['sig'][ib, ia], tb['mu'][ib, ia]
    rc = 3 * sig
    alpha = 2 * 3 ** (2 * mu) * (3 / (2 * (3 ** (2 * mu) - 1))) ** 3
    shape = np.where(r < rc, alpha * ((sig / r) ** (2 * mu) - 1) * ((rc / r) ** (2 * mu) - 1) ** 2, 0.0)
    ex = np.where(r < 3.5, np.exp(-k * r), 0.0)
    A = tb['A'][ib, ia]
    return shape * eps, A * ex / r, -dk * A * ex, TC.wf_attractive(r, sig, mu, 1, rc, alpha) * 0.5 * (fa + fb) * 4.184


def one(args):
    uid, d, want_basins = args
    t0 = time.time()
    try:
        import mpipi_recharged_fdt as R
        z = np.load(os.path.join(TC.RC_NPZ, f'{uid}.npz'))
        xyz, seq, info = R.load_run(d, uid)
        fi = info['production_start'] + info['stride'] * np.arange(len(xyz))
        if not np.array_equal(fi, z['frame_index']):
            raise ValueError('frames differ from the stored FDT frames')
        tb = MS.tables()
        k = MS.kappa_of(T0)
        dk = MS.kappa_of(T0 + 0.5) - MS.kappa_of(T0 - 0.5)
        k2, dk2 = k * np.sqrt(2.0), dk * np.sqrt(2.0)
        rg, G = z['rg'].astype(float), z['G_T'].astype(float)
        N, nf, L = len(seq), len(xyz), R.box_nm(len(seq))
        ix = np.array([tb['ix'][a] for a in seq])
        ft = {a: MS.G_TEA * TF.DL['charmm36'][a] + MS.G_HYD * TF.DL_CLASSES['hyd'][a] for a in MS.AA}
        fv = np.array([ft[a] for a in seq])
        i, j = np.triu_indices(N, k=2)
        npair = len(i)
        typ = pair_types(seq, i, j)
        half = sparse.csr_matrix((np.full(2 * npair, 0.5), (np.concatenate([i, j]), np.concatenate([np.arange(npair)] * 2))), shape=(N, npair))
        Et = {t: np.zeros(nf) for t in TYPES}
        Hres, dHres, Dres = np.zeros((nf, N)), np.zeros((nf, N)), np.zeros((nf, N))
        W, Y1, dY1, D, Y2, dY2 = (np.zeros(nf) for _ in range(6))
        step = max(1, int(TC.MAX_ELEMS // max(1, npair)))
        for a in range(0, nf, step):
            dv = xyz[a:a + step, i].astype(float) - xyz[a:a + step, j].astype(float)
            dv -= L * np.round(dv / L)
            r = np.linalg.norm(dv, axis=2)
            w_, y_, dy_, d_ = per_pair(r, ix[i], ix[j], tb, k, dk, fv[i], fv[j])
            _, y2_, dy2_, _ = per_pair(r, ix[i], ix[j], tb, k2, dk2, fv[i], fv[j])
            W[a:a + step], Y1[a:a + step], dY1[a:a + step], D[a:a + step] = w_.sum(1), y_.sum(1), dy_.sum(1), d_.sum(1)
            Y2[a:a + step], dY2[a:a + step] = y2_.sum(1), dy2_.sum(1)
            for t in TYPES:
                Et[t][a:a + step] = w_[:, typ == t].sum(1)
            Hres[a:a + step] = (half @ (w_ + y_).T).T
            dHres[a:a + step] = (half @ dy_.T).T
            Dres[a:a + step] = (half @ d_.T).T
        if np.abs(W - z['V_AH']).max() > 1e-3 * max(1.0, np.abs(z['V_AH']).max()):
            raise ValueError('WF energy differs from the npz')
        Gtot = G + D / RT
        m = rg.mean()
        cv = lambda x: float(((rg - m) * (x - x.mean())).sum() / (nf - 1))  # noqa: E731
        row = {'UniProt_ID': uid, 'N': N, 'chi_rel': -cv(G) / m * 100, 'total_rel': -cv(Gtot) / m * 100,
               'chi_WF_rel': cv(W) / (RT * T0) / m * 100}
        # P
        cs = sum(cv(Et[t]) for t in TYPES)
        pairs = {'UniProt_ID': uid}
        for t in TYPES:
            nt = int((typ == t).sum())
            pairs[f'share_chi_{t}'] = cv(Et[t]) / cs if cs else np.nan
            pairs[f'share_pairs_{t}'] = nt / npair
            pairs[f'n_pairs_{t}'] = nt
            pairs[f'chi_per_pair_{t}'] = (cv(Et[t]) / nt) / (cs / npair) if nt and cs else np.nan
        # H
        res = pd.DataFrame({'UniProt_ID': uid, 'pos': np.arange(1, N + 1), 'residue': list(seq),
                            'boltz_share': [cv(Hres[:, q]) / (RT * T0) / m * 100 for q in range(N)],
                            'screen_share': [-cv(dHres[:, q]) / RT / m * 100 for q in range(N)],
                            'solv_share': [-cv(Dres[:, q]) / RT / m * 100 for q in range(N)]})
        res['total_share'] = res[['boltz_share', 'screen_share', 'solv_share']].sum(1)

        def perturbed(dH, dG):
            x = -dH / RT
            w = np.exp(x - x.max())
            f = lambda ixs: (-HC.wcov(rg[ixs], Gtot[ixs] + dG[ixs], w[ixs]) / ((w[ixs] * rg[ixs]).sum() / w[ixs].sum())  # noqa: E731
                             + HC.fdt._cov(rg[ixs], Gtot[ixs]) * (len(ixs) - 1) / len(ixs) / rg[ixs].mean()) * 100
            return f(np.arange(nf)), HC.jk(f, nf), float(w.sum() ** 2 / (w ** 2).sum() / nf), float(((w * rg).sum() / w.sum() / m - 1) * 100)
        # S
        dHs = Y2 - Y1
        row['salt_dtotal'], row['salt_dtotal_se'], row['salt_ess'], row['salt_dRg_pct'] = perturbed(
            dHs, -dHs / (R_KJ * T0 ** 2) + (dY2 - dY1) / RT)
        # F
        fm = []
        for p in range(N):
            if seq[p] not in 'ST':
                continue
            other = np.array([q for q in range(N) if abs(q - p) >= 2])
            dv = xyz[:, [p]].astype(float) - xyz[:, other].astype(float)
            dv -= L * np.round(dv / L)
            r = np.linalg.norm(dv, axis=2)
            base = [x.sum(1) for x in per_pair(r, np.full(len(other), ix[p]), ix[other], tb, k, dk, np.full(len(other), fv[p]), fv[other])]
            for b in 'ED':
                new = [x.sum(1) for x in per_pair(r, np.full(len(other), tb['ix'][b]), ix[other], tb, k, dk, np.full(len(other), ft[b]), fv[other])]
                dH = (new[0] - base[0]) + (new[1] - base[1])
                dG = -dH / (R_KJ * T0 ** 2) + (new[2] - base[2]) / RT + (new[3] - base[3]) / RT
                dt, se, ess, drg = perturbed(dH, dG)
                fm.append({'UniProt_ID': uid, 'pos': p + 1, 'residue': seq[p], 'mimic': b, 'dtotal': dt, 'dtotal_se': se,
                           'ess': ess, 'dRg_pct': drg})
        # B
        if want_basins:
            import basin_susceptibility as BS
            ls = BS.Landscape(BS.seg_features(xyz))
            Yp = ls.project(BS.seg_features(xyz))
            lab, _, _ = ls.basins(Yp)
            kb = int(lab.max()) + 1
            h = len(Yp) // 2
            k1, k2b = int(ls.basins(Yp[:h])[0].max()) + 1, int(ls.basins(Yp[h:])[0].max()) + 1
            row.update({'n_basins': kb, 'multi_basin': bool(kb > 1 and k1 > 1 and k2b > 1)})
            if kb > 1:
                for name, force in (('plain', G), ('total', Gtot)):
                    tot, btw, _, _ = BS.split(rg, force, lab, kb)
                    row[f'between_share_{name}'] = btw / tot if tot else np.nan
        row['seconds'] = round(time.time() - t0, 1)
        return row, pairs, res, pd.DataFrame(fm), None
    except Exception as ex:  # noqa: BLE001
        return {'UniProt_ID': uid}, None, None, None, f'{uid}: {ex!r} | {traceback.format_exc()[-300:]}'.replace('\n', ' ')


def sample():
    hx = pd.read_csv(os.path.join(TLA, 'Analysis_FDT', 'total_response', 'heat_crosses', 'sample.csv')).UniProt_ID
    pt = pd.read_csv(os.path.join(TLA, 'Analysis_FDT_MPIPI', 'compact', 'pair_type_split', 'per_protein.csv')).UniProt_ID
    per = pd.read_csv(os.path.join(TLA, 'Analysis_FDT', 'energy_mechanism_split', 'per_protein.csv'), index_col=0).dropna(subset=['chi_C', 'chi_M'])
    basins = set(per.index)
    return sorted(set(hx) | set(pt) | basins), basins


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', required=True)
    ap.add_argument('--workers', type=int, default=6)
    ap.add_argument('--limit', type=int, default=None)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    ids, basins = sample()
    pc = os.path.join(a.out, 'chains.csv')
    done = set(pd.read_csv(pc, usecols=['UniProt_ID']).UniProt_ID) if os.path.exists(pc) else set()
    dirs = TC.rc_dirs()
    todo = [(u, dirs[u], u in basins) for u in ids if u not in done and u in dirs and u != 'Q983T3'][:a.limit]
    print(f'{len(ids)} chains, {len(done)} done, {len(todo)} to run', flush=True)
    nf = 0
    with mp.get_context('spawn').Pool(a.workers) as pool:
        for kk, (row, pairs, res, fm, err) in enumerate(pool.imap_unordered(one, todo), 1):
            if err:
                nf += 1
                print('failed', err[:300], flush=True)
                continue
            cols = ['UniProt_ID', 'N', 'chi_rel', 'total_rel', 'chi_WF_rel', 'salt_dtotal', 'salt_dtotal_se', 'salt_ess',
                    'salt_dRg_pct', 'n_basins', 'multi_basin', 'between_share_plain', 'between_share_total', 'seconds']
            # fixed columns (2026-10-10 fix); pairs/residues restored 2026-10-10 after the comment had swallowed them
            for df, name in ((pd.DataFrame([row]).reindex(columns=cols), 'chains.csv'), (pd.DataFrame([pairs]), 'pairs.csv'),
                             (res, 'residues.csv'), (fm, 'phosphomimetic.csv')):
                p = os.path.join(a.out, name)
                if df is not None and len(df):
                    df.to_csv(p, mode='a', header=not os.path.exists(p), index=False)
            if kk % 50 == 0:
                print(f'{kk}/{len(todo)} done, {nf} failed', flush=True)
    print('finished', flush=True)


if __name__ == '__main__':
    main()
