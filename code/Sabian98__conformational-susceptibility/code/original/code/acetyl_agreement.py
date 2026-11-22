"""
acetyl_agreement.py — why CALVADOS-2 and Mpipi-Recharged agree on the change of Rg on acetylating the same Lys (rank
correlation 0.90) but not on the change of the total response (0.26). Exploratory, after the fact.

Same reweighting as acetyl_heat.py, with the change of the total response split exactly into three parts:
  shift      -Cov_w(Rg, G_tot)/<Rg>_w + Cov(Rg, G_tot)/<Rg>      (the ensemble moves, thermal force unchanged)
  boltzmann  -Cov_w(Rg, dG_B)/<Rg>_w,  dG_B = -dH/(R T^2)          (the Boltzmann part of the force change)
  electro    -Cov_w(Rg, dG - dG_B)/<Rg>_w                          (temperature dependence of the electrostatics:
             dielectric strength and screening in CALVADOS-2, screening only in Mpipi-Recharged)
shift + boltzmann + electro = the stored change (checked). Chains and charge model as acetyl_heat.py.
    python acetyl_agreement.py --model C|R --out DIR [--workers N]
Output: DIR/parts_{C,R}.csv (append; resumable)
"""
import sys

sys.dont_write_bytecode = True

import argparse  # noqa: E402
import multiprocessing as mp  # noqa: E402
import os  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

HERE = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, os.environ.get('BENDER_REANALYSIS_DIR', HERE))
sys.path.insert(0, HERE)
import acetyl_heat as AH  # noqa: E402

HC = AH.HC
T0 = 300.0


def sites_parts(uid, seq, rg, Gtot, phi, dG_of, RT, Gfix=None):
    nf, m = len(rg), rg.mean()
    base = np.cov(rg, Gtot, ddof=0)[0, 1] / m

    def resp(w, G):
        return -HC.wcov(rg, G, w) / ((w * rg).sum() / w.sum())
    rows = []
    for p_, c in enumerate(seq):
        if c != 'K':
            continue
        dH, dG = dG_of(p_)
        dGB = -dH / (RT * T0)
        x = -dH / RT
        w = np.exp(x - x.max())
        shift = (resp(w, Gtot) + base) * 100
        boltz = -HC.wcov(rg, dGB, w) / ((w * rg).sum() / w.sum()) * 100
        elec = -HC.wcov(rg, dG - dGB, w) / ((w * rg).sum() / w.sum()) * 100
        rows.append({'UniProt_ID': uid, 'pos': p_ + 1, 'dtotal': shift + boltz + elec, 'shift': shift, 'boltzmann': boltz,
                     'electro': elec, 'ess': float(w.sum() ** 2 / (w ** 2).sum() / nf),
                     'dRg_pct': float(((w * rg).sum() / w.sum() / m - 1) * 100)})
    row = {'UniProt_ID': uid, 'N': len(seq), 'total_rel': -np.cov(rg, Gtot)[0, 1] / m * 100}
    return row, pd.DataFrame(rows)


def work_C(args):
    AH.sites_from = sites_parts
    return AH.one_calvados(args)


def work_R(args):
    AH.sites_from = sites_parts
    return AH.one_recharged(args)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--model', choices=['C', 'R'], required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--workers', type=int, default=6)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    ids = sorted(pd.read_csv(AH.IDS[a.model]).UniProt_ID)
    ps = os.path.join(a.out, f'parts_{a.model}.csv')
    done = set(pd.read_csv(ps, usecols=['UniProt_ID']).UniProt_ID) if os.path.exists(ps) else set()
    if a.model == 'C':
        todo, fn = list(HC.TF.index_dirs([u for u in ids if u not in done]).items()), work_C
    else:
        import tea_crossff as TC
        dirs = TC.rc_dirs()
        todo, fn = [(u, dirs[u]) for u in ids if u not in done and u in dirs], work_R
    print(f'{a.model}: {len(ids)} chains, {len(done)} done, {len(todo)} to run', flush=True)
    nfail = 0
    with mp.get_context('spawn').Pool(a.workers) as pool:
        for k, (row, S, err) in enumerate(pool.imap_unordered(fn, todo), 1):
            if err:
                nfail += 1
                print('failed', err[:300], flush=True)
                continue
            if S is not None and len(S):
                S.assign(total_rel=row['total_rel']).to_csv(ps, mode='a', header=not os.path.exists(ps), index=False)
            if k % 100 == 0:
                print(f'{k}/{len(todo)} done, {nfail} failed', flush=True)
    print('finished', flush=True)


if __name__ == '__main__':
    main()
