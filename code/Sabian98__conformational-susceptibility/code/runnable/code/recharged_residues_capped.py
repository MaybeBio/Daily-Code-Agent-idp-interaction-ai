"""
recharged_residues_capped.py — per-residue thermal hotspots of Mpipi-Recharged (recharged_traj_pass.py part H) under the
calibrated model with saturating forcing (from 2026-10-25), for the random chains of heat_crosses' sample.

The Boltzmann and screening shares of each residue contain no forcing and are taken from the stored pass
(Analysis_FDT_RECHARGED/traj_pass/residues.csv). The solvation share becomes each residue's own term of the saturating
forcing, -Cov(Rg, D_i)/RT / <Rg> x 100 with
    D_i = -1/2 x 4.184 x (g_tea f_tea(a_i) + g_h f_chg(a_i)) h(c_i)      (kJ/mol/K; capped_proteome.py, capped_frames strengths)
and c_i the per-frame contact numbers of capped_proteome.recharged_one (its finish() intercepted). Check per chain: the
residue solvation shares add up to the chain's capped forcing term in capped_dchi_recharged_300K.csv to 1e-8 relative.
    BENDER_REANALYSIS_DIR=... python recharged_residues_capped.py [--workers 6]      (openmpipi env)
Output: Analysis_FDT_RECHARGED/traj_pass_capped/{residues.csv, analysis_residue_types.csv, summary.txt}
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
import capped_frames as CF  # noqa: E402
import capped_proteome as CP  # noqa: E402

TLA = os.path.join(CP.TF.L.BENDER_MT, 'Traj Level Analysis')
SRC = os.path.join(TLA, 'Analysis_FDT_RECHARGED', 'traj_pass')
OUT = os.path.join(TLA, 'Analysis_FDT_RECHARGED', 'traj_pass_capped')
HX = os.path.join(TLA, 'Analysis_FDT', 'total_response', 'heat_crosses')
_C = {}


def _catch(uid, seq, obs, D, c, t0, Rgas):
    _C[uid] = (np.asarray(obs['rg'], float), np.asarray(c, float), seq, Rgas * CP.T0)
    return {'UniProt_ID': uid}


CP.finish = _catch


def one(task):
    uid, d = task
    try:
        _, err = CP.recharged_one((uid, d))
        if err:
            raise RuntimeError(err)
        rg, c, seq, RT = _C.pop(uid)
        gh, gt = CF.strengths()['R']
        f = 4.184 * np.array([gt * CP.F_TEA[a] + gh * CP.F_CHG[a] for a in seq])
        Di = -0.5 * CP.h(c) * f                                   # frames x residues
        rc = rg - rg.mean()
        share = -100 * (rc @ (Di - Di.mean(0))) / (len(rg) - 1) / RT / rg.mean()
        return pd.DataFrame({'UniProt_ID': uid, 'pos': np.arange(1, len(seq) + 1), 'residue': list(seq),
                             'solv_share_cap': share}), None
    except Exception as ex:  # noqa: BLE001
        return None, f'{uid}: {ex!r} | {traceback.format_exc()[-300:]}'.replace('\n', ' ')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--workers', type=int, default=6)
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    smp = pd.read_csv(os.path.join(HX, 'sample.csv'))
    rnd = sorted(set(smp.loc[smp.random.astype(bool), 'UniProt_ID']))
    R = pd.read_csv(os.path.join(SRC, 'residues.csv')).drop_duplicates(['UniProt_ID', 'pos'])
    R = R[R.UniProt_ID.isin(rnd)]
    dirs = CP.dirs('recharged')
    tasks = [(u, dirs[u]) for u in sorted(set(R.UniProt_ID)) if u in dirs]
    rows, fails = [], []
    with mp.get_context('spawn').Pool(a.workers) as pool:
        for df, err in pool.imap_unordered(one, tasks):
            (fails.append(err) if err else rows.append(df))
    S = pd.concat(rows)
    M = R.merge(S, on=['UniProt_ID', 'pos', 'residue'], how='inner')
    cap = pd.read_csv(os.path.join(CP.OUT, 'capped_dchi_recharged_300K.csv')).set_index('UniProt_ID')
    gh, gt = CF.strengths()['R']
    term = (gt * cap.dchi_rg_rel_tea_cap + gh * cap.dchi_rg_rel_chg_cap)
    sums = M.groupby('UniProt_ID').solv_share_cap.sum()
    dev = ((sums - term.reindex(sums.index)).abs() / term.reindex(sums.index).abs()).max()
    if not dev < 1e-8:
        raise SystemExit(f'residue solvation shares do not add up to the chain terms ({dev:.1e})')
    M['solv_share'] = M.solv_share_cap
    M['total_share'] = M.boltz_share + M.screen_share + M.solv_share
    M = M.drop(columns='solv_share_cap')
    M.to_csv(os.path.join(OUT, 'residues.csv'), index=False)
    by = M.groupby('residue')[['solv_share', 'boltz_share', 'screen_share', 'total_share']].median().sort_values('total_share')
    by.to_csv(os.path.join(OUT, 'analysis_residue_types.csv'))
    Rc = pd.read_csv(os.path.join(HX, 'residues.csv'))
    byc = Rc[Rc.UniProt_ID.isin(rnd)].groupby('residue').total_share.median()
    from scipy import stats
    L = [f'{M.UniProt_ID.nunique()} random chains, {len(M):,} residues, {len(fails)} failures; residue solvation shares add up '
         f'to the chain capped forcing terms to {dev:.1e}', by.round(4).to_string(),
         f'residue-type ranking vs CALVADOS-2: Spearman {stats.spearmanr(by.total_share, byc.reindex(by.index))[0]:.2f}']
    open(os.path.join(OUT, 'summary.txt'), 'w').write('\n'.join(L) + '\n')
    print('\n'.join(L))
    for f_ in fails[:5]:
        print(f_)


if __name__ == '__main__':
    main()
