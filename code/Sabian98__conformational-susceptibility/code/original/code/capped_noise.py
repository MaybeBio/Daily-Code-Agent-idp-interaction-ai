"""
capped_noise.py — sampling noise of the total heating response under the capped forcing (the main model from 2026-10-24),
for the noise ceilings of total_response_evidence.py (section F), as total_response_noise.py does for the linear A0 forcing.

Same 400-chain subsets per model as total_response_noise.py (the IDs of total_response/noise/total_response_se_<k>.csv).
Per frame the capped forcings D_chg_cap and D_tea_cap come from capped_proteome.py (its finish() intercepted), G as in
total_response_noise.py (CALVADOS-2 and Mpipi-Recharged npz 'G_T', Mpipi-GG -PotEng/(R T^2)); strengths from
total_response_evidence.load() with BENDER_FORCING=capped (capped_model.py fit, Mpipi-GG by the SD-ratio rule). Then
    total = -Cov(Rg, G + (g_h D_chg_cap + g_tea D_tea_cap) / RT)
with the same 10-block delete-one jackknife SE. Check (the script stops otherwise): the capped forcing terms must equal the
loader's to 1e-8 relative, the total to 1e-4 (the stored fixed-solvent chi reproduces to about 1.5e-5 relative).
    BENDER_REANALYSIS_DIR=... python capped_noise.py [--workers 6]      (openmpipi env)
Output: Traj Level Analysis/Analysis_FDT/total_response/noise_capped/total_response_se_<model>.csv (relative values)
"""
import sys

sys.dont_write_bytecode = True

import argparse  # noqa: E402
import multiprocessing as mp  # noqa: E402
import os  # noqa: E402
import traceback  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

os.environ['BENDER_FORCING'] = 'capped'
HERE = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, os.environ.get('BENDER_REANALYSIS_DIR', HERE))
sys.path.insert(0, HERE)
import capped_proteome as CP  # noqa: E402
import total_response_noise as TN  # noqa: E402

KEY = {'C': 'calvados', 'M': 'gg', 'R': 'recharged'}
NOISE_A0 = os.path.join(CP.TF.L.BENDER_MT, 'Traj Level Analysis', 'Analysis_FDT', 'total_response', 'noise')
OUT = os.path.join(os.path.dirname(NOISE_A0), 'noise_capped')

_CAUGHT = {}


def _catch(uid, seq, obs, D, c, t0, Rgas):
    _CAUGHT[uid] = (np.asarray(obs['rg'], float), D['chg_cap'], D['tea_cap'], Rgas * CP.T0)
    return {'UniProt_ID': uid}


CP.finish = _catch   # module level, so spawned workers intercept too


def one(task):
    uid, d, model, g_h, g_tea = task
    try:
        row, err = CP.RUNNERS[KEY[model]]((uid, d))
        if err:
            raise RuntimeError(err)
        rg, D, Dtea, RT = _CAUGHT.pop(uid)
        if model == 'M':
            G = -np.load(os.path.join(CP.X.GG_NPZ, f'{uid}.npz'))['PotEng_kcal'].astype(float) / (CP.X.R_KCAL * CP.T0 ** 2)
        else:
            G = np.load(os.path.join(CP.TF.NPZ if model == 'C' else CP.X.RC_NPZ, f'{uid}.npz'))['G_T'].astype(float)
        Gt = G + (g_h * D + g_tea * Dtea) / RT

        def cv(x):
            return lambda idx: -float(np.cov(rg[idx], x[idx], ddof=1)[0, 1])

        n, m = len(rg), rg.mean()
        return {'UniProt_ID': uid, 'model': model, 'gamma_hyd': g_h, 'gamma_tea': g_tea, 'n_frames': n,
                'chi_rel': 100 * cv(G)(np.arange(n)) / m, 'chi_rel_se': 100 * TN.jk_se(cv(G), n) / m,
                'hyd_rel': 100 * cv(D / RT)(np.arange(n)) / m, 'hyd_rel_se': 100 * TN.jk_se(cv(D / RT), n) / m,
                'tea_rel': 100 * cv(Dtea / RT)(np.arange(n)) / m,
                'total_rel': 100 * cv(Gt)(np.arange(n)) / m, 'total_rel_se': 100 * TN.jk_se(cv(Gt), n) / m}, None
    except Exception as ex:
        return None, f'{uid}: {ex!r} | {traceback.format_exc()[-300:]}'.replace('\n', ' ')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--workers', type=int, default=6)
    ap.add_argument('--models', default='CRM')
    a = ap.parse_args()
    import total_response_evidence as E
    assert E.FORCING == 'capped'
    Dl, g = E.load()
    os.makedirs(OUT, exist_ok=True)
    for model in a.models:
        ids = pd.read_csv(os.path.join(NOISE_A0, f'total_response_se_{model}.csv')).UniProt_ID.tolist()
        dirs = CP.dirs(KEY[model])
        tasks = [(u, dirs[u], model, g[model], g[f'tea_{model}']) for u in ids if u in dirs]
        rows, fails = [], []
        with mp.get_context('spawn').Pool(a.workers) as pool:
            for r, err in pool.imap_unordered(one, tasks):
                (fails.append(err) if err else rows.append(r))
        t = pd.DataFrame(rows).set_index('UniProt_ID').loc[[u for u in ids if u in {r['UniProt_ID'] for r in rows}]]
        ref = Dl[model].reindex(t.index).dropna(subset=['total'])     # Mpipi-GG subset: 62 chains lie outside the
        t_ = t.loc[ref.index]                                          # loader's main set and are kept but not checked
        forc = g[model] * t_.hyd_rel + g[f'tea_{model}'] * t_.tea_rel
        dev_f = ((forc - (ref.H + ref.TEA)).abs() / (ref.H + ref.TEA).abs().clip(lower=1e-6)).max()
        scale = np.maximum(ref.total.abs(), ref.chi.abs()).clip(lower=1e-6)   # the stored chi's error scales with chi
        dev = ((t_.total_rel - ref.total).abs() / scale).max()
        if t_[['total_rel', 'total_rel_se', 'chi_rel_se']].isna().any().any():
            raise SystemExit(f'{model}: non-finite values in the noise table')
        print(f'{model}: {len(t)} of {len(ids)} chains ({len(t_)} in the loader set, checked), {len(fails)} failures, '
              f'max relative deviation from loader: '
              f'forcing terms {dev_f:.2e}, total {dev:.2e}; median SE/|total| '
              f'{(t.total_rel_se / t.total_rel.abs()).median():.3f}', flush=True)
        for f in fails[:5]:
            print('  ', f)
        # the capped terms must be exact; the total carries the stored chi's precision (about 1.5e-5 relative,
        # total_response_noise.py docstring)
        if dev_f > 1e-8 or dev > 1e-4 or len(t_) < 0.8 * len(t):
            raise SystemExit(f'{model}: totals do not reproduce the loader')
        t.reset_index().to_csv(os.path.join(OUT, f'total_response_se_{model}.csv'), index=False)


if __name__ == '__main__':
    main()
