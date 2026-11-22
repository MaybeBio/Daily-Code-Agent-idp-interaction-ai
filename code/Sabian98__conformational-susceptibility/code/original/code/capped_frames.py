"""
capped_frames.py — per-frame capped hydration forcing of a chain at 300 K (the main model from 2026-10-24), for the
reweighting pipelines whose perturbation leaves every residue's forcing coefficient unchanged (phosphate charge, Lys
acetylation, ionic strength): there only the baseline thermal force changes,
    G_tot = G + (g_h D_chg_cap + g_tea D_tea_cap) / RT
with D per unit strength from capped_proteome.py (its finish() intercepted) and the strengths of
total_response_evidence.load() under the capped forcing (CALVADOS-2 and Mpipi-Recharged: capped_model.py; Mpipi-GG: the
SD-ratio rule), stored once in capped_300K/strengths.csv by `python capped_frames.py --strengths`.
Per-frame forcings are cached in capped_300K/frames/<model>/<uid>.npz (D_chg_cap, D_tea_cap, rg, RT). term() checks
that the caller's Rg series equals the cached one (same frames, same order) before returning.
    from capped_frames import term;  Gtot = G + term('C' | 'M' | 'R', uid, rg)
"""
import sys

sys.dont_write_bytecode = True

import os  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

HERE = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, os.environ.get('BENDER_REANALYSIS_DIR', HERE))
sys.path.insert(0, HERE)
import capped_proteome as CP  # noqa: E402

KEY = {'C': 'calvados', 'M': 'gg', 'R': 'recharged'}
STRENGTHS = os.path.join(CP.OUT, 'strengths.csv')
CACHE = os.path.join(CP.OUT, 'frames')
_CAUGHT = {}
_G = None
_DIRS = {}


def _catch(uid, seq, obs, D, c, t0, Rgas):
    _CAUGHT[uid] = {'D_chg_cap': D['chg_cap'], 'D_tea_cap': D['tea_cap'], 'rg': np.asarray(obs['rg'], float),
                    'RT': Rgas * CP.T0}
    return {'UniProt_ID': uid}


def strengths():
    """Strengths per model; CALVADOS-2 and Mpipi-Recharged must equal capped_model.params() (the source the scan
    scripts use), so a refit without rerunning --strengths stops here instead of mixing two sets."""
    global _G
    if _G is None:
        import capped_model as CM
        t = pd.read_csv(STRENGTHS).set_index('model')
        _G = {k: (float(t.loc[k, 'gamma_h']), float(t.loc[k, 'gamma_tea'])) for k in KEY}
        gh, gt = CM.params()
        for k in ('C', 'R'):
            if abs(_G[k][0] - gh[k][0]) > 1e-9 or abs(_G[k][1] - gt[k]) > 1e-9:
                raise RuntimeError(f'strengths.csv ({_G[k]}) differs from capped_model params for {k}; rerun --strengths')
    return _G


def params_key():
    """What the cached forcings depend on: c_sat and the per-residue coefficients of capped_proteome."""
    import hashlib
    import json
    blob = json.dumps({'csat': CP.CSAT, 'tea': CP.F_TEA, 'chg': CP.F_CHG}, sort_keys=True)
    return hashlib.sha256(blob.encode()).hexdigest()[:16]


_KEY_OK = False


def check_cache_key():
    """The cache folder records the parameters its files were made with; a mismatch stops instead of reusing stale
    forcings (move the folder aside to recompute)."""
    global _KEY_OK
    if _KEY_OK:
        return
    os.makedirs(CACHE, exist_ok=True)
    p = os.path.join(CACHE, 'params_key.txt')
    key = params_key()
    if os.path.exists(p):
        old = open(p).read().strip()
        if old != key:
            raise RuntimeError(f'capped forcing cache {CACHE} was made with other parameters ({old} vs {key})')
    else:
        open(p, 'w').write(key + '\n')
    _KEY_OK = True


def _dir(model, uid):
    if model == 'C':
        return CP.TF.index_dirs([uid])[uid]
    if model not in _DIRS:
        _DIRS[model] = CP.dirs(KEY[model])
    return _DIRS[model][uid]


def frames(model, uid, d=None):
    check_cache_key()
    p = os.path.join(CACHE, KEY[model], f'{uid}.npz')
    if os.path.exists(p):
        z = np.load(p)
        return {k: z[k] for k in z.files}
    orig = CP.finish
    CP.finish = _catch
    try:
        _, err = CP.RUNNERS[KEY[model]]((uid, d or _dir(model, uid)))    # d: the caller's run directory
    finally:
        CP.finish = orig
    if err:
        raise RuntimeError(err)
    out = _CAUGHT.pop(uid)
    os.makedirs(os.path.dirname(p), exist_ok=True)
    tmp = p[:-4] + f'.tmp{os.getpid()}.npz'
    np.savez(tmp, **out)
    try:
        os.replace(tmp, p)
    except OSError:                               # another worker wrote it first (SMB may refuse the rename): keep theirs
        os.remove(tmp)                            # (only this worker's own temporary copy)
        if not os.path.exists(p):
            raise
        z = np.load(p)
        return {k: z[k] for k in z.files}
    return out


def term(model, uid, rg, d=None):
    """(g_h D_chg_cap + g_tea D_tea_cap) / RT per frame, after checking the caller's Rg against the cached frames."""
    f = frames(model, uid, d)
    rg = np.asarray(rg, float)
    if (f['rg'].shape != rg.shape or not np.isfinite(rg).all() or not np.isfinite(f['rg']).all()
            or not np.allclose(f['rg'], rg, rtol=0, atol=1e-6 * rg.mean())):
        raise ValueError(f'{uid}: Rg series differs from the capped forcing frames (or is not finite)')
    gh, gt = strengths()[model]
    return (gh * f['D_chg_cap'] + gt * f['D_tea_cap']) / float(f['RT'])


def write_strengths():
    os.environ['BENDER_FORCING'] = 'capped'
    import capped_model as CM
    import total_response_evidence as E
    assert E.FORCING == 'capped'
    _, g = E.load()
    gh, gt = CM.params()
    for k in ('C', 'R'):
        assert abs(g[k] - gh[k][0]) < 1e-12 and abs(g[f'tea_{k}'] - gt[k]) < 1e-12
    t = pd.DataFrame([{'model': k, 'gamma_h': g[k], 'gamma_tea': g[f'tea_{k}']} for k in KEY])
    t.to_csv(STRENGTHS, index=False)
    print(t.to_string(index=False))


if __name__ == '__main__':
    if '--strengths' in sys.argv:
        write_strengths()
