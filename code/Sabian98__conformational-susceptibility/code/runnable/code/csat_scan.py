"""
csat_scan.py — the capped hydration forcing at several saturation values c_sat for every chain at 300 K (2026-10-26), for
the sensitivity check of the calibrated model. Same frames, pair shapes and per-residue coefficients as capped_proteome.py,
whose per-chain loaders are reused; only the saturation function changes. One pass over each trajectory gives all values.
Per chain, for each c_sat in CSATS: dchi_rg_rel_chg_cap<c>, dchi_rg_rel_tea_cap<c> (% per K per unit strength) and the
D spreads, as in capped_proteome.finish.
Check: the c_sat = 2 columns must equal the stored capped_proteome terms (dchi_rg_rel_chg_cap, dchi_rg_rel_tea_cap) to
1e-6 relative.
    BENDER_REANALYSIS_DIR=... python csat_scan.py --model calvados|recharged|gg [--workers 6] [--limit N]
Output: Traj Level Analysis/Analysis_FDT/capped_300K/csat_scan/csat_dchi_<model>.csv (one row per chain, resumable)
"""
import sys

sys.dont_write_bytecode = True

import argparse  # noqa: E402
import multiprocessing as mp  # noqa: E402
import os  # noqa: E402
import time  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

HERE = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, os.environ.get('BENDER_REANALYSIS_DIR', HERE))
sys.path.insert(0, HERE)
import capped_proteome as CP  # noqa: E402

CSATS = [1.25, 1.5, 2.0, 3.0]
OUT = os.path.join(CP.OUT, 'csat_scan')


def forcings_scan(seq, beta, s0, unit):
    """Per-frame capped forcings of both terms at every c_sat in CSATS (replaces capped_proteome.forcings here)."""
    c = -beta / s0
    out = {}
    for k in ('chg', 'tea'):
        f = np.array([CP.TERMS[k][a] for a in seq]) * unit
        for cs in CSATS:
            out[f'{k}_cap{cs:g}'] = -0.5 * s0 * ((cs * (1.0 - np.exp(-c / cs))) @ f)
    return out, c


# top level, so that spawned workers (which re-import this module) use the scan version too
CP.forcings = forcings_scan


def check(model, path):
    new = pd.read_csv(path).set_index('UniProt_ID')
    old = pd.read_csv(os.path.join(CP.OUT, f'capped_dchi_{model}_300K.csv')).set_index('UniProt_ID')
    ids = new.index.intersection(old.index)
    rel = {}
    for k in ('chg', 'tea'):
        a, b = new.loc[ids, f'dchi_rg_rel_{k}_cap2'], old.loc[ids, f'dchi_rg_rel_{k}_cap']
        rel[k] = float((a - b).abs().max() / max(b.abs().max(), 1e-12))
    return len(ids), rel


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--model', choices=list(CP.RUNNERS), required=True)
    ap.add_argument('--workers', type=int, default=6)
    ap.add_argument('--limit', type=int, default=None)
    ap.add_argument('--check-only', action='store_true')
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    path = os.path.join(OUT, f'csat_dchi_{a.model}.csv')
    if not a.check_only:
        done = set(pd.read_csv(path, usecols=['UniProt_ID']).UniProt_ID) if os.path.exists(path) else set()
        todo = [(u, d) for u, d in CP.dirs(a.model).items() if u not in done][:a.limit]
        print(f'{a.model}: {len(done)} done, {len(todo)} to run, c_sat {CSATS}', flush=True)
        t0, nf = time.time(), 0
        with mp.get_context('spawn').Pool(a.workers) as pool:
            for k, (row, err) in enumerate(pool.imap_unordered(CP.RUNNERS[a.model], todo), 1):
                if err:
                    nf += 1
                    print(f"failed {row['UniProt_ID']}: {err[:200]}", flush=True)
                    continue
                pd.DataFrame([row]).to_csv(path, mode='a', header=not os.path.exists(path), index=False)
                if k % 500 == 0:
                    el = (time.time() - t0) / 3600
                    print(f'{k}/{len(todo)} {el:.2f} h, ~{el / k * (len(todo) - k):.1f} h left, {nf} failed', flush=True)
    n, rel = check(a.model, path)
    print(f'check on {n} chains: c_sat 2 matches the stored capped terms to {rel["chg"]:.1e} (hydrophilic), '
          f'{rel["tea"]:.1e} (TEA), relative', flush=True)


if __name__ == '__main__':
    main()
