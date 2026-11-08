"""
phospho_function_sites.py — the per-site ensemble table of Part F, phospho_charge_function_plan.md (hash 28befc80) with
Amendment 1 (hash a78bf482). Written and hashed before it was run, and run before the Ochoa et al. functional scores
were available; the output is hashed before any merge with the scores.

For every Ser/Thr of every human BENDER entry (taxonomy 9606, CALVADOS-2, 300 K):
  log ESS     charge-model weights of phospho_network (dq = -2, the stored 1,000 frames, salt_ph_phospho potential)
  dRg_rw_pct  reweighted Rg change (%)
  seqpred     frozen sequence charge model (phospho_charge/seq_model.json) prediction of log ESS
  dSCD        change in sequence charge decoration when q_s goes from 0 to -2 (Amendment 1), and |dSCD|
  KR5, DE5    K/R and D/E counts within |i-j| <= 5
  R_hub_rw, R_deg_rw   reweighted network change (sum over residues of |change| in hub occupancy and in degree), chains
                       of at most 300 residues, with phospho_network.networks (the paper's rule)

  BENDER_REANALYSIS_DIR=/Users/taseef/Downloads/Scripts/bender_reanalysis ~/.venvs/finches/bin/python phospho_function_sites.py
      ->  phospho_function/sites_human.csv, phospho_function/failed.txt
"""
import json
import os
import sys

os.environ.setdefault('BENDER_REANALYSIS_DIR', '/Users/taseef/Downloads/Scripts/bender_reanalysis')
DYN = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, DYN)

import multiprocessing as mp  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import phospho_network as pn  # noqa: E402  (networks(), sp = salt_ph_phospho, DQ = -2)
from phospho_charge_test import features, local_counts  # noqa: E402
from phospho_direct_test import check_frozen  # noqa: E402

BIO = '/Volumes/User Homes/Nu Project/BENDER-BIO/BENDER_BIO.csv'
MODEL = os.path.join(DYN, 'phospho_charge', 'seq_model.json')
OUT = os.path.join(DYN, 'phospho_function')
NMAX_NET = 300
WORKERS = 8
Q = {'K': 1, 'R': 1, 'D': -1, 'E': -1}


def dscd(seq, pos):
    i, n = pos - 1, len(seq)
    return sum(-2.0 * Q.get(c, 0) * abs(i - j) ** 0.5 for j, c in enumerate(seq) if j != i) / n


def one(task):
    uid, hint, bio_seq, coef = task
    try:
        seq, rg, _, _, phi, beta = pn.sp.calvados_terms(uid, hint)
        n = len(seq)
        sites = [k + 1 for k, c in enumerate(seq) if c in 'ST']
        if not sites:
            return [], None
        net = n <= NMAX_NET
        if net:
            xyz = pn.sp.bs.frames_calvados(uid, 300, hint)[0]
            D, H = pn.networks(xyz, n, uid)
            deg0, occ0, Hf = D.mean(0), H.mean(0), H.astype(float)
        rows = []
        for pos in sites:
            i = pos - 1
            x = -beta * pn.DQ * phi[:, i]
            w = np.exp(x - x.max())
            w /= w.sum()
            ess = float(1.0 / (w ** 2).sum() / len(w))
            f = features(seq, pos)
            d = dscd(seq, pos)
            kr, de = local_counts(seq, pos)
            rows.append({'UniProt_ID': uid, 'N': n, 'pos': pos, 'residue': seq[i], 'seq_matches_bio': seq == bio_seq,
                         'ess': ess, 'log_ess': float(np.log(ess)),
                         'dRg_rw_pct': (float(w @ rg) / rg.mean() - 1) * 100,
                         'seqpred': float(coef[0] + np.dot(f, coef[1:])), 'dSCD': d, 'abs_dSCD': abs(d),
                         'KR5': kr, 'DE5': de,
                         'R_hub_rw': float(np.abs(w @ Hf - occ0).sum()) if net else np.nan,
                         'R_deg_rw': float(np.abs(w @ D - deg0).sum()) if net else np.nan})
        return rows, None
    except Exception as e:  # noqa: BLE001
        return None, f'{uid}: {e!r}'


def main():
    check_frozen(MODEL)
    os.makedirs(OUT, exist_ok=True)
    m = json.load(open(MODEL))
    coef = np.array([m['intercept']] + [m['coef'][k] for k in m['features']])
    bio = pd.read_csv(BIO, usecols=['UniProt_ID', 'taxonomy_id', 'sequence'], low_memory=False)
    hum = bio[bio.taxonomy_id == 9606]
    kingdom = pd.read_csv(pn.sp.bc.BIO_CSV, usecols=['UniProt_ID', 'kingdom']).set_index('UniProt_ID')['kingdom']
    tasks = [(u, kingdom.get(u), s, coef) for u, s in zip(hum.UniProt_ID, hum.sequence)]
    print(f'{len(tasks)} human entries', flush=True)
    rows, fails = [], []
    with mp.Pool(WORKERS) as pool:
        for k, (r, err) in enumerate(pool.imap_unordered(one, tasks), 1):
            if err:
                fails.append(err)
            else:
                rows += r
            if k % 50 == 0:
                print(f'  {k}/{len(tasks)}, {len(fails)} failed', flush=True)
    df = pd.DataFrame(rows).sort_values(['UniProt_ID', 'pos'])
    df.to_csv(os.path.join(OUT, 'sites_human.csv'), index=False)
    open(os.path.join(OUT, 'failed.txt'), 'w').write('\n'.join(fails) + '\n')
    print(f'{len(df)} Ser/Thr in {df.UniProt_ID.nunique()} proteins; {len(fails)} proteins failed', flush=True)


if __name__ == '__main__':
    main()
