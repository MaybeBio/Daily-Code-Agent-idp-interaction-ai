"""
pair_type_split.py — which residue pairs carry the sticker susceptibility? The sticker covariance splits exactly
over pair types,
    chi_s = Cov(Rg, E_s) / (R T^2) = sum over types t of Cov(Rg, E_t) / (R T^2),
with E_t the sticker energy (Wang-Frenkel for Mpipi-GG, Ashbaugh-Hatch for CALVADOS-2) summed over the nonbonded
pairs of type t. Types, from the residues of the pair (aro = F, W, Y):
    aro-aro, R-aro, K-aro, R-other, K-other, aro-other, other-other      (R-K, R-R, K-K count as R-other / K-other
    with R taking precedence)
Per protein at 300 K, on 200 of the 1,000 FDT frames (contact_count_check.py frame selection and pair energies):
  share_chi_<t>     share of Cov(Rg, E_s) from type t (sums to 1)
  share_E_<t>       share of the mean sticker energy
  share_pairs_<t>   share of the nonbonded pairs (the composition baseline)
  chi_per_pair_<t>  Cov(Rg, E_t) / n_pairs_t, relative to the protein's all-pair mean (1 = average pair)
Within-protein ratios of chi per pair: R-aro / K-aro and R-other / K-other (proteins with >= 20 pairs of both).
Sample: up to 40 proteins per own-nu bin (all six bins; random, seed 0) from chi_factor_split/per_protein.csv.

Output: Traj Level Analysis/Analysis_FDT_MPIPI/compact/pair_type_split/ (per_protein.csv, summary.csv, bins.csv)
Usage: python3 pair_type_split.py [--workers 6]
"""
import sys

sys.dont_write_bytecode = True  # keep the NAS code folder free of __pycache__

import argparse  # noqa: E402
import multiprocessing as mp  # noqa: E402
import os  # noqa: E402
import traceback  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

HERE = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, HERE)
import contact_count_check as cc  # noqa: E402

bc = cc.bc
MP = cc.MP
OUT_DIR = f'{MP}/compact/pair_type_split'
TYPES = ['aro-aro', 'R-aro', 'K-aro', 'R-other', 'K-other', 'aro-other', 'other-other']
PER_BIN, SEED, MIN_PAIRS = 40, 0, 20
ARO = set('FWY')


def pair_types(seq, i, j):
    a, b = np.array(list(seq))[i], np.array(list(seq))[j]
    ia, ib = np.isin(a, list(ARO)), np.isin(b, list(ARO))
    ra, rb, ka, kb = a == 'R', b == 'R', a == 'K', b == 'K'
    t = np.full(len(i), 'other-other', dtype=object)
    t[ia | ib] = 'aro-other'
    t[(ka | kb)] = 'K-other'
    t[(ra | rb)] = 'R-other'
    t[(ka & ib) | (kb & ia)] = 'K-aro'
    t[(ra & ib) | (rb & ia)] = 'R-aro'
    t[ia & ib] = 'aro-aro'
    return t


def npz_path(uid, model):
    if model == 'CALVADOS-2':
        return f'{cc.MT}/Analysis_FDT/npz/300K/{uid}.npz'
    p = f'{MP}/npz_300K/{uid}.npz'
    return p if os.path.exists(p) else f'{MP}/compact/npz_300K/{uid}.npz'


def analyse(args):
    uid, model, hint = args
    try:
        xyz, i, j, ufun, _ = cc.pair_energies_mpipi(uid) if model == 'Mpipi-GG' else cc.pair_energies_calvados(uid, hint)
        seq = str(np.load(npz_path(uid, model))['sequence'])
        if len(seq) != xyz.shape[1]:
            raise ValueError(f'sequence length {len(seq)} != {xyz.shape[1]} beads')
        typ = pair_types(seq, i, j)
        nf = xyz.shape[0]
        rg = np.sqrt(((xyz - xyz.mean(1, keepdims=True)) ** 2).sum(2).mean(1))
        Et = {t: np.zeros(nf) for t in TYPES}
        masks = {t: typ == t for t in TYPES}
        step = max(1, int(cc.MAX_ELEMS // len(i)))
        for a in range(0, nf, step):
            u = ufun(np.linalg.norm(xyz[a:a + step, i] - xyz[a:a + step, j], axis=2))
            for t in TYPES:
                if masks[t].any():
                    Et[t][a:a + step] = u[:, masks[t]].sum(1)
        E = sum(Et.values())
        cov = {t: np.cov(rg, Et[t])[0, 1] for t in TYPES}
        ctot = np.cov(rg, E)[0, 1]
        npairs = {t: int(masks[t].sum()) for t in TYPES}
        per_pair_all = ctot / len(i)
        row = {'UniProt_ID': uid, 'model': model, 'N': len(seq), 'status': 'ok',
               'f_R': seq.count('R') / len(seq), 'f_K': seq.count('K') / len(seq),
               'f_aro': sum(seq.count(c) for c in ARO) / len(seq)}
        for t in TYPES:
            row[f'share_chi_{t}'] = cov[t] / ctot
            row[f'share_E_{t}'] = Et[t].mean() / E.mean()
            row[f'share_pairs_{t}'] = npairs[t] / len(i)
            row[f'n_pairs_{t}'] = npairs[t]
            row[f'chi_per_pair_{t}'] = cov[t] / npairs[t] / per_pair_all if npairs[t] else np.nan
        for x in ('aro', 'other'):
            ok = npairs[f'R-{x}'] >= MIN_PAIRS and npairs[f'K-{x}'] >= MIN_PAIRS
            row[f'ratio_R_vs_K_{x}'] = row[f'chi_per_pair_R-{x}'] / row[f'chi_per_pair_K-{x}'] if ok else np.nan
        return row
    except Exception as e:
        return {'UniProt_ID': uid, 'model': model, 'status': f'failed: {e!r}'[:200], 'trace': traceback.format_exc()[-300:]}


def main(workers):
    os.makedirs(OUT_DIR, exist_ok=True)
    src = pd.read_csv(cc.SRC, index_col=0)
    pick = pd.concat([g.sample(min(len(g), PER_BIN), random_state=SEED) for _, g in src.groupby(['model', 'b'])])
    kingdom = pd.read_csv(bc.BIO_CSV, usecols=['UniProt_ID', 'kingdom']).set_index('UniProt_ID')['kingdom']
    tasks = [(u, m, kingdom.get(u)) for u, m in zip(pick.index, pick['model'])]
    print(f'{len(tasks)} proteins, {workers} workers', flush=True)
    res = []
    with mp.Pool(workers) as pool:
        for k, r in enumerate(pool.imap_unordered(analyse, tasks), 1):
            res.append(r)
            if r['status'] != 'ok' or k % 60 == 0:
                print(f"[{k}/{len(tasks)}] {r['model']} {r['UniProt_ID']} {r['status']}", flush=True)
    per = pd.DataFrame(res).merge(pick[['model', 'nu', 'b']].reset_index(), on=['UniProt_ID', 'model'])
    per.to_csv(f'{OUT_DIR}/per_protein.csv', index=False)
    ok = per[per['status'] == 'ok']

    rows = []
    for model, g in ok.groupby('model'):
        for t in TYPES:
            has = g[f'n_pairs_{t}'] >= MIN_PAIRS
            rows.append({'model': model, 'type': t, 'n_proteins': int(has.sum()),
                         'median_share_chi': g[f'share_chi_{t}'].median(), 'median_share_E': g[f'share_E_{t}'].median(),
                         'median_share_pairs': g[f'share_pairs_{t}'].median(),
                         'median_chi_per_pair': g.loc[has, f'chi_per_pair_{t}'].median()})
        for x in ('aro', 'other'):
            v = g[f'ratio_R_vs_K_{x}'].dropna()
            rows.append({'model': model, 'type': f'ratio R-{x} / K-{x} per pair', 'n_proteins': len(v),
                         'median_chi_per_pair': v.median(), 'frac_gt_1': (v > 1).mean()})
    summ = pd.DataFrame(rows)
    summ.to_csv(f'{OUT_DIR}/summary.csv', index=False)
    bins = ok.groupby(['model', 'b'])[[f'share_chi_{t}' for t in TYPES]].median()
    bins.insert(0, 'n', ok.groupby(['model', 'b']).size())
    bins.to_csv(f'{OUT_DIR}/bins.csv')
    bc.write_provenance(f'{OUT_DIR}/summary.csv', __file__, {'types': TYPES, 'frames': f'every {cc.STRIDE}th of the 1,000 FDT frames',
                                                           'per_bin': PER_BIN, 'seed': SEED, 'min_pairs': MIN_PAIRS}, [cc.SRC])
    with pd.option_context('display.width', 250, 'display.max_columns', 20):
        print(summ.round(3).to_string(index=False))
        print('\nmedian share of the sticker covariance by own-nu bin')
        print(bins.round(3).to_string())
    print(f"{(per['status'] != 'ok').sum()} failed; saved {OUT_DIR}")


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--workers', type=int, default=6)
    main(ap.parse_args().workers)
