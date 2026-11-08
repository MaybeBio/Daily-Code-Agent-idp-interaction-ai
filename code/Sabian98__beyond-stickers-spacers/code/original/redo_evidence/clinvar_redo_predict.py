"""
clinvar_redo_predict.py — step 1 of clinvar_redo_plan.md (hash 3ee66781). From sequence only, no class is read:
degree predictions (frozen CALVADOS-2 degree model, wild type and mutant, isolated IDR), AlphaMissense lookup,
conservation and IDRome's Delta_nu for every variant of clinvar_hub/variants.csv.

  OMP_NUM_THREADS=1 ~/.venvs/bender/bin/python clinvar_redo_predict.py  ->  clinvar_redo/predictions.csv
"""
import csv
import gzip
import os
import pickle
import subprocess
import sys
from collections import defaultdict
from multiprocessing import Pool

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, 'hub_predictor'))
sys.path.insert(0, HERE)
import train_hub_model as T  # noqa: E402
import clinvar_hub as CH  # noqa: E402

DEG = os.path.join(HERE, 'hub_vs_degree', 'degree_model_cal.pkl')
AM = os.path.join(os.path.dirname(HERE), 'external_data', 'alphamissense', 'AlphaMissense_aa_substitutions.tsv.gz')
OUT = os.path.join(HERE, 'clinvar_redo')
_M = None


def _init():
    global _M
    _M = pickle.load(open(DEG, 'rb'))


def _job(job):
    idr, seq, muts = job                          # muts: [(variant, k, wt, mut)]
    w = _M.predict(T.features(seq, 10))
    out = []
    for v, k, wt, mt in muts:
        ms = seq[:k] + mt + seq[k + 1:]
        m = _M.predict(T.features(ms, 10))
        out.append((v, float(m[k] - w[k]), float(np.abs(m - w).sum()), float(w[k])))
    return out


def main():
    os.makedirs(OUT, exist_ok=True)
    D = CH.idrome()
    V = list(csv.DictReader(open(os.path.join(HERE, 'clinvar_hub', 'variants.csv'))))
    by_idr = defaultdict(list)
    for r in V:
        k = int(r['k'])
        assert D[r['idr']]['seq'][k] == r['wt'], r['variant']
        by_idr[r['idr']].append((r['variant'], k, r['wt'], r['mut']))
    jobs = sorted(((i, D[i]['seq'], m) for i, m in by_idr.items()), key=lambda j: -len(j[1]) * len(j[2]))
    pred = {}
    with Pool(8, initializer=_init) as pool:
        for res in pool.imap_unordered(_job, jobs, chunksize=1):
            for v, ds, dg, ws in res:
                pred[v] = (ds, dg, ws)
    print('degree predictions for', len(pred), 'variants', flush=True)

    cons = {}
    for r in csv.DictReader(open(os.path.join(HERE, 'clinvar_hub', 'conservation.csv'))):
        cons[(r['idr'], int(r['pos']))] = float(r['identity'])
    dnu = {}
    rd = csv.reader(gzip.open(os.path.join(CH.EXT, 'missense_variants.csv.gz'), 'rt'))
    h = next(rd)
    ix = {c: i for i, c in enumerate(h)}
    for row in rd:
        try:
            dnu[row[0]] = float(row[ix['Delta_nu']])
        except ValueError:
            pass

    need = {}
    for r in V:
        d = D[r['idr']]
        need[(d['uniprot'], f'{r["wt"]}{d["first"] + int(r["k"])}{r["mut"]}')] = r['variant']
    accs = {a for a, _ in need}
    am = {}
    p = subprocess.Popen(['gzip', '-dc', AM], stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
    n_lines = 0
    for line in p.stdout:
        n_lines += 1
        t = line.find('\t')
        if t < 0 or line[:t] not in accs:
            continue
        f = line.rstrip('\n').split('\t')
        key = (f[0], f[1])
        if key in need:
            am[need[key]] = float(f[2])
    p.wait()
    print(f'AlphaMissense lines read {n_lines}; variants matched {len(am)} of {len(V)}', flush=True)

    with open(os.path.join(OUT, 'predictions.csv'), 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['set', 'variant', 'variation_id', 'idr', 'gene', 'k', 'wt', 'mut', 'N', 'dist_end', 'dsite', 'dglob',
                    'wsite', 'cons', 'am', 'delta_nu'])
        for r in V:
            k, n = int(r['k']), len(D[r['idr']]['seq'])
            ds, dg, ws = pred[r['variant']]
            c = cons.get((r['idr'], k))
            w.writerow([r['set'], r['variant'], r['variation_id'], r['idr'], r['gene'], k, r['wt'], r['mut'], n,
                        min(k, n - 1 - k), f'{ds:.6f}', f'{dg:.6f}', f'{ws:.6f}', '' if c is None else f'{c:.4f}',
                        '' if r['variant'] not in am else f'{am[r["variant"]]:.4f}',
                        '' if r['variant'] not in dnu else f'{dnu[r["variant"]]:.6f}'])
    print('written', len(V), flush=True)


if __name__ == '__main__':
    main()
