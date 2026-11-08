"""
freeze_v3_models.py — train the full sequence model and the sticker-and-spacer context model (with proline) on all v3
labels of each force field, pickle them, and write blind predictions for a design set before any of its simulations
exist. Hashes go to ../gaivi/FROZEN.sha256; nothing is ever overwritten.
  Models per force field (CALVADOS-2, Mpipi-GG, Mpipi-Recharged), same booster as cv_full.py:
    full      the 487 sequence features of train_hub_model.features (float32)
    context   chain position and length + sticker_spacer features with the force field's own stickiness scale + proline
              (the 'context own + proline' model of context_proline_full.py)
  ~/.venvs/bender/bin/python freeze_v3_models.py --seqs ../gaivi/hub_design_glycine_seqs.csv --tag glycine_2026-10-04
      ->  models_v3_2026-10-04.pkl, ../gaivi/predictions_v3_2026-10-04_<tag>.csv
"""
import argparse
import csv
import datetime
import hashlib
import os
import pickle

import numpy as np

import context_model_by_ff as C
import context_proline_full as CP
import sticker_spacer_test as SS
import train_hub_model as T

HERE = os.path.dirname(os.path.abspath(__file__))
GAIVI = os.path.join(os.path.dirname(HERE), 'gaivi')
SHORT = {'CALVADOS-2': 'cal', 'MPIPI-GG': 'gg', 'Mpipi-Recharged': 'rc'}


def register(path, note):
    h = hashlib.sha256(open(path, 'rb').read()).hexdigest()
    with open(os.path.join(GAIVI, 'FROZEN.sha256'), 'a') as fh:
        fh.write(f'{h}  {os.path.relpath(path, GAIVI)}'          # every entry relative to gaivi/ (review fix)
                 f'   # {datetime.datetime.now().isoformat(timespec="seconds")}; {note}\n')
    return h


def ctx_features(seq, ff, sc, chg):
    P = T.features(seq, 10).astype(np.float32)[:, C.POS_COLS]
    return np.hstack([P, SS.sticker_features(seq, sc[C.OWN[ff]], chg).astype(np.float32), CP.proline(seq)])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--seqs', required=True)
    ap.add_argument('--tag', required=True)
    ap.add_argument('--date', default='2026-10-04')
    a = ap.parse_args()
    mp = os.path.join(HERE, f'models_v3_{a.date}.pkl')
    out = os.path.join(GAIVI, f'predictions_v3_{a.date}_{a.tag}.csv')
    for f in (out,):
        if os.path.exists(f):
            raise SystemExit(f'{f} exists; frozen predictions are never overwritten')
    sc, chg = C.scales()
    if os.path.exists(mp):
        h = hashlib.sha256(open(mp, 'rb').read()).hexdigest()
        if h not in open(os.path.join(GAIVI, 'FROZEN.sha256')).read():
            raise SystemExit(f'{mp} exists but its hash is not in FROZEN.sha256')
        models = pickle.load(open(mp, 'rb'))
        print('loaded', mp, flush=True)
    else:
        prots, y = T.load(os.path.join(HERE, 'data_v3'), 'occupancy')
        models, sizes = {}, {}
        for ff in C.FFS:
            keys, seqs, Y = C.label_set(prots, y, ff)
            Yc = np.concatenate(Y)
            Xc = np.vstack([ctx_features(s, ff, sc, chg) for s in seqs])
            models[(ff, 'context')] = SS.hgb().fit(Xc, Yc)
            del Xc
            Xf = np.vstack([T.features(s, 10).astype(np.float32) for s in seqs])
            models[(ff, 'full')] = SS.hgb().fit(Xf, Yc)
            del Xf
            sizes[ff] = len(keys)
            print(f'{ff}: trained on {len(keys)} proteins, {len(Yc)} residues', flush=True)
        pickle.dump(models, open(mp, 'wb'))
        register(mp, 'v3 labels (data_v3): full and context+proline HGB per force field, trained on all proteins '
                     + ', '.join(f'{SHORT[f]} {n}' for f, n in sizes.items()))
    seqs = [(r['id'], r['sequence']) for r in csv.DictReader(open(a.seqs))]
    cols = [f'{SHORT[ff]}_{kind}' for ff in C.FFS for kind in ('full', 'context')]
    with open(out, 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['id', 'pos', 'aa'] + cols)
        for uid, seq in seqs:
            F = T.features(seq, 10).astype(np.float32)
            pr = {}
            for ff in C.FFS:
                pr[f'{SHORT[ff]}_full'] = models[(ff, 'full')].predict(F)
                pr[f'{SHORT[ff]}_context'] = models[(ff, 'context')].predict(ctx_features(seq, ff, sc, chg))
            w.writerows([uid, i, aa] + [f'{pr[c][i]:.5f}' for c in cols] for i, aa in enumerate(seq))
    h = register(out, f'models models_v3_{a.date}.pkl; {len(seqs)} sequences from {os.path.basename(a.seqs)}')
    print(f'-> {out}\nsha256 {h}')


if __name__ == '__main__':
    main()
