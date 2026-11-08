"""
freeze_predictions.py — blind, frozen hub-occupancy predictions for the GAIVI design batch (../gaivi/hub_design_table.csv),
made before any of those sequences has been simulated, so the design batch is a prospective test.

Models (gradient boosting, same features and settings as train_hub_model.py / train_hub_multiff.py; target occupancy):
  cal        CALVADOS-2 only, all proteins in DATA
  pooled     all force fields, residues weighted so that each force field contributes equally
  consensus  proteins with labels in >= 2 force fields; target = mean occupancy across force fields
Writes ../gaivi/predictions_frozen_<date>.csv (id, pos, aa, cal, pooled, consensus) and appends its SHA-256 and the
training-set sizes to ../gaivi/FROZEN.sha256. Never overwrites an existing frozen file: if it exists, the retrained
models must reproduce it exactly (the run stops otherwise), which shows the models are the frozen ones.
The fitted models are pickled to frozen_models_<date>.pkl (hash recorded too), so later sequence sets can be scored
with exactly the same models.
  python3 freeze_predictions.py [--data data_full] [--extra SEQS.csv --tag NAME]   # SEQS.csv: id, sequence
"""
import argparse
import csv
import datetime
import hashlib
import os
import pickle
from collections import defaultdict

import numpy as np

import train_hub_model as T
from train_hub_multiff import fit_w

HERE = os.path.dirname(os.path.abspath(__file__))
GAIVI = os.path.join(os.path.dirname(HERE), 'gaivi')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--data', default=os.path.join(HERE, 'data_full'))
    ap.add_argument('--window', type=int, default=10)
    ap.add_argument('--extra', help='CSV with columns id, sequence: another set to freeze predictions for')
    ap.add_argument('--tag', default='extra')
    ap.add_argument('--float32', action='store_true', help='predict from float32 features, as the models were trained')
    ap.add_argument('--date', default=datetime.date.today().isoformat(), help='date of the frozen design file to check')
    a = ap.parse_args()
    design = list(csv.DictReader(open(os.path.join(GAIVI, 'hub_design_table.csv'))))
    out = os.path.join(GAIVI, f'predictions_frozen_{a.date}.csv')
    out_extra = os.path.join(GAIVI, f'predictions_frozen_{a.date}_{a.tag}.csv') if a.extra else None
    if out_extra and os.path.exists(out_extra):
        raise SystemExit(f'{out_extra} exists; frozen predictions are never overwritten')

    prots, y = T.load(a.data, 'occupancy')
    ffs = sorted({k[0] for k in prots})
    feats = {}

    def F(seq):
        if seq not in feats:
            feats[seq] = T.features(seq, a.window).astype(np.float32)
        return feats[seq]

    rows = defaultdict(list)                     # ff -> [(uid, seq, labels)]
    for (ff, uid), r in prots.items():
        lab = y.get((ff, uid), {})
        if len(lab) == len(r['sequence']):
            rows[ff].append((uid, r['sequence'], np.array([lab[i] for i in range(len(lab))], np.float32)))
    sizes = {ff: len(v) for ff, v in rows.items()}
    print('training proteins:', sizes, flush=True)
    models = {}

    X = np.vstack([F(s) for _, s, _ in rows['CALVADOS-2']])
    Y = np.concatenate([l for _, _, l in rows['CALVADOS-2']])
    models['cal'] = fit_w('hgb', X, Y, np.ones(len(Y)))
    print('cal model:', len(Y), 'residues', flush=True)
    del X, Y

    Xs, Ys, W = [], [], []
    for ff in ffs:
        n_ff = sum(len(l) for _, _, l in rows[ff])
        for _, s, l in rows[ff]:
            Xs.append(F(s))
            Ys.append(l)
            W.append(np.full(len(l), 1.0 / n_ff, np.float32))
    X, Y, w = np.vstack(Xs), np.concatenate(Ys), np.concatenate(W)
    models['pooled'] = fit_w('hgb', X, Y, w * len(Y) / len(ffs))
    print('pooled model:', len(Y), 'residues', flush=True)
    del X, Y, w, Xs, Ys, W

    by = defaultdict(dict)
    for ff in ffs:
        for uid, s, l in rows[ff]:
            by[(uid, s)][ff] = l
    multi = [k for k, d in by.items() if len(d) >= 2]
    X = np.vstack([F(s) for _, s in multi])
    Y = np.concatenate([np.mean(list(by[k].values()), axis=0) for k in multi])
    models['consensus'] = fit_w('hgb', X, Y, np.ones(len(Y)))
    print('consensus model:', len(multi), 'proteins,', len(Y), 'residues', flush=True)
    del X, Y

    def rows_for(items):
        rows = [['id', 'pos', 'aa', 'cal', 'pooled', 'consensus']]
        for uid, seq in items:
            Xd = T.features(seq, a.window)
            if a.float32:                          # same dtype as training (the 2026-10-02 files used float64 here)
                Xd = Xd.astype(np.float32)
            preds = {k: m.predict(Xd) for k, m in models.items()}
            rows += [[uid, str(i), aa] + [f'{preds[k][i]:.5f}' for k in ('cal', 'pooled', 'consensus')]
                     for i, aa in enumerate(seq)]
        return rows

    def freeze(path, rows, note):
        with open(path, 'w', newline='') as fh:
            csv.writer(fh).writerows(rows)
        h = hashlib.sha256(open(path, 'rb').read()).hexdigest()
        with open(os.path.join(GAIVI, 'FROZEN.sha256'), 'a') as fh:
            fh.write(f'{h}  {os.path.basename(path)}   # {datetime.datetime.now().isoformat(timespec="seconds")}; {note}\n')
        print(f'-> {path}\nsha256 {h}')

    if os.path.exists(out):
        old = list(csv.reader(open(out)))
        frozen_ids = list(dict.fromkeys(r[0] for r in old[1:]))      # the design table may have grown since
        seq_of = {d['id']: d['sequence'] for d in design}
        design_rows = rows_for([(u, seq_of[u]) for u in frozen_ids])
        same = sum(x == y for x, y in zip(old, design_rows))
        print(f'design predictions reproduced: {same}/{len(old)} rows identical to {os.path.basename(out)}')
        if same != len(old) or len(old) != len(design_rows):
            raise SystemExit('retrained models differ from the frozen ones; nothing written')
    else:
        design_rows = rows_for([(d['id'], d['sequence']) for d in design])
        freeze(out, design_rows, f'training proteins {sizes}; consensus {len(multi)}; design sequences {len(design)}')
    mp = os.path.join(HERE, f'frozen_models_{a.date}.pkl')
    if not os.path.exists(mp):
        pickle.dump(models, open(mp, 'wb'))
        h = hashlib.sha256(open(mp, 'rb').read()).hexdigest()
        with open(os.path.join(GAIVI, 'FROZEN.sha256'), 'a') as fh:
            fh.write(f'{h}  hub_predictor/{os.path.basename(mp)}   # the models behind predictions_frozen_{a.date}*.csv\n')
        print(f'-> {mp}')
    if out_extra:
        extra = [(r['id'], r['sequence']) for r in csv.DictReader(open(a.extra))]
        freeze(out_extra, rows_for(extra), f'same models; {len(extra)} sequences from {os.path.basename(a.extra)}')


if __name__ == '__main__':
    main()
