"""
predict_frozen.py — score a new set of sequences with the pickled frozen models (frozen_models_<date>.pkl, written by
freeze_predictions.py), so later test sets get blind predictions from exactly the same models without retraining.
Checks the pickle against its SHA-256 in ../gaivi/FROZEN.sha256 first, never overwrites, and records the new file's hash.
  python3 predict_frozen.py --seqs SEQS.csv --tag NAME [--date 2026-10-02]      # SEQS.csv: id, sequence
"""
import argparse
import csv
import datetime
import hashlib
import os
import pickle

import numpy as np

import train_hub_model as T

HERE = os.path.dirname(os.path.abspath(__file__))
GAIVI = os.path.join(os.path.dirname(HERE), 'gaivi')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--seqs', required=True)
    ap.add_argument('--tag', required=True)
    ap.add_argument('--date', default='2026-10-02')
    ap.add_argument('--window', type=int, default=10)
    ap.add_argument('--float32', action='store_true', help='predict from float32 features, as the models were trained')
    a = ap.parse_args()
    mp = os.path.join(HERE, f'frozen_models_{a.date}.pkl')
    h = hashlib.sha256(open(mp, 'rb').read()).hexdigest()
    reg = open(os.path.join(GAIVI, 'FROZEN.sha256')).read()
    if h not in reg:
        raise SystemExit(f'{mp} does not match any hash in FROZEN.sha256')
    models = pickle.load(open(mp, 'rb'))
    out = os.path.join(GAIVI, f'predictions_frozen_{a.date}_{a.tag}.csv')
    if os.path.exists(out):
        raise SystemExit(f'{out} exists; frozen predictions are never overwritten')
    seqs = [(r['id'], r['sequence']) for r in csv.DictReader(open(a.seqs))]
    with open(out, 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['id', 'pos', 'aa', 'cal', 'pooled', 'consensus'])
        for uid, seq in seqs:
            X = T.features(seq, a.window)
            if a.float32:
                X = X.astype(np.float32)
            pr = {k: m.predict(X) for k, m in models.items()}
            w.writerows([uid, i, aa] + [f'{pr[k][i]:.5f}' for k in ('cal', 'pooled', 'consensus')]
                        for i, aa in enumerate(seq))
    ho = hashlib.sha256(open(out, 'rb').read()).hexdigest()
    with open(os.path.join(GAIVI, 'FROZEN.sha256'), 'a') as fh:
        fh.write(f'{ho}  {os.path.basename(out)}   # {datetime.datetime.now().isoformat(timespec="seconds")}; '
                 f'models {os.path.basename(mp)}; {len(seqs)} sequences from {os.path.basename(a.seqs)}\n')
    print(f'-> {out}\nsha256 {ho}')


if __name__ == '__main__':
    main()
