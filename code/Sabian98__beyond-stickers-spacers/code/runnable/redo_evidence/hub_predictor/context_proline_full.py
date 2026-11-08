"""
context_proline_full.py — the context-own model plus proline (the residue's own proline flag and the proline fraction in
windows of +-2, +-5 and +-10) on the full label sets, same folds and booster as context_model_by_ff.py.
  ~/.venvs/bender/bin/python context_proline_full.py  ->  context_proline_full.txt, context_proline_full_per_protein.csv
"""
import csv
import os

import numpy as np

import context_model_by_ff as C
import sticker_spacer_test as SS
import train_hub_model as T


def proline(seq):
    x = np.array([a == 'P' for a in seq], float)
    return np.vstack([x] + [SS.window_mean(x, w) for w in (2, 5, 10)]).T.astype(np.float32)


def main():
    sc, chg = C.scales()
    prots, y = T.load(os.path.join(C.HERE, 'data_v3'), 'occupancy')
    out, rows = [], []
    for ff in C.FFS:
        keys, seqs, Y = C.label_set(prots, y, ff)
        P = np.vstack([T.features(s, 10).astype(np.float32)[:, C.POS_COLS] for s in seqs])
        M = np.hstack([P, np.vstack([SS.sticker_features(s, sc[C.OWN[ff]], chg).astype(np.float32) for s in seqs]),
                       np.vstack([proline(s) for s in seqs])])
        Yc = np.concatenate(Y)
        G = np.concatenate([[0], np.cumsum([len(s) for s in seqs])])
        FO = np.concatenate([np.full(len(s), int(prots[k]['fold'])) for k, s in zip(keys, seqs)])
        pred = np.zeros(len(Yc), np.float32)
        for f in range(5):
            tr, te = FO != f, FO == f
            pred[te] = SS.hgb().fit(M[tr], Yc[tr]).predict(M[te])
        r = C.per_protein(pred, Yc, G)
        rows += [{'protein_id': k[1], 'force_field': C.NICE[ff], 'condition': 'context own + proline', 'spearman': v}
                 for k, v in zip(keys, r)]
        out.append(f'{C.NICE[ff]:16s} {len(keys):6d} proteins   context own + proline {np.nanmedian(r):.3f}')
        print(out[-1], flush=True)
    open(os.path.join(C.HERE, 'context_proline_full.txt'), 'w').write('\n'.join(out) + '\n')
    with open(os.path.join(C.HERE, 'context_proline_full_per_protein.csv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


if __name__ == '__main__':
    main()
