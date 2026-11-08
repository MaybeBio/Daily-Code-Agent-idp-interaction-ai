"""
context_size_control.py — is the gap between the context model and the full model in Mpipi-GG a property of the force
field, or of its small, particular protein set (908 proteins)?

Same folds, booster and features as context_model_by_ff.py, on size-matched label sets:
  Mpipi-GG                         its 908 proteins
  CALVADOS-2 on the GG proteins    the same 908 sequences with CALVADOS-2 labels
  CALVADOS-2, 908 random proteins  a seeded random draw of 908 CALVADOS-2 proteins (size without the GG selection)
  Mpipi-Recharged on GG proteins   the GG proteins that also have Mpipi-Recharged labels
Conditions: identity, context own, context + spacer types, full.
  context + spacer types   context own plus the fractions of P, D/E, K/R, N/Q/T and A/V/I/L/M in windows of +-2, +-5 and
                           +-10, and whether the residue itself is a proline (spacers are not inert in the framework)
  ~/.venvs/bender/bin/python context_size_control.py  ->  context_size_control_summary.txt
"""
import os
import time

import numpy as np

import context_model_by_ff as C
import sticker_spacer_test as SS
import train_hub_model as T

SPACER_TYPES = ['P', 'DE', 'KR', 'NQT', 'AVILM']
out = []


def p(s=''):
    out.append(s)
    print(s, flush=True)


def spacer_features(seq):
    cols = [np.array([a == 'P' for a in seq], float)]
    for w in (2, 5, 10):
        for letters in SPACER_TYPES:
            cols.append(SS.window_mean(np.array([a in letters for a in seq], float), w))
    return np.vstack(cols).T.astype(np.float32)


def run(name, keys, prots, y, scale, sc, chg):
    seqs = [prots[k]['sequence'] for k in keys]
    Y = np.concatenate([np.array([y[k][i] for i in range(len(s))], np.float32) for k, s in zip(keys, seqs)])
    F = [T.features(s, 10).astype(np.float32) for s in seqs]
    P = np.vstack([f[:, C.POS_COLS] for f in F])
    S = np.vstack([SS.sticker_features(s, sc[scale], chg).astype(np.float32) for s in seqs])
    X = {'identity': np.hstack([P, np.vstack([C.onehot(s) for s in seqs])]),
         'context own': np.hstack([P, S]),
         'context + spacer types': np.hstack([P, S, np.vstack([spacer_features(s) for s in seqs])]),
         'full': np.vstack(F)}
    G = np.concatenate([[0], np.cumsum([len(s) for s in seqs])])
    FO = np.concatenate([np.full(len(s), int(prots[k]['fold'])) for k, s in zip(keys, seqs)])
    pred = {c: np.zeros(len(Y), np.float32) for c in X}
    for f in range(5):
        tr, te = FO != f, FO == f
        for c, M in X.items():
            pred[c][te] = SS.hgb().fit(M[tr], Y[tr]).predict(M[te])
    med = {c: np.nanmedian(C.per_protein(pred[c], Y, G)) for c in X}
    p(f'  {name:34s} {len(keys):4d} proteins  ' + '   '.join(f'{c} {v:.3f}' for c, v in med.items()))


def main():
    t0 = time.time()
    sc, chg = C.scales()
    prots, y = T.load(os.path.join(C.HERE, 'data_v3'), 'occupancy')

    def complete(ff):
        return {k[1]: k for k in prots if k[0] == ff and len(y.get(k, {})) == len(prots[k]['sequence'])}
    gg, cal, rc = complete('MPIPI-GG'), complete('CALVADOS-2'), complete('Mpipi-Recharged')
    ids = sorted(gg)
    rng = np.random.default_rng(7)
    rand = sorted(rng.choice(sorted(cal), size=len(ids), replace=False))
    p('Context model versus full model on size-matched label sets (5-fold CV by protein, median per-protein Spearman)')
    run('Mpipi-GG', [gg[u] for u in ids], prots, y, 'eps GG', sc, chg)
    run('CALVADOS-2 on the GG proteins', [cal[u] for u in ids if u in cal], prots, y, 'lambda', sc, chg)
    run('CALVADOS-2, 908 random proteins', [cal[u] for u in rand], prots, y, 'lambda', sc, chg)
    run('Mpipi-Recharged on the GG proteins', [rc[u] for u in ids if u in rc], prots, y, 'eps RC', sc, chg)
    p(f'({time.time() - t0:.0f} s)')
    open(os.path.join(C.HERE, 'context_size_control_summary.txt'), 'w').write('\n'.join(out) + '\n')


if __name__ == '__main__':
    main()
