"""
fresh_ped_test.py — the untouched external test: hub occupancy in 8 PED ensembles of disordered proteins that were
never used before (ped_fresh_set.csv), scored against predictions frozen BEFORE these labels were computed.

Protocol (enforced): the frozen predictions file must exist and its SHA-256 must be in gaivi/FROZEN.sha256, otherwise
nothing is computed. Labels are computed exactly as for the original 18 ensembles (exp_occupancy.py): top-10%
betweenness hubs per model, fixed 0.8 nm rule, |i-j| >= 3, at most 1,000 evenly spaced models; baselines are the
Rg-matched residue-neutral null (ped_null.py) and distance from the nearest chain end.

Caveat built into the report: PED00533 (Aβ40) and PED00536 (ACTR) are time-ordered a99SB-disp MD frames whose
maximum-entropy weights are not in the deposited PDB, so as used here they are an independent all-atom force field,
not experiment. They are reported separately from the six experiment-selected ensembles.
  ~/.venvs/bender/bin/python fresh_ped_test.py  ->  hub_predictor/fresh_ped_occupancy.csv, fresh_ped_test_summary.txt
"""
import csv
import hashlib
import os
import sys
from multiprocessing import Pool

import numpy as np
from scipy.stats import spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import exp_occupancy as eo  # noqa: E402

FROZEN = os.path.join(HERE, 'gaivi', 'predictions_frozen_2026-10-02_ped_fresh.csv')
RAW_MD = {'PED00533', 'PED00536'}
out = []


def p(s=''):
    out.append(s)
    print(s, flush=True)


def main():
    h = hashlib.sha256(open(FROZEN, 'rb').read()).hexdigest()
    if h not in open(os.path.join(HERE, 'gaivi', 'FROZEN.sha256')).read():
        raise SystemExit('frozen predictions not registered; refusing to compute labels')
    pred = {}
    for r in csv.DictReader(open(FROZEN)):
        pred.setdefault(r['id'], []).append((float(r['cal']), float(r['pooled']), float(r['consensus'])))
    fresh = list(csv.DictReader(open(os.path.join(HERE, 'ped_fresh_set.csv'))))
    jobs = [(100 + k, (r['id'], r['file'], r['protein'], r['restraints'])) for k, r in enumerate(fresh)]
    with Pool(8) as pool:
        res = sorted(pool.imap_unordered(eo.one, jobs))
    rows, table = [], []
    for entry, prot, restraints, seq, oe, on, used in res:
        n = len(seq)
        pr = np.array(pred[entry])
        assert len(pr) == n, (entry, len(pr), n)
        oe_, on_ = np.array(oe), np.array(on)
        end = np.minimum(np.arange(n), n - 1 - np.arange(n))
        rho = {k: spearmanr(pr[:, j], oe_).correlation for j, k in enumerate(('cal', 'pooled', 'consensus'))}
        rho['null'] = spearmanr(on_, oe_).correlation
        rho['position'] = spearmanr(end, oe_).correlation
        table.append((entry, prot, restraints, n, used, rho))
        for i, a in enumerate(seq):
            rows.append([entry, prot, restraints, i, a, f'{oe[i]:.4f}', f'{on[i]:.4f}', used]
                        + [f'{x:.5f}' for x in pr[i]])
    with open(os.path.join(HERE, 'hub_predictor', 'fresh_ped_occupancy.csv'), 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['entry', 'protein', 'restraints', 'pos', 'aa', 'occ_exp', 'occ_null', 'models', 'pred_cal',
                    'pred_pooled', 'pred_consensus'])
        w.writerows(rows)

    keys = ('consensus', 'pooled', 'cal', 'null', 'position')
    p('Fresh external test: 8 never-used PED ensembles of disordered proteins; per-protein Spearman between predicted '
      'and ensemble hub occupancy')
    p(f'(predictions frozen before labels: {os.path.basename(FROZEN)}, sha256 {h[:16]}...)\n')
    p(f'  {"entry":9s} {"protein":38s} {"N":>4s} {"models":>6s}' + ''.join(f'{k:>11s}' for k in keys))
    for entry, prot, restraints, n, used, rho in table:
        tag = ' [raw MD]' if entry in RAW_MD else ''
        p(f'  {entry:9s} {(prot + tag)[:38]:38s} {n:4d} {used:6d}' + ''.join(f'{rho[k]:+11.2f}' for k in keys))
    for label, sel in (('experiment-selected (6)', [t for t in table if t[0] not in RAW_MD]),
                       ('raw a99SB-disp MD (2)', [t for t in table if t[0] in RAW_MD]),
                       ('all 8', table)):
        med = {k: float(np.nanmedian([t[5][k] for t in sel])) for k in keys}
        beats = {b: sum(t[5]['consensus'] > t[5][b] for t in sel) for b in ('null', 'position')}
        p(f'\n  {label}: median ' + ', '.join(f'{k} {med[k]:+.2f}' for k in keys)
          + f'; consensus beats the null in {beats["null"]}/{len(sel)} and position-only in {beats["position"]}/{len(sel)}')
    p('\n  Reference (the 18 used ensembles, not blind; frozen models): consensus +0.46, null +0.30, position +0.31 (single-pool null).')
    open(os.path.join(HERE, 'fresh_ped_test_summary.txt'), 'w').write('\n'.join(out) + '\n')


if __name__ == '__main__':
    main()
