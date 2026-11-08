"""
design_response.py — design batch: is the response of the hub network to a perturbation predictable?
Implements design_response_plan.md (hashed 2026-10-09 before any design label existed).
  ~/.venvs/bender/bin/python design_response.py --labels Mpipi-Recharged=labels_design_recharged:hub_frames_v3_recharged
      [--labels CALVADOS-2=labels_design_calvados:hub_frames_v3_cal]  ->  design_response/design_response_summary.txt
Each --labels gives FORCE_FIELD=DESIGN_LABEL_DIR:ORIGINAL_BENDER_LABEL_DIR (the latter for the replica noise).
"""
import argparse
import csv
import math
import os
from collections import defaultdict

import numpy as np
from scipy.stats import binomtest, spearmanr, wilcoxon

HERE = os.path.dirname(os.path.abspath(__file__))
GAIVI = os.path.join(HERE, 'gaivi')
OUT = os.path.join(HERE, 'design_response')
CAL_TABLE = os.path.expanduser((__import__('os').environ['BSS_DATA'] + '/Scripts/CALVADOS IDP/input/residues_CALVADOS2.csv'))


def read_labels(d, keep=None):
    occ = defaultdict(dict)
    for r in csv.DictReader(open(os.path.join(HERE, d, 'per_residue.csv'))):
        if keep is None or r['UniProt_ID'] in keep:
            occ[r['UniProt_ID']][int(r['pos'])] = float(r['occupancy'])
    return {k: np.array([v[i] for i in range(len(v))]) for k, v in occ.items()}


def window(x, w=5):
    return np.array([x[max(0, i - w):i + w + 1].mean() for i in range(len(x))])


def centre(o):
    return float((np.arange(len(o)) * o).sum() / o.sum() / (len(o) - 1))


def entropy(o):
    q = o / o.sum()
    return float(-(q * np.log(q + 1e-12)).sum() / math.log(len(o)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--labels', action='append', required=True)
    ap.add_argument('--pred', default=os.path.join(GAIVI, 'predictions_frozen_2026-10-02.csv'))
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    lines = []

    def p(s=''):
        lines.append(s)
        print(s, flush=True)

    table = {r['id']: r for r in csv.DictReader(open(os.path.join(GAIVI, 'hub_design_table.csv')))}
    pred = defaultdict(list)
    for r in csv.DictReader(open(a.pred)):
        pred[r['id']].append(float(r['consensus']))
    pred = {k: np.array(v) for k, v in pred.items()}
    lam = {r['one']: float(r['lambdas']) for r in csv.DictReader(open(CAL_TABLE))}
    charge = {'K': 1, 'R': 1, 'D': -1, 'E': -1}
    subs = [i for i, r in table.items() if r['category'] == 'variant' and 'shuf' not in i]
    shufs = [i for i, r in table.items() if r['category'] == 'variant' and 'shuf' in i]
    reps = [i for i, r in table.items() if r['category'] == 'replica']
    p('Design batch: response of the hub network (design_response_plan.md, pre-registered 2026-10-09)')
    for spec in a.labels:
        ff, dirs = spec.split('=')
        ddir, odir = dirs.split(':')
        L = read_labels(ddir)
        orig = read_labels(odir, {i[4:] for i in reps})
        p(f'\n=== {ff}: {len(L)} design labels ===')
        # replica noise: d = run2 - run1, s2(m) = a + b m fitted to d^2 / 2
        D, M = [], []
        for i in reps:
            u = i[4:]
            if i in L and u in orig and len(L[i]) == len(orig[u]):
                D.append(L[i] - orig[u])
                M.append((L[i] + orig[u]) / 2)
        D, M = np.concatenate(D), np.concatenate(M)
        A = np.vstack([np.ones_like(M), M]).T
        coef = np.linalg.lstsq(A, D ** 2 / 2, rcond=None)[0]
        s2 = lambda m: np.maximum(coef[0] + coef[1] * m, 1e-8)  # noqa: E731
        p(f'replica noise from {len(reps)} pairs ({len(M)} residues): s2(m) = {coef[0]:.2e} + {coef[1]:.2e} m; '
          f'replica profile agreement median {np.median([spearmanr(L[i], orig[i[4:]]).correlation for i in reps if i in L and i[4:] in orig]):+.2f}')
        rows = []
        for i in sorted(subs):
            par = table[i]['parent']
            if i not in L or par not in L or len(L[i]) != len(L[par]):
                continue
            sim = L[i] - L[par]
            prd = pred[i] - pred[par]
            seq_v, seq_p = table[i]['sequence'], table[par]['sequence']
            noise = (s2(L[i]) + s2(L[par])).mean()
            R = max(0.0, 1 - noise / sim.var()) if sim.var() > 0 else 0.0
            dl = window(np.array([lam.get(c, 0.5) for c in seq_v])) - window(np.array([lam.get(c, 0.5) for c in seq_p]))
            dq = window(np.array([charge.get(c, 0) for c in seq_v], float)) - window(np.array([charge.get(c, 0) for c in seq_p], float))
            z = lambda x: (x - x.mean()) / x.std() if x.std() > 0 else x * 0  # noqa: E731
            base = z(dl) + z(dq)
            sites = np.array([x != y for x, y in zip(seq_v, seq_p)])
            dist = np.array([np.min(np.abs(np.nonzero(sites)[0] - k)) for k in range(len(sites))]) if sites.any() else None
            qwt = window(np.array([charge.get(c, 0) for c in seq_p], float))
            rows.append(dict(id=i, kind=i.split('_')[-1], rho=spearmanr(prd, sim).correlation, ceil=math.sqrt(R),
                             base=spearmanr(base, sim).correlation, dist=dist, sim=sim, prd=prd, qwt=qwt,
                             cs=(centre(L[i]) - centre(L[par]), centre(pred[i]) - centre(pred[par])),
                             es=(entropy(L[i]) - entropy(L[par]), entropy(pred[i]) - entropy(pred[par]))))
        r = np.array([x['rho'] for x in rows])
        k = int((r > 0).sum())
        pv = binomtest(k, len(r), 0.5, 'greater').pvalue
        cr = np.array([x['rho'] / x['ceil'] for x in rows if x['ceil'] > 0.3])
        p(f'PRIMARY: {len(r)} substitution variants; Spearman(predicted, simulated response) median {np.median(r):+.2f}, '
          f'positive in {k}/{len(r)} (one-sided sign test P = {pv:.2g}); median correlation / ceiling '
          f'{np.median(cr) if len(cr) else float("nan"):.2f} over {len(cr)} variants with ceiling > 0.3')
        p(f'  decision (P < 0.01 and median > 0.2): {"predictable" if pv < 0.01 and np.median(r) > 0.2 else "not established"}')
        b = np.array([x['base'] for x in rows])
        p(f'physics baseline (local stickiness + charge change): median {np.median(b):+.2f}; model vs baseline paired '
          f'Wilcoxon P = {wilcoxon(r, b).pvalue:.2g} (model better in {int((r > b).sum())}/{len(r)})')
        p('by variant type: ' + '; '.join(f'{t} {np.median([x["rho"] for x in rows if x["kind"] == t]):+.2f} '
                                         f'(n = {sum(x["kind"] == t for x in rows)})' for t in sorted({x["kind"] for x in rows})))
        for nm, lo, hi in (('sites', 0, 0), ('1-3', 1, 3), ('4-10', 4, 10), ('>10', 11, 10 ** 6)):
            s_ = [np.abs(x['sim'][(x['dist'] >= lo) & (x['dist'] <= hi)]).mean() for x in rows if x['dist'] is not None
                  and ((x['dist'] >= lo) & (x['dist'] <= hi)).any()]
            q_ = [np.abs(x['prd'][(x['dist'] >= lo) & (x['dist'] <= hi)]).mean() for x in rows if x['dist'] is not None
                  and ((x['dist'] >= lo) & (x['dist'] <= hi)).any()]
            p(f'  range {nm:5s}: mean |simulated| {np.mean(s_):.4f}, mean |predicted| {np.mean(q_):.4f}')
        ch = [spearmanr(x['qwt'], x['sim']).correlation for x in rows if x['kind'] in ('DE2NQ', 'K2R', 'R2K')]
        p(f'charge context (charge-altering variants, n = {len(ch)}): Spearman(simulated response, local net charge of '
          f'the wild type) median {np.median(ch):+.2f}')
        p(f'relocation: centre shift simulated vs predicted Spearman {spearmanr([x["cs"][0] for x in rows], [x["cs"][1] for x in rows]).correlation:+.2f}; '
          f'entropy change {spearmanr([x["es"][0] for x in rows], [x["es"][1] for x in rows]).correlation:+.2f}')
        sh = [spearmanr(pred[i], L[i]).correlation for i in shufs if i in L]
        p(f'shuffles: predicted vs simulated profile, median Spearman {np.median(sh):+.2f} over {len(sh)}')
        with open(os.path.join(OUT, f'design_response_{ff}.csv'), 'w', newline='') as fh:
            w = csv.writer(fh)
            w.writerow(['id', 'kind', 'rho_pred_sim', 'ceiling', 'rho_baseline_sim'])
            w.writerows([x['id'], x['kind'], f'{x["rho"]:.4f}', f'{x["ceil"]:.4f}', f'{x["base"]:.4f}'] for x in rows)
    open(os.path.join(OUT, 'design_response_summary.txt'), 'w').write('\n'.join(lines) + '\n')


if __name__ == '__main__':
    main()
