"""
analyze_glycine.py — the pre-registered analysis of the glycine sub-batch (glycine_subbatch_plan.md). Written and hashed
before any of these sequences was simulated; do not edit after the labels exist (record amendments separately).

Inputs: hub_design_glycine_table.csv (design key); per_residue.csv label folders from hub_labels_gaivi.py for the variants
in each force field; the natural (WT) references, primarily same-batch WT runs (hub_design_glycine_wt_300K.fasta,
Amendment 1) and otherwise the existing v3 labels (CALVADOS-2 hub_frames_v3_cal, Mpipi-Recharged hub_frames_v3_recharged,
Mpipi-GG hub_frames_v3_mpipi); frozen predictions (predictions_frozen_2026-10-02_glycine_2026-10-04.csv,
predictions_v3_2026-10-04_glycine_2026-10-04.csv).
  ~/.venvs/bender/bin/python analyze_glycine.py --cal DIR --rc DIR [--gg DIR] [--cal-wt DIR --rc-wt DIR --gg-wt DIR]
      ->  glycine_results.txt, glycine_deltas.csv
Amendment 1 (2026-10-04, before any of these sequences was simulated): same-batch WT references as primary, the Mpipi-GG
fallback is hub_frames_v3_mpipi (labels after the Oct 3 tie-breaking fix), and H5 is reported for the matched
model / force-field pairs of the plan only.
"""
import argparse
import csv
import os
import re
from collections import defaultdict

import numpy as np
from scipy.stats import spearmanr, wilcoxon

HERE = os.path.dirname(os.path.abspath(__file__))
EV = os.path.dirname(HERE)
FF = {'cal': 'CALVADOS-2', 'gg': 'Mpipi-GG', 'rc': 'Mpipi-Recharged'}
WT_DIRS = {'cal': ['hub_frames_v3_cal'], 'rc': ['hub_frames_v3_recharged'], 'gg': ['hub_frames_v3_mpipi']}
MATCHED = {'cal': ('cal', 'consensus', 'cal_full', 'cal_context'), 'rc': ('consensus', 'rc_full', 'rc_context'),
           'gg': ('consensus', 'gg_full', 'gg_context')}      # H5 model / force-field pairs of the plan
out = []


def p(s=''):
    out.append(s)
    print(s, flush=True)


def read_labels(path):
    """{id: {pos: occupancy}} for proteins whose rows are complete (per_protein N matches)."""
    if not os.path.isabs(path):                         # label folders next to this script, references in redo_evidence
        path = next((os.path.join(b, path) for b in (HERE, EV) if os.path.isdir(os.path.join(b, path))),
                    os.path.join(EV, path))
    done = {r['UniProt_ID']: int(r['N']) for r in csv.DictReader(open(os.path.join(path, 'per_protein.csv')))}
    d = defaultdict(dict)
    with open(os.path.join(path, 'per_residue.csv')) as fh:
        for r in csv.DictReader(fh):
            if r['UniProt_ID'] in done and r.get('occupancy', '') != '':
                d[r['UniProt_ID']][int(r['pos'])] = float(r['occupancy'])
    return {u: v for u, v in d.items() if len(v) == done[u]}


def key():
    info = {}
    for r in csv.DictReader(open(os.path.join(HERE, 'hub_design_glycine_table.csv'))):
        b = info.setdefault(r['parent'], {'n': int(r['n'])})
        b['h'] = int(re.search(r'hub [A-Z](\d+)', r['notes']).group(1))
        c = re.search(r'centred at (\d+)', r['notes'])
        if c:
            b['c'] = int(c.group(1))
        b.setdefault('variants', []).append(r['guest'])
    return info


def window(prof, centre):
    return [i for i in range(centre - 2, centre + 3) if i in prof]


def delta(var, wt, centre):
    idx = window(wt, centre)
    return float(np.mean([var[i] - wt[i] for i in idx]))


def predictions():
    pr = defaultdict(lambda: defaultdict(dict))
    for f, cols in (('predictions_frozen_2026-10-02_glycine_2026-10-04.csv', ('cal', 'consensus')),
                    ('predictions_v3_2026-10-04_glycine_2026-10-04.csv',
                     ('cal_full', 'cal_context', 'gg_full', 'gg_context', 'rc_full', 'rc_context'))):
        path = os.path.join(HERE, f)
        if not os.path.exists(path):
            continue
        for r in csv.DictReader(open(path)):
            for c in cols:
                pr[c][r['id']][int(r['pos'])] = float(r[c])
    return pr


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--cal', required=True)
    ap.add_argument('--rc', required=True)
    ap.add_argument('--gg')
    ap.add_argument('--cal-wt')
    ap.add_argument('--rc-wt')
    ap.add_argument('--gg-wt')
    a = ap.parse_args()
    same_batch = {'cal': a.cal_wt, 'rc': a.rc_wt, 'gg': a.gg_wt}
    info = key()
    pr = predictions()
    rows = []
    obs = {}
    for ff, vdir in (('cal', a.cal), ('rc', a.rc), ('gg', a.gg)):
        if not vdir:
            continue
        lab = read_labels(vdir)
        wt, n_same = {}, 0
        for d in WT_DIRS[ff] + ([same_batch[ff]] if same_batch[ff] else []):     # same-batch runs, read last, take precedence
            for u, v in read_labels(d).items():
                if u.startswith('GLY_') and u.endswith('_WT'):
                    n_same += 1
                wt[u.replace('GLY_', '').replace('_WT', '')] = v
        p(f'{FF[ff]}: {n_same} same-batch WT references, the rest from {WT_DIRS[ff][0]}')
        for u, b in info.items():
            if u not in wt or len(wt[u]) != b['n']:
                continue
            for var in b['variants']:
                vid = f'GLY_{u}_{var}'
                if vid not in lab or len(lab[vid]) != b['n']:
                    continue
                dh = delta(lab[vid], wt[u], b['h'])
                dc = delta(lab[vid], wt[u], b['c']) if 'c' in b else np.nan
                obs[(ff, u, var)] = (dh, dc)
                row = {'force_field': FF[ff], 'base': u, 'variant': var, 'delta_hub': dh, 'delta_far': dc}
                for m in pr:
                    wtid = f'GLY_{u}_WT'
                    if vid in pr[m] and wtid in pr[m]:
                        row[f'pred_{m}'] = delta(pr[m][vid], pr[m][wtid], b['h'])
                rows.append(row)
    p('Glycine sub-batch, pre-registered analysis (glycine_subbatch_plan.md)')
    for ff in [f for f in ('cal', 'rc', 'gg') if any(k[0] == f for k in obs)]:
        get = lambda var, i=0: {u: obs[(ff, u, var)][i] for (f, u, v) in obs if f == ff and v == var}
        near, far = get('G2S_near'), get('G2S_far')
        both = sorted(set(near) & set(far))
        p(f'\n{FF[ff]}  ({len(both)} bases with both G2S variants)')
        if both:
            d = np.array([near[u] - far[u] for u in both])
            p(f'  H1  near - far: median {np.median(d):+.4f}, negative in {np.mean(d < 0):.0%}, '
              f'Wilcoxon one-sided p = {wilcoxon(d, alternative="less").pvalue:.4g}')
        mv_h, mv_c = get('Gmove'), get('Gmove', 1)
        if mv_h:
            vh, vc = np.array(list(mv_h.values())), np.array(list(mv_c.values()))
            p(f'  H2  Gmove hub window: median {np.median(vh):+.4f}, p(<0) = {wilcoxon(vh, alternative="less").pvalue:.4g};'
              f' far window: median {np.median(vc):+.4f}, p(>0) = {wilcoxon(vc, alternative="greater").pvalue:.4g}')
        ga = np.array(list(get('G2A_near').values()))
        if len(ga):
            p(f'  H3  G2A_near: median {np.median(ga):+.4f}, p(<0) = {wilcoxon(ga, alternative="less").pvalue:.4g}')
        for m in [m for m in MATCHED[ff] if m in pr]:
            pairs = [(r[f'pred_{m}'], r['delta_hub']) for r in rows if r['force_field'] == FF[ff] and f'pred_{m}' in r]
            if len(pairs) > 5:
                x, y = zip(*pairs)
                rho, pv = spearmanr(x, y)
                pv1 = pv / 2 if rho > 0 else 1 - pv / 2
                agree = np.mean([np.sign(i) == np.sign(j) for i, j in pairs])
                p(f'  H5  {m:12s} predicted vs observed delta_hub: Spearman {rho:+.2f} (one-sided p = {pv1:.3g}), '
                  f'sign agreement {agree:.0%}, n = {len(pairs)}')
    if any(k[0] == 'gg' for k in obs):
        pg = {u: v for (f, u, var), (v, _) in obs.items() if f == 'gg' and var == 'P_near'}
        for other in ('cal', 'rc'):
            po = {u: v for (f, u, var), (v, _) in obs.items() if f == other and var == 'P_near'}
            common = sorted(set(pg) & set(po))
            if len(common) >= 5:
                d = np.array([pg[u] - po[u] for u in common])
                p(f'\n  H4  P_near Mpipi-GG - {FF[other]}: median {np.median(d):+.4f}, n = {len(common)}, '
                  f'p(<0) = {wilcoxon(d, alternative="less").pvalue:.4g}')
    with open(os.path.join(HERE, 'glycine_deltas.csv'), 'w', newline='') as fh:
        cols = sorted({k for r in rows for k in r}, key=lambda k: (k.startswith('pred_'), k))
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)
    open(os.path.join(HERE, 'glycine_results.txt'), 'w').write('\n'.join(out) + '\n')


if __name__ == '__main__':
    main()
