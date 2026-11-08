"""
degree_wetlab.py — degree_wetlab_plan.md (hash 0e53a33b). Contact degree of the PED ensembles against the frozen degree
model (hub_vs_degree/degree_model_cal.pkl, 3a77ea2a), the residue-neutral null and chain position, on the paper's sets.
Each set is labelled exactly as its hub occupancy was (same model files, subsampling to at most 1,000 models and null
seeds), with degree (mean contacts per model, Calpha < 0.8 nm, |i - j| >= 3) in place of hub occupancy:
  development 18 + fresh 8   exp_occupancy.py / null_multi.py   null seeds 50000 + 97k + j (j over the 18 then the 8)
  extra 3                    ped_extra_test.py                  null seeds 70000 + 97k + sum(ord(entry)), files in ped_local
  v2 (sets A, B, C)          ped_v2_labels.py                   pooled files, null seeds 70000 + 97k + j (j over the rows)
K = 5 null pools everywhere. Scores as score_ped_frozen.py / ped_v2_score.py.

  ~/.venvs/bender/bin/python degree_wetlab.py   ->  degree_wetlab/{per_residue.csv, scores.csv, summary.txt}
"""
import csv
import math
import os
import pickle
import sys
from collections import defaultdict
from multiprocessing import Pool

import numpy as np
from scipy.stats import binomtest, spearmanr, wilcoxon

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [HERE, os.path.join(HERE, 'hub_predictor')]
os.environ['HUB_MIN_SEP'] = '3'
import hub_frames as hf  # noqa: E402
import ped_null as pn  # noqa: E402
import train_hub_model as T  # noqa: E402
from score_ped_frozen import partial  # noqa: E402

OUT = os.path.join(HERE, 'degree_wetlab')
LOCAL = os.path.join(HERE, 'ped_local')
MODEL = os.path.join(HERE, 'hub_vs_degree', 'degree_model_cal.pkl')
K = 5
NOVEL = {'Sendai NTAIL', 'KISS-1', 'MKK7', 'E1A 36-146', 'N-WASP domain V', 'Tif2 NRID', 'AR Tau-5 R2-R3', 'Tir 388-550',
         'SARS-CoV-2 N linker'}
EXTRA = [('PED00158', 'PED00158/PED00158e006.pdb'), ('PED00001', 'PED00001/PED00001e003.pdb'),
         ('PED00489', 'PED00489/pdbfile.pdb')]


def degree(seq, models):
    n = len(seq)
    deg, used = np.zeros(n), 0
    for m in models:
        if len(m[0]) != n:
            continue
        adj = hf.contact_lists([list(m[0]), list(m[1]), list(m[2])], n)
        deg += np.array([len(a) for a in adj])
        used += 1
    return deg / used, used


def sub(models):
    return models[::math.ceil(len(models) / 1000)] if len(models) > 1000 else models


def job(spec):
    kind, entry, files, seeds = spec
    if kind == 'pooled':
        seq, models = None, []
        for f in files.split(';'):
            s, m = pn.load(f)
            seq = seq or s
            if s == seq:
                models += m
    else:
        seq, models = pn.load(files)
    models = sub(models)
    de, used = degree(seq, models)
    dn = np.mean([degree(seq, pn.null_models(seq, models, sd)[0])[0] for sd in seeds], axis=0)
    return entry, seq, de, dn, used


def main():
    os.makedirs(OUT, exist_ok=True)
    model = pickle.load(open(MODEL, 'rb'))
    specs, meta = [], {}
    ens = [(e, f) for e, f, _, _ in pn.ENSEMBLES]
    fresh = [(r['id'], r['file']) for r in csv.DictReader(open(os.path.join(HERE, 'ped_fresh_set.csv')))]
    for j, (e, f) in enumerate(ens + fresh):
        specs.append(('single', e, f, [50000 + 97 * k + j for k in range(K)]))
        meta[e] = 'used18' if j < len(ens) else 'fresh8'
    for e, f in EXTRA:
        specs.append(('single', e, os.path.join(LOCAL, f), [70000 + 97 * k + sum(map(ord, e)) for k in range(K)]))
        meta[e] = 'extra3'
    v2 = list(csv.DictReader(open(os.path.join(HERE, 'ped_v2_set.csv'))))
    v2meta = {r['id']: r for r in v2}
    for j, r in enumerate(v2):
        specs.append(('pooled', 'v2:' + r['id'], r['files'], [70000 + 97 * k + j for k in range(K)]))
        meta['v2:' + r['id']] = 'v2' + r['set']
    specs.sort(key=lambda s: -1 if s[0] == 'pooled' else 0)
    with Pool(8) as p:
        res = list(p.imap_unordered(job, specs))
    rows, scores = [], []
    for entry, seq, de, dn, used in sorted(res):
        pr = model.predict(T.features(seq, 10).astype(np.float32))
        n = len(seq)
        end = np.minimum(np.arange(n), n - 1 - np.arange(n)).astype(float)
        for i in range(n):
            rows.append([meta[entry], entry, i, seq[i], f'{de[i]:.4f}', f'{dn[i]:.4f}', f'{pr[i]:.4f}'])
        e0 = entry.replace('v2:', '')
        scores.append({'set': meta[entry], 'entry': e0, 'protein': v2meta[e0]['protein'] if entry.startswith('v2:') else '',
                       'N': n, 'models': used, 'rho_model': spearmanr(pr, de).correlation,
                       'rho_null': spearmanr(dn, de).correlation, 'rho_position': spearmanr(end, de).correlation,
                       'partial_model': partial(pr, de, end), 'partial_null': partial(dn, de, end)})
    with open(os.path.join(OUT, 'per_residue.csv'), 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['set', 'entry', 'pos', 'aa', 'deg_exp', 'deg_null_mean', 'pred_degree'])
        w.writerows(rows)
    keys = list(scores[0])
    with open(os.path.join(OUT, 'scores.csv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(scores)

    S = {s['set'] + ':' + s['entry']: s for s in scores}
    occ = {r['entry']: r for r in csv.DictReader(open(os.path.join(HERE, 'hub_predictor', 'ped_scores_all.csv')))}
    occ_v2 = {r['entry']: r for r in csv.DictReader(open(os.path.join(HERE, 'ped_v2_scores.csv')))}
    L = ['Contact degree against ensembles built from experimental data (degree_wetlab_plan.md, 0e53a33b)',
         'frozen degree model (CALVADOS-2 degree, 3a77ea2a); partial = Spearman given distance from the chain ends', '']

    def block(name, sel, occ_lookup, occ_key):
        A = [S[k] for k in sel]
        pm = np.array([a['partial_model'] for a in A])
        pnl = np.array([a['partial_null'] for a in A])
        L.append(f'{name}: {len(A)} ensembles')
        L.append(f'  median Spearman: degree model {np.median([a["rho_model"] for a in A]):+.2f}, null '
                 f'{np.median([a["rho_null"] for a in A]):+.2f}, position {np.median([a["rho_position"] for a in A]):+.2f}')
        L.append(f'  median partial: degree model {np.median(pm):+.2f}, null {np.median(pnl):+.2f}; partial > 0 in '
                 f'{int((pm > 0).sum())}/{len(pm)} (sign P {binomtest(int((pm > 0).sum()), len(pm), 0.5, alternative="greater").pvalue:.4g}); '
                 f'model > null in {int((pm > pnl).sum())}/{len(pm)} (sign P {binomtest(int((pm > pnl).sum()), len(pm), 0.5, alternative="greater").pvalue:.4g})')
        if occ_lookup is not None:
            po = np.array([float(occ_lookup[a['entry']][occ_key]) for a in A])
            d = pm - po
            L.append(f'  occupancy model partial (paper): median {np.median(po):+.2f}; degree minus occupancy: median {np.median(d):+.2f}, '
                     f'degree higher in {int((d > 0).sum())}/{len(d)} (two-sided Wilcoxon P {wilcoxon(d).pvalue:.3g})')
        return pm, pnl

    block('Development set (not blind)', [k for k in S if k.startswith('used18:')], occ, 'partial_consensus')
    pm, pnl = block('D1  Blind set (8 fresh + 3 added)', [k for k in S if k.startswith(('fresh8:', 'extra3:'))], occ,
                    'partial_consensus')
    d1 = binomtest(int((pm > 0).sum()), len(pm), 0.5, alternative='greater').pvalue < 0.05 and \
        binomtest(int((pm > pnl).sum()), len(pm), 0.5, alternative='greater').pvalue < 0.05
    L.append(f'  D1 decision (partial > 0 and model > null, both sign P < 0.05): {"MET" if d1 else "NOT MET"}')
    L.append('')
    # second blind set: set A without IAPP, one value per protein (mean over entries)
    byprot = defaultdict(list)
    for k, s in S.items():
        if s['set'] == 'v2A' and s['protein'] != 'IAPP':
            byprot[s['protein']].append(s)
    prots = sorted(byprot)
    pm = np.array([np.mean([s['partial_model'] for s in byprot[p]]) for p in prots])
    pnl = np.array([np.mean([s['partial_null'] for s in byprot[p]]) for p in prots])
    po = np.array([np.mean([float(occ_v2[s['entry']]['partial_model']) for s in byprot[p]]) for p in prots])
    sp = binomtest(int((pm > 0).sum()), len(pm), 0.5, alternative='greater').pvalue
    sn = binomtest(int((pm > pnl).sum()), len(pm), 0.5, alternative='greater').pvalue
    L.append(f'D2  Second blind set (set A without IAPP): {len(prots)} proteins')
    L.append(f'  median partial: degree model {np.median(pm):+.2f}, null {np.median(pnl):+.2f}; partial > 0 in {int((pm > 0).sum())}/{len(pm)} '
             f'(P {sp:.4g}); model > null in {int((pm > pnl).sum())}/{len(pm)} (P {sn:.4g})')
    L.append(f'  occupancy model partial (paper): median {np.median(po):+.2f} (model > null was 10/14, criterion not met); degree minus '
             f'occupancy median {np.median(pm - po):+.2f}, degree higher in {int((pm > po).sum())}/{len(pm)}')
    L.append(f'  D2 decision (both sign P < 0.05): {"MET" if (sp < 0.05 and sn < 0.05) else "NOT MET"}')
    nov = [i for i, p in enumerate(prots) if p in NOVEL]
    L.append(f'  sensitivity, {len(nov)} novel proteins: partial > 0 in {int((pm[nov] > 0).sum())}/{len(nov)}, model > null in '
             f'{int((pm[nov] > pnl[nov]).sum())}/{len(nov)}, median partial {np.median(pm[nov]):+.2f}')
    L.append('')
    B = [s for s in scores if s['set'] == 'v2B']
    L.append('Set B (per entry): ' + ', '.join(f'{s["entry"]} {s["partial_model"]:+.2f} (null {s["partial_null"]:+.2f})' for s in B))
    fus = S.get('used18:PED00494')
    if fus:
        L.append(f'FUS NTD (PED00494): degree model rho {fus["rho_model"]:+.2f}, partial {fus["partial_model"]:+.2f}; occupancy model '
                 f'rho {float(occ["PED00494"]["rho_consensus"]):+.2f}, partial {float(occ["PED00494"]["partial_consensus"]):+.2f}')
    txt = '\n'.join(L) + '\n'
    open(os.path.join(OUT, 'summary.txt'), 'w').write(txt)
    print(txt)


if __name__ == '__main__':
    main()
