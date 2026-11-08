"""
llps_csat.py — implements llps_csat_plan.md (hashed 2026-10-09, Amendment 1 hashed before any feature was computed):
do predicted hub features of the frozen consensus model explain the measured saturation concentration (csat) of
hnRNPA1 LCD variants beyond composition?

Inputs (external_data):
  llps_pappulab_plcd/Phase_Separation_Data.xlsx   sheet "non averaged data", csat (uM) at 4 C, replicates
  llps_bremer2022_si/bremer2022_table_s1.fasta    sequences from Supplementary Table 1 of Bremer et al. 2022
Features (normalised predicted occupancy q = p / sum p):
  F1  sticker hub share     sum of q over F, W, Y, R                                   (primary)
  F2  hub spread            entropy of q / log N
  F3  hub spacing           coefficient of variation of the gaps between consecutive top-10% residues
Baseline: log10 csat ~ 1 + n(F, W, Y, R) + NCPR (least squares). Primary test: Spearman(F1, baseline residual),
two-sided permutation P (10,000 permutations of the residuals), P < 0.05.
  OMP_NUM_THREADS=1 ~/.venvs/bender/bin/python llps_csat.py  ->  llps_csat/llps_csat_summary.txt, llps_features.csv
"""
import csv
import math
import os
import pickle
import re
import sys
import xml.etree.ElementTree as ET
import zipfile
from collections import Counter

import numpy as np
from scipy.stats import spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, 'hub_predictor'))
import train_hub_model as T  # noqa: E402

EXT = os.path.join(HERE, '..', 'external_data')
XLSX = os.path.join(EXT, 'llps_pappulab_plcd', 'Phase_Separation_Data.xlsx')
FASTA = os.path.join(EXT, 'llps_bremer2022_si', 'bremer2022_table_s1.fasta')
OUT = os.path.join(HERE, 'llps_csat')
NS = {'m': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
STICK = set('FWYR')


def xlsx_rows(path, sheet):
    """Rows of one sheet as {column letter: value} (stdlib only)."""
    z = zipfile.ZipFile(path)
    ss = [''.join(t.text or '' for t in si.iter('{%s}t' % NS['m']))
          for si in ET.fromstring(z.read('xl/sharedStrings.xml')).findall('m:si', NS)]
    wb = ET.fromstring(z.read('xl/workbook.xml'))
    rels = {r.get('Id'): r.get('Target') for r in ET.fromstring(z.read('xl/_rels/workbook.xml.rels'))}
    rid = [s.get('{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id')
           for s in wb.find('m:sheets', NS) if s.get('name').strip() == sheet][0]
    t = rels[rid].lstrip('/')
    t = t if t.startswith('xl/') else 'xl/' + t
    rows = []
    for row in ET.fromstring(z.read(t)).find('m:sheetData', NS).findall('m:row', NS):
        d = {}
        for c in row.findall('m:c', NS):
            v = c.find('m:v', NS)
            if v is None:
                continue
            val = ss[int(v.text)] if c.get('t') == 's' else v.text
            d[re.match(r'[A-Z]+', c.get('r')).group()] = val
        rows.append(d)
    return rows


def csat_4c():
    """Mean csat (M) at 4 C per construct, from the first block of the sheet (rows up to the first empty name)."""
    out, started = {}, False
    for d in xlsx_rows(XLSX, 'non averaged data'):
        name = (d.get('B') or '').strip()
        if name == 'Construct':
            if started:
                break                      # second block (20 C)
            started = True
            continue
        if started and name:
            vals = [float(d[c]) for c in 'CDEFGH' if c in d and d[c] not in ('', None)]
            if vals:
                out[name] = float(np.mean(vals)) * 1e-6
            else:
                out[name] = None
        elif started and not name and out:
            break
    return out


def fasta(path):
    seqs, cur = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith('>'):
            cur = line[1:]
            seqs[cur] = ''
        elif cur:
            seqs[cur] += line
    return seqs


def features(q, seq):
    n = len(seq)
    f1 = float(sum(q[i] for i, c in enumerate(seq) if c in STICK))
    f2 = float(-(q * np.log(q + 1e-12)).sum() / math.log(n))
    top = np.sort(np.argsort(-q, kind='stable')[:math.ceil(0.1 * n)])
    gaps = np.diff(top)
    f3 = float(gaps.std() / gaps.mean())
    return f1, f2, f3


def perm_p(x, y, obs, n=10000, seed=0):
    rng = np.random.default_rng(seed)
    k = sum(abs(spearmanr(x, rng.permutation(y)).correlation) >= abs(obs) for _ in range(n))
    return (k + 1) / (n + 1)


def holm(ps):
    order = np.argsort(ps)
    adj, m, run = np.empty(len(ps)), len(ps), 0.0
    for r, i in enumerate(order):
        run = max(run, min(1.0, (m - r) * ps[i]))
        adj[i] = run
    return adj


def main():
    os.makedirs(OUT, exist_ok=True)
    lines = []

    def p(s=''):
        lines.append(s)
        print(s, flush=True)
    cs = csat_4c()
    seqs = fasta(FASTA)
    model = pickle.load(open(os.path.join(HERE, 'hub_predictor', 'frozen_models_2026-10-02.pkl'), 'rb'))['consensus']
    rows = []
    for name, c in cs.items():
        if name not in seqs:
            p(f'  {name}: no sequence in Table S1')
            continue
        s = seqs[name]
        pr = np.clip(model.predict(T.features(s, 10)), 0, None)
        q = pr / pr.sum()
        f1, f2, f3 = features(q, s)
        cnt = Counter(s)
        rows.append(dict(name=name, seq=s, csat=c, y=math.log10(c) if c else None, f1=f1, f2=f2, f3=f3,
                         nstick=sum(cnt[a] for a in STICK),
                         ncpr=(cnt['K'] + cnt['R'] - cnt['D'] - cnt['E']) / len(s),
                         comp=tuple(sorted(cnt.items()))))
    test = [r for r in rows if r['y'] is not None]
    p('LLPS: predicted hub features and measured csat of A1-LCD variants (llps_csat_plan.md, Amendment 1)')
    p(f'{len(rows)} constructs in the 4 C block with a sequence; {len(test)} with a csat value (tested)')
    for r in rows:
        if r['y'] is None:
            p(f'  not tested (no 4 C value): {r["name"]}')
    y = np.array([r['y'] for r in test])
    X = np.vstack([np.ones(len(test)), [r['nstick'] for r in test], [r['ncpr'] for r in test]]).T
    beta = np.linalg.lstsq(X, y, rcond=None)[0]
    res = y - X @ beta
    r2 = 1 - res.var() / y.var()
    p(f'\nBaseline log10 csat = {beta[0]:+.3f} {beta[1]:+.4f} n(FWYR) {beta[2]:+.3f} NCPR;  R^2 = {r2:.2f}')
    F = {k: np.array([r[k] for r in test]) for k in ('f1', 'f2', 'f3')}
    obs = spearmanr(F['f1'], res).correlation
    pp = perm_p(F['f1'], res, obs)
    p(f'\nPRIMARY: Spearman(F1 sticker hub share, baseline residual) = {obs:+.2f}, two-sided permutation P = {pp:.3g} '
      f'(n = {len(test)})')
    p(f'  decision (P < 0.05): {"hub arrangement adds information beyond composition" if pp < 0.05 else "not established"}')
    p('\nRaw Spearman with log10 csat (no decision):')
    for k, lab in (('f1', 'F1 sticker hub share'), ('f2', 'F2 hub spread'), ('f3', 'F3 hub spacing'),
                   ('nstick', 'n(FWYR) (composition)'), ('ncpr', 'NCPR (composition)')):
        v = np.array([r[k] for r in test], float)
        p(f'  {lab:24s} {spearmanr(v, y).correlation:+.2f}')
    p('\nSecondary: F2 and F3 against the baseline residual (Holm-corrected)')
    sec = [(k, spearmanr(F[k], res).correlation) for k in ('f2', 'f3')]
    ps = np.array([perm_p(F[k], res, r, seed=1) for k, r in sec])
    for (k, r), q_, a in zip(sec, ps, holm(ps)):
        p(f'  {k}: Spearman {r:+.2f}, permutation P {q_:.3g}, Holm {a:.3g}')
    err0, err1 = [], []
    for i in range(len(test)):
        m = np.arange(len(test)) != i
        b0 = np.linalg.lstsq(X[m], y[m], rcond=None)[0]
        X1 = np.column_stack([X, F['f1']])
        b1 = np.linalg.lstsq(X1[m], y[m], rcond=None)[0]
        err0.append(abs(y[i] - X[i] @ b0))
        err1.append(abs(y[i] - X1[i] @ b1))
    p(f'  leave-one-out median absolute error of log10 csat: baseline {np.median(err0):.3f}, baseline + F1 '
      f'{np.median(err1):.3f}')
    comps = Counter(r['comp'] for r in test)
    p(f'  same-composition pairs among tested constructs: {sum(v * (v - 1) // 2 for v in comps.values())}')
    p('\nPer construct')
    p(f'  {"construct":24s} {"log10 csat":>10s} {"n(FWYR)":>7s} {"NCPR":>6s} {"F1":>6s} {"F2":>6s} {"F3":>6s} {"resid":>6s}')
    for r, e in zip(test, res):
        p(f'  {r["name"]:24s} {r["y"]:10.2f} {r["nstick"]:7d} {r["ncpr"]:+6.3f} {r["f1"]:6.3f} {r["f2"]:6.3f} '
          f'{r["f3"]:6.2f} {e:+6.2f}')
    with open(os.path.join(OUT, 'llps_features.csv'), 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['construct', 'csat_M_4C', 'n_FWYR', 'NCPR', 'F1_sticker_hub_share', 'F2_hub_spread', 'F3_hub_spacing'])
        for r in rows:
            w.writerow([r['name'], r['csat'], r['nstick'], f'{r["ncpr"]:.4f}', f'{r["f1"]:.4f}', f'{r["f2"]:.4f}',
                        f'{r["f3"]:.3f}'])
    open(os.path.join(OUT, 'llps_csat_summary.txt'), 'w').write('\n'.join(lines) + '\n')


if __name__ == '__main__':
    main()
