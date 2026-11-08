"""
ped_blind11_report.py — one place for every number the paper quotes about the PED ensembles.
  Development set: 18 ensembles used while the method was built (not blind).
  Blind set: 11 ensembles whose predictions were frozen and hashed before their labels were computed, namely the 8 fresh
    ones of 2026-10-02 and the 3 added on 2026-10-03 (PED00158, PED00001, PED00489), which were scored only after the
    partial correlation had been fixed as a co-primary metric.
  Ceiling: agreement in hub occupancy between independent PED ensembles of the same protein (drkN x 4, Sic1/pSic1 x 5).
    If each ensemble is the same true profile plus independent error, a perfect predictor would correlate with a single
    ensemble at about the square root of that agreement.
Inputs: hub_predictor/ped_scores_frozen.csv (score_ped_frozen.py; 18 + 8), hub_predictor/ped_extra_occupancy.csv and
  gaivi/predictions_frozen_2026-10-02_ped_extra.csv (ped_extra_test.py; the 3 new ensembles and the ceiling groups).
Outputs: hub_predictor/ped_scores_all.csv, hub_predictor/ped_ceiling_pairs.csv, ped_blind11_summary.txt and
  ../methods_paper/ped_table_s3_rows.tex (rows for Table S3).
  ~/.venvs/bender/bin/python ped_blind11_report.py
"""
import csv
import hashlib
import itertools
import math
import os
from collections import defaultdict

import numpy as np
from scipy.stats import binom, rankdata, spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
HP = os.path.join(HERE, 'hub_predictor')
FROZEN_EXTRA = os.path.join(HERE, 'gaivi', 'predictions_frozen_2026-10-02_ped_extra.csv')
NEW = ('PED00158', 'PED00001', 'PED00489')
RAW_MD = {'PED00533', 'PED00536'}
GROUPS = {'drkN': ['PED00022', 'PED00156', 'PED00157', 'PED00158'],
          'Sic1': ['PED00023', 'PED00159', 'PED00160', 'PED00161', 'PED00001']}
# names and ensemble methods as printed in Table S3. The method is the conformer generator and the selection or
# refinement step that the depositors report in PED (generation-method annotations and entry titles; for MKK4 and
# CBP ID4 also the depositors' description and the source paper). Data types are named only where the entry title states
# them. The earlier 'restraints' labels were PED's list of everything measured (e.g. heteronuclear relaxation NOE), not
# the data the ensembles were fitted to (2026-10-03).
# FM, flexible-meccano; MD, molecular dynamics; BME, Bayesian maximum entropy.
LATEX = {
    'PED00003': (r'$\beta$-synuclein', 'MD, PRE-restrained'), 'PED00004': ('p27 KID', 'MD'),
    'PED00006': (r'$\alpha$/$\beta$-synuclein hybrid', 'MD, PRE-restrained'), 'PED00013': ('PaaA2', 'ARIA, EOM'),
    'PED00022': ('drkN', 'TraDES, ENSEMBLE'), 'PED00023': ('Sic1', 'TraDES, ENSEMBLE'),
    'PED00024': (r'$\alpha$-synuclein', 'MD, PRE-restrained'),
    'PED00156': ('drkN', 'X-EISD, TraDES pool'), 'PED00157': ('drkN', 'X-EISD, ENSEMBLE pool'),
    'PED00159': ('Sic1', 'ENSEMBLE, PRE and SAXS'), 'PED00160': ('Sic1', 'ENSEMBLE, PRE, SAXS and CS'),
    'PED00161': ('pSic1', 'ENSEMBLE, PRE and SAXS'),
    'PED00225': ('eIF4G1 (1--249)', 'CYANA, EOM'), 'PED00494': ('FUS NTD', 'DEER-restrained'),
    'PED00498': ('V2R C-tail', 'FM, RDC-refined'), 'PED00499': ('GHSR C-tail', 'FM, RDC-refined'),
    'PED00500': (r'$\beta_2$AR C-tail', 'FM, RDC-refined'), 'PED00511': (r'PDE$\gamma$', 'X-PLOR'),
    'PED00017': ('tau fragment', 'FM, ASTEROIDS'), 'PED00154': ('MKK4 N-terminal region', 'FM, ASTEROIDS'),
    'PED00216': ('CBP ID4', 'FM, no selection'), 'PED00501': ('V2R C-tail, phosphomimetic', 'FM, RDC-refined'),
    'PED00502': ('GHSR C-tail, phosphomimetic', 'FM, RDC-refined'),
    'PED00503': (r'$\beta_2$AR C-tail, phosphomimetic', 'FM, RDC-refined'),
    'PED00533': (r'A$\beta$40', 'MD, a99SB-disp'), 'PED00536': ('ACTR', 'MD, a99SB-disp'),
    'PED00158': ('drkN, unfolded state', 'X-EISD, mixed pool'), 'PED00001': ('pSic1, 11 models', 'TraDES, ENSEMBLE'),
    'PED00489': ('4E-BP2, 5-phosphorylated', 'FastFloppyTail, BME')}
NOT_SELECTED = {'PED00216'}   # CBP ID4: flexible-meccano with chemical-shift propensities, no selection (Piai 2016)
out = []


def p(s=''):
    out.append(s)
    print(s, flush=True)


def partial(x, y):
    n = len(y)
    z = np.minimum(np.arange(n), n - 1 - np.arange(n)).astype(float)
    rx, ry, rz = (rankdata(v) for v in (x, y, z))
    A = np.vstack([rz, np.ones_like(rz)]).T
    ex = rx - A @ np.linalg.lstsq(A, rx, rcond=None)[0]
    ey = ry - A @ np.linalg.lstsq(A, ry, rcond=None)[0]
    return float(np.corrcoef(ex, ey)[0, 1])


def fmt(v):
    """Two decimals, with negative values in math mode so that the minus sign is typeset correctly."""
    return f'${v:.2f}$' if v < 0 else f'{v:.2f}'


def sign_p(k, n):
    """One-sided sign test, P(at least k successes of n) under p = 0.5."""
    return float(binom.sf(k - 1, n, 0.5))


def main():
    h = hashlib.sha256(open(FROZEN_EXTRA, 'rb').read()).hexdigest()
    if h not in open(os.path.join(HERE, 'gaivi', 'FROZEN.sha256')).read():
        raise SystemExit('frozen predictions for the extra ensembles are not registered')
    rows = []
    for r in csv.DictReader(open(os.path.join(HP, 'ped_scores_frozen.csv'))):
        rows.append({k: (float(v) if k not in ('set', 'entry', 'protein', 'restraints') and v != '' else v)
                     for k, v in r.items()})
    occ = defaultdict(list)
    for r in csv.DictReader(open(os.path.join(HP, 'ped_extra_occupancy.csv'))):
        occ[r['entry']].append(r)
    pred = defaultdict(list)
    for r in csv.DictReader(open(FROZEN_EXTRA)):
        pred[r['id']].append((float(r['cal']), float(r['pooled']), float(r['consensus'])))
    for e in NEW:
        oe = np.array([float(x['occ_exp']) for x in occ[e]])
        on = np.array([float(x['occ_null_mean5']) for x in occ[e]])
        pr = np.array(pred[e])
        n = len(oe)
        assert pr.shape == (n, 3), (e, pr.shape, n)
        end = np.minimum(np.arange(n), n - 1 - np.arange(n)).astype(float)
        row = {'set': 'extra3', 'entry': e, 'protein': LATEX[e][0], 'restraints': LATEX[e][1], 'N': n}
        for j, k in enumerate(('cal', 'pooled', 'consensus')):
            row[f'rho_{k}'] = spearmanr(pr[:, j], oe).correlation
            row[f'partial_{k}'] = partial(pr[:, j], oe)
        row.update(rho_null=spearmanr(on, oe).correlation, partial_null=partial(on, oe),
                   rho_position=spearmanr(end, oe).correlation, partial_position='')
        rows.append(row)
    cols = list(rows[0].keys())
    with open(os.path.join(HP, 'ped_scores_all.csv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow({k: (f'{v:.4f}' if isinstance(v, float) else v) for k, v in r.items()})

    # ceiling
    pairs = []
    for g, ents in GROUPS.items():
        for a, b in itertools.combinations(ents, 2):
            ra = np.array([float(x['occ_exp']) for x in occ[a]])
            rb = np.array([float(x['occ_exp']) for x in occ[b]])
            pairs.append((g, a, b, spearmanr(ra, rb).correlation))
    with open(os.path.join(HP, 'ped_ceiling_pairs.csv'), 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['protein', 'entry_a', 'entry_b', 'spearman'])
        w.writerows([(g, a, b, f'{r:.4f}') for g, a, b, r in pairs])
    allr = [r for *_, r in pairs]
    med = float(np.median(allr))
    p('Experimental ceiling (hub occupancy agreement between independent PED ensembles of the same protein)')
    for g in GROUPS:
        v = [r for gg, *_, r in pairs if gg == g]
        p(f'  {g}: {len(v)} pairs, median {np.median(v):+.2f} (range {min(v):+.2f} to {max(v):+.2f})')
    p(f'  all: {len(allr)} pairs, median {med:+.2f}; a perfect predictor would reach about sqrt({med:.2f}) = '
      f'{math.sqrt(max(med, 0)):.2f} against a single ensemble')

    def summ(label, S):
        md = lambda k: float(np.median([s[k] for s in S]))
        n = len(S)
        k0 = sum(s['partial_consensus'] > 0 for s in S)
        kn = sum(s['partial_consensus'] > s['partial_null'] for s in S)
        p(f'\n{label} (n = {n})')
        p(f'  raw Spearman, median: model {md("rho_consensus"):+.2f}, null {md("rho_null"):+.2f}, '
          f'position {md("rho_position"):+.2f}; model > null in {sum(s["rho_consensus"] > s["rho_null"] for s in S)}/{n}, '
          f'model > position in {sum(s["rho_consensus"] > s["rho_position"] for s in S)}/{n}')
        p(f'  partial on end distance, median: model {md("partial_consensus"):+.2f}, null {md("partial_null"):+.2f}; '
          f'model > 0 in {k0}/{n} (sign test p = {sign_p(k0, n):.4f}), model > null in {kn}/{n} (p = {sign_p(kn, n):.4f})')
        return md

    dev = [r for r in rows if r['set'] == 'used18']
    fresh = [r for r in rows if r['set'] == 'fresh8']
    extra = [r for r in rows if r['set'] == 'extra3']
    blind = fresh + extra
    summ('Development set, not blind', dev)
    summ('Blind, fresh 8 (2026-10-02)', fresh)
    summ('Blind, 3 scored after the partial metric was fixed (2026-10-03)', extra)
    summ('Blind, all 11', blind)
    summ('Blind, the 9 built from experimental data (without the two unweighted MD ensembles)',
         [r for r in blind if r['entry'] not in RAW_MD])
    summ('Blind, the 8 selected or refined against experimental data (also without CBP ID4, flexible-meccano only)',
         [r for r in blind if r['entry'] not in RAW_MD | NOT_SELECTED])
    by = {r['entry']: r for r in rows}
    for g, ents in GROUPS.items():
        v = [by[e]['rho_consensus'] for e in ents]
        p(f'\n{g} family ({len(ents)} ensembles): model median {np.median(v):+.2f}, null '
          f'{np.median([by[e]["rho_null"] for e in ents]):+.2f}, position {np.median([by[e]["rho_position"] for e in ents]):+.2f}, '
          f'agreement between its ensembles {np.median([r for gg, *_, r in pairs if gg == g]):+.2f}')
    open(os.path.join(HERE, 'ped_blind11_summary.txt'), 'w').write('\n'.join(out) + '\n')

    # Table S3 rows
    tex = [r'\multicolumn{9}{l}{\emph{Development set}} \\']
    tex += [f'{r["entry"]} & {LATEX[r["entry"]][0]} & {LATEX[r["entry"]][1]} & {int(r["N"])} & '
            f'{fmt(r["rho_consensus"])} & {fmt(r["rho_null"])} & {fmt(r["rho_position"])} & '
            f'{fmt(r["partial_consensus"])} & {fmt(r["partial_null"])} \\\\' for r in dev]
    tex += [r'\midrule', r'\multicolumn{9}{l}{\emph{Blind set, predictions frozen before labels}} \\']
    for r in blind:
        mark = r'$^{a}$' if r['entry'] in NEW else ''
        tex.append(f'{r["entry"]}{mark} & {LATEX[r["entry"]][0]} & {LATEX[r["entry"]][1]} & {int(r["N"])} & '
                   f'{fmt(r["rho_consensus"])} & {fmt(r["rho_null"])} & {fmt(r["rho_position"])} & '
                   f'{fmt(r["partial_consensus"])} & {fmt(r["partial_null"])} \\\\')
    path = os.path.join(os.path.dirname(HERE), 'methods_paper', 'ped_table_s3_rows.tex')
    open(path, 'w').write('\n'.join(tex) + '\n')
    print(f'\nTable S3 rows -> {path}')
    # ceiling table rows: protein, the two ensembles with their restraints, agreement
    tex = []
    for g in GROUPS:
        tex.append(r'\multicolumn{3}{l}{\emph{' + g + r'}} \\')
        for gg, a, b, r in pairs:
            if gg == g:
                tex.append(f'{a} ({LATEX[a][1]}) & {b} ({LATEX[b][1]}) & {fmt(r)} \\\\')
    path = os.path.join(os.path.dirname(HERE), 'methods_paper', 'ped_table_ceiling_rows.tex')
    open(path, 'w').write('\n'.join(tex) + '\n')
    print(f'ceiling table rows -> {path}')


if __name__ == '__main__':
    main()
