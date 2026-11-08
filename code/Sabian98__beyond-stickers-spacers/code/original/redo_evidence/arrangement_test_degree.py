"""
arrangement_test_degree.py — the arrangement test of arrangement_test.py (plan c4ca3409) with contact degree in place of
hub occupancy; plan arrangement_test_degree_plan.md (hashed with this script before any degree value of the shuffled or
Y-to-F labels was read). D1-D3 repeat H1-H3 with the same functions and decision rule; D4 and D5 are descriptive.
  ~/.venvs/bender/bin/python -I arrangement_test_degree.py  ->  arrangement_test_degree_summary.txt,
                                                                arrangement_test_degree_values.csv (per-protein agreements)
"""
import csv
import json
import os
import pickle
import sys
from collections import defaultdict

import numpy as np
from scipy.stats import pearsonr, spearmanr, wilcoxon

HERE = '/Volumes/User Homes/Nu Project/BENDER-DYN/redo_evidence'
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, 'hub_predictor'))
import arrangement_test as at  # noqa: E402
import train_hub_model as T  # noqa: E402

CAL_TABLE = os.path.expanduser('~/Downloads/Scripts/CALVADOS IDP/input/residues_CALVADOS2.csv')
AAS = 'ACDEFGHIKLMNPQRSTVWY'
out = []


def p(s=''):
    out.append(s)
    print(s, flush=True)


def read(dirname, col, keep=None):
    """{id: (sequence, values)} for proteins whose residue rows are complete (as arrangement_test.read, any column)."""
    path = os.path.join(HERE, dirname, 'per_residue.csv')
    done = {r['UniProt_ID']: int(r['N']) for r in csv.DictReader(open(os.path.join(HERE, dirname, 'per_protein.csv')))}
    d, aa = defaultdict(dict), defaultdict(dict)
    with open(path) as fh:
        rd = csv.reader(fh)
        h = next(rd)
        iu, ip, ia, io = h.index('UniProt_ID'), h.index('pos'), h.index('aa'), h.index(col)
        for row in rd:
            if len(row) <= io or row[iu] not in done or (keep is not None and row[iu] not in keep):
                continue
            d[row[iu]][int(row[ip])] = float(row[io])
            aa[row[iu]][int(row[ip])] = row[ia]
    res = {}
    for u, v in d.items():
        if len(v) == done[u] and max(v) == len(v) - 1:
            res[u] = (''.join(aa[u][i] for i in range(len(v))), np.array([v[i] for i in range(len(v))]))
    return res


def load(col):
    SH = read('hub_frames_v3_shuffle', col)
    YF = read('hub_frames_v3_yf', col)
    wt_sh = defaultdict(list)
    for u in SH:
        wt_sh[u.rsplit('_global_', 1)[0]].append(u)
    wt_yf = {u[:-3]: u for u in YF}
    need = set(wt_sh) | set(wt_yf)
    CAL = read('hub_frames_v3_cal', col, need)
    RCH = read('hub_frames_v3_recharged', col, need)
    GG = read('hub_frames_v3_mpipi', col, need)
    return SH, YF, wt_sh, wt_yf, CAL, RCH, GG


def d1(col, rows_out):
    SH, YF, wt_sh, wt_yf, CAL, RCH, GG = load(col)
    tm = {name: at.type_means(D) for name, D in (('cal', CAL), ('sh', SH), ('yf', YF), ('rch', RCH), ('gg', GG))}
    cohort = sorted(u for u in wt_sh if u in CAL and len(wt_sh[u]) == 3)
    res = {}
    r3 = {k: [] for k in ('raw', 'par', 'res')}
    for u in cohort:
        a, b, c = at.agreement(CAL, SH, [(u, v) for v in wt_sh[u]], tm['cal'], tm['sh'])
        r3['raw'].append(np.mean(a))
        r3['par'].append(np.mean(b))
        r3['res'].append(np.mean(c))
        rows_out.append((col, 'shuffled', u, np.mean(b)))
    res['shuffled'] = tuple(np.array(r3[k]) for k in ('raw', 'par', 'res'))
    pr = [(u, u) for u in cohort if u in RCH]
    res['Mpipi-Recharged'] = at.agreement(CAL, RCH, pr, tm['cal'], tm['rch'])
    pg = [(u, u) for u in cohort if u in GG]
    res['Mpipi-GG'] = at.agreement(CAL, GG, pg, tm['cal'], tm['gg'])
    pgr = [(u, u) for u in cohort if u in GG and u in RCH]
    res['Mpipi-GG vs Recharged'] = at.agreement(GG, RCH, pgr, tm['gg'], tm['rch'])
    py = [(u, wt_yf[u]) for u in wt_yf if u in CAL]
    res['Y-to-F'] = at.agreement(CAL, YF, py, tm['cal'], tm['yf'])
    for k in ('Mpipi-Recharged', 'Mpipi-GG', 'Mpipi-GG vs Recharged', 'Y-to-F'):
        for v in res[k][1]:
            rows_out.append((col, k, '', v))
    return res, cohort, (SH, YF, wt_sh, wt_yf, CAL, RCH, GG)


def main():
    rows_out = []
    p('Arrangement test with contact degree (arrangement_test_degree_plan.md; amendment to arrangement_test.py c4ca3409)')
    resd, cohort, (SH, YF, wt_sh, wt_yf, CAL, RCH, GG) = d1('deg_mean', rows_out)
    reso, _, _ = d1('occupancy', rows_out)
    p(f'natural CALVADOS-2 sequences with all 3 shuffles labelled: {len(cohort)}; Y-to-F mutants with a labelled wild type: '
      f'{sum(u in CAL for u in wt_yf)}')
    p('\nD1  agreement at matched positions, median over proteins (degree; hub occupancy for comparison)')
    p(f'  {"comparison":28s} {"n":>5s}  {"degree raw":>10s} {"pos. removed":>12s} {"pos.+identity":>13s}   {"occupancy pos. removed":>22s}')
    for k in ('shuffled', 'Y-to-F', 'Mpipi-Recharged', 'Mpipi-GG', 'Mpipi-GG vs Recharged'):
        raw, par, rs = resd[k]
        p(f'  {k:28s} {len(raw):5d}  {np.nanmedian(raw):10.2f} {np.nanmedian(par):12.2f} {np.nanmedian(rs):13.2f}   '
          f'{np.nanmedian(reso[k][1]):22.2f}')
    ok = (np.nanmedian(resd['shuffled'][1]) < 0.2 and np.nanmedian(resd['Mpipi-Recharged'][1]) >= 0.4
          and np.nanmedian(resd['Y-to-F'][1]) >= 0.8)
    p(f'  D1 decision rule (shuffle < 0.2, Mpipi-Recharged >= 0.4, Y-to-F >= 0.8, position removed): '
      f'{"MET" if ok else "NOT MET"}')

    p('\nD2  share of the variance of position-adjusted degree explained by residue identity (same proteins)')
    nat = at.identity_share(CAL, cohort)
    shu = at.identity_share(SH, [v for u in cohort for v in wt_sh[u]])
    s_nat, s_shu = at.share(*nat), at.share(*shu)
    rng = np.random.default_rng(0)
    diffs = []
    for _ in range(2000):
        pick = rng.integers(0, len(cohort), len(cohort))
        rn = ([nat[0][k] for k in pick], [nat[1][k] for k in pick])
        rs = ([shu[0][3 * k + j] for k in pick for j in range(3)], [shu[1][3 * k + j] for k in pick for j in range(3)])
        diffs.append(at.share(*rs) - at.share(*rn))
    lo, hi = np.percentile(diffs, [2.5, 97.5])
    p(f'  natural {100 * s_nat:.1f}%   shuffled {100 * s_shu:.1f}%   difference {100 * (s_shu - s_nat):+.1f} points '
      f'(95% bootstrap interval {100 * lo:+.1f} to {100 * hi:+.1f})')

    p('\nD3  natural versus shuffled degree profiles (two-sided Wilcoxon; share with natural above all 3 shuffles)')
    for name, f in (('SD of position-adjusted degree', lambda o: float(at.posadj(o).std())),
                    ('lag-1 autocorrelation along the chain', at.lag1), ('maximum degree', lambda o: float(o.max()))):
        nv = np.array([f(CAL[u][1]) for u in cohort])
        sv = np.array([[f(SH[v][1]) for v in wt_sh[u]] for u in cohort])
        above = np.mean([nv[k] > sv[k].max() for k in range(len(cohort))])
        st = wilcoxon(nv, sv.mean(1))
        p(f'  {name:40s} natural median {np.median(nv):.3f}, shuffled {np.median(sv):.3f}; natural above all 3 shuffles '
          f'in {100 * above:.0f}% (25% expected); P = {st.pvalue:.1e}')

    p('\nD4  residue-type part (mean within-protein z of the position-adjusted readout, no null) against interaction scales,')
    p('    proteins labelled in all three force fields; Pearson r over the 20 residue types (descriptive)')
    lam = {r['one']: float(r['lambdas']) for r in csv.DictReader(open(CAL_TABLE))}
    eps = json.load(open(os.path.join(HERE, 'mpipi_epsilon_table.json')))
    epsm = {a: np.mean([eps[a + b] for b in AAS if a + b in eps]) for a in AAS}
    rc = json.load(open(os.path.join(HERE, 'mpipi_recharged_parameters.json')))['epsilon_kJ_per_mol']
    epsr = {a: np.mean([rc[a + b] for b in AAS]) for a in AAS}
    scales = {'CALVADOS-2 lambda': lam, 'Mpipi-GG epsilon': epsm, 'Mpipi-Recharged epsilon': epsr}
    own = {'CALVADOS-2': 'CALVADOS-2 lambda', 'MPIPI-GG': 'Mpipi-GG epsilon', 'Mpipi-Recharged': 'Mpipi-Recharged epsilon'}
    dirs = {'CALVADOS-2': 'hub_frames_v3_cal', 'MPIPI-GG': 'hub_frames_v3_mpipi', 'Mpipi-Recharged': 'hub_frames_v3_recharged'}
    tm_rows = []
    for col in ('deg_mean', 'occupancy'):
        L = {ff: read(d, col) for ff, d in dirs.items()}
        common = set(L['CALVADOS-2']) & set(L['MPIPI-GG']) & set(L['Mpipi-Recharged'])
        common = {u for u in common if len({len(L[ff][u][1]) for ff in L}) == 1}
        p(f'  {col}: {len(common)} proteins')
        for ff in dirs:
            acc = defaultdict(list)
            for u in common:
                s, o = L[ff][u]
                r = at.posadj(o)
                z = (r - r.mean()) / r.std() if r.std() > 0 else r * 0
                for a, v in zip(s, z):
                    acc[a].append(v)
            m = {a: float(np.mean(acc[a])) for a in AAS}
            for a in AAS:
                tm_rows.append((col, ff, a, m[a]))
            xs = np.array([m[a] for a in AAS])
            rs = '  '.join(f'{k} r = {pearsonr([v[a] for a in AAS], xs)[0]:+.2f}' for k, v in scales.items())
            p(f'    {ff:16s} own scale r = {pearsonr([scales[own[ff]][a] for a in AAS], xs)[0]:+.2f}   ({rs})')

    p('\nD5  frozen degree models on arrangements they never saw (median per-protein Spearman, raw / position-partial)')
    ids_v3 = {r['protein_id'] for r in csv.DictReader(open(os.path.join(HERE, 'hub_predictor', 'data_v3', 'proteins.csv')))}
    unseen = [u for u in list(SH) + list(YF) if u in ids_v3]
    p(f'  shuffle or Y-to-F identifiers present in data_v3: {len(unseen)}')
    models = {'CALVADOS-2': pickle.load(open(os.path.join(HERE, 'hub_vs_degree', 'degree_model_cal.pkl'), 'rb')),
              'consensus': pickle.load(open(os.path.join(HERE, 'hub_predictor', 'degree_consensus_2026-10-13.pkl'), 'rb'))}
    for label, keys, D in (('shuffles', [v for u in cohort for v in wt_sh[u]], SH),
                           ('Y-to-F mutants', [wt_yf[u] for u in wt_yf if u in CAL], YF),
                           ('their natural wild types (seen in training)', cohort, CAL)):
        for name, m in models.items():
            raw, par = [], []
            for u in keys:
                s, o = D[u]
                x = m.predict(T.features(s, 10))
                if o.std() == 0:
                    continue
                raw.append(spearmanr(x, o).correlation)
                par.append(at.partial(x, o))
                rows_out.append(('deg_mean', f'D5 {label} {name}', u, par[-1]))
            p(f'  {label:44s} {name:11s} n={len(raw):5d}  {np.nanmedian(raw):.2f} / {np.nanmedian(par):.2f}')
    p('  reference: held-out CV degree on natural CALVADOS-2 proteins 0.934 raw (hub_vs_degree Part C)')
    open(os.path.join(HERE, 'arrangement_test_degree_summary.txt'), 'w').write('\n'.join(out) + '\n')
    with open(os.path.join(HERE, 'arrangement_test_degree_values.csv'), 'w') as fh:
        w = csv.writer(fh)
        w.writerow(['readout', 'comparison', 'id', 'partial'])
        w.writerows(rows_out)
    with open(os.path.join(HERE, 'arrangement_test_degree_typemeans.csv'), 'w') as fh:
        w = csv.writer(fh)
        w.writerow(['readout', 'force_field', 'aa', 'mean_z'])
        w.writerows(tm_rows)


if __name__ == '__main__':
    main()
