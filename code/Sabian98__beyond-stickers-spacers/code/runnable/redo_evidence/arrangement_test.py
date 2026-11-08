"""
arrangement_test.py — does the arrangement of residues, rather than their chemistry, set hub occupancy?

Plan fixed on 2026-10-04 at 16:20, before any label of the shuffled or Y-to-F trajectories was looked at, and after
frozen predictions for them were written (gaivi/predictions_frozen_2026-10-02_arrangement_2026-10-04.csv). Data:
  arrangement axis   597 natural IDRs, each with 3 composition-preserving global shuffles simulated in CALVADOS-2 with the
                     BENDER protocol (earlier BENDER-DYN version), labels hub_frames_v3_shuffle; the natural sequences are
                     the CALVADOS-2 v3 labels (hub_frames_v3_cal)
  chemistry axis     (i) the same sequences in Mpipi-Recharged and Mpipi-GG where simulated (hub_frames_v3_recharged,
                     hub_frames_v3_mpipi); (ii) 800 IDRs with every tyrosine replaced by phenylalanine in CALVADOS-2
                     (lambda 0.977 -> 0.867, sigma almost unchanged), labels hub_frames_v3_yf
Agreement between two occupancy profiles of the same length is the per-protein Spearman correlation at matched chain
positions, raw and as a partial correlation that removes distance from the nearest chain end (ranks residualised on the
rank of end distance, as for the PED tests), and also after removing each label set's residue-type means of
position-adjusted occupancy (identity removed). Medians over proteins; a protein's three shuffles are averaged first.

H1  Arrangement versus chemistry. Changing the arrangement at fixed composition changes the occupancy profile far more
    than changing the chemistry at fixed arrangement. Supported if the position-removed agreement between a natural
    sequence and its shuffles is below 0.2, while it is at least 0.4 between force fields for the same sequences and at
    least 0.8 between wild type and Y-to-F.
H2  Identity share. Residue identity explains a larger share of the within-protein variance of position-adjusted occupancy
    in the shuffled sequences than in their natural counterparts (same proteins), because natural arrangements add
    context-driven variance. Paired bootstrap over proteins (2,000 resamples) for the difference.
H3  Are natural hub profiles atypical for their composition? (descriptive, two-sided) Natural versus the mean of its
    shuffles for (a) the standard deviation of position-adjusted occupancy, (b) its lag-1 autocorrelation along the chain
    and (c) the maximum occupancy; Wilcoxon signed-rank over proteins, and the share of proteins in which the natural value
    exceeds all three shuffles (25% expected if the natural arrangement were exchangeable with random ones). A difference
    shows that natural arrangements are non-random with respect to hubs, not that hubs were selected.
H4  Does the model generalise to arrangements it never saw? Frozen consensus and CALVADOS-2 models on the shuffles
    (Spearman, raw and position-partial, medians), next to the held-out CV median on natural sequences (0.89).
H5  Context-dependent substitution effects. For the Y-to-F mutants, the simulated change in occupancy (mutant - wild type)
    against the frozen model's predicted change, at the substituted positions and separately at all other positions.
  ~/.venvs/bender/bin/python arrangement_test.py  ->  arrangement_test_summary.txt
"""
import csv
import os
from collections import defaultdict

import numpy as np
from scipy.stats import rankdata, spearmanr, wilcoxon

HERE = os.path.dirname(os.path.abspath(__file__))
PRED = os.path.join(HERE, 'gaivi', 'predictions_frozen_2026-10-02_arrangement_2026-10-04.csv')
out = []


def p(s=''):
    out.append(s)
    print(s, flush=True)


def read(dirname, keep=None):
    """{id: (sequence, occupancy array)} for proteins whose residue rows are complete."""
    path = os.path.join(HERE, dirname, 'per_residue.csv')
    done = {r['UniProt_ID']: int(r['N']) for r in csv.DictReader(open(os.path.join(HERE, dirname, 'per_protein.csv')))}
    d, aa = defaultdict(dict), defaultdict(dict)
    with open(path) as fh:
        rd = csv.reader(fh)
        h = next(rd)
        iu, ip, ia, io = h.index('UniProt_ID'), h.index('pos'), h.index('aa'), h.index('occupancy')
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


def posadj(o):
    n = len(o)
    x = np.arange(n) / max(n - 1, 1)
    A = np.vstack([np.ones(n), x, x ** 2]).T
    return o - A @ np.linalg.lstsq(A, o, rcond=None)[0]


def partial(a, b):
    n = len(a)
    e = np.minimum(np.arange(n), n - 1 - np.arange(n)).astype(float)
    ra, rb, re = rankdata(a), rankdata(b), rankdata(e)
    X = np.vstack([re, np.ones(n)]).T
    ea = ra - X @ np.linalg.lstsq(X, ra, rcond=None)[0]
    eb = rb - X @ np.linalg.lstsq(X, rb, rcond=None)[0]
    return float(np.corrcoef(ea, eb)[0, 1]) if ea.std() > 0 and eb.std() > 0 else np.nan


def type_means(D):
    acc = defaultdict(list)
    for s, o in D.values():
        for a, v in zip(s, posadj(o)):
            acc[a].append(v)
    return {a: float(np.mean(v)) for a, v in acc.items()}


def noid(tm, s, o):
    return posadj(o) - np.array([tm[a] for a in s])


def agreement(A, B, pairs, tmA, tmB):
    raw, par, res = [], [], []
    for ua, ub in pairs:
        (sa, oa), (sb, ob) = A[ua], B[ub]
        if len(oa) != len(ob) or oa.std() == 0 or ob.std() == 0:
            continue
        raw.append(spearmanr(oa, ob).correlation)
        par.append(partial(oa, ob))
        res.append(spearmanr(noid(tmA, sa, oa), noid(tmB, sb, ob)).correlation)
    return np.array(raw), np.array(par), np.array(res)


def identity_share(D, keys):
    resid, aas = [], []
    for u in keys:
        s, o = D[u]
        r = posadj(o)
        resid.append(r / r.std() if r.std() > 0 else r * 0.0)     # keep one entry per protein so the bootstrap stays aligned
        aas.append(np.array(list(s)))
    return resid, aas


def share(resid, aas):
    r, a = np.concatenate(resid), np.concatenate(aas)
    between = sum((a == t).sum() * r[a == t].mean() ** 2 for t in set(a))
    return between / (r ** 2).sum()


def lag1(o):
    r = posadj(o)
    return float(np.corrcoef(r[:-1], r[1:])[0, 1]) if r.std() > 0 else np.nan


def main():
    SH = read('hub_frames_v3_shuffle')
    YF = read('hub_frames_v3_yf')
    wt_sh = defaultdict(list)
    for u in SH:
        wt_sh[u.rsplit('_global_', 1)[0]].append(u)
    wt_yf = {u[:-3]: u for u in YF}
    need = set(wt_sh) | set(wt_yf)
    CAL = read('hub_frames_v3_cal', need)
    RCH = read('hub_frames_v3_recharged', need)
    GG = read('hub_frames_v3_mpipi', need)
    tm = {name: type_means(D) for name, D in (('cal', CAL), ('sh', SH), ('yf', YF), ('rch', RCH), ('gg', GG))}
    cohort = sorted(u for u in wt_sh if u in CAL and len(wt_sh[u]) == 3)
    p(f'Arrangement test (plan fixed before the shuffled and Y-to-F labels were examined). Natural CALVADOS-2 sequences with '
      f'all 3 shuffles labelled: {len(cohort)}; Y-to-F mutants with a labelled wild type: {sum(u in CAL for u in wt_yf)}')

    p('\nH1  agreement with the natural CALVADOS-2 profile at matched positions (median over proteins)')
    p(f'  {"comparison":52s} {"n":>5s} {"raw":>6s} {"position removed":>17s} {"position + identity removed":>28s}')

    def row(label, raw, par, res):
        p(f'  {label:52s} {len(raw):5d} {np.nanmedian(raw):6.2f} {np.nanmedian(par):17.2f} {np.nanmedian(res):28.2f}')
    # shuffles: average the three shuffles of a protein first
    r3 = {k: [] for k in ('raw', 'par', 'res')}
    for u in cohort:
        a, b, c = agreement(CAL, SH, [(u, v) for v in wt_sh[u]], tm['cal'], tm['sh'])
        r3['raw'].append(np.mean(a))
        r3['par'].append(np.mean(b))
        r3['res'].append(np.mean(c))
    row('shuffled arrangement, same chemistry (CALVADOS-2)', np.array(r3['raw']), np.array(r3['par']), np.array(r3['res']))
    sh_par = np.array(r3['par'])
    pr = [(u, u) for u in cohort if u in RCH]
    rr = agreement(CAL, RCH, pr, tm['cal'], tm['rch'])
    row('same arrangement, Mpipi-Recharged chemistry', *rr)
    pg = [(u, u) for u in cohort if u in GG]
    if pg:
        rg = agreement(CAL, GG, pg, tm['cal'], tm['gg'])
        row('same arrangement, Mpipi-GG chemistry', *rg)
    py = [(u, wt_yf[u]) for u in wt_yf if u in CAL]
    ry = agreement(CAL, YF, py, tm['cal'], tm['yf'])
    row('same arrangement, tyrosines to phenylalanine (CALVADOS-2)', *ry)
    ok = (np.nanmedian(sh_par) < 0.2 and np.nanmedian(rr[1]) >= 0.4 and np.nanmedian(ry[1]) >= 0.8)
    p(f'  H1 decision rule (shuffle < 0.2, other force field >= 0.4, Y-to-F >= 0.8, position removed): '
      f'{"met" if ok else "not met"}')

    p('\nH2  share of the variance of position-adjusted occupancy explained by residue identity (same proteins)')
    nat = identity_share(CAL, cohort)
    shu = identity_share(SH, [v for u in cohort for v in wt_sh[u]])
    s_nat, s_shu = share(*nat), share(*shu)
    rng = np.random.default_rng(0)
    idx_nat = {u: k for k, u in enumerate(cohort)}
    diffs = []
    for _ in range(2000):
        pick = rng.integers(0, len(cohort), len(cohort))
        rn = ([nat[0][k] for k in pick], [nat[1][k] for k in pick])
        rs = ([shu[0][3 * k + j] for k in pick for j in range(3)], [shu[1][3 * k + j] for k in pick for j in range(3)])
        diffs.append(share(*rs) - share(*rn))
    lo, hi = np.percentile(diffs, [2.5, 97.5])
    p(f'  natural {100 * s_nat:.1f}%   shuffled {100 * s_shu:.1f}%   difference {100 * (s_shu - s_nat):+.1f} points '
      f'(95% bootstrap interval {100 * lo:+.1f} to {100 * hi:+.1f})')

    p('\nH3  natural versus shuffled hub profiles (two-sided Wilcoxon over proteins; share with natural above all 3 shuffles)')
    for name, f in (('SD of position-adjusted occupancy', lambda o: float(posadj(o).std())),
                    ('lag-1 autocorrelation along the chain', lag1), ('maximum occupancy', lambda o: float(o.max()))):
        nv = np.array([f(CAL[u][1]) for u in cohort])
        sv = np.array([[f(SH[v][1]) for v in wt_sh[u]] for u in cohort])
        above = np.mean([nv[k] > sv[k].max() for k in range(len(cohort))])
        st = wilcoxon(nv, sv.mean(1))
        p(f'  {name:40s} natural median {np.median(nv):.3f}, shuffled {np.median(sv):.3f}; natural above all 3 shuffles '
          f'in {100 * above:.0f}% (25% expected); P = {st.pvalue:.1e}')

    p('\nH4  frozen models on arrangements they never saw (median per-protein Spearman, raw / position-partial)')
    pred = defaultdict(list)
    for r in csv.DictReader(open(PRED)):
        pred[r['id']].append((float(r['consensus']), float(r['cal'])))
    for label, keys, D in (('shuffles', [v for u in cohort for v in wt_sh[u]], SH),
                           ('Y-to-F mutants', [wt_yf[u] for u in wt_yf if u in CAL], YF),
                           ('their natural wild types (seen in training)', cohort, CAL)):
        for j, model in enumerate(('consensus', 'CALVADOS-2')):
            raw, par = [], []
            for u in keys:
                if u not in pred or len(pred[u]) != len(D[u][1]):
                    continue
                x = np.array([q[j] for q in pred[u]])
                raw.append(spearmanr(x, D[u][1]).correlation)
                par.append(partial(x, D[u][1]))
            p(f'  {label:44s} {model:11s} n={len(raw):5d}  {np.nanmedian(raw):.2f} / {np.nanmedian(par):.2f}')
    p('  reference: held-out CV on natural CALVADOS-2 sequences 0.89 (cv_full.py)')

    p('\nH5  Y-to-F: simulated versus predicted change in occupancy (mutant - wild type), pooled Spearman')
    for label, at_site in (('substituted positions', True), ('all other positions', False)):
        sim_d, pred_d = [], []
        for u, m in wt_yf.items():
            if u not in CAL or u not in pred or m not in pred:
                continue
            s_wt, o_wt = CAL[u]
            _, o_m = YF[m]
            if len(o_wt) != len(o_m):
                continue
            pw = np.array([q[0] for q in pred[u]])
            pm = np.array([q[0] for q in pred[m]])
            mask = np.array([a == 'Y' for a in s_wt]) if at_site else np.array([a != 'Y' for a in s_wt])
            sim_d += list((o_m - o_wt)[mask])
            pred_d += list((pm - pw)[mask])
        sim_d, pred_d = np.array(sim_d), np.array(pred_d)
        p(f'  {label:24s} residues {len(sim_d):6d}: mean simulated change {sim_d.mean():+.4f}, predicted {pred_d.mean():+.4f}; '
          f'Spearman {spearmanr(sim_d, pred_d).correlation:+.2f}')
    open(os.path.join(HERE, 'arrangement_test_summary.txt'), 'w').write('\n'.join(out) + '\n')


if __name__ == '__main__':
    main()
