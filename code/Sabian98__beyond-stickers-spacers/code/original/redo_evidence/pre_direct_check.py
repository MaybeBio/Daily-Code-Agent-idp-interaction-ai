"""
pre_direct_check.py — a quick, model-free look at the PeptoneDB-Integrative PRE data against frozen hub-occupancy
predictions. This is NOT the planned test (that one reweights our simulations against the PRE); it only asks whether
residues predicted to be hubs show more long-range PRE than their sequence separation from the spin labels implies.

Per spin label s (file PRE-<s>.dat = Gamma2 rates, PRE-I-<s>.dat = intensity ratios I_para/I_dia; residues 1-based):
  contact signal c_i = 1 - I_i (ratios, clipped to [0, 1]) or Gamma2_i / max Gamma2 (rates)
  residues with |i - s| < 10 are dropped; c is regressed on log|i - s| separately on each side of the label, and the
  residual is the excess over what sequence separation predicts
Per residue: excess = mean residual over the labels that cover it (at least 2). Compared with the frozen predictions
(consensus / pooled / CALVADOS-2-only) and with distance from the nearest chain end, by Spearman per protein.
Disordered proteins only (FLN5A3A3 and d131d are folded or partly folded and are left out).
  ~/.venvs/bender/bin/python pre_direct_check.py  ->  pre_direct_check_summary.txt
"""
import csv
import hashlib
import os
import re

import numpy as np
from scipy.stats import spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
PEPTONE = os.path.expanduser('~/Downloads/PeptoneDBs/PeptoneDB-Integrative')
GAIVI = os.path.join(HERE, 'gaivi')
# PeptoneDB folder -> (frozen predictions file, id in that file)
PROTEINS = {'OPN': ('predictions_frozen_2026-10-02_peptone_pre.csv', 'PRE_OPN'),
            'Tau': ('predictions_frozen_2026-10-02_peptone_pre.csv', 'PRE_Tau441'),
            'ERNTD': ('predictions_frozen_2026-10-02_peptone_pre.csv', 'PRE_ERNTD'),
            'ANAC046': ('predictions_frozen_2026-10-02_peptone_pre.csv', 'PRE_ANAC046'),
            'ACTR': ('predictions_frozen_2026-10-02_ped_fresh.csv', 'PED00536'),
            'Sic1': ('predictions_frozen_2026-10-02.csv', 'PED_Sic1'),
            'aSyn': ('predictions_frozen_2026-10-02.csv', 'PED_aSyn'),
            'drkN': ('predictions_frozen_2026-10-02.csv', 'PED_drkN')}
NEW = {'OPN', 'Tau', 'ERNTD', 'ANAC046'}
MIN_SEP = 10
out = []


def p(s=''):
    out.append(s)
    print(s, flush=True)


def frozen(fname):
    path = os.path.join(GAIVI, fname)
    h = hashlib.sha256(open(path, 'rb').read()).hexdigest()
    if h not in open(os.path.join(GAIVI, 'FROZEN.sha256')).read():
        raise SystemExit(f'{fname} is not registered in FROZEN.sha256')
    d = {}
    for r in csv.DictReader(open(path)):
        d.setdefault(r['id'], []).append((r['aa'], float(r['consensus']), float(r['pooled']), float(r['cal'])))
    return d


def excess_profile(folder, n):
    resid = [[] for _ in range(n)]
    n_labels = 0
    for f in sorted(os.listdir(folder)):
        m = re.match(r'PRE-(I-)?(\d+)\.dat$', f)
        if not m:
            continue
        ratio, s = bool(m.group(1)), int(m.group(2)) - 1
        vals = {}
        for line in open(os.path.join(folder, f)):
            t = line.split()
            if len(t) >= 2 and t[1].lower() != 'nan':
                i = int(float(t[0])) - 1
                if 0 <= i < n:
                    vals[i] = float(t[1])
        if not vals:
            continue
        top = max(vals.values()) or 1.0
        c = {i: (min(1.0, max(0.0, 1.0 - v)) if ratio else v / top) for i, v in vals.items() if abs(i - s) >= MIN_SEP}
        n_labels += 1
        for side in (-1, 1):
            idx = [i for i in c if (i - s) * side > 0]
            if len(idx) < 5:
                continue
            x = np.log([abs(i - s) for i in idx])
            y = np.array([c[i] for i in idx])
            A = np.vstack([x, np.ones_like(x)]).T
            coef = np.linalg.lstsq(A, y, rcond=None)[0]
            for i, r in zip(idx, y - A @ coef):
                resid[i].append(r)
    prof = np.array([np.mean(v) if len(v) >= 2 else np.nan for v in resid])
    return prof, n_labels


def main():
    cache = {}
    p('Model-free PRE check (PeptoneDB-Integrative): excess long-range PRE beyond sequence separation vs frozen '
      'hub-occupancy predictions; per-protein Spearman over residues covered by >= 2 labels\n')
    keys = ('consensus', 'pooled', 'cal', 'position')
    p(f'  {"protein":9s} {"N":>4s} {"labels":>6s} {"residues":>8s}' + ''.join(f'{k:>11s}' for k in keys))
    res = {}
    for prot, (fname, pid) in PROTEINS.items():
        if prot == 'Sic1':   # PeptoneDB values include the diamagnetic baseline, so they are not Gamma2 (code review 2026-10-03)
            p(f'  {prot:9s} left out: values include the diamagnetic baseline (not Gamma2)')
            continue
        if fname not in cache:
            cache[fname] = frozen(fname)
        pr = cache[fname][pid]
        folder = os.path.join(PEPTONE, prot)
        seq = ''.join(line.strip() for line in open(os.path.join(folder, 'seq.fasta')) if not line.startswith('>'))
        assert ''.join(a for a, *_ in pr) == seq, f'{prot}: sequence differs from the frozen predictions'
        n = len(seq)
        prof, n_labels = excess_profile(folder, n)
        ok = ~np.isnan(prof)
        P = np.array([x[1:] for x in pr])
        end = np.minimum(np.arange(n), n - 1 - np.arange(n))
        rho = {k: spearmanr(P[ok, j], prof[ok]).correlation for j, k in enumerate(('consensus', 'pooled', 'cal'))}
        rho['position'] = spearmanr(end[ok], prof[ok]).correlation
        res[prot] = rho
        p(f'  {prot + (" *" if prot in NEW else ""):9s} {n:4d} {n_labels:6d} {int(ok.sum()):8d}'
          + ''.join(f'{rho[k]:+11.2f}' for k in keys))
    for label, sel in (('new proteins (*)', [k for k in res if k in NEW]), ('all with usable PRE', list(res))):
        p(f'\n  {label}: median ' + ', '.join(f'{k} {np.nanmedian([res[x][k] for x in sel]):+.2f}' for k in keys)
          + f'; consensus > 0 in {sum(res[x]["consensus"] > 0 for x in sel)}/{len(sel)}')
    p('\n  Caveats: PRE is measured at different temperatures, pH and ionic strengths (277-310 K) and reports contacts '
      'with only a few label sites; excess PRE is closer to contact degree than to betweenness. The planned test '
      'reweights our simulations against these data instead.')
    open(os.path.join(HERE, 'pre_direct_check_summary.txt'), 'w').write('\n'.join(out) + '\n')


if __name__ == '__main__':
    main()
