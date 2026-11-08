"""
positional_null_check.py — are the hub POSITIONS of experimental ensembles predicted by geometry alone?
For each PRE ensemble: Spearman between the experimental per-residue centrality (position-adjusted z) and that of the
Rg-matched residue-neutral null (ped_null.py), and, where simulations of the same construct exist (α-syn SASDQJ7,
eIF4G1 P39935), of CALVADOS-2 and MPIPI-GG; plus the share of each predictor's top 10% that lies in the experimental
top 20% (chance 20%). Both contact rules.  python3 positional_null_check.py -> positional_null_check.txt
"""
import difflib
import json
import os
import sys
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import exp_anchor as ea  # noqa: E402
import hub_frames_tests as ht  # noqa: E402
import ped_null as pn  # noqa: E402

SIMS = {'PED00024': 'SASDQJ7', 'PED00225': 'P39935'}
IDS = ['PED00024', 'PED00006', 'PED00003', 'PED00157', 'PED00160', 'PED00225', 'PED00499', 'PED00500', 'PED00511']


def top_hits(pred, exp, f_pred=0.10, f_exp=0.20):
    n = len(exp)
    tp = set(sorted(range(n), key=lambda i: -pred[i])[:max(1, round(f_pred * n))])
    te = set(sorted(range(n), key=lambda i: -exp[i])[:max(1, round(f_exp * n))])
    return len(tp & te), len(tp)


def one(job):
    entry, c = job
    path = dict((e, f) for e, f, *_ in pn.ENSEMBLES)[entry]
    seq, models = pn.load(path)
    if len(models) > 1000:
        models = models[::-(-len(models) // 1000)]
    nm = pn.null_models(seq, models, 1000 + IDS.index(entry))[0]
    out = {}
    for rule in ('fixed', 'size'):
        ze = pn.centrality(seq, models, rule, c)
        zn = pn.centrality(seq, nm, rule, c)
        row = {'null': (ea.spearman(zn, ze), top_hits(zn, ze))}
        if entry in SIMS:
            for ff in ('CALVADOS-2', 'MPIPI-GG'):
                r = ea.sim_results(ea.SIM[(ff, rule)], SIMS[entry])
                if not r:
                    continue
                sm = difflib.SequenceMatcher(None, r[0], seq, autojunk=False)
                pairs = [(b.a + k, b.b + k) for b in sm.get_matching_blocks() for k in range(b.size)]
                zs_ = [r[1][i] for i, j in pairs]
                ze_ = [ze[j] for i, j in pairs]
                zn_ = [zn[j] for i, j in pairs]
                row[ff] = (ea.spearman(zs_, ze_), top_hits(zs_, ze_), ea.spearman(zs_, zn_))
        out[rule] = row
    return entry, out


if __name__ == '__main__':
    c = json.load(open(ht.C_FILE))['calvados']
    names = {e: prot for e, f, prot, r in pn.ENSEMBLES}
    with Pool(9) as pool:
        res = dict(pool.imap_unordered(one, [(e, c) for e in IDS]))
    lines = ['Do residue-neutral null chains predict WHERE the experimental hubs are? (Spearman of per-residue centrality; '
             'hits = predictor top 10% within experimental top 20%, chance 20%)']
    for rule in ('fixed', 'size'):
        lines.append(f'\nRule: {rule}')
        for e in IDS:
            row = res[e][rule]
            s, (h, n) = row['null']
            line = f'  {e} {names[e][:22]:22s} null: rho {s:+.2f}, hits {h}/{n}'
            for ff in ('CALVADOS-2', 'MPIPI-GG'):
                if ff in row:
                    s2, (h2, n2), s3 = row[ff]
                    line += f' | {ff}: rho {s2:+.2f}, hits {h2}/{n2} (rho with null {s3:+.2f})'
            lines.append(line)
    open(os.path.join(HERE, 'positional_null_check.txt'), 'w').write('\n'.join(lines) + '\n')
    print('\n'.join(lines))
