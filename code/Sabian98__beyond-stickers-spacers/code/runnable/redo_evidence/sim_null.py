"""
sim_null.py — how much of each force field's residue-type hub ranking is beyond geometry?

For each simulated protein: a pool of self-avoiding Calpha chains with the force field's own bead sizes as excluded
volume and no attraction (ped_null.grow_chain; three chain stiffnesses), Rg-quantile-matched to the simulation's own
radius-of-gyration distribution (computed from the same trajectory frames). Both ensembles analysed identically
(fixed 0.8 nm rule = the near-neutral rule; bead-size rule for completeness), exact betweenness, mean over frames,
position-adjusted z. Excess = simulation - null, per residue type and class.
Sets: CALVADOS-2 and MPIPI-GG on the same N_SAMPLE proteins of the 908 paired set; Mpipi-Recharged on its 289.
Simulation centralities are the 1,000-frame values already computed (hub_frames_*/per_residue.csv).
Standard library only.  python3 sim_null.py  ->  sim_null_summary.txt, sim_null_by_type.csv
"""
import csv
import json
import math
import os
import random
import sys
import bisect
from collections import defaultdict
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import exp_anchor as ea  # noqa: E402
import hub_frames as hf  # noqa: E402
import hub_frames_tests as ht  # noqa: E402
import ped_null as pn  # noqa: E402

N_SAMPLE = int(os.environ.get('SIM_NULL_N', '300'))        # proteins of the 908 paired set (908 = all)
TAG = os.environ.get('SIM_NULL_TAG', '')                     # output suffix, e.g. _all
N_POOL, N_NULL, SEED = 2000, 500, 7
SIM_DIRS = {('CALVADOS-2', 'fixed'): 'hub_frames_paired_calvados', ('CALVADOS-2', 'size'): 'hub_frames_sizeaware_calvados',
            ('MPIPI-GG', 'fixed'): 'hub_frames_mpipi', ('MPIPI-GG', 'size'): 'hub_frames_sizeaware_mpipi',
            ('Mpipi-Recharged', 'fixed'): 'hub_frames_recharged', ('Mpipi-Recharged', 'size'): 'hub_frames_sizeaware_recharged'}
AAS = 'GSACYFWNTPQHDRELIVMK'
CLASS = {**dict.fromkeys('FWY', 'FWY'), **dict.fromkeys('DE', 'DE'), **dict.fromkeys('KR', 'KR'),
         **dict.fromkeys('GSA', 'GSA'), **dict.fromkeys('ILVM', 'ILVM')}
out = []


def p(s=''):
    out.append(s)
    print(s, flush=True)


def sim_z(dirname, ids):
    want = set(ids)
    z = defaultdict(dict)
    for r in csv.DictReader(open(os.path.join(HERE, dirname, 'per_residue.csv'))):
        if r['UniProt_ID'] in want:
            z[r['UniProt_ID']][int(r['pos'])] = (r['aa'], float(r['z_posadj']))
    return {u: [d[i] for i in range(len(d))] for u, d in z.items()}


def rg3(X, Y, Z):
    n = len(X)
    cx, cy, cz = sum(X) / n, sum(Y) / n, sum(Z) / n
    return math.sqrt(sum((x - cx) ** 2 + (y - cy) ** 2 + (z - cz) ** 2 for x, y, z in zip(X, Y, Z)) / n)


def sim_rgs(ff, uid, group):
    if ff == 'CALVADOS-2':
        seq, get_frame, n_frames, idx, close = ht.cal_frames(uid, group)
    elif ff == 'MPIPI-GG':
        seq, get_frame, n_frames, idx, close = ht.mpi_frames(uid)
    else:
        seq, get_frame, n_frames, idx, close = ht.cal_frames(uid, group, base=ht.RECHARGED)
    rgs = []
    for k in idx[:: max(1, len(idx) // N_NULL)][:N_NULL]:
        xyz, box = get_frame(k)
        hf.unwrap(xyz, box)
        rgs.append(rg3(*xyz))
    close()
    return seq, rgs


def ff_sigma(ff, seq, tab):
    if ff == 'CALVADOS-2':
        s = {r['one']: float(r['sigmas']) * 10.0 for r in csv.DictReader(open(ht.CAL_PARAMS))}
        return [s[a] for a in seq], ht.cal_sigma_matrix(seq)
    return [tab[a + a] for a in seq], ht.table_sigma_matrix(seq, tab)


def one(job):
    ff, uid, group, zfix, zsize, c, seed = job
    try:
        tab = json.load(open(ht.SIGMA_TABLE))
        seq, rgs = sim_rgs(ff, uid, group)
        if len(seq) != len(zfix) or (zsize and len(seq) != len(zsize)):
            return None
        sig, S = ff_sigma(ff, seq, tab)
        rng = random.Random(seed)
        pool = []
        for angles in ((85.0, 150.0), (100.0, 150.0), (115.0, 155.0)):
            target = len(pool) + N_POOL // 3
            while len(pool) < target:
                ch = pn.grow_chain(sig, rng, angles=angles)
                if ch is not None:
                    pool.append(ch)
        prg = sorted((rg3(*ch), k) for k, ch in enumerate(pool))
        vals = [v for v, _ in prg]
        chosen = []
        for r in rgs:
            k = bisect.bisect_left(vals, r)
            k = min(range(max(0, k - 1), min(len(vals), k + 1)), key=lambda t: abs(vals[t] - r))
            chosen.append(pool[prg[k][1]])
        res = {'ff': ff, 'uid': uid, 'seq': seq, 'rg_sim': sum(rgs) / len(rgs),
               'rg_null': sum(rg3(*ch) for ch in chosen) / len(chosen), 'rg_pool': sum(vals) / len(vals)}
        n = len(seq)
        for rule in ('fixed', 'size'):
            C2, maxcut = (None, 8.0) if rule == 'fixed' else ht.c2_matrix(S, c)
            tot = [0.0] * n
            for X, Y, Z in chosen:
                adj = hf.contact_lists([X, Y, Z], n, C2, maxcut)
                b = hf.betweenness(adj, n)
                tot = [a + v for a, v in zip(tot, b)]
            zn = ea.zs([v / len(chosen) for v in tot]) if hasattr(ea, 'zs') else None
            if zn is None:
                from hub_z_v2_stdlib import zs as _zs
                zn = _zs([v / len(chosen) for v in tot])
            from hub_z_v2_stdlib import quad_resid
            res[f'null_{rule}'] = quad_resid(zn) if zn else None
        res['sim_fixed'] = zfix
        res['sim_size'] = zsize
        return res
    except Exception as e:
        return {'error': repr(e), 'ff': ff, 'uid': uid}


def main():
    cs = json.load(open(ht.C_FILE))
    paired = [(r['UniProt_ID'], r['group']) for r in csv.DictReader(open(os.path.join(HERE, 'hub_frames_mpipi', 'ids.csv')))]
    paired = random.Random(SEED).sample(paired, min(N_SAMPLE, len(paired)))
    rech = [(r['UniProt_ID'], r['group']) for r in csv.DictReader(open(os.path.join(HERE, 'recharged_ids.csv')))]
    sets = {'CALVADOS-2': paired, 'MPIPI-GG': paired, 'Mpipi-Recharged': rech}
    ckey = {'CALVADOS-2': 'calvados', 'MPIPI-GG': 'mpipi', 'Mpipi-Recharged': 'recharged'}
    jobs = []
    for ff, ids in sets.items():
        zf = sim_z(SIM_DIRS[(ff, 'fixed')], [u for u, g in ids])
        zs_ = sim_z(SIM_DIRS[(ff, 'size')], [u for u, g in ids])
        for k, (u, g) in enumerate(ids):
            if u in zf:
                jobs.append((ff, u, g, [z for a, z in zf[u]], [z for a, z in zs_[u]] if u in zs_ else None, cs[ckey[ff]], SEED * 1000 + k))
    print(f'{len(jobs)} jobs', flush=True)
    acc = {(ff, rule, w): defaultdict(list) for ff in sets for rule in ('fixed', 'size') for w in ('sim', 'null')}
    prot = defaultdict(list)          # (ff, rule) -> per-protein class excess dicts
    rg_rows = defaultdict(list)
    fails = 0
    per_res = open(os.path.join(HERE, f'sim_null_per_residue{TAG}.csv'), 'w', newline='')
    pw = csv.writer(per_res)
    pw.writerow(['force_field', 'UniProt_ID', 'pos', 'aa', 'N', 'z_sim_fixed', 'z_null_fixed', 'z_sim_size', 'z_null_size'])
    with Pool(9) as pool:
        for done, r in enumerate(pool.imap_unordered(one, jobs, chunksize=1), 1):
            if done % 100 == 0:
                print(f'{done}/{len(jobs)}', flush=True)
            if not r or 'error' in r:
                fails += 1
                continue
            ff, seq = r['ff'], r['seq']
            for i, a in enumerate(seq):
                pw.writerow([ff, r['uid'], i, a, len(seq), f"{r['sim_fixed'][i]:.4f}", f"{r['null_fixed'][i]:.4f}" if r['null_fixed'] else '',
                             f"{r['sim_size'][i]:.4f}" if r['sim_size'] else '', f"{r['null_size'][i]:.4f}" if r['null_size'] else ''])
            rg_rows[ff].append((r['rg_sim'], r['rg_null'], r['rg_pool']))
            for rule in ('fixed', 'size'):
                zs_sim, zs_null = r[f'sim_{rule}'], r[f'null_{rule}']
                if not zs_sim or not zs_null:
                    continue
                cls_s, cls_n = defaultdict(list), defaultdict(list)
                for a, x, y in zip(seq, zs_sim, zs_null):
                    acc[(ff, rule, 'sim')][a].append(x)
                    acc[(ff, rule, 'null')][a].append(y)
                    if a in CLASS:
                        cls_s[CLASS[a]].append(x)
                        cls_n[CLASS[a]].append(y)
                prot[(ff, rule)].append({k: sum(cls_s[k]) / len(cls_s[k]) - sum(cls_n[k]) / len(cls_n[k]) for k in cls_s})
    p(f'Simulation vs residue-neutral null (self-avoiding chains with each force field\'s bead sizes, no attraction, '
      f'Rg-matched to each simulation). {len(jobs) - fails} protein runs, {fails} failed.')
    for ff, rows in rg_rows.items():
        p(f'  {ff}: mean Rg sim {sum(r[0] for r in rows) / len(rows):.1f} Å, matched null {sum(r[1] for r in rows) / len(rows):.1f} Å, '
          f'unmatched pool {sum(r[2] for r in rows) / len(rows):.1f} Å ({len(rows)} proteins)')
    csv_rows = []
    for rule in ('fixed', 'size'):
        p(f'\nRule: {rule}   (mean position-adjusted z: simulation / null / excess)')
        p('  aa ' + ''.join(f'{ff:>26s}' for ff in sets))
        prof = {ff: {} for ff in sets}
        for aa in AAS:
            line = f'  {aa}  '
            for ff in sets:
                s, n_ = acc[(ff, rule, 'sim')][aa], acc[(ff, rule, 'null')][aa]
                if s:
                    ms, mn = sum(s) / len(s), sum(n_) / len(n_)
                    prof[ff][aa] = ms - mn
                    line += f'   {ms:+.2f} / {mn:+.2f} / {ms - mn:+.2f}  '
                    csv_rows.append([rule, ff, aa, len(s), f'{ms:.4f}', f'{mn:.4f}', f'{ms - mn:.4f}'])
                else:
                    line += f'{"-":>26s}'
            p(line)
        ffs = list(sets)
        for i in range(len(ffs)):
            for j in range(i + 1, len(ffs)):
                ks = [aa for aa in AAS if aa in prof[ffs[i]] and aa in prof[ffs[j]]]
                p(f'  excess-profile agreement {ffs[i]} vs {ffs[j]}: r = {ea.pearson([prof[ffs[i]][a] for a in ks], [prof[ffs[j]][a] for a in ks]):+.2f}'
                  f' ({ea.pearson([prof[ffs[i]][a] for a in ks if a != "P"], [prof[ffs[j]][a] for a in ks if a != "P"]):+.2f} without P)')
        for ff in sets:
            pr = prot[(ff, rule)]
            both = [d['FWY'] - d['DE'] for d in pr if 'FWY' in d and 'DE' in d]
            k = sum(v > 0 for v in both)
            p(f'  {ff}: per-protein (FWY excess) - (DE excess) > 0 in {k}/{len(both)} ({k / max(1, len(both)):.0%}); '
              + '  '.join(f'{c} excess {sum(d[c] for d in pr if c in d) / max(1, sum(1 for d in pr if c in d)):+.2f}'
                          for c in ('FWY', 'DE', 'KR', 'GSA', 'ILVM')))
        if rule == 'fixed':
            lam = {r['one']: float(r['lambdas']) for r in csv.DictReader(open(ht.CAL_PARAMS))}
            ks = [aa for aa in AAS if aa in prof['CALVADOS-2']]
            p(f'  CALVADOS-2 excess vs its lambda across types: r = {ea.pearson([prof["CALVADOS-2"][a] for a in ks], [lam[a] for a in ks]):+.2f}')
    per_res.close()
    with open(os.path.join(HERE, f'sim_null_by_type{TAG}.csv'), 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['rule', 'force_field', 'aa', 'n_residues', 'z_sim', 'z_null', 'excess'])
        w.writerows(csv_rows)
    open(os.path.join(HERE, f'sim_null_summary{TAG}.txt'), 'w').write('\n'.join(out) + '\n')


if __name__ == '__main__':
    main()
