"""
sim_null_v3.py — the residue-neutral null of sim_null.py on the sequences simulated with all three force fields, with
the simulation side taken from the v3 labels (hub_frames_v3_cal, hub_frames_v3_mpipi, hub_frames_v3_recharged).

Set: the proteins of the 908 CALVADOS-2 / Mpipi-GG paired set that also have a labelled Mpipi-Recharged run
(507 on 2026-10-04). The CALVADOS-2 and Mpipi-GG nulls are the ones already computed by sim_null.py for the 908
(sim_null_per_residue_all.csv, column z_null_fixed; the null depends only on sequence, bead sizes and the Rg of the
simulation). The Mpipi-Recharged null is computed here exactly as in sim_null.py (2,000 self-avoiding C-alpha chains
with the Mpipi-Recharged bead sizes (OpenMpipi) in three stiffnesses, Rg-quantile matching to 500 frames of the simulation, fixed 0.8 nm rule,
mean betweenness over the matched chains, position-adjusted z), written incrementally so that a run can resume.
Excess = simulation z_posadj - null z, per residue type and class, for all three force fields on the same proteins.
  ~/.venvs/bender/bin/python sim_null_v3.py run [workers]   ->  sim_null_v3_recharged_null.csv (resumable)
  ~/.venvs/bender/bin/python sim_null_v3.py summarise       ->  sim_null_v3_summary.txt, sim_null_v3_by_type.csv,
                                                                sim_null_v3_per_residue.csv
"""
import bisect
import csv
import json
import os
import random
import sys
import time
import zlib
from collections import defaultdict
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import exp_anchor as ea  # noqa: E402
import hub_frames as hf  # noqa: E402
import hub_frames_tests as ht  # noqa: E402
import ped_null as pn  # noqa: E402
import sim_null as sn  # noqa: E402  (rg3, sim_rgs, ff_sigma, N_POOL, N_NULL, SEED, AAS, CLASS)
from hub_z_v2_stdlib import quad_resid, zs  # noqa: E402

V3 = {'CALVADOS-2': 'hub_frames_v3_cal', 'MPIPI-GG': 'hub_frames_v3_mpipi', 'Mpipi-Recharged': 'hub_frames_v3_recharged'}
NULL_OUT = os.path.join(HERE, 'sim_null_v3_recharged_null.csv')


def shared_set(labelled=True):
    """Paired-set proteins with a Mpipi-Recharged run: labelled ones (summary) or those with a trajectory on the NAS
    (labels_input_recharged_present.fasta, for running the null before the labels are finished)."""
    paired = [(r['UniProt_ID'], r['group']) for r in csv.DictReader(open(os.path.join(HERE, 'hub_frames_mpipi', 'ids.csv')))]
    if labelled:
        rech = {r['UniProt_ID'] for r in csv.DictReader(open(os.path.join(HERE, V3['Mpipi-Recharged'], 'per_protein.csv')))}
    else:
        rech = {l[1:].split('|')[0].strip() for l in open(os.path.join(HERE, 'labels_input_recharged_present.fasta'))
                if l.startswith('>')}
    return [(u, g) for u, g in paired if u in rech]


def rc_sigma_table():
    """Mpipi-Recharged bead sizes in angstrom (as the Mpipi-GG table), from OpenMpipi recharged_params.txt."""
    sig = json.load(open(os.path.join(HERE, 'mpipi_recharged_parameters.json')))['sigma_nm']
    return {k: 10.0 * v for k, v in sig.items()}


def one(job):
    uid, group, seed = job
    try:
        tab = rc_sigma_table()   # Recharged's own bead sizes; the Mpipi-GG table made proline 33% too large (fixed 2026-10-04)
        seq, rgs = sn.sim_rgs('Mpipi-Recharged', uid, group)
        sig, _ = sn.ff_sigma('Mpipi-Recharged', seq, tab)
        rng = random.Random(seed)
        pool = []
        for angles in ((85.0, 150.0), (100.0, 150.0), (115.0, 155.0)):
            target = len(pool) + sn.N_POOL // 3
            while len(pool) < target:
                ch = pn.grow_chain(sig, rng, angles=angles)
                if ch is not None:
                    pool.append(ch)
        prg = sorted((sn.rg3(*ch), k) for k, ch in enumerate(pool))
        vals = [v for v, _ in prg]
        chosen = []
        for r in rgs:
            k = bisect.bisect_left(vals, r)
            k = min(range(max(0, k - 1), min(len(vals), k + 1)), key=lambda t: abs(vals[t] - r))
            chosen.append(pool[prg[k][1]])
        n = len(seq)
        tot = [0.0] * n
        for X, Y, Z in chosen:
            b = hf.betweenness(hf.contact_lists([X, Y, Z], n, None, 8.0), n)
            tot = [a + v for a, v in zip(tot, b)]
        zn = quad_resid(zs([v / len(chosen) for v in tot]))
        return {'uid': uid, 'seq': seq, 'z': zn, 'rg_sim': sum(rgs) / len(rgs),
                'rg_null': sum(sn.rg3(*ch) for ch in chosen) / len(chosen)}
    except Exception as e:
        return {'uid': uid, 'error': repr(e)}


def run(workers):
    ids = shared_set(labelled=False)
    done = set()
    if os.path.exists(NULL_OUT):
        data = open(NULL_OUT, 'rb').read()
        if data and not data.endswith(b'\n'):
            open(NULL_OUT, 'wb').write(data[:data.rfind(b'\n') + 1])      # drop a cut-off last line
        count, length = defaultdict(int), {}
        for r in csv.DictReader(open(NULL_OUT)):
            count[r['UniProt_ID']] += 1
            length[r['UniProt_ID']] = int(r['N'])
        done = {u for u in count if count[u] == length[u]}               # a protein counts once all its rows are in
    jobs = [(u, g, sn.SEED * 1000 + zlib.crc32(u.encode()) % 1000003) for u, g in ids if u not in done]
    print(f'{len(ids)} shared proteins, {len(done)} done, {len(jobs)} to run with {workers} workers', flush=True)
    new = not os.path.exists(NULL_OUT)
    t0, fails = time.time(), []
    with open(NULL_OUT, 'a', newline='') as fh, Pool(workers) as pool:
        w = csv.writer(fh)
        if new:
            w.writerow(['UniProt_ID', 'pos', 'aa', 'N', 'z_null_fixed', 'rg_sim', 'rg_null'])
        for k, r in enumerate(pool.imap_unordered(one, jobs, chunksize=1), 1):
            if 'error' in r:
                fails.append((r['uid'], r['error']))
            else:
                n = len(r['seq'])
                w.writerows([[r['uid'], i, a, n, f"{r['z'][i]:.4f}", f"{r['rg_sim']:.2f}", f"{r['rg_null']:.2f}"]
                             for i, a in enumerate(r['seq'])])
                fh.flush()
            if k % 25 == 0:
                print(f'{k}/{len(jobs)}  {time.time() - t0:.0f}s', flush=True)
    print(f'done: {len(jobs) - len(fails)} new, {len(fails)} failed {fails[:5]}, {time.time() - t0:.0f}s', flush=True)


def summarise():
    out = []

    def p(s=''):
        out.append(s)
        print(s, flush=True)
    ids = [u for u, g in shared_set()]
    want = set(ids)
    null = {'Mpipi-Recharged': defaultdict(dict), 'CALVADOS-2': defaultdict(dict), 'MPIPI-GG': defaultdict(dict)}
    for r in csv.DictReader(open(NULL_OUT)):
        null['Mpipi-Recharged'][r['UniProt_ID']][int(r['pos'])] = float(r['z_null_fixed'])
    for r in csv.DictReader(open(os.path.join(HERE, 'sim_null_per_residue_all.csv'))):
        if r['force_field'] in ('CALVADOS-2', 'MPIPI-GG') and r['UniProt_ID'] in want and r['z_null_fixed']:
            null[r['force_field']][r['UniProt_ID']][int(r['pos'])] = float(r['z_null_fixed'])
    # Mpipi-GG null with its own (non-additive) pair sizes, which enlarge proline (sim_null_gg_pairs.py, 2026-10-04)
    gg_pairs = os.path.join(HERE, 'sim_null_gg_pairsigma_per_residue.csv')
    if os.path.exists(gg_pairs):
        null['MPIPI-GG'] = defaultdict(dict)
        for r in csv.DictReader(open(gg_pairs)):
            if r['UniProt_ID'] in want:
                null['MPIPI-GG'][r['UniProt_ID']][int(r['pos'])] = float(r['z_null_fixed'])
    sim = {ff: defaultdict(dict) for ff in V3}
    for ff, d in V3.items():
        for r in csv.DictReader(open(os.path.join(HERE, d, 'per_residue.csv'))):
            if r['UniProt_ID'] in want:
                sim[ff][r['UniProt_ID']][int(r['pos'])] = (r['aa'], float(r['z_posadj']))
    # proteins with all six profiles complete and of one length
    keep = []
    for u in ids:
        lens = {len(sim[ff].get(u, {})) for ff in V3} | {len(null[ff].get(u, {})) for ff in V3}
        if len(lens) == 1 and 0 not in lens:
            keep.append(u)
    p(f'Simulation vs residue-neutral null on the {len(keep)} sequences simulated with all three force fields '
      f'(of {len(ids)} shared; v3 labels; fixed 0.8 nm rule). Mean position-adjusted z: simulation / null / excess')
    acc = {(ff, w): defaultdict(list) for ff in V3 for w in ('sim', 'null')}
    prot = defaultdict(list)
    with open(os.path.join(HERE, 'sim_null_v3_per_residue.csv'), 'w', newline='') as fh:
        wr = csv.writer(fh)
        wr.writerow(['force_field', 'UniProt_ID', 'pos', 'aa', 'N', 'z_sim_fixed', 'z_null_fixed'])
        for ff in V3:
            for u in keep:
                n = len(sim[ff][u])
                cls_s, cls_n = defaultdict(list), defaultdict(list)
                for i in range(n):
                    a, x = sim[ff][u][i]
                    y = null[ff][u][i]
                    wr.writerow([ff, u, i, a, n, f'{x:.4f}', f'{y:.4f}'])
                    acc[(ff, 'sim')][a].append(x)
                    acc[(ff, 'null')][a].append(y)
                    if a in sn.CLASS:
                        cls_s[sn.CLASS[a]].append(x)
                        cls_n[sn.CLASS[a]].append(y)
                prot[ff].append({k: sum(cls_s[k]) / len(cls_s[k]) - sum(cls_n[k]) / len(cls_n[k]) for k in cls_s})
    prof = {ff: {} for ff in V3}
    rows = []
    p('  aa ' + ''.join(f'{ff:>26s}' for ff in V3))
    for aa in sn.AAS:
        line = f'  {aa}  '
        for ff in V3:
            s, n_ = acc[(ff, 'sim')][aa], acc[(ff, 'null')][aa]
            ms, mn = sum(s) / len(s), sum(n_) / len(n_)
            prof[ff][aa] = ms - mn
            line += f'   {ms:+.2f} / {mn:+.2f} / {ms - mn:+.2f}  '
            rows.append(['fixed', ff, aa, len(s), f'{ms:.4f}', f'{mn:.4f}', f'{ms - mn:.4f}'])
        p(line)
    ffs = list(V3)
    for i in range(len(ffs)):
        for j in range(i + 1, len(ffs)):
            a, b = prof[ffs[i]], prof[ffs[j]]
            p(f'  excess-profile agreement {ffs[i]} vs {ffs[j]}: r = {ea.pearson([a[k] for k in sn.AAS], [b[k] for k in sn.AAS]):+.2f}'
              f' ({ea.pearson([a[k] for k in sn.AAS if k != "P"], [b[k] for k in sn.AAS if k != "P"]):+.2f} without P)')
    for ff in V3:
        pr = prot[ff]
        both = [d['FWY'] - d['DE'] for d in pr if 'FWY' in d and 'DE' in d]
        k = sum(v > 0 for v in both)
        p(f'  {ff}: per-protein (FWY excess) - (DE excess) > 0 in {k}/{len(both)} ({k / max(1, len(both)):.0%}); '
          + '  '.join(f'{c} excess {sum(d[c] for d in pr if c in d) / max(1, sum(1 for d in pr if c in d)):+.2f}'
                      for c in ('FWY', 'DE', 'KR', 'GSA', 'ILVM')))
    lam = {r['one']: float(r['lambdas']) for r in csv.DictReader(open(ht.CAL_PARAMS))}
    p(f'  CALVADOS-2 excess vs its lambda across types: r = '
      f'{ea.pearson([prof["CALVADOS-2"][a] for a in sn.AAS], [lam[a] for a in sn.AAS]):+.2f}')
    with open(os.path.join(HERE, 'sim_null_v3_by_type.csv'), 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['rule', 'force_field', 'aa', 'n_residues', 'z_sim', 'z_null', 'excess'])
        w.writerows(rows)
    open(os.path.join(HERE, 'sim_null_v3_summary.txt'), 'w').write('\n'.join(out) + '\n')


if __name__ == '__main__':
    if sys.argv[1:2] == ['run']:
        run(int(sys.argv[2]) if len(sys.argv) > 2 else 9)
    else:
        summarise()
