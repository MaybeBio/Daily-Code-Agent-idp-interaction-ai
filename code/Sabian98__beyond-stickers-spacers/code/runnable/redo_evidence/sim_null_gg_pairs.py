"""
sim_null_gg_pairs.py — the Mpipi-GG residue-neutral null with Mpipi-GG's own pair sizes (2026-10-04).

Mpipi-GG's sigma table is not additive for proline: sigma(P, x) is 0.8-1.2 A larger than (sigma_P + sigma_x) / 2, which
is how ALBATROSS enlarges proline. sim_null.py grew its null chains with ped_null.grow_chain, whose hard-core distance is
0.8 * (sigma_i + sigma_j) / 2 from the self sizes, so the Mpipi-GG null let residues come closer to proline than
Mpipi-GG does. This rerun grows the same chains (same pool construction, seeds SEED * 1000 + k of the 908-protein run,
Rg matching and fixed-rule analysis as sim_null.py) with the hard-core distance 0.8 * sigma_ij taken from the table, for
the proteins of the sim_null_v3 shared set.
  ~/.venvs/bender/bin/python sim_null_gg_pairs.py check     # additive mixing on 3 proteins must reproduce the old null
  ~/.venvs/bender/bin/python sim_null_gg_pairs.py run [workers]  ->  sim_null_gg_pairsigma_per_residue.csv
"""
import bisect
import csv
import json
import math
import os
import random
import sys
from multiprocessing import Pool

os.environ.setdefault('SIM_NULL_N', '908')          # the 908-protein run behind sim_null_per_residue_all.csv
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import hub_frames as hf  # noqa: E402
import hub_frames_tests as ht  # noqa: E402
import ped_null as pn  # noqa: E402
import sim_null as sn  # noqa: E402
import sim_null_v3 as s3  # noqa: E402
from hub_z_v2_stdlib import quad_resid, zs  # noqa: E402

OUT = os.path.join(HERE, 'sim_null_gg_pairsigma_per_residue.csv')


def grow_chain_pairs(S, rng, max_tries=60, back=5, angles=(85.0, 150.0)):
    """ped_null.grow_chain with the hard core 0.8 * S[i][j] (pair sizes) instead of 0.8 * (sigma_i + sigma_j) / 2."""
    n = len(S)
    X, Y, Z = [0.0], [0.0], [0.0]
    X.append(pn.BOND)
    Y.append(0.0)
    Z.append(0.0)
    fails = 0
    while len(X) < n:
        i = len(X)
        ok = False
        for _ in range(max_tries):
            if i >= 2:
                ax, ay, az = X[i - 1] - X[i - 2], Y[i - 1] - Y[i - 2], Z[i - 1] - Z[i - 2]
            else:
                ax, ay, az = 1.0, 0.0, 0.0
            na = math.sqrt(ax * ax + ay * ay + az * az)
            ax, ay, az = ax / na, ay / na, az / na
            tx, ty, tz = (0.0, 1.0, 0.0) if abs(ax) > 0.9 else (1.0, 0.0, 0.0)
            ux, uy, uz = ay * tz - az * ty, az * tx - ax * tz, ax * ty - ay * tx
            nu = math.sqrt(ux * ux + uy * uy + uz * uz)
            ux, uy, uz = ux / nu, uy / nu, uz / nu
            vx, vy, vz = ay * uz - az * uy, az * ux - ax * uz, ax * uy - ay * ux
            theta = math.radians(rng.uniform(*angles))
            phi = rng.uniform(-math.pi, math.pi)
            ca, sa = math.cos(math.pi - theta), math.sin(math.pi - theta)
            dx = ca * ax + sa * (math.cos(phi) * ux + math.sin(phi) * vx)
            dy = ca * ay + sa * (math.cos(phi) * uy + math.sin(phi) * vy)
            dz = ca * az + sa * (math.cos(phi) * uz + math.sin(phi) * vz)
            x, y, z = X[i - 1] + pn.BOND * dx, Y[i - 1] + pn.BOND * dy, Z[i - 1] + pn.BOND * dz
            Si = S[i]
            clash = False
            for j in range(i - 2):
                dmin = 0.8 * Si[j]
                ex = x - X[j]
                if ex > dmin or ex < -dmin:
                    continue
                ey = y - Y[j]
                if ey > dmin or ey < -dmin:
                    continue
                ez = z - Z[j]
                if ex * ex + ey * ey + ez * ez < dmin * dmin:
                    clash = True
                    break
            if not clash:
                X.append(x)
                Y.append(y)
                Z.append(z)
                ok = True
                break
        if not ok:
            fails += 1
            if fails > 200:
                return None
            k = max(2, len(X) - back)
            del X[k:], Y[k:], Z[k:]
    return X, Y, Z


def one(job):
    uid, group, seed, additive = job
    try:
        tab = json.load(open(ht.SIGMA_TABLE))
        seq, rgs = sn.sim_rgs('MPIPI-GG', uid, group)
        sig, S = sn.ff_sigma('MPIPI-GG', seq, tab)
        if additive:                                   # reproduces sim_null.py exactly (check mode)
            S = [[(a + b) / 2 for b in sig] for a in sig]
        rng = random.Random(seed)
        pool = []
        for angles in ((85.0, 150.0), (100.0, 150.0), (115.0, 155.0)):
            target = len(pool) + sn.N_POOL // 3
            while len(pool) < target:
                ch = grow_chain_pairs(S, rng, angles=angles)
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
            adj = hf.contact_lists([X, Y, Z], n, None, 8.0)
            b = hf.betweenness(adj, n)
            tot = [a + v for a, v in zip(tot, b)]
        return uid, seq, quad_resid(zs([v / len(chosen) for v in tot]))
    except Exception as e:
        return uid, None, repr(e)


def jobs(additive=False):
    paired = [(r['UniProt_ID'], r['group']) for r in csv.DictReader(open(os.path.join(HERE, 'hub_frames_mpipi', 'ids.csv')))]
    paired = random.Random(sn.SEED).sample(paired, min(sn.N_SAMPLE, len(paired)))
    want = {u for u, g in s3.shared_set()}
    return [(u, g, sn.SEED * 1000 + k, additive) for k, (u, g) in enumerate(paired) if u in want]


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else 'check'
    if mode == 'check':
        old = {}
        for r in csv.DictReader(open(os.path.join(HERE, 'sim_null_per_residue_all.csv'))):
            if r['force_field'] == 'MPIPI-GG':
                old.setdefault(r['UniProt_ID'], {})[int(r['pos'])] = float(r['z_null_fixed'])
        for job in jobs(additive=True)[:3]:
            uid, seq, z = one(job)
            d = max(abs(z[i] - old[uid][i]) for i in range(len(seq)))
            print(f'{uid}: additive rerun vs sim_null_per_residue_all.csv, max |diff| = {d:.4f}', flush=True)
        return
    workers = int(sys.argv[2]) if len(sys.argv) > 2 else 8
    js = jobs()
    print(f'{len(js)} proteins, {workers} workers', flush=True)
    with open(OUT, 'w', newline='') as fh, Pool(workers) as pool:
        w = csv.writer(fh)
        w.writerow(['UniProt_ID', 'pos', 'aa', 'z_null_fixed'])
        for done, (uid, seq, z) in enumerate(pool.imap_unordered(one, js, chunksize=1), 1):
            if seq is None:
                print('FAILED', uid, z, flush=True)
                continue
            w.writerows([uid, i, a, f'{z[i]:.4f}'] for i, a in enumerate(seq))
            if done % 50 == 0:
                print(f'{done}/{len(js)}', flush=True)
    print('done', flush=True)


if __name__ == '__main__':
    main()
