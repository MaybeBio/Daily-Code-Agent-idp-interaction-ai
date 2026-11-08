"""
hub_frames.py — do residues that are hubs on average stay hubs frame by frame, and what contacts make them hubs?
(2026-10-02, option 2 of the "how do hub and anti-hub contacts evolve" question.)

The persistence pilot (bender_reanalysis/exploratory/hub_persistence_pilot.py) found that contacts decorrelate
within one saved frame (10-22 ps), so frames are independent snapshots. This analysis therefore asks about the
distribution over snapshots, not kinetics.

Sample: analysis_ok proteins of at most N_MAX residues (98% of them), N_PER_GROUP per BENDER group
(random, seed 0), CALVADOS-2 300 K
(BENDER-MT/BENDER_300K/<group>/<ID>/<ID>.dcd), FRAMES evenly spaced frames from the production half.
Per frame: contact graph with Calpha-Calpha < 0.8 nm and |i-j| >= 3 (bender_common rule; chain rebuilt from
minimum-image bonds first), exact Brandes betweenness (normalised as networkx), degree.
A residue is a hub in a frame if it is among the frame's top 10% by betweenness (only residues with betweenness > 0).

Per residue:
  class         hub / middle / anti-hub = top 10% / middle / bottom 10% of the protein by time-averaged betweenness
                (also class_evenframes: from even frames only, to pair with occupancy_oddframes without circularity),
                position-adjusted (residual of the within-protein z from a quadratic in chain position)
  occupancy     fraction of frames in which the residue is a hub (10% if hub status were random)
  participation fraction of frames with betweenness > 0
  burstiness    share of the residue's summed betweenness that comes from its 10% highest frames
                (0.10 if betweenness were the same in every frame; 1 if it all comes from a few frames)
  contacts      degree, long-range share (|i-j| >= 6), mean |i-j| and partner classes, separately for frames in which
                the residue is a hub and frames in which it is not
Per protein: Jaccard overlap of the hub sets of successive sampled frames (vs random sets of the same size), and a
check of the time-averaged betweenness against the stored v2 npz_residue means (50 frames, k = 30 pivots).

Standard library only. Read-only on the NAS; writes per_residue.csv, per_protein.csv to OUT.
  python3 hub_frames.py              # full sample
  python3 hub_frames.py --test 4     # timing on a few proteins
  python3 hub_frames.py --ids FILE --out DIR --resume   # keep the proteins already in DIR, label the rest
If the NAS stalls (all workers blocked, no CPU, no result for HUB_STALL_SECONDS), the run stops with exit code 75 and
can be continued with --resume; results do not depend on where a run was resumed.
"""
import array
import csv
import math
import multiprocessing
import os
import random
import signal
import struct
import subprocess
import sys
import time
import zipfile
import zlib
from collections import deque
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from hub_z_v2_stdlib import read_npy, zs, quad_resid  # noqa: E402

BIO = '/Volumes/User Homes/Nu Project/BENDER-BIO/BENDER_BIO.csv'
SCALAR = '/Volumes/User Homes/Nu Project/BENDER-MT/Traj Level Analysis/Analysis_300K_v2_graph/scalar_features_v2.csv'
TRAJ = '/Volumes/User Homes/Nu Project/BENDER-MT/BENDER_300K'
NPZ = '/Volumes/User Homes/Nu Project/BENDER-MT/Traj Level Analysis/Analysis_300K_v2_graph/npz_residue'
OUT = os.path.join(HERE, 'hub_frames')
N_PER_GROUP, SEED, N_MAX = 200, 0, 300
FRAMES = int(os.environ.get('HUB_FRAMES', '1000'))        # frames per protein (default 1,000)
CUT2 = 8.0 ** 2          # 0.8 nm in Å (DCD units), squared
MIN_SEP = int(os.environ.get('HUB_MIN_SEP', '3'))          # |i-j| >= 3 (v2); 2 reproduces the original pipeline
LR_SEP = 6
TOP = 0.10
PCLASS = {**dict.fromkeys('FWY', 'aro'), **dict.fromkeys('KR', 'pos'), **dict.fromkeys('DE', 'neg'),
          **dict.fromkeys('GSA', 'GSA'), **dict.fromkeys('ILVM', 'ILVM'), 'P': 'P'}
CLASSES = ['aro', 'pos', 'neg', 'GSA', 'ILVM', 'P', 'other']
AA3 = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H',
       'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
       'TYR': 'Y', 'VAL': 'V'}


def frames_for(n):
    """Frames per protein, the rule of all labels: 1,000 up to 300 residues, 500 up to 600, 250 above.
    Setting HUB_FRAMES fixes the number instead (sensitivity analyses only)."""
    if os.environ.get('HUB_FRAMES'):
        return int(os.environ['HUB_FRAMES'])
    return 1000 if n <= 300 else 500 if n <= 600 else 250


def top_hubs(bc, n_top, rng):
    """The n_top residues with the highest non-zero betweenness. Ties (equal to 1e-12) are broken at random with the
    caller's seeded rng, so residues tied at the cutoff are not favoured by their index (a stable sort would hand
    them to the N-terminal residue; ties at the cutoff occur in about a fifth of frames). n random numbers are drawn
    per call whatever the ties, so the stream, and the result, are reproducible for a given seed."""
    tb = [rng.random() for _ in range(len(bc))]
    return sorted((i for i in range(len(bc)) if bc[i] > 0), key=lambda i: (-round(bc[i], 12), tb[i]))[:n_top]


def tie_rng(key):
    """Deterministic tie-breaking generator for one protein or ensemble (same on every machine)."""
    return random.Random(zlib.crc32(str(key).encode()))


def expected_frames(n):
    return 10000 if n <= 196 else 15000 if n <= 312 else 20000


# ── DCD (CHARMM/OpenMM layout, little-endian) ────────────────────────────────
def dcd_open(path):
    fh = open(path, 'rb')
    h = fh.read(4096)
    n1 = struct.unpack('<i', h[:4])[0]
    assert h[4:8] == b'CORD' and n1 == 84, 'not a little-endian CHARMM DCD'
    ic = struct.unpack('<9if10i', h[8:88])
    off = 4 + n1 + 4
    n2 = struct.unpack('<i', h[off:off + 4])[0]
    off += 4 + n2 + 4
    natom = struct.unpack('<i', h[off + 4:off + 8])[0]
    off += 4 + 4 + 4
    has_cell = ic[10] == 1
    frame = (56 if has_cell else 0) + 3 * (8 + 4 * natom)
    n_frames = (os.path.getsize(path) - off) // frame
    if ic[0] > n_frames:            # header lists more frames than the file holds: an incomplete copy
        fh.close()
        raise ValueError(f'truncated DCD: header {ic[0]} frames, file holds {n_frames}')
    return fh, off, frame, natom, has_cell, n_frames


def dcd_frame(fh, off, frame, natom, has_cell, k):
    fh.seek(off + k * frame)
    b = fh.read(frame)
    p = 0
    box = None
    if has_cell:
        c = struct.unpack('<6d', b[4:52])          # A, gamma, B, beta, alpha, C
        box = (c[0], c[2], c[5])
        p = 56
    xyz = []
    for _ in range(3):
        a = array.array('f')
        a.frombytes(b[p + 4:p + 4 + 4 * natom])
        xyz.append(a)
        p += 8 + 4 * natom
    return xyz, box


def unwrap(xyz, box):
    """Rebuild the chain from minimum-image bond vectors; returns number of bonds that changed."""
    if not box or min(box) <= 0:
        return 0
    fixed = 0
    for d in range(3):
        L, a = box[d], xyz[d]
        for i in range(1, len(a)):
            b = a[i] - a[i - 1]
            s = round(b / L)
            if s:
                a[i] -= s * L
                fixed += 1
    return fixed


# ── graph ────────────────────────────────────────────────────────────────────
def contact_lists(xyz, n, C2=None, maxcut=8.0):
    """Contacts with |i-j| >= MIN_SEP: d < 0.8 nm, or, if C2 is given, d^2 < C2[i][j] (pair-specific, Å^2)."""
    if C2 is not None:
        return _contact_lists_pair(xyz, n, C2, maxcut)
    X, Y, Z = xyz
    adj = [[] for _ in range(n)]
    for i in range(n - MIN_SEP):
        xi, yi, zi = X[i], Y[i], Z[i]
        ai = adj[i]
        for j in range(i + MIN_SEP, n):
            dx = X[j] - xi
            if dx > 8.0 or dx < -8.0:
                continue
            dy = Y[j] - yi
            if dy > 8.0 or dy < -8.0:
                continue
            dz = Z[j] - zi
            if dx * dx + dy * dy + dz * dz < CUT2:
                ai.append(j)
                adj[j].append(i)
    return adj


def _contact_lists_pair(xyz, n, C2, maxcut):
    X, Y, Z = xyz
    adj = [[] for _ in range(n)]
    for i in range(n - MIN_SEP):
        xi, yi, zi = X[i], Y[i], Z[i]
        ai, ci = adj[i], C2[i]
        for j in range(i + MIN_SEP, n):
            dx = X[j] - xi
            if dx > maxcut or dx < -maxcut:
                continue
            dy = Y[j] - yi
            if dy > maxcut or dy < -maxcut:
                continue
            dz = Z[j] - zi
            if dx * dx + dy * dy + dz * dz < ci[j]:
                ai.append(j)
                adj[j].append(i)
    return adj


def betweenness(adj, n):
    """Exact Brandes, unweighted undirected, normalised like networkx (raw / ((n-1)(n-2)))."""
    bc = [0.0] * n
    for s in range(n):
        if not adj[s]:
            continue
        S = []
        preds = {s: []}
        sigma = {s: 1.0}
        dist = {s: 0}
        Q = deque([s])
        while Q:
            v = Q.popleft()
            S.append(v)
            dv = dist[v] + 1
            for w in adj[v]:
                if w not in dist:
                    dist[w] = dv
                    sigma[w] = 0.0
                    preds[w] = []
                    Q.append(w)
                if dist[w] == dv:
                    sigma[w] += sigma[v]
                    preds[w].append(v)
        delta = dict.fromkeys(S, 0.0)
        while S:
            w = S.pop()
            dw = (1.0 + delta[w]) / sigma[w]
            for v in preds[w]:
                delta[v] += sigma[v] * dw
            if w != s:
                bc[w] += delta[w]
    scale = 1.0 / ((n - 1) * (n - 2)) if n > 2 else 0.0
    return [x * scale for x in bc]


# ── per protein ──────────────────────────────────────────────────────────────
def find_traj(uid, group):
    """BENDER_300K/<group>/<ID>; falls back to the other lineage folders (some IDs moved in the taxon fix)."""
    d = os.path.join(TRAJ, group, uid)
    if os.path.exists(os.path.join(d, f'{uid}.dcd')):
        return d
    for g in ('Bacteria', 'Fungi', 'Mammals', 'Plants', 'Protists', 'Viruses'):
        d2 = os.path.join(TRAJ, g, uid)
        if os.path.exists(os.path.join(d2, f'{uid}.dcd')):
            return d2
    return d


def one(job):
    uid, group = job
    d = find_traj(uid, group)
    path = os.path.join(d, f'{uid}.dcd')
    try:
        seq = ''.join(AA3.get(line[17:20], 'X') for line in open(os.path.join(d, 'top.pdb'))
                      if line.startswith('ATOM') and line[12:16].strip() == 'CA')
        fh, off, frame, natom, has_cell, n_frames = dcd_open(path)
    except Exception as e:
        return ('fail', uid, group, repr(e))
    n = natom
    if len(seq) != n:
        return ('fail', uid, group, f'top {len(seq)} vs dcd {n} atoms')
    if n_frames < expected_frames(n):
        return ('fail', uid, group, f'truncated {n_frames} < {expected_frames(n)}')
    start = n_frames // 2
    nf = frames_for(n)
    idx = [start + (k * (n_frames - start)) // nf for k in range(nf)]
    try:
        return analyse(uid, group, seq, lambda k: dcd_frame(fh, off, frame, natom, has_cell, k), n_frames, idx)
    finally:
        fh.close()


def analyse(uid, group, seq, get_frame, n_frames, idx, ref=True, C2=None, maxcut=8.0):
    """Force-field independent part. get_frame(k) -> ([X, Y, Z] in Å, (Lx, Ly, Lz) in Å or None).
    C2: optional pair-specific squared contact cutoffs (Å^2), e.g. for the bead-size test."""
    FRAMES = len(idx)
    n = len(seq)
    n_top = max(1, round(TOP * n))
    pcls = [PCLASS.get(a, 'other') for a in seq]

    btw_f = [[0.0] * FRAMES for _ in range(n)]
    hub_n = [0] * n
    pos_n = [0] * n
    acc = {key: [0.0] * n for key in ('deg_hub', 'deg_non', 'lr_hub', 'lr_non', 'sep_hub', 'sep_non')}
    part = {(h, c): [0] * n for h in ('hub', 'non') for c in CLASSES}
    fixed_bonds = 0
    jacc, jacc_null = [], []
    rng = tie_rng(uid)
    hub_sets = []
    prev = None
    for f, k in enumerate(idx):
        xyz, box = get_frame(k)
        fixed_bonds += unwrap(xyz, box)
        adj = contact_lists(xyz, n, C2, maxcut)
        bc = betweenness(adj, n)
        hubs = set(top_hubs(bc, n_top, rng))
        hub_sets.append(hubs)
        if prev is not None and (hubs or prev):
            jacc.append(len(hubs & prev) / len(hubs | prev))
            a, b2 = len(hubs), len(prev)          # random sets of these sizes: E|A&B| = ab/n, J ~ E / (a + b - E)
            e = a * b2 / n
            jacc_null.append(e / (a + b2 - e) if a + b2 - e > 0 else 0.0)
        prev = hubs
        for i in range(n):
            b = bc[i]
            btw_f[i][f] = b
            if b > 0:
                pos_n[i] += 1
            h = 'hub' if i in hubs else 'non'
            if h == 'hub':
                hub_n[i] += 1
            nb = adj[i]
            acc[f'deg_{h}'][i] += len(nb)
            for j in nb:
                s = j - i if j > i else i - j
                acc[f'sep_{h}'][i] += s
                if s >= LR_SEP:
                    acc[f'lr_{h}'][i] += 1
                part[(h, pcls[j])][i] += 1

    mean_b = [sum(v) / FRAMES for v in btw_f]
    z = zs(mean_b)
    zr = quad_resid(z) if z else None
    if zr is None:
        return ('fail', uid, group, 'flat betweenness')
    cls = classify(zr, n, n_top)
    # cross-validated version: class from even frames only, occupancy measured on odd frames (no circularity)
    ze = zs([sum(v[0::2]) for v in btw_f])
    zre = quad_resid(ze) if ze else None
    cls_even = classify(zre, n, n_top) if zre else [''] * n
    zo = zs([sum(v[1::2]) for v in btw_f])
    zro = quad_resid(zo) if zo else None
    sh_z = pearson(zre, zro) if zre and zro else ''     # split-half of the position-adjusted z (the z_btw label)

    # check against the stored v2 means (50 frames, k = 30 pivots)
    r_check = ''
    try:
        if ref:
            zf = zipfile.ZipFile(f'{NPZ}/{uid}_residue_features.npz')
            stored = read_npy(zf.read('betweenness_centrality_mean.npy'))[:n]
            if len(stored) == n:
                r_check = pearson(mean_b, stored)
    except Exception:
        pass
    # split-half reliability (even vs odd sampled frames) of mean betweenness and of hub occupancy
    sh_b = pearson([sum(v[0::2]) for v in btw_f], [sum(v[1::2]) for v in btw_f])
    occ_e = [sum(1 for f in range(0, FRAMES, 2) if i in hub_sets[f]) for i in range(n)]
    occ_o = [sum(1 for f in range(1, FRAMES, 2) if i in hub_sets[f]) for i in range(n)]
    sh_o = pearson(occ_e, occ_o)

    rows = []
    n_top_frames = max(1, round(TOP * FRAMES))
    for i in range(n):
        v = btw_f[i]
        tot = sum(v)
        burst = sum(sorted(v, reverse=True)[:n_top_frames]) / tot if tot > 0 else ''
        mu = tot / FRAMES
        sd = math.sqrt(sum((x - mu) ** 2 for x in v) / FRAMES)
        nh, nn = hub_n[i], FRAMES - hub_n[i]
        c_h = sum(part[('hub', c)][i] for c in CLASSES)
        c_n = sum(part[('non', c)][i] for c in CLASSES)
        row = [uid, group, n, i, seq[i], cls[i], f'{zr[i]:.4f}', f'{mu:.6g}', f'{sd / mu:.4f}' if mu > 0 else '',
               f'{nh / FRAMES:.4f}', f'{pos_n[i] / FRAMES:.4f}', f'{burst:.4f}' if burst != '' else '',
               f'{(acc["deg_hub"][i] + acc["deg_non"][i]) / FRAMES:.4f}',
               f'{acc["deg_hub"][i] / nh:.4f}' if nh else '', f'{acc["deg_non"][i] / nn:.4f}' if nn else '',
               f'{acc["lr_hub"][i] / c_h:.4f}' if c_h else '', f'{acc["lr_non"][i] / c_n:.4f}' if c_n else '',
               f'{acc["sep_hub"][i] / c_h:.3f}' if c_h else '', f'{acc["sep_non"][i] / c_n:.3f}' if c_n else '',
               nh, c_h, c_n, cls_even[i], f'{occ_o[i] / len(range(1, FRAMES, 2)):.4f}']
        row += [part[('hub', c)][i] for c in CLASSES] + [part[('non', c)][i] for c in CLASSES]
        rows.append(row)
    exp_j = sum(jacc_null) / len(jacc_null) if jacc_null else 0.0   # random sets of the observed hub-set sizes
    fmt = lambda x: f'{x:.4f}' if isinstance(x, float) else ''  # noqa: E731
    prot = [uid, group, n, n_frames, FRAMES, fixed_bonds, f'{sum(jacc) / len(jacc):.4f}' if jacc else '',
            f'{exp_j:.4f}', fmt(r_check), fmt(sh_b), fmt(2 * sh_b / (1 + sh_b)) if isinstance(sh_b, float) else '',
            fmt(sh_o), fmt(sh_z)]
    return ('ok', uid, group, (prot, rows))


def classify(zr, n, n_top):
    """hub / middle / anti = top / middle / bottom n_top residues by position-adjusted z."""
    rank = sorted(range(n), key=lambda i: -zr[i])
    cls = ['middle'] * n
    for i in rank[:n_top]:
        cls[i] = 'hub'
    for i in rank[-n_top:]:
        cls[i] = 'anti'
    return cls


def pearson(a, b):
    n = len(a)
    ma, mb = sum(a) / n, sum(b) / n
    sa = math.sqrt(sum((x - ma) ** 2 for x in a))
    sb = math.sqrt(sum((y - mb) ** 2 for y in b))
    if sa == 0 or sb == 0:
        return ''
    return sum((x - ma) * (y - mb) for x, y in zip(a, b)) / (sa * sb)


RES_HEADER = (['UniProt_ID', 'group', 'N', 'pos', 'aa', 'class', 'z_posadj', 'mean_btw', 'cv_btw', 'occupancy',
               'participation', 'burstiness', 'deg_mean', 'deg_hubframes', 'deg_nonhubframes', 'lr_share_hubframes',
               'lr_share_nonhubframes', 'mean_sep_hubframes', 'mean_sep_nonhubframes', 'n_hubframes',
               'contacts_hubframes', 'contacts_nonhubframes', 'class_evenframes', 'occupancy_oddframes']
              + [f'hub_{c}' for c in CLASSES] + [f'non_{c}' for c in CLASSES])
PROT_HEADER = ['UniProt_ID', 'group', 'N', 'n_frames_dcd', 'frames_used', 'bonds_unwrapped', 'jaccard_successive',
               'jaccard_random_approx', 'r_vs_v2_npz_mean', 'splithalf_r_btw', 'splithalf_spearman_brown_btw',
               'splithalf_r_occupancy', 'splithalf_r_zposadj']


def sample():
    ok = {r['UniProt_ID'] for r in csv.DictReader(open(BIO)) if r['analysis_ok'] == 'True'}
    by = {}
    for r in csv.DictReader(open(SCALAR)):
        if r['UniProt_ID'] in ok and r['status'] == 'ok' and (os.environ.get('HUB_ALL') == '1' or int(r['n_residues']) <= N_MAX):
            by.setdefault(r['folder'], []).append(r['UniProt_ID'])
    rng = random.Random(SEED)
    jobs = []
    for g in sorted(by):
        ids = sorted(by[g])
        rng.shuffle(ids)
        jobs += [(u, g) for u in (ids if os.environ.get('HUB_ALL') == '1' else ids[:N_PER_GROUP])]
    skip = set()
    for d in filter(None, os.environ.get('HUB_SKIP', '').split(',')):
        f = os.path.join(HERE, d, 'per_protein.csv')
        if os.path.exists(f):
            skip |= {r['UniProt_ID'] for r in csv.DictReader(open(f))}
    jobs = [j for j in jobs if j[0] not in skip]
    # longest first, so the slow chains do not end up running alone at the end
    nres = {r['UniProt_ID']: int(r['n_residues']) for r in csv.DictReader(open(SCALAR))}
    return sorted(jobs, key=lambda j: -nres.get(j[0], 0))


def arg(name, default=None):
    return sys.argv[sys.argv.index(name) + 1] if name in sys.argv else default


STALL = int(os.environ.get('HUB_STALL_SECONDS', '900'))
STALL_EXIT = 75


def pool_cpu(pool):
    """CPU seconds used so far by the pool's workers (from ps), to tell a stalled pool from a busy one."""
    out = subprocess.run(['ps', '-o', 'time=', '-p', ','.join(str(w.pid) for w in pool._pool)],
                         capture_output=True, text=True).stdout
    total = 0.0
    for t in out.split():
        days, _, t = t.rpartition('-')
        sec = 0.0
        for part in t.split(':'):
            sec = sec * 60 + float(part)
        total += sec + 86400 * int(days or 0)
    return total


def watched(pool, worker, jobs, stall=STALL, poll=60):
    """pool.imap_unordered, except that when no result arrives and the workers use no CPU for `stall` seconds (all of
    them blocked on a stalled network share), the workers are killed and the process exits with code 75 so that the
    run can be resumed. Callers must flush their output after every result."""
    it = pool.imap_unordered(worker, jobs, chunksize=1)
    cpu, quiet = -1.0, time.time()
    while True:
        try:
            res = it.next(timeout=poll)
        except StopIteration:
            return
        except multiprocessing.TimeoutError:
            now = pool_cpu(pool)
            if now > cpu + 1.0:
                cpu, quiet = now, time.time()
            elif time.time() - quiet > stall:
                print(f'stalled: no result and no worker CPU for {stall} s; exiting with code {STALL_EXIT} so the '
                      'run can be resumed', flush=True)
                for w in pool._pool:
                    try:
                        os.kill(w.pid, signal.SIGKILL)
                    except OSError:
                        pass
                os._exit(STALL_EXIT)
            continue
        quiet = time.time()
        yield res


def repair(pp, pr):
    """Before resuming: drop a cut-off last line and the residue rows of proteins without a protein row. Returns the
    IDs already labelled."""
    for path in (pp, pr):
        data = open(path, 'rb').read()
        if data and not data.endswith(b'\n'):
            open(path, 'wb').write(data[:data.rfind(b'\n') + 1])
    done = {r['UniProt_ID'] for r in csv.DictReader(open(pp))}
    if not os.path.getsize(pr):
        return done
    dropped = 0
    with open(pr, newline='') as src, open(pr + '.tmp', 'w', newline='') as dst:
        rd, wt = csv.reader(src), csv.writer(dst)
        wt.writerow(next(rd))
        for r in rd:
            if r and r[0] in done:
                wt.writerow(r)
            else:
                dropped += 1
    if dropped:
        os.replace(pr + '.tmp', pr)
        print(f'repaired {pr}: dropped {dropped} rows of unfinished proteins', flush=True)
    else:
        os.remove(pr + '.tmp')
    return done


def main(worker=None, jobs=None, out=None):
    """--ids FILE (UniProt_ID,group rows) replaces the random sample; --out DIR names the output folder; --resume keeps
    the proteins already labelled in DIR."""
    worker = worker or one
    if jobs is None:
        ids = arg('--ids')
        jobs = [(r['UniProt_ID'], r['group']) for r in csv.DictReader(open(ids))] if ids else sample()
    out = out or os.path.join(HERE, arg('--out', 'hub_frames'))
    if '--test' in sys.argv:
        k = int(arg('--test'))
        for job in jobs[::max(1, len(jobs) // k)][:k]:
            t = time.time()
            res = worker(job)
            print(job, res[0], f'{time.time() - t:.1f}s',
                  res[3][0] if res[0] == 'ok' else res[3])
        return
    OUT = out
    os.makedirs(OUT, exist_ok=True)
    pp, pr = f'{OUT}/per_protein.csv', f'{OUT}/per_residue.csv'
    done = set()
    if '--resume' in sys.argv and os.path.exists(pp) and os.path.exists(pr):
        done = repair(pp, pr)
        jobs = [j for j in jobs if j[0] not in done]
        print(f'resuming: {len(done)} already labelled, {len(jobs)} to do', flush=True)
    mode = 'a' if done else 'w'
    fails = []
    t0 = time.time()
    with open(pr, mode, newline='') as fr, open(pp, mode, newline='') as fp:
        wr, wp = csv.writer(fr), csv.writer(fp)
        if not done:
            wr.writerow(RES_HEADER)
            wp.writerow(PROT_HEADER)
        with Pool(int(os.environ.get('WORKERS', '8'))) as pool:
            for k, res in enumerate(watched(pool, worker, jobs), 1):
                if res[0] == 'ok':
                    prot, rows = res[3]
                    wr.writerows(rows)
                    wp.writerow(prot)              # protein row last: marked done only once its residues are written
                    fr.flush()
                    fp.flush()
                else:
                    fails.append(res[1:])
                if k % 50 == 0:
                    print(f'{k}/{len(jobs)}  {time.time() - t0:.0f}s', flush=True)
    with open(f'{OUT}/failed.csv', 'w', newline='') as fh:
        csv.writer(fh).writerows([('UniProt_ID', 'group', 'reason')] + fails)
    print(f'done: {len(jobs) - len(fails)} ok, {len(fails)} failed, {time.time() - t0:.0f}s -> {OUT}'
          + (f' ({len(done)} kept from before)' if done else ''))


if __name__ == '__main__':
    main()
