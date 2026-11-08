"""
hub_frames_tests.py — two checks on the consensus hubs of §3d (see consensus_hubs.py).

1. Bead-size test. Both CALVADOS-2 and MPIPI-GG give small residues (G, S, A) small beads, and a fixed
   Calpha-Calpha < 0.8 nm cutoff lets small beads touch more neighbours. Here the contact rule is d < c * sigma_ij
   (sigma_ij: CALVADOS-2 (sigma_i + sigma_j) / 2 from residues_CALVADOS2.csv; MPIPI-GG the wf/cut sigma of the pair
   from each run's run.in), with c calibrated per force field so that the total number of contacts equals that of the
   0.8 nm rule on a calibration sample. Only the size bias changes, not how dense the graphs are.
2. Third model. Mpipi-Recharged (OpenMpipi, BENDER protocol: 10,000 frames, first half equilibration), 301 proteins
   in Nu Project/MPIPI Recharged/300K/<ID>/<ID>.dcd, standard 0.8 nm rule; CALVADOS-2 is run on the same IDs.

  python3 hub_frames_tests.py calibrate        # prints c for each force field (writes sizeaware_c.json)
  python3 hub_frames_tests.py size             # -> hub_frames_sizeaware_calvados/, hub_frames_sizeaware_mpipi/
  python3 hub_frames_tests.py recharged        # -> hub_frames_recharged/ and recharged_ids.csv
  python3 hub_frames.py --ids recharged_ids.csv --out hub_frames_recharged_pair_calvados
  python3 hub_frames_tests.py size_recharged   # bead-size rule for Mpipi-Recharged (Mpipi sigma table) and CALVADOS-2
"""
import csv
import json
import os
import random
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import hub_frames as hf  # noqa: E402
import hub_frames_mpipi as hm  # noqa: E402
import xtc_stdlib as xtc  # noqa: E402

CAL_PARAMS = '/Users/taseef/Downloads/Scripts/CALVADOS-Slab/input/residues_CALVADOS2.csv'
RECHARGED = '/Volumes/User Homes/Nu Project/MPIPI Recharged/300K'
C_FILE = os.path.join(HERE, 'sizeaware_c.json')
PAIRED_IDS = os.path.join(HERE, 'hub_frames_mpipi', 'ids.csv')


# ── bead sizes ───────────────────────────────────────────────────────────────
def cal_sigma_matrix(seq):
    s = {r['one']: float(r['sigmas']) * 10.0 for r in csv.DictReader(open(CAL_PARAMS))}   # Å
    v = [s[a] for a in seq]
    return [[(v[i] + v[j]) / 2.0 for j in range(len(v))] for i in range(len(v))]


def mpi_sigma_matrix(uid):
    d = os.path.join(hm.RUNS, uid)
    data = open(os.path.join(d, 'seq.dat')).read()
    atoms = data.split('Atoms')[1].strip().splitlines()
    types = []
    for line in atoms:
        t = line.split()
        if len(t) >= 7 and t[0].isdigit():
            types.append((int(t[0]), int(t[2])))
        elif types:
            break
    types = [ty for _, ty in sorted(types)]
    pair = {}
    for m in re.finditer(r'^pair_coeff\s+(\d+)\s+(\d+)\s+wf/cut\s+[\d.eE+-]+\s+([\d.eE+-]+)',
                         open(os.path.join(d, 'run.in')).read(), re.M):
        i, j, sg = int(m.group(1)), int(m.group(2)), float(m.group(3))
        pair[(i, j)] = pair[(j, i)] = sg
    return [[pair[(ti, tj)] for tj in types] for ti in types]


def c2_matrix(S, c):
    return [[(c * x) ** 2 for x in row] for row in S], c * max(max(row) for row in S)


# ── frame providers ──────────────────────────────────────────────────────────
def cal_frames(uid, group, base=None):
    """(seq, get_frame, n_frames, idx, closer) for a CALVADOS-2-style OpenMM DCD (production half)."""
    if base is None:
        d = hf.find_traj(uid, group)
    else:
        d = os.path.join(base, uid)
    seq = ''.join(hf.AA3.get(line[17:20], 'X') for line in open(os.path.join(d, 'top.pdb'))
                  if line.startswith('ATOM') and line[12:16].strip() == 'CA')
    fh, off, frame, natom, has_cell, n_frames = hf.dcd_open(os.path.join(d, f'{uid}.dcd'))
    if len(seq) != natom:
        fh.close()
        raise ValueError(f'top {len(seq)} vs dcd {natom} atoms')
    if n_frames < hf.expected_frames(natom):
        fh.close()
        raise ValueError(f'truncated {n_frames} < {hf.expected_frames(natom)}')
    start = n_frames // 2
    nf = hf.frames_for(natom)
    idx = [start + (k * (n_frames - start)) // nf for k in range(nf)]
    return seq, (lambda k: hf.dcd_frame(fh, off, frame, natom, has_cell, k)), n_frames, idx, fh.close


def mpi_frames(uid):
    d = os.path.join(hm.RUNS, uid)
    seq = hm.read_sequence(os.path.join(d, 'seq.dat'))
    data = open(os.path.join(d, 'traj.xtc'), 'rb').read()
    offs = xtc.index(data)
    n_frames = len(offs)
    nf = hf.frames_for(len(seq))
    idx = [(k * n_frames) // nf for k in range(nf)]

    def get_frame(k):
        natoms, step, time, box, X, Y, Z = xtc.decode(data, offs[k])
        return [[v * 10.0 for v in X], [v * 10.0 for v in Y], [v * 10.0 for v in Z]], \
            (box[0] * 10.0, box[4] * 10.0, box[8] * 10.0)
    return seq, get_frame, n_frames, idx, (lambda: None)


# ── workers (top level, so the spawn start method can pickle them) ───────────
def one_cal_size(job):
    uid, group, c = job
    try:
        seq, get_frame, n_frames, idx, close = cal_frames(uid, group)
    except Exception as e:
        return ('fail', uid, group, repr(e))
    try:
        C2, maxcut = c2_matrix(cal_sigma_matrix(seq), c)
        return hf.analyse(uid, group, seq, get_frame, n_frames, idx, ref=False, C2=C2, maxcut=maxcut)
    except Exception as e:
        return ('fail', uid, group, repr(e))
    finally:
        close()


def one_mpi_size(job):
    uid, group, c = job
    try:
        seq, get_frame, n_frames, idx, close = mpi_frames(uid)
        S = mpi_sigma_matrix(uid)
        if len(S) != len(seq):
            raise ValueError(f'{len(S)} typed atoms vs sequence {len(seq)}')
        C2, maxcut = c2_matrix(S, c)
        return hf.analyse(uid, group, seq, get_frame, n_frames, idx, ref=False, C2=C2, maxcut=maxcut)
    except Exception as e:
        return ('fail', uid, group, repr(e))


def one_recharged(job):
    uid, group = job
    try:
        seq, get_frame, n_frames, idx, close = cal_frames(uid, group, base=RECHARGED)
    except Exception as e:
        return ('fail', uid, group, repr(e))
    try:
        return hf.analyse(uid, group, seq, get_frame, n_frames, idx, ref=False)
    except Exception as e:
        return ('fail', uid, group, repr(e))
    finally:
        close()


# ── calibration of c ─────────────────────────────────────────────────────────
def calibrate_one(job):
    """Histogram of d / sigma_ij over all pairs |i-j| >= 3 on 100 frames, plus the 0.8 nm contact count."""
    uid, group, ff = job
    try:
        if ff == 'calvados':
            seq, get_frame, n_frames, idx, close = cal_frames(uid, group)
            S = cal_sigma_matrix(seq)
        else:
            seq, get_frame, n_frames, idx, close = mpi_frames(uid)
            S = mpi_sigma_matrix(uid)
    except Exception as e:
        return None
    n = len(seq)
    hist = [0] * 400                   # d / sigma in [0, 4) in steps of 0.01
    n_std = 0
    for k in idx[::10]:
        xyz, box = get_frame(k)
        hf.unwrap(xyz, box)
        X, Y, Z = xyz
        for i in range(n - hf.MIN_SEP):
            for j in range(i + hf.MIN_SEP, n):
                d2 = (X[j] - X[i]) ** 2 + (Y[j] - Y[i]) ** 2 + (Z[j] - Z[i]) ** 2
                if d2 < hf.CUT2:
                    n_std += 1
                r = d2 ** 0.5 / S[i][j]
                if r < 4.0:
                    hist[int(r * 100)] += 1
    close()
    return ff, n_std, hist


def calibrate():
    jobs_ids = [(r['UniProt_ID'], r['group']) for r in csv.DictReader(open(PAIRED_IDS))]
    rng = random.Random(1)
    sample = rng.sample(jobs_ids, 40)
    from multiprocessing import Pool
    jobs = [(u, g, ff) for u, g in sample for ff in ('calvados', 'mpipi')]
    tot = {'calvados': [0, [0] * 400], 'mpipi': [0, [0] * 400]}
    with Pool(9) as pool:
        for res in pool.imap_unordered(calibrate_one, jobs):
            if res:
                ff, n_std, hist = res
                tot[ff][0] += n_std
                tot[ff][1] = [a + b for a, b in zip(tot[ff][1], hist)]
    cs = {}
    for ff, (n_std, hist) in tot.items():
        cum = 0
        for b, h in enumerate(hist):
            if cum + h >= n_std:
                frac = (n_std - cum) / h if h else 0.0
                cs[ff] = round((b + frac) / 100.0, 4)
                break
            cum += h
        print(f'{ff}: {n_std} contacts at 0.8 nm on the calibration sample -> c = {cs[ff]} (cutoff = c * sigma_ij)')
    json.dump(cs, open(C_FILE, 'w'), indent=1)
    return cs


# ── Mpipi-Recharged with the bead-size rule (assumes Mpipi's Wang-Frenkel sigma table; Recharged changed the
#    charge terms) ──────────────────────────────────────────────────────────────
SIGMA_TABLE = os.path.join(HERE, 'mpipi_sigma_table.json')


def build_mpipi_sigma_table(n_runs=400):
    """sigma (Å) for every residue-letter pair, from the wf/cut lines of MPIPI-GG run.in files."""
    tab = {}
    for u in sorted(os.listdir(hm.RUNS))[:n_runs]:
        d = os.path.join(hm.RUNS, u)
        try:
            masses = open(os.path.join(d, 'seq.dat')).read().split('Masses')[1].split('Atoms')[0]
            t2aa = {int(m.group(1)): m.group(2) for m in re.finditer(r'^\s*(\d+)\s+[\d.]+\s+#\s*([A-Z])', masses, re.M)}
            for m in re.finditer(r'^pair_coeff\s+(\d+)\s+(\d+)\s+wf/cut\s+[\d.eE+-]+\s+([\d.eE+-]+)',
                                 open(os.path.join(d, 'run.in')).read(), re.M):
                a, b = t2aa.get(int(m.group(1))), t2aa.get(int(m.group(2)))
                if a and b:
                    tab[a + b] = tab[b + a] = float(m.group(3))
        except Exception:
            continue
        if len(tab) >= 400:
            break
    json.dump(tab, open(SIGMA_TABLE, 'w'))
    return tab


def table_sigma_matrix(seq, tab):
    return [[tab[a + b] for b in seq] for a in seq]


def one_recharged_size(job):
    uid, group, c = job
    try:
        seq, get_frame, n_frames, idx, close = cal_frames(uid, group, base=RECHARGED)
    except Exception as e:
        return ('fail', uid, group, repr(e))
    try:
        tab = json.load(open(SIGMA_TABLE))
        C2, maxcut = c2_matrix(table_sigma_matrix(seq, tab), c)
        return hf.analyse(uid, group, seq, get_frame, n_frames, idx, ref=False, C2=C2, maxcut=maxcut)
    except Exception as e:
        return ('fail', uid, group, repr(e))
    finally:
        close()


def calibrate_recharged_one(job):
    uid, group = job
    try:
        seq, get_frame, n_frames, idx, close = cal_frames(uid, group, base=RECHARGED)
        S = table_sigma_matrix(seq, json.load(open(SIGMA_TABLE)))
    except Exception:
        return None
    n = len(seq)
    hist = [0] * 400
    n_std = 0
    for k in idx[::10]:
        xyz, box = get_frame(k)
        hf.unwrap(xyz, box)
        X, Y, Z = xyz
        for i in range(n - hf.MIN_SEP):
            for j in range(i + hf.MIN_SEP, n):
                d2 = (X[j] - X[i]) ** 2 + (Y[j] - Y[i]) ** 2 + (Z[j] - Z[i]) ** 2
                if d2 < hf.CUT2:
                    n_std += 1
                r = d2 ** 0.5 / S[i][j]
                if r < 4.0:
                    hist[int(r * 100)] += 1
    close()
    return n_std, hist


def size_recharged():
    from multiprocessing import Pool
    tab = build_mpipi_sigma_table()
    print(f'sigma table: {len(tab)} ordered residue pairs', flush=True)
    ids = [(r['UniProt_ID'], r['group']) for r in csv.DictReader(open(os.path.join(HERE, 'recharged_ids.csv')))]
    sample = random.Random(1).sample(ids, 40)
    n_std, hist = 0, [0] * 400
    with Pool(9) as pool:
        for res in pool.imap_unordered(calibrate_recharged_one, sample):
            if res:
                n_std += res[0]
                hist = [a + b for a, b in zip(hist, res[1])]
    cum, c = 0, None
    for b, h in enumerate(hist):
        if cum + h >= n_std:
            c = round((b + ((n_std - cum) / h if h else 0.0)) / 100.0, 4)
            break
        cum += h
    cs = json.load(open(C_FILE))
    cs['recharged'] = c
    json.dump(cs, open(C_FILE, 'w'), indent=1)
    print(f'recharged: c = {c}', flush=True)
    hf.main(worker=one_recharged_size, jobs=[(u, g, c) for u, g in ids],
            out=os.path.join(HERE, 'hub_frames_sizeaware_recharged'))
    hf.main(worker=one_cal_size, jobs=[(u, g, cs['calvados']) for u, g in ids],
            out=os.path.join(HERE, 'hub_frames_sizeaware_recharged_pair_calvados'))


def main():
    cmd = sys.argv[1] if len(sys.argv) > 1 else ''
    if cmd == 'calibrate':
        calibrate()
    elif cmd == 'size':
        cs = json.load(open(C_FILE))
        ids = [(r['UniProt_ID'], r['group']) for r in csv.DictReader(open(PAIRED_IDS))]
        hf.main(worker=one_cal_size, jobs=[(u, g, cs['calvados']) for u, g in ids],
                out=os.path.join(HERE, 'hub_frames_sizeaware_calvados'))
        hf.main(worker=one_mpi_size, jobs=[(u, g, cs['mpipi']) for u, g in ids],
                out=os.path.join(HERE, 'hub_frames_sizeaware_mpipi'))
    elif cmd == 'recharged':
        ok = {r['UniProt_ID']: r['kingdom'] for r in csv.DictReader(open(hf.BIO)) if r['analysis_ok'] == 'True'}
        ids = []
        for u in sorted(os.listdir(RECHARGED)):
            if u in ok and os.path.exists(os.path.join(RECHARGED, u, f'{u}.dcd')):
                ids.append((u, ok[u]))
        with open(os.path.join(HERE, 'recharged_ids.csv'), 'w', newline='') as fh:
            csv.writer(fh).writerows([('UniProt_ID', 'group')] + ids)
        print(f'{len(ids)} Mpipi-Recharged proteins', flush=True)
        hf.main(worker=one_recharged, jobs=ids, out=os.path.join(HERE, 'hub_frames_recharged'))
    elif cmd == 'size_recharged':
        size_recharged()
    else:
        print(__doc__)


if __name__ == '__main__':
    main()
