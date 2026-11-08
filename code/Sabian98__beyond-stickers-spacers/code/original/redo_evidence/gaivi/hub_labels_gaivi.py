"""
hub_labels_gaivi.py — per-frame hub labels (the hub_frames.py pipeline) for simulations that live on GAIVI, so the
trajectories never have to be copied off the cluster.

Layouts (--layout):
  openmm  <root>/<ID>/<ID>.dcd + top.pdb    CALVADOS-2 (CALVADOS IDP/prepare.py) and Mpipi-Recharged
                                            (OpenMPIPI/run_openmpipi_idp.py); production = second half of the frames.
                                            A run counts as finished once <root>/data/<ID>/conf_prop.csv exists.
  lammps  <root>/<ID>/traj.xtc + seq.dat    MPIPI-GG (mpipi_sim/scripts/run_mpipi_from_fasta.py); every saved frame
                                            is production. A run counts as finished once it has the expected number
                                            of frames (--mpipi-time-exponent, 2.0 as in batch_submit.sh).
IDs and groups come from a FASTA with headers >ID|group|... (the group field is optional).
Analysis is identical to the lab-Mac labels in hub_frames_all/: Calpha-Calpha < 0.8 nm, |i-j| >= 3, exact
betweenness, 1,000 / 500 / 250 evenly spaced production frames for N <= 300 / <= 600 / > 600.

Incremental: IDs already in OUT/per_protein.csv are skipped and new results are appended, so the script can be rerun
as simulation chunks finish. Unfinished or failed runs are listed in OUT/failed.csv (rewritten on every run) and are
tried again next time. If the file system stalls (no result and no worker CPU for HUB_STALL_SECONDS, default 900),
the script exits with code 75 and can simply be run again.

Needs hub_frames.py, hub_z_v2_stdlib.py and xtc_stdlib.py next to this file (or in the parent folder). Standard
library only.
  python3 hub_labels_gaivi.py --layout openmm --root /general/taseefr/300K --fasta bender_10282.fasta --out labels_recharged
"""
import argparse
import csv
import os
import sys
import time
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [HERE, os.path.dirname(HERE)]
os.environ['HUB_MIN_SEP'] = '3'                    # contacts need |i - j| >= 3 (fixed here, never from the environment)
os.environ.pop('HUB_FRAMES', None)                 # frames by length: 1,000 / 500 / 250 for N <= 300 / 600 / longer
import hashlib  # noqa: E402
import json  # noqa: E402

import hub_frames as hf  # noqa: E402
import xtc_stdlib as xtc  # noqa: E402

CFG = {}


def read_fasta(path):
    """[(ID, group, sequence)] for headers >ID|group|... (prepare.py / batch_submit.sh naming)."""
    out, head, seq = [], None, []
    for line in open(path):
        line = line.strip()
        if line.startswith('>'):
            if head is not None:
                out.append((*head, ''.join(seq)))
            f = [x.strip() for x in line[1:].split()[0].split('|')]
            uid = f[1] if f[0] in ('sp', 'tr') and len(f) > 1 else f[0]
            head, seq = (uid, f[1] if f[0] not in ('sp', 'tr') and len(f) > 1 else 'NA'), []
        elif line:
            seq.append(line)
    if head is not None:
        out.append((*head, ''.join(seq)))
    return out


def mpipi_production_frames(n, p):
    """Saved (= production) frames of a finished MPIPI-GG run, as in run_mpipi_from_fasta.py."""
    t = 200 if n <= 150 else int(200 * (n / 150) ** p)
    total = 10000 if t <= 300 else 15000 if t <= 600 else 20000
    return total // 2


def read_seq_dat(path):
    """As hub_frames_mpipi.read_sequence: the longest alphabetic token of more than 5 letters."""
    toks = [t for t in open(path).read().split() if t.isalpha() and len(t) > 5]
    return max(toks, key=len) if toks else ''


def one(job):
    uid, group, fasta_seq = job
    root, layout = CFG['root'], CFG['layout']
    d = os.path.join(root, uid)
    if layout == 'openmm':
        if CFG['require_done'] and not os.path.exists(os.path.join(root, 'data', uid, 'conf_prop.csv')):
            return ('fail', uid, group, 'not finished (no data/<ID>/conf_prop.csv)')
        try:
            seq = ''.join(hf.AA3.get(line[17:20], 'X') for line in open(os.path.join(d, 'top.pdb'))
                          if line.startswith('ATOM') and line[12:16].strip() == 'CA')
            fh, off, frame, natom, has_cell, n_frames = hf.dcd_open(os.path.join(d, f'{uid}.dcd'))
        except Exception as e:
            return ('fail', uid, group, repr(e))
        try:
            if len(seq) != natom:
                return ('fail', uid, group, f'top {len(seq)} vs dcd {natom} atoms')
            if n_frames < hf.expected_frames(natom):
                return ('fail', uid, group, f'truncated {n_frames} < {hf.expected_frames(natom)}')
            start = n_frames // 2
            nf = n_frames - start if CFG['frames'] == 'all' else hf.frames_for(natom)
            idx = [start + (k * (n_frames - start)) // nf for k in range(nf)]
            res = hf.analyse(uid, group, seq, lambda k: hf.dcd_frame(fh, off, frame, natom, has_cell, k),
                             n_frames, idx, ref=False)
        except Exception as e:
            return ('fail', uid, group, repr(e))
        finally:
            try:
                fh.close()
            except OSError:            # a network share can fail on close after every frame was read (2026-10-04)
                pass
    else:
        try:
            seq = read_seq_dat(os.path.join(d, 'seq.dat'))
            data = open(os.path.join(d, 'traj.xtc'), 'rb').read()
            offs = xtc.index(data)
        except Exception as e:
            return ('fail', uid, group, repr(e))
        n_frames, need = len(offs), mpipi_production_frames(len(seq), CFG['p'])
        if n_frames < need:
            return ('fail', uid, group, f'not finished ({n_frames} < {need} frames)')
        nf = n_frames if CFG['frames'] == 'all' else hf.frames_for(len(seq))
        idx = [(k * n_frames) // nf for k in range(nf)]

        def get_frame(k):
            natoms, step, t, box, X, Y, Z = xtc.decode(data, offs[k])
            if natoms != len(seq):
                raise ValueError(f'{natoms} atoms vs sequence {len(seq)}')
            return [[v * 10.0 for v in X], [v * 10.0 for v in Y], [v * 10.0 for v in Z]], \
                (box[0] * 10.0, box[4] * 10.0, box[8] * 10.0)
        try:
            res = hf.analyse(uid, group, seq, get_frame, n_frames, idx, ref=False)
        except Exception as e:
            return ('fail', uid, group, repr(e))
    if fasta_seq and len(seq) != len(fasta_seq):
        # the simulated chain is the reference (setup scripts map non-standard letters), but its length must match
        return ('fail', uid, group, f'sequence length {len(seq)} vs FASTA {len(fasta_seq)}')
    return res


def init(cfg):
    CFG.update(cfg)


def repair(pp, pr):
    """After an interrupted pass: drop a cut-off last line and residue rows of proteins without a protein row."""
    if not (os.path.exists(pp) and os.path.exists(pr)):
        return
    for path in (pp, pr):
        data = open(path, 'rb').read()
        if data and not data.endswith(b'\n'):
            open(path, 'wb').write(data[:data.rfind(b'\n') + 1])
    done = {r['UniProt_ID'] for r in csv.DictReader(open(pp))}
    rows = list(csv.reader(open(pr)))
    keep = [rows[0]] + [r for r in rows[1:] if r and r[0] in done]
    if len(keep) != len(rows):
        with open(pr, 'w', newline='') as fh:
            csv.writer(fh).writerows(keep)
        print(f'repaired {pr}: dropped {len(rows) - len(keep)} rows of unfinished proteins', flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--layout', choices=['openmm', 'lammps'], required=True)
    ap.add_argument('--root', required=True, help='folder holding one sub-folder per run, e.g. /general/taseefr/300K')
    ap.add_argument('--fasta', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--workers', type=int, default=int(os.environ.get('SLURM_CPUS_PER_TASK', os.cpu_count() or 4)))
    ap.add_argument('--no-require-done', action='store_true', help='openmm: do not require data/<ID>/conf_prop.csv')
    ap.add_argument('--mpipi-time-exponent', type=float, default=2.0)
    ap.add_argument('--limit', type=int, default=0, help='only the first N pending IDs (testing)')
    ap.add_argument('--frames', choices=['default', 'all'], default='default',
                    help='default: 1,000 / 500 / 250 snapshots as for the training labels; all: every production frame '
                         '(less within-run noise, same expected occupancy)')
    a = ap.parse_args()
    cfg = {'root': a.root, 'layout': a.layout, 'require_done': not a.no_require_done, 'p': a.mpipi_time_exponent,
           'frames': a.frames}
    os.makedirs(a.out, exist_ok=True)
    pp, pr = os.path.join(a.out, 'per_protein.csv'), os.path.join(a.out, 'per_residue.csv')
    repair(pp, pr)
    code = hashlib.sha256(open(hf.__file__, 'rb').read()).hexdigest()[:16]
    json.dump({'contact_cutoff_A': 8.0, 'min_sep': hf.MIN_SEP, 'top_fraction': hf.TOP, 'frames_rule': '1000/500/250 by length'
               if a.frames == 'default' else 'all production frames', 'ties': 'random, seeded by protein ID',
               'layout': a.layout, 'root': a.root, 'hub_frames_sha256': code},
              open(os.path.join(a.out, 'settings.json'), 'w'), indent=1)
    done = {r['UniProt_ID'] for r in csv.DictReader(open(pp))} if os.path.exists(pp) else set()
    jobs = [(u, g, s) for u, g, s in read_fasta(a.fasta) if u not in done]
    jobs.sort(key=lambda j: -len(j[2]))            # longest first, so slow chains do not finish alone at the end
    if a.limit:
        jobs = jobs[:a.limit]
    print(f'{len(done)} already labelled; {len(jobs)} to try with {a.workers} workers', flush=True)
    new_p, new_r = not os.path.exists(pp), not os.path.exists(pr)
    fails, ok, t0 = [], 0, time.time()
    with open(pp, 'a', newline='') as fp, open(pr, 'a', newline='') as fr:
        wp, wr = csv.writer(fp), csv.writer(fr)
        if new_p:
            wp.writerow(hf.PROT_HEADER)
        if new_r:
            wr.writerow(hf.RES_HEADER)
        with Pool(a.workers, initializer=init, initargs=(cfg,)) as pool:
            for k, res in enumerate(hf.watched(pool, one, jobs), 1):   # exits with code 75 if the share stalls
                if res[0] == 'ok':
                    prot, rows = res[3]
                    wr.writerows(rows)
                    wp.writerow(prot)              # protein row last: marked done only once its residues are written
                    fr.flush()
                    fp.flush()
                    ok += 1
                else:
                    fails.append(res[1:])
                if k % 100 == 0:
                    print(f'{k}/{len(jobs)}  ok {ok}  {time.time() - t0:.0f}s', flush=True)
    with open(os.path.join(a.out, 'failed.csv'), 'w', newline='') as fh:
        csv.writer(fh).writerows([('UniProt_ID', 'group', 'reason')] + fails)
    reasons = {}
    for f in fails:
        key = f[2].split('(')[0].strip()[:40]
        reasons[key] = reasons.get(key, 0) + 1
    print(f'done: {ok} new, {len(fails)} not done/failed {reasons}, {time.time() - t0:.0f}s; '
          f'{len(done) + ok} labelled in {a.out}', flush=True)


if __name__ == '__main__':
    main()
