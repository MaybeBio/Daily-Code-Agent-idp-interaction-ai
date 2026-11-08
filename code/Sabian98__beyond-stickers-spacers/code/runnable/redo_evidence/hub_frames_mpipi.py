"""
hub_frames_mpipi.py — the per-frame hub analysis of hub_frames.py on the MPIPI-GG 300 K trajectories
(BENDER-BIO/MPIPI/mpipi_sim/300k_runs/Mpipi 300k /<ID>/traj.xtc; LAMMPS, all saved frames are production,
coordinates wrapped in the 50 nm box and rebuilt from minimum-image bonds as in bender_common.load_mpipi).

Sample: the sequences simulated with both force fields (mpipi_results_v2.csv: status ok, in_paper_set), analysis_ok
in BENDER_BIO, at most hub_frames.N_MAX residues. FRAMES evenly spaced frames over the whole trajectory.
The same IDs are written to hub_frames_mpipi/ids.csv so CALVADOS-2 can be run on exactly the same set:
  python3 hub_frames_mpipi.py                                          # MPIPI-GG -> hub_frames_mpipi/
  python3 hub_frames.py --ids hub_frames_mpipi/ids.csv --out hub_frames_paired_calvados
XTC frames are read with xtc_stdlib.py (validated against LAMMPS Rg.out: |ΔRg| <= 0.001 Å).
"""
import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import hub_frames as hf  # noqa: E402
import xtc_stdlib as xtc  # noqa: E402

RUNS = (__import__('os').environ['BSS_DATA'] + '/BENDER-BIO/MPIPI/mpipi_sim/300k_runs/Mpipi 300k ')   # trailing space is in the name
RESULTS = (__import__('os').environ['BSS_DATA'] + '/BENDER-BIO/MPIPI/MPIPI Results/mpipi_results_v2.csv')
OUT = os.path.join(HERE, 'hub_frames_mpipi')


def read_sequence(seq_dat):
    """As bender_common.read_mpipi_sequence: the longest alphabetic token of more than 5 letters."""
    toks = open(seq_dat).read().split()
    seqs = [t for t in toks if t.isalpha() and len(t) > 5]
    return max(seqs, key=len).upper()


def one(job):
    uid, group = job
    d = os.path.join(RUNS, uid)
    try:
        seq = read_sequence(os.path.join(d, 'seq.dat'))
        data = open(os.path.join(d, 'traj.xtc'), 'rb').read()
        offs = xtc.index(data)
    except Exception as e:
        return ('fail', uid, group, repr(e))
    n_frames = len(offs)
    nf = hf.frames_for(len(seq))
    if n_frames < nf:
        return ('fail', uid, group, f'only {n_frames} frames')
    idx = [(k * n_frames) // nf for k in range(nf)]

    def get_frame(k):
        natoms, step, time, box, X, Y, Z = xtc.decode(data, offs[k])
        if natoms != len(seq):
            raise ValueError(f'{natoms} atoms vs sequence {len(seq)}')
        return [[v * 10.0 for v in X], [v * 10.0 for v in Y], [v * 10.0 for v in Z]], \
            (box[0] * 10.0, box[4] * 10.0, box[8] * 10.0)

    try:
        return hf.analyse(uid, group, seq, get_frame, n_frames, idx, ref=False)
    except Exception as e:
        return ('fail', uid, group, repr(e))


def sample():
    ok = {r['UniProt_ID']: r['kingdom'] for r in csv.DictReader(open(hf.BIO)) if r['analysis_ok'] == 'True'}
    jobs = []
    for r in csv.DictReader(open(RESULTS)):
        u = r['UniProt_ID']
        if (r['status'] == 'ok' and r['in_paper_set'] == 'True' and u in ok
                and float(r['n_residues']) <= hf.N_MAX and os.path.exists(os.path.join(RUNS, u, 'traj.xtc'))):
            jobs.append((u, ok[u]))
    return sorted(jobs)


if __name__ == '__main__':
    jobs = sample()
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, 'ids.csv'), 'w', newline='') as fh:
        csv.writer(fh).writerows([('UniProt_ID', 'group')] + jobs)
    print(f'{len(jobs)} sequences with both force fields', flush=True)
    hf.main(worker=one, jobs=jobs, out=OUT)
