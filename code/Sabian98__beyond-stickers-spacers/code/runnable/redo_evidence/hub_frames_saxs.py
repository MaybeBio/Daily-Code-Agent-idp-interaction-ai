"""
hub_frames_saxs.py — the per-frame hub analysis (hub_frames.py, 1,000 frames, v2 contact rule) on the 77 SAXS-benchmark
constructs simulated with both force fields:
  CALVADOS-2  $BSS_DATA/When PLMs fail- ICLR 2027/Experimental Validation Traj/<id>/<id>.dcd + top.pdb
              (production half)
  MPIPI-GG    $BSS_DATA/BENDER-BIO/MPIPI/mpipi_sim/300k_runs/Mpipi SAXS/<id>/traj.xtc (all frames)
Constructs and labels from ZOIDBERG .../experimental_data/bender_benchmark/traj_vs_saxs_rg.csv (traj_id, protein_name,
taxon, calvados2_overlap). No length cap.
  python3 hub_frames_saxs.py   ->  hub_frames_saxs_calvados/, hub_frames_saxs_mpipi/, saxs_ids.csv
"""
import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import hub_frames as hf  # noqa: E402
import hub_frames_mpipi as hm  # noqa: E402
import hub_frames_tests as ht  # noqa: E402
import xtc_stdlib as xtc  # noqa: E402

CAL_SAXS = (__import__('os').environ['BSS_DATA'] + '/When PLMs fail- ICLR 2027/Experimental Validation Traj')
MPI_SAXS = (__import__('os').environ['BSS_DATA'] + '/BENDER-BIO/MPIPI/mpipi_sim/300k_runs/Mpipi SAXS')
TABLE = (__import__('os').environ['BSS_DATA'] + '/ZOIDBERG-ICML 2027/experimental_data/bender_benchmark/traj_vs_saxs_rg.csv')


def one_cal(job):
    uid, group = job
    d = os.path.join(CAL_SAXS, uid)
    try:
        seq = ''.join(hf.AA3.get(line[17:20], 'X') for line in open(os.path.join(d, 'top.pdb'))
                      if line.startswith('ATOM') and line[12:16].strip() == 'CA')
        fh, off, frame, natom, has_cell, n_frames = hf.dcd_open(os.path.join(d, f'{uid}.dcd'))
    except Exception as e:
        return ('fail', uid, group, repr(e))
    try:
        if len(seq) != natom:
            return ('fail', uid, group, f'top {len(seq)} vs dcd {natom} atoms')
        start = n_frames // 2
        nf = hf.frames_for(natom)
        idx = [start + (k * (n_frames - start)) // nf for k in range(nf)]
        return hf.analyse(uid, group, seq, lambda k: hf.dcd_frame(fh, off, frame, natom, has_cell, k), n_frames, idx,
                          ref=False)
    except Exception as e:
        return ('fail', uid, group, repr(e))
    finally:
        fh.close()


def one_mpi(job):
    uid, group = job
    d = os.path.join(MPI_SAXS, uid)
    try:
        seq = hm.read_sequence(os.path.join(d, 'seq.dat'))
        data = open(os.path.join(d, 'traj.xtc'), 'rb').read()
        offs = xtc.index(data)
        n_frames = len(offs)
        nf = hf.frames_for(len(seq))
        idx = [(k * n_frames) // nf for k in range(nf)]

        def get_frame(k):
            natoms, step, time, box, X, Y, Z = xtc.decode(data, offs[k])
            if natoms != len(seq):
                raise ValueError(f'{natoms} atoms vs sequence {len(seq)}')
            return [[v * 10.0 for v in X], [v * 10.0 for v in Y], [v * 10.0 for v in Z]], \
                (box[0] * 10.0, box[4] * 10.0, box[8] * 10.0)
        return hf.analyse(uid, group, seq, get_frame, n_frames, idx, ref=False)
    except Exception as e:
        return ('fail', uid, group, repr(e))


if __name__ == '__main__':
    rows = list(csv.DictReader(open(TABLE)))
    jobs = [(r['traj_id'], r['taxon'] or 'NA') for r in rows
            if os.path.exists(os.path.join(CAL_SAXS, r['traj_id'])) and os.path.exists(os.path.join(MPI_SAXS, r['traj_id']))]
    with open(os.path.join(HERE, 'saxs_ids.csv'), 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['UniProt_ID', 'group', 'protein_name', 'calvados2_overlap'])
        names = {r['traj_id']: (r['protein_name'], r['calvados2_overlap']) for r in rows}
        w.writerows([(u, g) + names[u] for u, g in jobs])
    print(f'{len(jobs)} SAXS constructs with both force fields', flush=True)
    hf.main(worker=one_cal, jobs=jobs, out=os.path.join(HERE, 'hub_frames_saxs_calvados'))
    hf.main(worker=one_mpi, jobs=jobs, out=os.path.join(HERE, 'hub_frames_saxs_mpipi'))
