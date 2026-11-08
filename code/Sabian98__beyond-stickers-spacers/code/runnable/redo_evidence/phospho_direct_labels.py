"""
phospho_direct_labels.py — labels for phospho_direct_plan.md (hash 5a69b002) with Amendment 1 (hash 807887ae).
Written and hashed with phospho_direct_test.py before any run was labelled.

1. Degree and hub occupancy per residue for the 240 runs with the paper's labeller (gaivi/hub_labels_gaivi.py, which
   calls hub_frames.py): Calpha-Calpha < 0.8 nm, |i-j| >= 3, exact betweenness, top round(0.1 N) non-zero residues per
   frame, ties by hub_frames.tie_rng (seeded by the run name), 1,000 evenly spaced frames of the second half of each
   trajectory. --no-require-done because only the .dcd and top.pdb files were copied back; completeness was checked
   instead by gaivi/phospho_direct/qc_runs.py (240/240, hashed 8739663c).
2. Calpha radius of gyration (unweighted, as bender_common.geometry and the PAGE4 validation) on the same frames,
   after the same minimum-image unwrap, for the reweighted-against-direct Rg comparison.

  ~/.venvs/finches/bin/python phospho_direct_labels.py
      ->  phospho_direct/labels/{per_residue,per_protein,failed}.csv, settings.json, phospho_direct/rg.csv
"""
import csv
import math
import os
import subprocess
import sys
from multiprocessing import Pool

DYN = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, DYN)
os.environ['HUB_MIN_SEP'] = '3'
os.environ.pop('HUB_FRAMES', None)
import hub_frames as hf  # noqa: E402

ROOT = (__import__('os').environ['BSS_DATA'] + '/BENDER-DYN/phospho_direct/300K')
FASTA = os.path.join(DYN, 'gaivi', 'phospho_direct', 'CALVADOS_IDP_snapshot_amend1', 'input_phospho',
                     'phospho_direct_300K.fasta')                     # the hashed copy of the submitted FASTA
OUT = os.path.join(DYN, 'phospho_direct')
WORKERS = 8


def runs():
    return [line[1:].split('|')[0].strip() for line in open(FASTA) if line.startswith('>')]


def rg_one(run):
    """Mean Calpha Rg (nm) over the label frames: second half, hf.frames_for(N) evenly spaced, as hub_labels_gaivi."""
    fh, off, frame, natom, has_cell, n_frames = hf.dcd_open(os.path.join(ROOT, run, f'{run}.dcd'))
    try:
        start = n_frames // 2
        nf = hf.frames_for(natom)
        idx = [start + (k * (n_frames - start)) // nf for k in range(nf)]
        tot = 0.0
        for k in idx:
            xyz, box = hf.dcd_frame(fh, off, frame, natom, has_cell, k)
            hf.unwrap(xyz, box)
            c = [sum(a) / natom for a in xyz]
            tot += math.sqrt(sum((xyz[d][i] - c[d]) ** 2 for d in range(3) for i in range(natom)) / natom)
        return run, natom, len(idx), tot / len(idx) / 10.0
    finally:
        fh.close()


def main():
    os.makedirs(OUT, exist_ok=True)
    names = runs()
    assert len(names) == 240, len(names)
    subprocess.run([sys.executable, os.path.join(DYN, 'gaivi', 'hub_labels_gaivi.py'), '--layout', 'openmm',
                    '--root', ROOT, '--fasta', FASTA, '--out', os.path.join(OUT, 'labels'),
                    '--no-require-done', '--workers', str(WORKERS)], check=True)
    with Pool(WORKERS) as pool:
        rows = sorted(pool.map(rg_one, names))
    with open(os.path.join(OUT, 'rg.csv'), 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['run', 'N', 'frames_used', 'rg_mean_nm'])
        w.writerows([(r, n, f, f'{g:.5f}') for r, n, f, g in rows])
    print(f'rg.csv: {len(rows)} runs', flush=True)


if __name__ == '__main__':
    main()
