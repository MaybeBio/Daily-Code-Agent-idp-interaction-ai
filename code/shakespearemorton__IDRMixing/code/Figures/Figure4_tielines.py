#!/usr/bin/env python3
"""FUS-A1 tie lines for Figure 4A, from DDM.py.

For each FUS-A1 / FUS-A1p12D composition the phases come from DDM.py on the
total voxel map (DDM_Heterotypic/results/known): the dilute and dense
voxels are DDM.phase_masks of DDM.phase_concentrations.  Each species'
concentration in a phase is its time-averaged mean over those voxels, and
its error is the block-averaged SEM of the per-frame mean over the same
voxels -- the recipe DDM.py uses for cdil_sem and cden_sem.  The per-frame
series need one pass over each trajectory (original plus continuations,
over the same frames DDM.py analysed), so the result is cached:

    python Figure4_tielines.py      -> figure4_data/fus_a1_tielines.csv
"""

import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd
import MDAnalysis as mda

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import DDM

RESULTS_DIR = HERE.parent / "DDM_Heterotypic" / "results" / "known"
TRAJ_DIR = HERE.parent / "TRAJECTORIES"
OUT_CSV = HERE / "figure4_data" / "fus_a1_tielines.csv"
PAIRS = ('FUS_A1', 'FUS_A1p12D')
N_WORKERS = 6           # one per system; shared storage, never more than 8


def species_series(name, m_dil, m_den, start, n_frames, chunk_size):
    """Per-frame mean residue count of each species over the dilute and the
    dense voxels, binned exactly as DDM.analyze_3d_chunks bins them."""
    pdb = TRAJ_DIR / f"{name}.pdb"
    dcds = [str(d) for d in sorted(TRAJ_DIR.glob(f"{name}*.dcd"), key=lambda p: len(p.name))]
    u = mda.Universe(str(pdb), dcds)
    ag1, _, _ = DDM.detect_species(u)
    ag = u.atoms

    box = u.dimensions[:3]
    n_chunks = np.maximum(np.round(box / chunk_size).astype(int), 1)
    chunk_dims = box / n_chunks
    ny, nz = n_chunks[1], n_chunks[2]
    n_voxels = int(np.prod(n_chunks))
    lz = u.dimensions[2]
    edges = np.arange(0, lz + 1, 1)
    z = edges[:-1] + 0.5

    series = {key: np.empty(n_frames) for key in ('dil1', 'dil2', 'den1', 'den2')}
    avg_tot = np.zeros(n_voxels)
    for i in range(n_frames):
        DDM.centre_slab(ag, u.trajectory[i + start], z, edges, lz)
        idx = np.clip((ag.positions // chunk_dims).astype(int), 0, n_chunks - 1)
        flat_idx = idx[:, 0] * ny * nz + idx[:, 1] * nz + idx[:, 2]
        counts = np.bincount(flat_idx, minlength=n_voxels)
        counts1 = np.bincount(flat_idx[ag1.indices], minlength=n_voxels)
        counts2 = counts - counts1
        avg_tot += (counts - avg_tot) / (i + 1)
        series['dil1'][i], series['dil2'][i] = counts1[m_dil].mean(), counts2[m_dil].mean()
        series['den1'][i], series['den2'][i] = counts1[m_den].mean(), counts2[m_den].mean()
    return series, avg_tot


def tie_line(path):
    name = path.stem
    with np.load(path) as d:
        s = d['avg_tot'].astype(np.float64)
        start, n_frames = int(d['start_frame']), int(d['n_frames'])
        chunk_size = float(d['chunk_size'])
    dil, den = DDM.phase_concentrations(s, n_frames)
    m_dil, m_den = DDM.phase_masks(s, dil, den)
    series, avg_tot = species_series(name, m_dil, m_den, start, n_frames, chunk_size)
    if not np.allclose(avg_tot, s, atol=1e-4):
        raise RuntimeError(f"{name}: recomputed map does not match the DDM.py map")

    c1, c2 = (int(t) for t in name.rsplit('_', 2)[1:])
    row = dict(base_name=name, pair=name.rsplit('_', 2)[0], comp1=c1, comp2=c2,
               n_dil_voxels=int(m_dil.sum()), n_den_voxels=int(m_den.sum()))
    for key, v in series.items():
        row[key] = float(v.mean())
        row[f'{key}_sem'] = DDM.block_sem(v)[0]
    print(f"{name}: dilute FUS {row['dil1']:.4f} +/- {row['dil1_sem']:.4f}, "
          f"A1 {row['dil2']:.4f} +/- {row['dil2_sem']:.4f}")
    return row


if __name__ == '__main__':
    paths = [p for pair in PAIRS for p in sorted(RESULTS_DIR.glob(f"{pair}_*_*.npz"))
             if p.stem.rsplit('_', 2)[0] == pair]
    with Pool(N_WORKERS) as pool:
        rows = pool.map(tie_line, paths)
    OUT_CSV.parent.mkdir(exist_ok=True)
    pd.DataFrame(rows).sort_values(['pair', 'comp1']).to_csv(OUT_CSV, index=False)
    print(f"wrote {OUT_CSV}")
