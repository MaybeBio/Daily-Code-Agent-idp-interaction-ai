"""
Heterotypic slab analysis with ../DDM.py.

    python run_ddm_het.py <index>     DDM of one system, written to
                                      results/<dataset>/<base_name>.npz
    python run_ddm_het.py --count     number of systems over all datasets
                                      (for sub.sh)
    python run_ddm_het.py --collect   build heterotypic_<dataset>.csv from
                                      the npz files

Systems are the rows of each dataset's prep CSV, taken in the order of
DATASETS, so an index means the same system between submissions. The CSVs
hold only the columns the Figures/ scripts read, plus cdil_onephase,
cdil1_onephase, cdil2_onephase and frac_vox_ge_cden for the two-phase gate
(>= 5% of voxels at or above cden), which is applied in the figure scripts.
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import MDAnalysis as mda

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import DDM

DATA = HERE.parent / 'DATASETS'

# tag -> (prep CSV, trajectory directory, minimum raw frames)
DATASETS = {
    'main':   (DATA / 'gin_prep.csv', HERE.parent / 'TRAJECTORIES', 4000),
    'sabari': (DATA / 'sabari_gin_prep.csv',
               HERE.parent.parent / 'mix_affinity' / 'run_sabari' / 'data', 6000),
    'exp':    (DATA / 'exp_gin_prep.csv', HERE.parent / 'TRAJECTORIES', 4000),
    'known':  (DATA / 'known_gin_prep.csv', HERE.parent / 'TRAJECTORIES', 4000),
}

CHUNK_SIZE = 20.0       # voxel side length S, in Angstrom

AA_MAP = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q',
    'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
    'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
    'TYR': 'Y', 'VAL': 'V', 'RNA': 'r',
}

CSV_COLUMNS = ['base_name', 'seq_name1', 'seq_name2', 'comp1', 'comp2',
               'len1', 'len2', 'n1', 'n2', 'phi1', 'phi2',
               'cdil', 'cden', 'cdil1', 'cden1', 'cdil2', 'cden2',
               'dGij', 'dGij_err', 'D_norm', 'cdil_onephase', 'cdil1_onephase',
               'cdil2_onephase', 'frac_vox_ge_cden', 'n_frames']


# --- locating systems ----------------------------------------------------

def base_name(row):
    """As the simulations were named: fractions truncated, not rounded."""
    frac1 = int(row['fraction_seq1'] * 100)
    frac2 = int(row['fraction_seq2'] * 100)
    return f"{row['seq_name1']}_{row['seq_name2']}_{frac1}_{frac2}"


def systems():
    """(dataset, row) for every row of every prep CSV, in a fixed order."""
    out = []
    for tag, (prep, _, _) in DATASETS.items():
        out += [(tag, row) for _, row in pd.read_csv(prep).iterrows()]
    return out


def trajectory_files(traj_dir, name):
    """Topology plus every dcd piece, restarts (<name>_2, <name>_2_3, ...)
    after the original."""
    pdb = traj_dir / f"{name}.pdb"
    dcds = sorted(traj_dir.glob(f"{name}*.dcd"), key=lambda p: len(p.name))
    return pdb, [str(d) for d in dcds]


def species_in_prep_order(u, row):
    """True if the first species in the topology is seq1 of the prep row,
    so DDM.py's species 1 and 2 line up with seq_name1 and seq_name2."""
    ag1, _, n_species = DDM.detect_species(u)
    if n_species != 2:
        return False
    fasta1 = "".join(AA_MAP.get(a.name, 'X') for a in ag1.segments[0].atoms)
    return fasta1 == row['seq1']


# --- running and collecting ----------------------------------------------

def run(index):
    tag, row = systems()[index]
    _, traj_dir, min_frames = DATASETS[tag]
    name = base_name(row)
    out = HERE / 'results' / tag / f"{name}.npz"

    pdb, dcds = trajectory_files(traj_dir, name)
    if not pdb.is_file() or not dcds:
        print(f"{name}: trajectory missing in {traj_dir}, skipping")
        return
    u = mda.Universe(str(pdb), dcds)
    total_frames = len(u.trajectory)

    # redo a finished system only if its trajectory has grown since
    if out.exists():
        with np.load(out) as d:
            done_frames = int(d['total_frames'])
        if total_frames <= done_frames:
            print(f"{name}: already done ({done_frames} frames), skipping")
            return
        print(f"{name}: {total_frames} frames now (was {done_frames}), redoing")

    if total_frames < min_frames:
        print(f"{name}: only {total_frames} frames (< {min_frames}), skipping")
        return
    if not species_in_prep_order(u, row):
        print(f"{name}: first species in the topology is not seq1, skipping")
        return

    DDM.process_sample(pdb, dcds, name=name, chunk_size=CHUNK_SIZE, out=out,
                       names=[row['seq_name1'], row['seq_name2']])


def csv_row(npz_path):
    comp1, comp2 = (int(c) for c in npz_path.stem.rsplit('_', 2)[1:])
    with np.load(npz_path) as d:
        return dict(base_name=npz_path.stem,
                    seq_name1=str(d['seq_name1']), seq_name2=str(d['seq_name2']),
                    comp1=comp1, comp2=comp2,
                    len1=int(d['len1']), len2=int(d['len2']),
                    n1=int(d['n1']), n2=int(d['n2']),
                    phi1=float(d['phi1']), phi2=float(d['phi2']),
                    cdil=float(d['cdil']), cden=float(d['cden']),
                    cdil1=float(d['cdil1']), cden1=float(d['cden1']),
                    cdil2=float(d['cdil2']), cden2=float(d['cden2']),
                    dGij=float(d['dGij']), dGij_err=float(d['dG_err']),
                    D_norm=float(d['D_norm']),
                    cdil_onephase=float(d['cdil_onephase']),
                    cdil1_onephase=float(d['avg_s1'].mean()),
                    cdil2_onephase=float(d['avg_s2'].mean()),
                    frac_vox_ge_cden=float(d['frac_vox_ge_cden']),
                    n_frames=int(d['n_frames']))


def collect():
    for tag in DATASETS:
        files = sorted((HERE / 'results' / tag).glob("*.npz"))
        if not files:
            continue
        df = pd.DataFrame([csv_row(f) for f in files], columns=CSV_COLUMNS)
        csv_path = HERE / f"heterotypic_{tag}.csv"
        df.to_csv(csv_path, index=False)
        print(f"{csv_path.name}: {len(df)} rows")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run DDM.py over the heterotypic slab simulations.")
    parser.add_argument("index", nargs="?", type=int,
                        help="Index of the system to analyse.")
    parser.add_argument("--count", action="store_true",
                        help="Print the number of systems over all datasets.")
    parser.add_argument("--collect", action="store_true",
                        help="Build heterotypic_<dataset>.csv from the npz files.")
    args = parser.parse_args()

    if args.count:
        print(len(systems()))
    elif args.collect:
        collect()
    elif args.index is not None:
        run(args.index)
    else:
        parser.error("give an index, --count or --collect")
