"""
Homotypic slab analysis with ../DDM.py.

    python run_ddm.py <index>      DDM of one sequence at every voxel size,
                                   written to results_<S>/<name>.npz
    python run_ddm.py --count      number of sequence directories (for sub.sh)
    python run_ddm.py --collect    build slab_results_<S>.csv from the npz
                                   files, plus frac_two_peak and
                                   frac_two_peak_occ from peaks_<S>/
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
import MDAnalysis as mda

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import DDM

SIM_ROOT = HERE.parent.parent / 'calv_comp'
SIM_DIRS = ['simulations_training',
            'simulations_validation_IDRome10',
            'simulations_validation_IDRome90']

S_LIST = [10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60]  # voxel side lengths in A
SHORT_RUN_START = 100   # start frame for runs shorter than the 1 us cutoff

CSV_COLUMNS = ['seq_name', 'cdil', 'cdil_sem', 'cden', 'cden_sem',
               'dG', 'dG_err', 'cdil_onephase', 'frac_vox_ge_cden', 'n_voxels']


# --- locating simulations ------------------------------------------------

def sequence_directories(sim_root=SIM_ROOT, sim_dirs=SIM_DIRS):
    dirs = []
    for d in sim_dirs:
        base = Path(sim_root) / d
        if not base.is_dir():
            continue
        for p in sorted(base.iterdir()):
            if not p.is_dir():
                continue
            for cand in (p / "293", p):
                if (cand / "top.pdb").is_file() and list(cand.glob("*.dcd")):
                    dirs.append(cand)
                    break
    return dirs


def sequence_name(direc):
    return direc.parent.name if direc.name == "293" else direc.name


def start_frame_from_config(direc, default=DDM.DEFAULT_START):
    """Equilibration cutoff from config.yaml, i.e. the frame at 1 us."""
    config_path = direc / "config.yaml"
    if not config_path.is_file():
        return default
    try:
        with open(config_path) as f:
            config = yaml.safe_load(f)
        wfreq = config['runtime_settings']['wfreq']
        return int(1e9 / (10 * wfreq))
    except (KeyError, TypeError, ValueError, ZeroDivisionError):
        return default


def sample(index):
    """(name, pdb, dcd, start) for one sequence directory, or None if the
    trajectory is too short to analyse."""
    direc = sequence_directories()[index]
    name = sequence_name(direc)
    pdb = direc / "top.pdb"
    dcd = direc / f"{name}.dcd"
    if not dcd.is_file():
        dcd = sorted(direc.glob("*.dcd"))[0]

    n_frames = len(mda.Universe(str(pdb), str(dcd)).trajectory)
    start = start_frame_from_config(direc)
    if n_frames <= start:
        start = SHORT_RUN_START
    if n_frames <= start:
        print(f"{name}: only {n_frames} frames, skipping")
        return None
    return name, pdb, dcd, start


# --- running and collecting ----------------------------------------------

def run(index, sizes=S_LIST):
    s = sample(index)
    if s is None:
        return
    name, pdb, dcd, start = s
    for S in sizes:
        out = HERE / f"results_{S}" / f"{name}.npz"
        if out.exists():
            print(f"{name}, S={S}: already done, skipping")
            continue
        DDM.process_sample(pdb, dcd, name=name, start=start,
                           chunk_size=S, out=out)


def csv_row(npz_path):
    """cdil_onephase and frac_vox_ge_cden are carried so the figure scripts
    can apply the two-phase gate (>= 5% of voxels at or above cden)."""
    with np.load(npz_path) as d:
        return dict(seq_name=str(d['base_name']),
                    cdil=float(d['cdil']), cdil_sem=float(d['cdil_sem']),
                    cden=float(d['cden']), cden_sem=float(d['cden_sem']),
                    dG=float(d['dG']), dG_err=float(d['dG_err']),
                    cdil_onephase=float(d['cdil_onephase']),
                    frac_vox_ge_cden=float(d['frac_vox_ge_cden']),
                    n_voxels=int(d['avg_tot'].size))


def peak_row(npz_path):
    with np.load(npz_path) as d:
        return dict(seq_name=npz_path.stem,
                    frac_two_peak=float(d['frac_two_peak']),
                    frac_two_peak_occ=float(d['frac_two_peak_occ']))


def collect(sizes=S_LIST):
    for S in sizes:
        files = sorted((HERE / f"results_{S}").glob("*.npz"))
        if not files:
            continue
        df = pd.DataFrame([csv_row(f) for f in files], columns=CSV_COLUMNS)

        peak_files = sorted((HERE / f"peaks_{S}").glob("*.npz"))
        if peak_files:
            peaks = pd.DataFrame([peak_row(f) for f in peak_files])
            df = df.merge(peaks, on='seq_name', how='left')

        csv_path = HERE / f"slab_results_{S}.csv"
        df.to_csv(csv_path, index=False)
        print(f"{csv_path.name}: {len(df)} rows")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run DDM.py over the homotypic slab simulations.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("index", nargs="?", type=int,
                        help="Index of the sequence directory to analyse.")
    parser.add_argument("--count", action="store_true",
                        help="Print the number of sequence directories.")
    parser.add_argument("--collect", action="store_true",
                        help="Build slab_results_<S>.csv from the npz files.")
    parser.add_argument("--sizes", type=int, nargs="+", default=S_LIST,
                        help="Voxel side lengths S, in Angstrom.")
    args = parser.parse_args()

    if args.count:
        print(len(sequence_directories()))
    elif args.collect:
        collect(args.sizes)
    elif args.index is not None:
        run(args.index, args.sizes)
    else:
        parser.error("give an index, --count or --collect")
