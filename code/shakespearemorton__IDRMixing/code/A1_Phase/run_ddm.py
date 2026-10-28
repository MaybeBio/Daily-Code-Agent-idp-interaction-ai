"""
hnRNPA1 phase diagram (Figure 1G) with ../DDM.py.

    python run_ddm.py              DDM of every run in jobs.csv, written to
                                   results/<run>.npz
    python run_ddm.py <index>      DDM of one run (row <index> of jobs.csv)
    python run_ddm.py --collect    build A1_phase.csv from the npz files

Each run is the directory A1_<sim_hash> beside this script: top.pdb plus its
trajectory pieces (backup_*.dcd in time order, then the live dcd). The
first 1 us is discarded, as in DDM_Homotypic/run_ddm.py. A1_phase.csv holds
only what Figure1.py reads; a run is two-phase when the DDM resolves a
dense phase, as for the other homotypic systems in Figure 1.
"""

import argparse
import re
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import DDM

CHUNK_SIZE = 20.0       # voxel side length S, in Angstrom
OUT_CSV = HERE / 'A1_phase.csv'
# backup_<name>_<YYYYDDMM>_<HHh>MM<m>SS<s>.dcd
BACKUP_RE = re.compile(r'^backup_.*_(\d{8})_(\d{2})h(\d{2})m(\d{2})s\.dcd$')

CSV_COLUMNS = ['condition', 'temp_K', 'box_x_nm', 'box_z_nm', 'two_phase',
               'c_dil_mM', 'c_dil_err', 'c_den_mM', 'c_den_err', 'c_mean_mM', 'dG_kT']


def jobs():
    return pd.read_csv(HERE / 'jobs.csv')


def run_dir(job):
    return HERE / f"A1_{job.sim_hash}"


def dcd_sort_key(path):
    """Backups first, oldest to newest by the timestamp in their name, then
    the live trajectory, which is still being appended to."""
    m = BACKUP_RE.match(path.name)
    if not m:
        return (1, path.name)
    date_str, hh, mm, ss = m.groups()
    return (0, datetime.strptime(f"{date_str}{hh}{mm}{ss}", "%Y%d%m%H%M%S"))


def start_frame(direc, default=DDM.DEFAULT_START):
    """Frame at 1 us, from wfreq in config.yaml."""
    try:
        with open(direc / "config.yaml") as f:
            config = yaml.safe_load(f)
        return int(1e9 / (10 * float(config['wfreq'])))
    except (OSError, KeyError, TypeError, ValueError, ZeroDivisionError):
        return default


def run(job):
    direc = run_dir(job)
    dcds = sorted(direc.glob("*.dcd"), key=dcd_sort_key)
    if not (direc / "top.pdb").is_file() or not dcds:
        print(f"{direc.name}: top.pdb or trajectory missing, skipping")
        return
    out = HERE / 'results' / f"{job.notes}.npz"
    DDM.process_sample(direc / "top.pdb", [str(d) for d in dcds], name=job.notes,
                       start=start_frame(direc), chunk_size=CHUNK_SIZE, out=out)


def csv_row(job):
    with np.load(HERE / 'results' / f"{job.notes}.npz") as d:
        box = d['chunk_dims'] * d['n_chunks'] / 10.           # nm
        return dict(condition=f"{box[0]:.0f}x{box[2]:.0f} nm", temp_K=float(job.temp),
                    box_x_nm=float(box[0]), box_z_nm=float(box[2]),
                    two_phase=bool(np.isfinite(d['cden'])),
                    c_dil_mM=float(d['cdil']), c_dil_err=float(d['cdil_sem']),
                    c_den_mM=float(d['cden']), c_den_err=float(d['cden_sem']),
                    c_mean_mM=float(d['cdil_onephase']), dG_kT=float(d['dG']))


def collect():
    rows = [csv_row(job) for job in jobs().itertuples()
            if (HERE / 'results' / f"{job.notes}.npz").exists()]
    df = pd.DataFrame(rows, columns=CSV_COLUMNS).sort_values(['condition', 'temp_K'])
    df.to_csv(OUT_CSV, index=False)
    print(f"{OUT_CSV.name}: {len(df)} rows")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="hnRNPA1 phase diagram with DDM.py.")
    parser.add_argument("index", nargs="?", type=int, help="Row of jobs.csv to analyse.")
    parser.add_argument("--collect", action="store_true", help="Build A1_phase.csv.")
    args = parser.parse_args()

    if args.collect:
        collect()
    elif args.index is not None:
        run(jobs().iloc[args.index])
    else:
        for job in jobs().itertuples():
            run(job)
