"""
Finite-size check for Figure 1F: DDM.py on the same sequences in smaller
and larger boxes, at the voxel size used for the main dataset.

    python finite.py               analyse every run beside this script
    python finite.py --force       re-analyse even when a result is cached

One npz per run is written to results_<S>/<name>.npz by DDM.process_sample,
and reused while the trajectory has not grown. ../DATASETS/finite_size.csv
holds only what Figure1.py reads, plus the concentrations behind dG.
"""

import argparse
import re
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import MDAnalysis as mda
import yaml

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import DDM

CHUNK_SIZE = 20.0       # voxel edge, A; matches CHOSEN_VOXEL_SIZE in Figure1.py
HASH_RE = re.compile(r'_[0-9a-f]{6,}$')
# backup_<name>_<YYYYDDMM>_<HHh>MM<m>SS<s>.dcd, e.g.
# backup_Q8IUB9_1_90_a68a4f49_20261308_10h55m52s.dcd -> 2026-08-13 10:55:52
BACKUP_RE = re.compile(r'^backup_.*_(\d{8})_(\d{2})h(\d{2})m(\d{2})s\.dcd$')
ROOT = HERE                                  # run directories live beside this script
OUT_CSV = HERE.parent / 'DATASETS' / 'finite_size.csv'

CSV_COLUMNS = ['name', 'base', 'box_x', 'cdil', 'cdil_sem',
               'cden', 'cden_sem', 'dG', 'dG_err']


# ==========================================================================
# run discovery
# ==========================================================================

def run_directories(root):
    """Every directory at or below root (one level) with a top.pdb and a dcd."""
    root = Path(root)
    if not root.is_dir():
        return []
    dirs = []
    for cand in [root] + sorted(p for p in root.iterdir() if p.is_dir()):
        if (cand / "top.pdb").is_file() and list(cand.glob("*.dcd")):
            dirs.append(cand)
            continue
        # calv_comp layout: <seq>/293/
        sub = cand / "293"
        if (sub / "top.pdb").is_file() and list(sub.glob("*.dcd")):
            dirs.append(sub)
    return dirs


def dcd_sort_key(path):
    """Backups first (chronological, oldest to newest, by the timestamp in
    their filename), then the live trajectory file last -- it's the one
    still being appended to, so its frames come after every backup."""
    m = BACKUP_RE.match(path.name)
    if not m:
        return (1, path.name)
    date_str, hh, mm, ss = m.groups()
    stamp = datetime.strptime(f"{date_str}{hh}{mm}{ss}", "%Y%d%m%H%M%S")
    return (0, stamp)


def read_config(direc):
    path = direc / "config.yaml"
    if not path.is_file():
        return {}
    with open(path) as f:
        return yaml.safe_load(f) or {}


def run_name(direc, config):
    name = (config or {}).get('sysname')
    if name:
        return name
    return direc.parent.name if direc.name == "293" else direc.name


def base_name(name):
    """Strip the bookkeeping hash: Q3LI60_1_44_0a470d24 -> Q3LI60_1_44."""
    return HASH_RE.sub('', name)


def equilibration_start(config, n_frames_total, min_frames=200):
    """First frame to analyse. The 1e9-step rule, with a 10% fallback for
    runs that have not reached full length yet -- which is every run while
    the simulations are still going."""
    start = None
    wfreq = (config or {}).get('runtime_settings', {}).get('wfreq')
    if wfreq is None:
        wfreq = (config or {}).get('wfreq')
    if wfreq:
        try:
            start = int(1e9 / (10 * float(wfreq)))
        except (TypeError, ValueError, ZeroDivisionError):
            start = None
    if start is None:
        start = n_frames_total // 10
    if n_frames_total - start < min_frames:
        fallback = n_frames_total // 10
        if fallback != start:
            print(f"  ! start={start} leaves only {n_frames_total - start} "
                  f"frames; using {fallback}")
        start = fallback
    return max(int(start), 0)


# ==========================================================================
# per-run analysis
# ==========================================================================

def result_path(chunk_size, name):
    return HERE / f"results_{int(chunk_size)}" / f"{name}.npz"


def process_run(direc, chunk_size=CHUNK_SIZE, force=False):
    """Analyse one run with DDM.py, or reuse its npz if the dcd has not grown.
    Returns the CSV row."""
    pdb = direc / "top.pdb"
    dcds = sorted(direc.glob("*.dcd"), key=dcd_sort_key)
    if not pdb.is_file() or not dcds:
        return None

    config = read_config(direc)
    name = run_name(direc, config)
    out = result_path(chunk_size, name)
    n_frames_total = len(mda.Universe(str(pdb), [str(d) for d in dcds]).trajectory)

    # reuse the result only while the trajectory is the same length -- the
    # runs are still writing frames, so a longer dcd must be re-analysed
    cached = False
    if out.exists() and not force:
        with np.load(out) as d:
            cached_frames = int(d['total_frames'])
        if cached_frames == n_frames_total:
            print(f"  {name}: cached ({cached_frames} frames)")
            cached = True
        else:
            print(f"  {name}: dcd grew {cached_frames} -> {n_frames_total} "
                  f"frames, re-analysing")

    if not cached:
        start = equilibration_start(config, n_frames_total)
        print(f"\n=== {name}  ({direc}) ===")
        DDM.process_sample(pdb, dcds, name=name, start=start,
                           chunk_size=chunk_size, out=out)

    with np.load(out) as d:
        box = d['chunk_dims'] * d['n_chunks']
        return dict(name=name, base=base_name(name), box_x=float(box[0]),
                    cdil=float(d['cdil']), cdil_sem=float(d['cdil_sem']),
                    cden=float(d['cden']), cden_sem=float(d['cden_sem']),
                    dG=float(d['dG']), dG_err=float(d['dG_err']))


# ==========================================================================

def main(force=False):
    dirs = run_directories(ROOT)
    if not dirs:
        sys.exit(f"no runs (top.pdb + *.dcd) found under {ROOT}")

    print(f"{len(dirs)} run(s) to process at S = {CHUNK_SIZE:g} A")
    rows = []
    for d in dirs:
        try:
            r = process_run(d, chunk_size=CHUNK_SIZE, force=force)
        except Exception as e:
            print(f"  ! {d}: {type(e).__name__}: {e}")
            continue
        if r:
            rows.append(r)

    df = pd.DataFrame(rows, columns=CSV_COLUMNS).sort_values(['base', 'box_x'])
    print()
    print(df.to_string(index=False))
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT_CSV, index=False)
    print(f"wrote {OUT_CSV}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Finite-size check: DDM.py on runs in other box sizes.")
    parser.add_argument("--force", action="store_true",
                        help="Re-analyse even when a cached result exists.")
    main(parser.parse_args().force)
