"""md_loader.py: definitions copied verbatim from the unpublished TODO_T2_trajectory_loader.py (BENDER-MT) by tools/extract_helpers.py.
Do not edit by hand."""
import os
import sys

import numpy as np
import pandas as pd
from scipy import sparse

sys.path.insert(0, os.environ.get('BENDER_REANALYSIS_DIR', os.path.dirname(os.path.abspath(__file__))))
import bender_common as bc  # noqa: E402


BENDER_MT = (__import__('os').environ['BSS_DATA'] + '/BENDER-MT')

TRAJ_ROOT_FMT = BENDER_MT + '/BENDER_{T}K'

OUT_ROOT = BENDER_MT + '/Traj Level Analysis/Analysis_FDT'

RESIDUES_CSV_PATH = (__import__('os').environ['BSS_DATA'] + '/Scripts/CALVADOS IDP/input/residues_CALVADOS2.csv')

if not os.path.exists(RESIDUES_CSV_PATH):   # this Mac's home folder (2026-10-10)
    RESIDUES_CSV_PATH = os.path.expanduser((__import__('os').environ['BSS_DATA'] + '/Scripts/CALVADOS IDP/input/residues_CALVADOS2.csv'))

K_FRAMES = 1000

def load_residue_params(csv_path: str = RESIDUES_CSV_PATH) -> pd.DataFrame:
    """residues_CALVADOS2.csv indexed by one-letter code (columns lambdas, sigmas, q, bondlength)."""
    return pd.read_csv(csv_path).set_index('one')

def is_flat_layout(base: str, ids: list) -> bool:
    """BENDER_<T>/<ID>/ (278, 288, 310, 320 K) rather than BENDER_300K/<group>/<ID>/ (as in the v2 runner)."""
    return any(os.path.isdir(os.path.join(base, u)) for u in ids[:50])

def find_traj_dir(uid: str, base: str, hint: str = None, flat: bool = None):
    """
    Folder holding <uid>.dcd and top.pdb: BENDER_<T>/<uid> (flat) or BENDER_300K/<group>/<uid>.
    Group folders are listed when the protein is processed, so a refiled 300 K folder is still found:
    the BENDER_BIO kingdom first, 'Other' last (it duplicates the lineage folders). In a flat layout
    a missing ID is reported at once (never scan ~10k protein folders as if they were groups).
    """
    has = lambda d: os.path.isfile(os.path.join(d, f'{uid}.dcd')) and os.path.isfile(os.path.join(d, 'top.pdb'))  # noqa: E731
    if has(os.path.join(base, uid)):
        return os.path.join(base, uid)
    if flat is None:
        flat = is_flat_layout(base, [uid])
    if flat:
        return None
    groups = sorted(g for g in os.listdir(base) if os.path.isdir(os.path.join(base, g)))
    groups = ([hint] if hint in groups else []) + [g for g in groups if g not in (hint, 'Other')] + \
             (['Other'] if 'Other' in groups else [])
    for g in groups:
        if has(os.path.join(base, g, uid)):
            return os.path.join(base, g, uid)
    return None
