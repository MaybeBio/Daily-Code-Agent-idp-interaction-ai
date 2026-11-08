"""basin_frames.py: definitions copied verbatim from the unpublished basin_susceptibility.py (BENDER-MT) by tools/extract_helpers.py.
Do not edit by hand."""
import os
import sys

import numpy as np
import pandas as pd
from scipy import sparse

import md_loader as L  # noqa: E402

bc = L.bc

def _rg(xyz):
    return np.sqrt(((xyz - xyz.mean(1, keepdims=True)) ** 2).sum(2).mean(1))

def _check_rg(xyz, rg, tag):
    r = _rg(xyz)
    if len(r) != len(rg) or np.corrcoef(r, rg)[0, 1] < 0.999 or np.abs(r - rg).max() > 0.02 * rg.mean():
        raise ValueError(f'{tag}: trajectory frames do not match the npz (Rg check)')

def frames_calvados(uid, T, hint):
    base = L.TRAJ_ROOT_FMT.format(T=T)
    d = L.find_traj_dir(uid, base, hint, flat=None if T == 300 else True)
    if d is None:
        return None
    xyz, _, _ = bc.load_calvados(os.path.join(d, f'{uid}.dcd'), os.path.join(d, 'top.pdb'), k=L.K_FRAMES)
    p = os.path.join(L.OUT_ROOT, 'npz', f'{T}K', f'{uid}.npz')
    if not os.path.exists(p):
        return None
    z = np.load(p)
    G, rg = z['G_T'].astype(float), z['rg'].astype(float)
    _check_rg(xyz, rg, f'CALVADOS-2 {uid} {T} K')
    return xyz, G, rg
