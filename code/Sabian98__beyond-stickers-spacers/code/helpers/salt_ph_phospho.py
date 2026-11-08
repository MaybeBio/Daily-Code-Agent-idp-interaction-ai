"""salt_ph_phospho.py: definitions copied verbatim from the unpublished salt_ph_phospho.py (BENDER-MT) by tools/extract_helpers.py.
Do not edit by hand."""
import os
import sys

import numpy as np
import pandas as pd
from scipy import sparse

import md_loader as L  # noqa: E402
import basin_frames as bs  # noqa: E402
import fdt_params as fdt  # noqa: E402

bc = L.bc

T = 300.0

CAL = fdt.get_dielectric_params(T)

BETA_C = 1.0 / (fdt.R_GAS * T)

MAX_ELEMS = 4_000_000

def site_matrix(i, j, q, n):
    """Sparse (n x pairs): phi = W @ M.T gives each residue's potential from per-pair weights W (no self term)."""
    p = np.arange(len(i))
    return sparse.csr_matrix((np.concatenate([q[j], q[i]]), (np.concatenate([i, j]), np.concatenate([p, p]))), shape=(n, len(i)))

def calvados_terms(uid, hint):
    f = bs.frames_calvados(uid, 300, hint)
    if f is None:
        raise FileNotFoundError('no 300 K CALVADOS-2 frames')
    xyz, _, rg = f
    z = np.load(os.path.join(L.OUT_ROOT, 'npz', '300K', f'{uid}.npz'))
    seq = str(z['sequence'])
    q, lam, sig, _ = fdt.calvados_residue_params(seq, L.load_residue_params())
    i, j, *_ = fdt.pair_setup(q, lam, sig)
    M = site_matrix(i, j, q, len(seq))
    k, eps, rc = CAL['kappa'], CAL['eps_yu'], fdt.RC_YU
    ec = np.exp(-k * rc)
    nf = len(xyz)
    V, dVdk, phi = np.zeros(nf), np.zeros(nf), np.zeros((nf, len(seq)))
    step = max(1, MAX_ELEMS // max(1, len(i)))
    for a in range(0, nf, step):
        r = np.linalg.norm(xyz[a:a + step, i] - xyz[a:a + step, j], axis=2)
        inside = r < rc
        er = np.exp(-k * r)
        Y = np.where(inside, er / r - ec / rc, 0.0)
        V[a:a + step] = (q[i] * q[j] * eps * Y).sum(1)
        dVdk[a:a + step] = (q[i] * q[j] * eps * np.where(inside, ec - er, 0.0)).sum(1)
        phi[a:a + step] = (M @ (eps * Y).T).T
    vz = z['V_DH'].astype(float)
    if np.corrcoef(V, vz)[0, 1] < 0.999 or abs(V.mean() - vz.mean()) > 1e-3 * max(1.0, abs(vz.mean())) + 1e-3:
        raise ValueError('recomputed V_DH does not match the npz')
    dVdk_npz = (z['dVdT_DH'].astype(float) - CAL['deps_yu_dT'] / CAL['eps_yu'] * vz) / CAL['dkappa_dT']
    return seq, rg, CAL['kappa'] / 2 * dVdk, CAL['kappa'] / 2 * dVdk_npz, phi, BETA_C
