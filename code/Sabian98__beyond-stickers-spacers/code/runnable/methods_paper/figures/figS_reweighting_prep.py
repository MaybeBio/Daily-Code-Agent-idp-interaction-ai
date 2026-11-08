"""
figS_reweighting_prep.py — cache for figS_reweighting.py (visual companion of the charge-model reweighting): the stored
BENDER-MT CALVADOS-2 frames of one direct-set protein, the electrostatic potential at its three direct sites in every
frame, and the radius of gyration. Same functions as phospho_network.py (no new analysis).
  BENDER_REANALYSIS_DIR=$BSS_REPO/bender_reanalysis ~/.venvs/finches/bin/python figS_reweighting_prep.py
"""
import os
import sys

import numpy as np
import pandas as pd

EV = (__import__('os').environ['BSS_REPO'] + '/redo_evidence')
os.environ.setdefault('BENDER_REANALYSIS_DIR', (__import__('os').environ['BSS_REPO'] + '/bender_reanalysis'))
sys.path[:0] = [(__import__('os').environ['BSS_DATA'] + '/BENDER-MT/code'), EV]
import phospho_network as pn  # noqa: E402

UID = 'P28667'
HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    sd = pd.read_csv(os.path.join(EV, 'phospho_direct', 'sites_direct.csv'))
    s = sd[sd.UniProt_ID == UID].set_index('role')
    kingdom = pd.read_csv(pn.sp.bc.BIO_CSV, usecols=['UniProt_ID', 'kingdom']).set_index('UniProt_ID')['kingdom']
    seq, rg, _, _, phi, beta = pn.sp.calvados_terms(UID, kingdom.get(UID))
    xyz = pn.sp.bs.frames_calvados(UID, 300, kingdom.get(UID))[0]
    roles = ['low', 'partner', 'annotated']
    pos = np.array([int(s.loc[r, 'pos']) for r in roles])
    np.savez_compressed(os.path.join(HERE, 'cache', f'reweighting_{UID}.npz'), seq=seq, xyz=np.asarray(xyz, np.float32),
                        rg=np.asarray(rg), phi_sites=np.asarray(phi[:, pos - 1]), beta=beta, pos=pos, roles=np.array(roles),
                        log_ess=np.array([s.loc[r, 'log_ess'] for r in roles]), E=np.array([s.loc[r, 'E_occ'] for r in roles]))
    w = np.exp(-beta * pn.DQ * phi[:, pos - 1])
    w /= w.sum(0)
    print(UID, len(seq), 'residues', xyz.shape, 'beta', beta, 'sites', dict(zip(roles, pos)),
          'ESS frac', np.round(1 / (w ** 2).sum(0) / len(w), 3), 'stored log ESS', s.log_ess.round(3).to_dict())


if __name__ == '__main__':
    main()
