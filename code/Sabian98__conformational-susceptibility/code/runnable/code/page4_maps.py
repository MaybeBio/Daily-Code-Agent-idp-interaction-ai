"""
page4_maps.py — mean residue-residue distance maps of PAGE4 from the direct runs (page4_direct.load_frames: 3,000
evenly spaced production frames, chains unwrapped), at 300 K and 150 mM, for WT, HIPK1 (occupancy) and CLK2 (occupancy)
in CALVADOS-2 and Mpipi-Recharged. Used for the WT-to-CLK2 distance-change panel of the PAGE4 figure.
    BENDER_REANALYSIS_DIR=... python page4_maps.py       (openmpipi env)
Output: Traj Level Analysis/Analysis_FDT/total_response/page4/maps/<model>_<form>_300K.npz (mean_r, contact_p, n_frames)
"""
import sys

sys.dont_write_bytecode = True

import os  # noqa: E402

import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, HERE)
import page4_direct as PD  # noqa: E402

OUT = os.path.join(PD.BM, 'Traj Level Analysis/Analysis_FDT/total_response/page4/maps')
CUT = 1.2        # nm, contact cutoff between bead centres


def main():
    os.makedirs(OUT, exist_ok=True)
    for model in ('calvados', 'recharged'):
        r = PD.runs_table(model)
        r = r[(r.temp == 300) & (r.ionic == 0.15) & r.form.isin(['WT', 'HIPK1occ', 'CLK2occ']) & r.present]
        for _, row in r.iterrows():
            out = os.path.join(OUT, f'{model}_{row.form}_300K.npz')
            if os.path.exists(out):
                print('exists', out)
                continue
            _, xyz, idx, mb, seq = PD.load_frames(row)
            if PD.HC.bc.max_bond(xyz) > 1.0:
                raise ValueError(f'{model} {row.form}: chain still broken after unwrapping')
            n = xyz.shape[1]
            s, c = np.zeros((n, n)), np.zeros((n, n))
            for a in range(0, len(xyz), 250):
                d = np.linalg.norm(xyz[a:a + 250, :, None, :] - xyz[a:a + 250, None, :, :], axis=3)
                s += d.sum(0)
                c += (d < CUT).sum(0)
            np.savez_compressed(out, mean_r=s / len(xyz), contact_p=c / len(xyz), n_frames=len(xyz), sequence=seq)
            print('written', out, len(xyz))


if __name__ == '__main__':
    main()
