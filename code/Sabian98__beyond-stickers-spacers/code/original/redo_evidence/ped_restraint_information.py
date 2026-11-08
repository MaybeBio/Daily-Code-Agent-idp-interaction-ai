"""
ped_restraint_information.py — how much do PRE restraints move hub occupancy away from the unrestrained prior?
For the three PRE-restrained PED ensembles whose restraint-free IDPConformerGenerator starting pools are also in PED
(α-synuclein PED00024 / pool PED00425, Sic1 PED00160 / PED00423, drkN PED00157 / PED00427), compare per-residue hub
occupancy (exp_occupancy.py definition) between restrained ensemble and pool. High agreement means the hubs are mostly
inherited from the prior; only where they differ does the ensemble test experiment rather than the conformer generator.
  ~/.venvs/bender/bin/python ped_restraint_information.py  ->  ped_restraint_information.txt
"""
import math
import os
import sys
from multiprocessing import Pool

import numpy as np
from scipy.stats import spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import exp_occupancy as eo  # noqa: E402
import ped_null as pn  # noqa: E402

PAIRS = [('α-synuclein', 'PED00024/PED00024e001.pdb', 'PED00425/pdbfile.pdb'),
         ('Sic1', 'PED00160/PED00160e002.pdb', 'PED00423/pdbfile.pdb'),
         ('drkN', 'PED00157/PED00157e002.pdb', 'PED00427/pdbfile.pdb')]


def occ(path):
    seq, models = pn.load(path)
    if len(models) > 1000:
        models = models[::math.ceil(len(models) / 1000)]
    o, used = eo.occupancy(seq, models)
    return path, seq, np.array(o), used


def main():
    paths = [x for _, a, b in PAIRS for x in (a, b)]
    with Pool(6) as pool:
        res = {p: (s, o, u) for p, s, o, u in pool.imap_unordered(occ, paths)}
    lines = ['Hub occupancy: PRE-restrained PED ensemble vs its own restraint-free starting pool (IDPConformerGenerator)']
    for name, a, b in PAIRS:
        sa, oa, ua = res[a]
        sb, ob, ub = res[b]
        assert sa == sb, name
        n = len(sa)
        end = np.minimum(np.arange(n), n - 1 - np.arange(n))
        lines.append(f'  {name:12s} N={n:3d}  restrained ({ua} models) vs pool ({ub}): Spearman {spearmanr(oa, ob).correlation:+.2f}'
                     f' | position vs restrained {spearmanr(end, oa).correlation:+.2f}, vs pool {spearmanr(end, ob).correlation:+.2f}'
                     f' | mean |Δ occupancy| {np.mean(np.abs(oa - ob)):.3f} (mean occupancy {oa.mean():.3f})')
    print('\n'.join(lines))
    open(os.path.join(HERE, 'ped_restraint_information.txt'), 'w').write('\n'.join(lines) + '\n')


if __name__ == '__main__':
    main()
