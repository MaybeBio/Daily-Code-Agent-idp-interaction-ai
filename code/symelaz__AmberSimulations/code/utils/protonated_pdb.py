#!/usr/bin/env python3

import sys
from openbabel import pybel

if len(sys.argv) != 4:
    print("Usage: python protonate_pdb.py input.pdb output.pdb pH")
    sys.exit(1)

in_pdb  = sys.argv[1]
out_pdb = sys.argv[2]
pH      = float(sys.argv[3])

# Load input PDB
mol = next(pybel.readfile("pdb", in_pdb))

# Protonate at given pH using positional arguments
mol.OBMol.AddHydrogens(False, True, pH)

# Write protonated PDB
mol.write("pdb", out_pdb, overwrite=True)

print(f"Protonated PDB saved to {out_pdb} at pH={pH}")
