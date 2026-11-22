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

# Open Babel encodes bond orders as REPEATED CONECT entries (a double bond
# C=O is written "CONECT c o o"). antechamber reads those as extra bonds, so
# e.g. every carboxylate carbon looks 4-connected and is typed c3 (sp3)
# instead of c, the carboxylate oxygens become os, and the prepin gets a
# wrong bond tree (tleap: "There is a bond of 3.06 angstroms ..."). Drop the
# CONECT records and the formal-charge column (Open Babel writes "O1-") so
# antechamber perceives connectivity and bond orders itself, from geometry
# and the -nc net charge, as it did before this protonation step existed.
lines = []
for line in open(out_pdb):
    if line.startswith("CONECT"):
        continue
    if line.startswith(("ATOM", "HETATM")):
        line = line[:76].ljust(76) + line[76:78] + "\n"  # element only
    lines.append(line)
open(out_pdb, "w").writelines(lines)

print(f"Protonated PDB saved to {out_pdb} at pH={pH}")
