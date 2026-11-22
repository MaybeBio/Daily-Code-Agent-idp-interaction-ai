#!/usr/bin/env python
from Bio.PDB import MMCIFParser, PDBIO
import sys


def cif_to_pdb(cif_file, pdb_file):
   parser = MMCIFParser()
   structure = parser.get_structure("structure", cif_file)
   io = PDBIO()
   io.set_structure(structure)
   io.save(pdb_file)
   

cif_to_pdb(sys.argv[1], sys.argv[2])
