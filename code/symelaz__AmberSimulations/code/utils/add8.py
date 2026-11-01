"""
add8.py <in.parm7> <out.parm7>

Adds the extra a99SB-disp backbone dihedral terms to an AMBER topology.
Same terms and parameters as the original pipeline (ff/add8_original_authors.py):
  psi: N(i)-CA(i)-C(i)-N(i+1)  for residues 2..n-1 of a chain (6 terms)
  phi: C(i-1)-N(i)-CA(i)-C(i)  for residues 2..n   of a chain (1 term)
Difference: chains are detected from the C(i)-N(i+1) bonds instead of assuming
one chain = residues 1..N, so multi-chain proteins and ligands/ions after the
protein are handled correctly. For a single chain the result is identical.
"""
import sys
from parmed.amber import AmberParm
from parmed.tools import addDihedral

PSI = [(-0.197, 1, 180.0), (-0.255, 1, 90.0), (0.075, 2, 90.0),
       (0.057, 3, 90.0), (0.025, 4, 90.0), (-0.011, 5, 90.0)]
PHI = [(0.40, 2, 90.0)]

parm = AmberParm(sys.argv[1])

def atom(res, name):
    for a in res.atoms:
        if a.name == name: return a
    return None

# protein residues = residues with N, CA, C; chains = runs connected by C(i)-N(i+1)
chains, current = [], []
for res in parm.residues:
    if not all(atom(res, n) for n in ('N', 'CA', 'C')):
        continue
    if current and atom(res, 'N') in atom(current[-1], 'C').bond_partners:
        current.append(res)
    else:
        if current: chains.append(current)
        current = [res]
if current: chains.append(current)

n_terms = 0
for chain in chains:
    first, last = chain[0].idx + 1, chain[-1].idx + 1    # 1-based residue numbers
    if last - first < 2:
        continue
    psi_i, psi_j = f":{first+1}-{last-1}", f":{first+2}-{last}"
    phi_p, phi_i = f":{first}-{last-1}", f":{first+1}-{last}"
    for k, per, ph in PSI:
        addDihedral(parm, psi_i+"@N", psi_i+"@CA", psi_i+"@C", psi_j+"@N", k, per, ph).execute()
        n_terms += last - first - 1
    for k, per, ph in PHI:
        addDihedral(parm, phi_p+"@C", phi_i+"@N", phi_i+"@CA", phi_i+"@C", k, per, ph).execute()
        n_terms += last - first
    print(f"chain residues {first}-{last}: extra backbone dihedrals added")

print(f"{len(chains)} chain(s), {n_terms} dihedral terms added")
parm.save(sys.argv[2], overwrite=True)
