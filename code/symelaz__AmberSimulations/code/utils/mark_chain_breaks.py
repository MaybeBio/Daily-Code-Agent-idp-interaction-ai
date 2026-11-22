"""
mark_chain_breaks.py <protein_pdb> [max_peptide_CN_angstrom]

Re-inserts TER records into the PDB that build.sh hands to tleap
(build/protein_noh.pdb), rewriting it in place.

Why: extract_protein.tcl writes that file with VMD, and VMD's writepdb drops
every TER of the input. tleap only starts a new chain at a TER, so without
them it joins *every* break with a peptide bond - including intentional
gaps (unresolved loops left open) and loops that loop-closure failed to
close. tleap only prints "There is a bond of X angstroms between C and N
atoms" and carries on; the stretched bond is then pulled shut during
minimisation/equilibration, distorting the structure around it.

A TER is written between consecutive residues when
  - the chain ID changes, or
  - the C(i) - N(i+1) distance exceeds the cutoff (default 2.0 A; a real
    peptide bond is ~1.33 A), i.e. there is a physical gap in the backbone.
Geometry rather than the input's own TER lines is used on purpose: it is
independent of whatever tool wrote the input, and it also catches breaks the
input never marked.

TER does not change tleap's sequential residue numbering, so the indices
find_disulfides.py writes stay valid.
"""
import sys

PDB = sys.argv[1]
CUTOFF = float(sys.argv[2]) if len(sys.argv) > 2 else 2.0

lines = [l for l in open(PDB).read().splitlines(keepends=True) if not l.startswith("TER")]

# residues in file order: key -> {'first': line index, 'last': line index, 'C': xyz, 'N': xyz}
order, res = [], {}
for i, line in enumerate(lines):
    if not line.startswith(("ATOM", "HETATM")):
        continue
    key = (line[21], line[22:26], line[26])
    if not order or order[-1] != key:
        order.append(key)
        res[key] = {"first": i, "atoms": {}, "name": line[17:20].strip()}
    res[key]["last"] = i
    res[key]["atoms"][line[12:16].strip()] = (
        float(line[30:38]), float(line[38:46]), float(line[46:54]))


def dist(a, b):
    return sum((a[k] - b[k]) ** 2 for k in range(3)) ** 0.5


label = lambda k: "%s/%s%s" % (k[0], res[k]["name"], k[1].strip())  # noqa: E731
ter_after = set()  # line indices after which a TER goes
n_chain, n_gap = 0, 0
for a, b in zip(order, order[1:]):
    if a[0] != b[0]:
        ter_after.add(res[a]["last"]); n_chain += 1
        continue
    c, n = res[a]["atoms"].get("C"), res[b]["atoms"].get("N")
    if c is None or n is None:
        continue
    d = dist(c, n)
    if d > CUTOFF:
        ter_after.add(res[a]["last"]); n_gap += 1
        print(f"mark_chain_breaks: backbone gap {label(a)} - {label(b)}  "
              f"C-N {d:.2f} A -> TER (left open, not bonded)")

out = []
for i, line in enumerate(lines):
    if line.startswith("END") and out and not out[-1].startswith("TER"):
        out.append("TER\n")
    out.append(line)
    if i in ter_after:
        out.append("TER\n")
open(PDB, "w").writelines(out)
print(f"mark_chain_breaks: {len(order)} residues, {n_chain} chain boundary TER(s), "
      f"{n_gap} backbone gap TER(s) (C-N > {CUTOFF} A)")
