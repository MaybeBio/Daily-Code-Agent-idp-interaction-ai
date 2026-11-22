"""
find_disulfides.py <protein_pdb> <output_dir> [cutoff_angstrom]

Scans the PDB that build.sh is about to hand to tleap (build/protein_noh.pdb)
for cysteine pairs whose SG atoms are closer than <cutoff_angstrom> (default
2.5 A, the standard SG-SG bond length tolerance for a real disulfide vs. two
free cysteines that happen to sit nearby).

Residue numbers in the PDB are NOT what tleap uses internally: a single
`loadpdb` unit numbers its residues 1, 2, 3, ... in the order they are READ,
regardless of the PDB's own chain/resSeq (this is also why this pipeline's
README warns that multi-chain systems come out of tleap with consecutive,
not original, residue numbers). This script enumerates residues in the same
order tleap will, so the indices it writes are the ones a `bond receptor.N.SG
receptor.M.SG` command must use.

Side effects:
  <protein_pdb>  rewritten in place: CYS residues found in a disulfide are
                 renamed to CYX (tleap would otherwise add an HG to them,
                 leaving the pair unbonded with a dangling thiol).

Writes to <output_dir>:
  disulfides.leap   one `bond receptor.N.SG receptor.M.SG` line per pair,
                     spliced into the tleap script by build.sh's fill().
  disulfides.json   the same pairs with original (chain, resSeq) labels and
                     the measured SG-SG distance, for check_topology.py to
                     verify the bond actually made it into the topology.

Prints one line per pair found, and a one-line summary. Exits 0 either way
(no disulfides is a normal outcome, not an error) so build.sh always
continues to tleap.
"""
import json
import sys

PDB = sys.argv[1]
OUTDIR = sys.argv[2]
CUTOFF = float(sys.argv[3]) if len(sys.argv) > 3 else 2.5


def read_residues(path):
    """[(chain, resseq, icode), ...] in file order (first appearance only),
    plus {key: {'resname': str, 'atoms': {name: (x, y, z)}}}.

    This is the same order tleap assigns 1, 2, 3, ... to when it loads the
    file as a single unit, TER lines included (TER breaks the backbone bond
    between consecutive residues; it does not reset or skip the numbering).
    """
    order, residues = [], {}
    for line in open(path):
        if not line.startswith(("ATOM", "HETATM")):
            continue
        key = (line[21], line[22:26], line[26])
        if key not in residues:
            order.append(key)
            residues[key] = {"resname": line[17:20].strip(), "atoms": {}}
        name = line[12:16].strip()
        residues[key]["atoms"][name] = (
            float(line[30:38]), float(line[38:46]), float(line[46:54]))
    return order, residues


def dist(a, b):
    return sum((a[i] - b[i]) ** 2 for i in range(3)) ** 0.5


order, residues = read_residues(PDB)
leap_index = {key: i + 1 for i, key in enumerate(order)}  # 1-based, tleap's own numbering

cys = [k for k in order if residues[k]["resname"] in ("CYS", "CYX")
       and "SG" in residues[k]["atoms"]]

pairs = []
for i, a in enumerate(cys):
    for b in cys[i + 1:]:
        d = dist(residues[a]["atoms"]["SG"], residues[b]["atoms"]["SG"])
        if d < CUTOFF:
            pairs.append((a, b, d))

label = lambda k: "%s/%s" % (k[0], k[1].strip())  # noqa: E731

# Near misses: SG-SG between the cutoff and 6 A is almost always a native
# disulfide whose geometry was damaged upstream (loop modelling, grafting,
# a predicted structure), not two genuinely free cysteines. It is NOT bonded,
# but warn loudly so it is checked by hand.
NEAR = 6.0
for i, a in enumerate(cys):
    for b in cys[i + 1:]:
        d = dist(residues[a]["atoms"]["SG"], residues[b]["atoms"]["SG"])
        if CUTOFF <= d < NEAR:
            print(f"find_disulfides: WARNING near-miss CYS {label(a)} - CYS {label(b)} "
                  f"SG-SG {d:.2f} A (>= cutoff {CUTOFF} A, NOT bonded) - "
                  f"check whether this is a native disulfide")

if not pairs:
    print(f"find_disulfides: {len(cys)} CYS residue(s), no SG-SG pair within "
          f"{CUTOFF} A -> 0 disulfides built")
    open(f"{OUTDIR}/disulfides.leap", "w").close()
    json.dump([], open(f"{OUTDIR}/disulfides.json", "w"), indent=2)
    sys.exit(0)

# Rename the bonded cysteines to CYX in the file tleap is about to load.
bonded = {k for pair in pairs for k in pair[:2]}
lines = open(PDB).read().splitlines(keepends=True)
out = []
for line in lines:
    if line.startswith(("ATOM", "HETATM")):
        key = (line[21], line[22:26], line[26])
        if key in bonded and line[17:20].strip() == "CYS":
            line = line[:17] + "CYX" + line[20:]
    out.append(line)
open(PDB, "w").writelines(out)

leap_lines, json_pairs = [], []
for a, b, d in pairs:
    leap_lines.append(f"bond receptor.{leap_index[a]}.SG receptor.{leap_index[b]}.SG\n")
    json_pairs.append({
        "leap_index": [leap_index[a], leap_index[b]],
        "chain_resnum": [label(a), label(b)],
        "distance_A": round(d, 3),
    })
    print(f"find_disulfides: CYX {label(a)} - CYX {label(b)}  "
          f"(SG-SG {d:.2f} A, tleap residue {leap_index[a]}-{leap_index[b]})")

open(f"{OUTDIR}/disulfides.leap", "w").writelines(leap_lines)
json.dump(json_pairs, open(f"{OUTDIR}/disulfides.json", "w"), indent=2)
print(f"find_disulfides: {len(pairs)} disulfide bond(s) detected and will be built "
      f"({len(cys) - 2 * len(pairs)} free CYS left as-is)")
