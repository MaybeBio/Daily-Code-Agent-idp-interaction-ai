"""Compare two GROMACS topologies of the same polymer atom by atom (charges and atom types).

    python devtools/compare_topology.py REFERENCE.top OURS.top

The molecules are matched by their bond graph, so atom names and order may differ.
"""

import collections
import re
import sys

from rdkit import Chem

ELEMENTS = {12.011: 6, 1.008: 1, 15.9994: 8, 14.0067: 7, 32.06: 16, 30.97376: 15, 18.9984: 9}


def read(path):
    """First molecule of a .top/.itp: {atoms: {nr: {...}}, bonds: [(i, j)]}."""
    mol, sec = None, None
    for line in open(path, encoding="utf-8"):
        m = re.match(r"\s*\[\s*(\w+)\s*\]", line)
        if m:
            sec = m.group(1)
            continue
        body = line.split(";")[0].split()
        if not body or body[0].startswith("#"):
            continue
        if sec == "moleculetype":
            if mol is not None:
                break
            mol = {"atoms": {}, "bonds": []}
        elif mol is not None and sec == "atoms":
            mol["atoms"][int(body[0])] = {"type": body[1], "res": body[3], "name": body[4], "q": float(body[6]),
                                          "mass": float(body[7])}
        elif mol is not None and sec == "bonds":
            mol["bonds"].append((int(body[0]), int(body[1])))
    if mol is None:
        sys.exit(f"{path}: no [ moleculetype ]")
    return mol


def graph(mol):
    rw, pos = Chem.RWMol(), {}
    for i, a in mol["atoms"].items():
        atom = Chem.Atom(ELEMENTS.get(round(a["mass"], 4), 6))
        atom.SetNoImplicit(True)
        pos[i] = rw.AddAtom(atom)
    for i, j in mol["bonds"]:
        rw.AddBond(pos[i], pos[j], Chem.BondType.SINGLE)
    return rw.GetMol(), {v: k for k, v in pos.items()}


def main(ref_path, our_path, tol=1.5e-4):
    ref, ours = read(ref_path), read(our_path)
    (rg, rinv), (og, oinv) = graph(ref), graph(ours)
    match = og.GetSubstructMatch(rg)
    if not (len(match) == len(ref["atoms"]) == len(ours["atoms"])):
        print(f"different molecules: {len(ref['atoms'])} vs {len(ours['atoms'])} atoms")
        return 1
    diffs, types = collections.defaultdict(list), collections.Counter()
    for k, o in enumerate(match):
        ra, oa = ref["atoms"][rinv[k]], ours["atoms"][oinv[o]]
        if abs(ra["q"] - oa["q"]) > tol:
            diffs[(ra["name"], ra["type"], round(ra["q"], 6), oa["type"], round(oa["q"], 6))].append(ra["res"])
        if ra["type"] != oa["type"]:
            types[(ra["name"], ra["type"], oa["type"])] += 1
    n = sum(map(len, diffs.values()))
    print(f"{len(match)} atoms, same molecule; {len(match) - n} charges agree within {tol:g} e, {n} differ")
    for (name, rt, rq, ot, oq), res in diffs.items():
        where = ", ".join(sorted(set(res)))
        print(f"  {name:5s} reference {rt} {rq:+.4f}  ours {ot} {oq:+.4f}  ({len(res)} atoms in {where})")
    for (name, rt, ot), c in types.items():
        print(f"  atom type: {name} reference {rt}, ours {ot} ({c} atoms)")
    print(f"  total charge: reference {sum(a['q'] for a in ref['atoms'].values()):+.6f}, "
          f"ours {sum(a['q'] for a in ours['atoms'].values()):+.6f}")
    return 0


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    sys.exit(main(sys.argv[1], sys.argv[2]))
