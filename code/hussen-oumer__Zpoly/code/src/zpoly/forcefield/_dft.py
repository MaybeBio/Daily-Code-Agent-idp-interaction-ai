"""DFT charges"""

# rq-9a042027

from __future__ import annotations

import os
import re
from dataclasses import dataclass, field
from pathlib import Path

import yaml

from ..structure._build import _make_template
from ..zwitterions import Initiator, Library, Zwitterion, default_library
from ._errors import ForceFieldError
from ._ions import fragment, ionic_groups, linker_group

MANIFEST = "dft_manifest.yaml"
SUM_TOLERANCE = 0.01  # e
_SAFE = re.compile(r"[^A-Za-z0-9._-]")
_ELEMENT = re.compile(r"[A-Za-z]+")


class DftChargeError(ForceFieldError):
    """A DFT charge file does not fit its exported molecule"""

    def __init__(self, message: str, path: str | None = None, atom: int | None = None) -> None:
        self.path = path
        self.atom = atom
        super().__init__(f"{path}: {message}" if path else message)


# rq-d125cfec
@dataclass(frozen=True)
class DftCharges:
    by_template: dict = field(default_factory=dict)  # template SMILES -> charges (file order)
    sources: dict = field(default_factory=dict)  # molecule name -> file used
    missing: tuple[str, ...] = ()
    names: dict = field(default_factory=dict)  # template SMILES -> molecule name
    folder: str = ""
    ions_by_template: dict = field(default_factory=dict)  # template SMILES -> {template atom: charge}
    ion_names: dict = field(default_factory=dict)  # template SMILES -> fragment molecule names


# Export
def _molecule_files(folder: Path, name: str, mol, xyz, net_charge: int) -> list[str]:
    lines = [str(mol.GetNumAtoms()), f"{name}  charge={net_charge}  multiplicity=1  (zpoly DFT input)"]
    for atom in mol.GetAtoms():
        x, y, z = xyz[atom.GetIdx()]
        lines.append(f"{atom.GetSymbol():<3s}{x:14.6f}{y:14.6f}{z:14.6f}")
    (folder / f"{name}.xyz").write_text("\n".join(lines) + "\n", encoding="ascii")
    counters: dict[str, int] = {}
    pdb = [f"REMARK   {name}: charge {net_charge}, multiplicity 1 (zpoly DFT input)"]
    for atom in mol.GetAtoms():
        sym = atom.GetSymbol()
        counters[sym] = counters.get(sym, 0) + 1
        aname = f"{sym}{counters[sym]}"
        aname = f" {aname:<3s}" if len(sym) == 1 and len(aname) < 4 else f"{aname:<4s}"
        x, y, z = xyz[atom.GetIdx()]
        q = atom.GetFormalCharge()
        qtxt = "  " if q == 0 else f"{abs(q)}{'+' if q > 0 else '-'}"
        pdb.append(f"HETATM{atom.GetIdx() + 1:5d} {aname} DFT A   1    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00"
                   f"          {sym.upper():>2s}{qtxt}")
    for atom in mol.GetAtoms():
        nbrs = sorted(n.GetIdx() + 1 for n in atom.GetNeighbors())
        for k in range(0, len(nbrs), 4):
            pdb.append(f"CONECT{atom.GetIdx() + 1:5d}" + "".join(f"{x:5d}" for x in nbrs[k:k + 4]))
    pdb.append("END")
    (folder / f"{name}.pdb").write_text("\n".join(pdb) + "\n", encoding="ascii")
    return [f"{name}.xyz", f"{name}.pdb"]


def _entry(folder: Path, name: str, templates: list[str], parts: str = "whole") -> list[dict]:
    """Manifest entries: the molecule, plus its fragments with parts=ions."""
    tpl = _make_template(templates[0], seed=0)
    caps = sorted({c for c in (tpl.cap_in, tpl.cap_out) if c is not None}
                  | {n.GetIdx() for c in (tpl.cap_in, tpl.cap_out) if c is not None
                     for n in tpl.mol.GetAtomWithIdx(c).GetNeighbors() if n.GetAtomicNum() == 1})
    net = sum(a.GetFormalCharge() for a in tpl.mol.GetAtoms())
    files = _molecule_files(folder, name, tpl.mol, tpl.positions[0], net)
    entries = [{"name": name, "files": files, "templates": templates, "net_charge": net, "multiplicity": 1,
                "atoms": [a.GetSymbol() for a in tpl.mol.GetAtoms()], "caps": [c + 1 for c in caps]}]
    if parts == "ions":
        entries += _fragment_entries(folder, name, templates, tpl)
    return entries


def fragment_names(parent: str, groups) -> list[str]:
    """Fragment names: <parent>_cation, _anion, _ion."""
    kinds = ["cation" if g.charge > 0 else "anion" if g.charge < 0 else "ion" for g in groups]
    seen: dict[str, int] = {}
    names = []
    for kind in kinds:
        seen[kind] = seen.get(kind, 0) + 1
        names.append(f"{parent}_{kind}" + (str(seen[kind]) if kinds.count(kind) > 1 else ""))
    return names


# rq-d4bd3142
def _fragment_entries(folder: Path, parent: str, templates: list[str], tpl) -> list[dict]:
    groups = ionic_groups(tpl.mol)
    entries = []
    for name, group in zip(fragment_names(parent, groups), groups, strict=True):
        frag = fragment(tpl.mol, tpl.positions[0], group)
        files = _molecule_files(folder, name, frag.mol, frag.mol.GetConformer().GetPositions(), frag.charge)
        entries.append({"name": name, "kind": "ion", "parent": parent, "files": files, "templates": list(templates),
                        "net_charge": frag.charge, "multiplicity": 1,
                        "atoms": [a.GetSymbol() for a in frag.mol.GetAtoms()],
                        "parent_atoms": [p + 1 for p in frag.parent_atoms],
                        "caps": [[c + 1, a + 1] for c, a in frag.caps]})
    linker = linker_group(tpl.mol, tpl.anchor_out)
    if linker is not None:
        group, apply, kind = linker
        frag = fragment(tpl.mol, tpl.positions[0], group)
        name = f"{parent}_{kind}"
        files = _molecule_files(folder, name, frag.mol, frag.mol.GetConformer().GetPositions(), 0)
        entries.append({"name": name, "kind": "ion", "parent": parent, "files": files, "templates": list(templates),
                        "net_charge": 0, "multiplicity": 1,
                        "atoms": [a.GetSymbol() for a in frag.mol.GetAtoms()],
                        "parent_atoms": [p + 1 for p in frag.parent_atoms],
                        "caps": [[c + 1, a + 1] for c, a in frag.caps],
                        "apply": [k + 1 for k, p in enumerate(frag.parent_atoms) if p in set(apply)]})
    return entries


# rq-396afb23 rq-17608b5f
def export_dft_inputs(monomer: Zwitterion | str, out_dir: str | os.PathLike, *,
                      initiator: Initiator | str | None = None, library: Library | None = None,
                      parts: str = "whole") -> Path:
    """Write methyl-capped DFT molecules (.xyz/.pdb) and dft_manifest.yaml"""
    lib = library if library is not None else default_library()
    z = lib.get(monomer) if isinstance(monomer, str) else monomer
    if not z.polymerizable:
        raise ForceFieldError(f"'{z.name}' is not polymerizable; there is no repeat unit to export.")
    _check_parts(parts)
    folder = Path(out_dir)
    folder.mkdir(parents=True, exist_ok=True)
    entries = _entry(folder, _SAFE.sub("_", z.name), [z.repeat_unit_smiles()], parts)
    if initiator is not None:
        init = lib.initiator(initiator) if isinstance(initiator, str) else initiator
        base = _SAFE.sub("_", init.name)
        head = init.head_smiles.replace("[*:1]", "[*:2]") if init.head_smiles else None
        tail = init.tail_smiles
        if head and tail and init.head_smiles == init.tail_smiles:
            entries += _entry(folder, f"{base}_end", [head, tail])
        else:
            if head:
                entries += _entry(folder, f"{base}_head", [head])
            if tail:
                entries += _entry(folder, f"{base}_tail", [tail])
    return _write_manifest(folder, entries)


def _check_parts(parts: str) -> None:
    if parts not in ("whole", "ions"):
        raise ForceFieldError(f"parts must be 'whole' or 'ions', got {parts!r}")


# rq-d4bd3142
def export_dft_molecule(molecule: Zwitterion | str, out_dir: str | os.PathLike, *, parts: str = "whole",
                        library: Library | None = None) -> Path:
    """A small molecule (e.g."""
    lib = library if library is not None else default_library()
    z = lib.get(molecule) if isinstance(molecule, str) else molecule
    _check_parts(parts)
    folder = Path(out_dir)
    folder.mkdir(parents=True, exist_ok=True)
    return _write_manifest(folder, _entry(folder, _SAFE.sub("_", z.name), [z.smiles], parts))


def _write_manifest(folder: Path, entries: list[dict]) -> Path:
    from .. import __version__

    manifest = folder / MANIFEST
    old = []
    if manifest.is_file():
        data = yaml.safe_load(manifest.read_text(encoding="utf-8")) or {}
        old = [m for m in data.get("molecules", []) if m.get("name") not in {e["name"] for e in entries}]
    doc = {"generator": f"zpoly {__version__}", "molecules": old + entries}
    manifest.write_text("# DFT inputs written by Zpoly. Put each molecule's charges next to its .xyz as\n"
                        "# <name>.charges (plain list), <name>.mol2, <name>.log (Gaussian) or <name>.out (ORCA).\n"
                        "# Keep the atom order of the .xyz file.\n"
                        + yaml.safe_dump(doc, sort_keys=False, width=120), encoding="utf-8")
    return manifest


# Charge readers: [(element or None, charge)]
def _read_plain(path: Path) -> list:
    out = []
    for n, raw in enumerate(path.read_text(encoding="utf-8", errors="replace").splitlines(), 1):
        tok = raw.split("#", 1)[0].split()
        if not tok:
            continue
        try:
            if len(tok) == 1:
                out.append((None, float(tok[0])))
            elif len(tok) == 2:
                out.append((tok[0], float(tok[1])))
            elif len(tok) == 3:
                if int(tok[0]) != len(out) + 1:
                    raise DftChargeError(f"line {n}: index {tok[0]} is out of order", str(path), len(out) + 1)
                out.append((tok[1], float(tok[2])))
            else:
                raise ValueError
        except ValueError:
            raise DftChargeError(f"line {n}: expected 'q', 'element q' or 'index element q'", str(path)) from None
    return out


def _read_mol2(path: Path) -> list:
    out, in_atoms = [], False
    for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if raw.startswith("@<TRIPOS>"):
            in_atoms = raw.strip() == "@<TRIPOS>ATOM"
            continue
        if in_atoms and raw.strip():
            tok = raw.split()
            m = _ELEMENT.match(tok[1])
            out.append((m.group(0) if m else None, float(tok[8])))
    return out


def _blocks(lines: list[str], start: re.Pattern, row: re.Pattern) -> list[list]:
    """All blocks that begin at a `start` line and continue with matching `row` lines"""
    blocks, i = [], 0
    while i < len(lines):
        if start.search(lines[i]):
            j, rows = i + 1, []
            while j < len(lines) and not row.match(lines[j]) and len(rows) == 0 and j - i < 4:
                j += 1  # header / dashes lines
            while j < len(lines) and (m := row.match(lines[j])):
                rows.append((m.group(2), float(m.group(3))))
                j += 1
            if rows:
                blocks.append(rows)
            i = j
        else:
            i += 1
    return blocks


_G_ROW = re.compile(r"^\s*(\d+)\s+([A-Za-z]{1,2})\s+(-?\d+\.\d+)")
_O_ROW = re.compile(r"^\s*(\d+)\s+([A-Za-z]{1,2})\s*:\s*(-?\d+\.\d+)")
# exact headers only (not the H-summed blocks)
GAUSSIAN = {"esp": re.compile(r"^\s*ESP charges:\s*$"),
            "mulliken": re.compile(r"^\s*Mulliken charges( and spin densities)?:\s*$")}
ORCA = {"chelpg": re.compile(r"CHELPG Charges"), "mulliken": re.compile(r"MULLIKEN ATOMIC CHARGES")}


def _read_qm(path: Path, table: dict, row: re.Pattern, preference: list[str], charge_type: str | None) -> list:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    kinds = [charge_type] if charge_type else preference
    for kind in kinds:
        if kind not in table:
            raise DftChargeError(f"charge_type '{kind}' is not available for this file type "
                                 f"(choose from {', '.join(table)})", str(path))
        found = _blocks(lines, table[kind], row)
        if found:
            return found[-1]
    raise DftChargeError(f"no {' or '.join(kinds)} charge block found", str(path))


def read_charge_file(path: Path, charge_type: str | None = None) -> list:
    suffix = path.suffix.lower()
    if suffix == ".mol2":
        return _read_mol2(path)
    if suffix == ".log":
        if charge_type == "chelpg":
            charge_type = "esp"  # Gaussian prints ChelpG fits as "ESP charges"
        return _read_qm(path, GAUSSIAN, _G_ROW, ["esp", "mulliken"], charge_type)
    if suffix == ".out":
        return _read_qm(path, ORCA, _O_ROW, ["chelpg", "mulliken"], charge_type)
    return _read_plain(path)


def _same_element(given: str, expected: str) -> bool:
    """Element from a charge file vs"""
    g = given.strip()
    if g.capitalize() == expected or g[:2].capitalize() == expected:
        return True
    return len(expected) == 1 and g[:1].upper() == expected and g[:2].capitalize() not in _TWO_LETTER


_TWO_LETTER = {"Cl", "Br", "Na", "Li", "Mg", "Ca", "Si", "Zn", "Fe", "Cu", "Se"}


# Loading a folder
# rq-31ab6f96 rq-c0666681
def load_dft_charges(folder: str | os.PathLike) -> DftCharges:
    """Read dft_manifest.yaml and every available charge file"""
    base = Path(folder)
    manifest = base / MANIFEST
    try:
        data = yaml.safe_load(manifest.read_text(encoding="utf-8")) or {}
    except OSError as exc:
        raise DftChargeError(f"cannot read the manifest ({exc.strerror or exc}); run dft-export first",
                             str(manifest)) from exc
    except yaml.YAMLError as exc:
        raise DftChargeError(f"not valid YAML ({exc})", str(manifest)) from exc
    by_template, sources, names, missing = {}, {}, {}, []
    ions: dict = {}
    ion_names: dict = {}
    complete: dict = {}
    for mol in data.get("molecules", []):
        if mol.get("kind") == "ion":
            for template in mol["templates"]:
                ion_names.setdefault(template, []).append(mol["name"])
                complete.setdefault(template, True)
        name = mol["name"]
        candidates = [base / mol["charges"]] if mol.get("charges") else [
            base / f"{name}{ext}" for ext in (".charges", ".mol2", ".log", ".out")]
        path = next((c for c in candidates if c.is_file()), None)
        if path is None:
            missing.append(name)
            if mol.get("kind") == "ion":
                for template in mol["templates"]:
                    complete[template] = False
            continue
        rows = read_charge_file(path, mol.get("charge_type"))
        atoms = mol["atoms"]
        if len(rows) != len(atoms):
            raise DftChargeError(f"{len(rows)} charges for {name}, which has {len(atoms)} atoms", str(path))
        for k, ((element, _), expected) in enumerate(zip(rows, atoms, strict=True), 1):
            if element is not None and not _same_element(element, expected):
                raise DftChargeError(f"atom {k} is {element} in the charge file but {expected} in {name}.xyz "
                                     "(keep the exported atom order)", str(path), k)
        charges = tuple(q for _, q in rows)
        if abs(sum(charges) - mol["net_charge"]) > SUM_TOLERANCE:
            raise DftChargeError(f"charges sum to {sum(charges):.4f}, but {name} has net charge "
                                 f"{mol['net_charge']}", str(path))
        sources[name] = path.name
        if mol.get("kind") == "ion":
            folded = _fold_caps(charges, mol, str(path))
            for template in mol["templates"]:
                ions.setdefault(template, {}).update(folded)
            continue
        for template in mol["templates"]:
            by_template[template] = charges
            names[template] = name
    ions = {t: q for t, q in ions.items() if complete.get(t)}
    return DftCharges(by_template, sources, tuple(missing), names, str(base), ions,
                      {t: tuple(v) for t, v in ion_names.items()})


def _fold_caps(charges, mol: dict, path: str) -> dict[int, float]:
    """Fragment charges -> {parent template atom (0-based)"""
    q = list(charges)
    parents = mol.get("parent_atoms")
    if not isinstance(parents, list) or len(parents) != len(q):
        raise DftChargeError(f"manifest entry {mol['name']} has no valid parent_atoms", path)
    for cap, atom in mol.get("caps", []):
        q[atom - 1] += q[cap - 1]
    use = set(mol["apply"]) if "apply" in mol else set(range(1, len(q) + 1))
    return {p - 1: q[k] for k, p in enumerate(parents) if p > 0 and k + 1 in use}


# rq-82b577fa
def dft_charges_from_values(values, units, initiator: Initiator | None = None,
                            aliases: dict | None = None) -> DftCharges:
    """Charges pasted as {name"""
    aliases = aliases or {}
    targets: dict[str, list[str]] = {}
    for z in units:
        targets.setdefault(z.name, [z.repeat_unit_smiles()])
    for alias, real in aliases.items():
        if real in targets:
            targets.setdefault(alias, targets[real])
    if initiator is not None:
        head = initiator.head_smiles.replace("[*:1]", "[*:2]") if initiator.head_smiles else None
        tail = initiator.tail_smiles
        ends = [t for t in (head, tail) if t]
        if head and tail and initiator.head_smiles == initiator.tail_smiles:
            targets[initiator.name] = ends
            targets[f"{initiator.name}_end"] = ends
        if head:
            targets[f"{initiator.name}_head"] = [head]
        if tail:
            targets[f"{initiator.name}_tail"] = [tail]
    by_template, names = {}, {}
    for key, charges in values.items():
        if key not in targets:
            raise DftChargeError(f"'{key}' is not a monomer or end group of this chain "
                                 f"(use one of: {', '.join(sorted(targets))})")
        tpl = _make_template(targets[key][0], seed=0).mol
        q = tuple(float(x) for x in charges)
        if len(q) != tpl.GetNumAtoms():
            raise DftChargeError(f"{len(q)} charges for {key}, whose DFT molecule has {tpl.GetNumAtoms()} atoms "
                                 "(see dft_inputs for the atom order)")
        net = sum(a.GetFormalCharge() for a in tpl.GetAtoms())
        if abs(sum(q) - net) > SUM_TOLERANCE:
            raise DftChargeError(f"charges for {key} sum to {sum(q):.4f}, but its net charge is {net}")
        for template in targets[key]:
            by_template[template] = q
            names[template] = key
    return DftCharges(by_template, {k: "values" for k in values}, (), names, "")
