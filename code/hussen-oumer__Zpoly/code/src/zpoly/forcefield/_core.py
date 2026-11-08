"""Parametrize a built chain"""

from __future__ import annotations

import dataclasses
import functools
import logging
import os
import re
import tempfile
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from importlib import resources
from pathlib import Path

import yaml
from rdkit import Chem

from ..structure import PolymerStructure, build_polymer
from ._charges import apply_user_charges, check_charge_mapping, expand_entries, load_charge_file, round_to_target
from ._dft import DftCharges, fragment_names, load_dft_charges
from ._errors import ForceFieldError, MissingParameterError
from ._ions import ionic_groups
from ._model import FragmentParameters, ParameterSet, ResidueBlock, Topology
from ._params import apply_user_parameters, check_parameter_mapping, load_parameter_file, parameter_template
from ._registry import get_backend

log = logging.getLogger(__name__)
WILDCARD = "X"


# Parameter matching
def match_exact(table: Mapping, key: tuple[str, ...]):
    """Bond/angle lookup: exact key or its reverse"""
    if key in table:
        return table[key]
    return table.get(tuple(reversed(key)))


def match_wildcard(table: Mapping, key: tuple[str, ...]):
    """Dihedral/improper lookup"""
    best, best_x = None, None
    for cand, value in table.items():
        for oriented in (key, tuple(reversed(key))):
            if all(c == WILDCARD or c == k for c, k in zip(cand, oriented, strict=True)):
                n_x = sum(c == WILDCARD for c in cand)
                if best_x is None or n_x < best_x:
                    best, best_x = value, n_x
                break
    return best


# Bonded terms of a molecule
def bonded_terms(mol: Chem.Mol):
    """Return (bonds, angles, propers, pairs) as sorted tuples of 0-based atom indices"""
    nbrs = [sorted(n.GetIdx() for n in a.GetNeighbors()) for a in mol.GetAtoms()]
    bonds = sorted((min(b.GetBeginAtomIdx(), b.GetEndAtomIdx()), max(b.GetBeginAtomIdx(), b.GetEndAtomIdx()))
                   for b in mol.GetBonds())
    angles = sorted((i, j, k) for j in range(len(nbrs)) for i in nbrs[j] for k in nbrs[j] if i < k)
    propers = []
    for j, k in bonds:
        for i in nbrs[j]:
            if i == k:
                continue
            for el in nbrs[k]:
                if el in (j, i):
                    continue
                propers.append((i, j, k, el))
    propers.sort()
    close = {tuple(sorted(b)) for b in bonds} | {(min(i, k), max(i, k)) for i, _, k in angles}
    pairs = sorted({(min(i, el), max(i, el)) for i, _, _, el in propers} - close)
    return tuple(bonds), tuple(angles), tuple(propers), tuple(pairs)


def bonded_classes(types: Sequence[str], params: ParameterSet) -> list[str]:
    """Bond type of each atom."""
    return [params.atomtypes[t].bonded_key if t in params.atomtypes else t for t in types]


def improper_terms(mol: Chem.Mol, types: Sequence[str], params: ParameterSet):
    """Impropers of the sp2 atoms."""
    types = bonded_classes(types, params)
    out = []
    for atom in mol.GetAtoms():
        if atom.GetDegree() != 3 or atom.GetHybridization() != Chem.HybridizationType.SP2:
            continue
        c = atom.GetIdx()
        n1, n2, n3 = sorted((n.GetIdx() for n in atom.GetNeighbors()), key=lambda i: (types[i], i))
        if match_wildcard(params.impropertypes, (types[n1], types[n2], types[c], types[n3])) is not None:
            out.append((n1, n2, c, n3))
    return tuple(out)


def find_missing(params: ParameterSet, types, bonds, angles, propers) -> dict[str, list]:
    """Bonded type combinations without parameters."""
    types = bonded_classes(types, params)
    missing: dict[str, list] = {}
    for kind, terms, table, matcher in (
        ("bond", bonds, params.bondtypes, match_exact),
        ("angle", angles, params.angletypes, match_exact),
        ("dihedral", propers, params.dihedraltypes, match_wildcard),
    ):
        keys = set()
        for t in terms:
            key = tuple(types[i] for i in t)
            if key not in keys and tuple(reversed(key)) not in keys and matcher(table, key) is None:
                keys.add(min(key, tuple(reversed(key))))
        if keys:
            missing[kind] = sorted(keys)
    return missing


# Merging parameter sets from several contexts
def merge_parameters(sets: Sequence[ParameterSet]) -> ParameterSet:
    first = sets[0]
    atomtypes: dict = {}
    tables = {"bondtypes": {}, "angletypes": {}, "dihedraltypes": {}, "impropertypes": {}}
    for ps in sets:
        if (ps.defaults != first.defaults or tuple(ps.rtp_bondedtypes) != tuple(first.rtp_bondedtypes)
                or ps.dihedral_funct != first.dihedral_funct):
            raise ForceFieldError("Contexts returned different [ defaults ], [ bondedtypes ] or dihedral styles.")
        for name, at in ps.atomtypes.items():
            old = atomtypes.get(name)
            if old is None:
                atomtypes[name] = at
            elif (old.atomic_number != at.atomic_number or old.bonded_key != at.bonded_key
                  or abs(old.mass - at.mass) > 1e-6
                  or abs(old.sigma - at.sigma) > 1e-9 or abs(old.epsilon - at.epsilon) > 1e-9):
                raise ForceFieldError(f"Atom type '{name}' is defined differently by two contexts.")
        for attr, merged in tables.items():
            for key, value in getattr(ps, attr).items():
                if key not in merged:
                    merged[key] = value
                elif merged[key] != value:
                    log.warning("Keeping the first %s parameters for %s; a later context differs.", attr, key)
    return ParameterSet(first.defaults, tuple(first.rtp_bondedtypes), atomtypes, **tables,
                        dihedral_funct=first.dihedral_funct)


# Charge normalization
def normalize_charges(charges: Sequence[float], target: int) -> tuple[float, ...]:
    """Shift to the integer target, then round to 4 decimals so the sum is exact"""
    n = len(charges)
    shifted = [q + (target - sum(charges)) / n for q in charges]
    return round_to_target(shifted, target, range(n))


# Context oligomers
@dataclass
class _Context:
    structure: PolymerStructure
    # residue name -> resnum in the context
    residues: dict[str, int]


def _residue_atoms(mol: Chem.Mol) -> dict[int, list[int]]:
    out: dict[int, list[int]] = {}
    for a in mol.GetAtoms():
        out.setdefault(a.GetIntProp("resnum"), []).append(a.GetIdx())
    return out


def _residue_names(mol: Chem.Mol) -> dict[int, str]:
    return {a.GetIntProp("resnum"): a.GetProp("res") for a in mol.GetAtoms()}


def _contexts(structure: PolymerStructure) -> list[_Context]:
    names = _residue_names(structure.mol)
    order = [names[k] for k in sorted(names)]
    units = list(structure.units)
    chain = structure.chain
    head = 1 if chain.head_group is not None else 0
    n = len(units)
    if n == 1:
        return [_Context(structure, {name: i + 1 for i, name in enumerate(order)})]

    def build(seq, ends):
        return build_polymer(seq, initiator=structure.initiator if ends != "none" else None, ends=ends,
                             conformation="coil", seed=structure.seed, relax_iterations=0, cap=chain.cap)

    contexts = []
    first_name = order[head]
    head_ends = "head" if chain.head_group is not None else "none"
    f_map = {first_name: 1 + head}
    if head:
        f_map[order[0]] = 1
    contexts.append(_Context(build([units[0]] * 3, head_ends), f_map))

    seen_r: dict[str, str] = {}  # repeat residue name -> monomer
    for i in range(1, n - 1):
        res = order[head + i]
        if res not in seen_r:
            seen_r[res] = units[i]
            contexts.append(_Context(build([units[i]] * 3, "none"), {res: 2}))

    last_name = order[head + n - 1]
    tail_ends = "tail" if chain.tail_group is not None else "none"
    e_map = {last_name: 3}
    if chain.tail_group is not None:
        e_map[order[-1]] = 4
    contexts.append(_Context(build([units[-1]] * 3, tail_ends), e_map))
    return contexts


# Public entry point
# rq-1229e2e7 rq-2dccba6a rq-f928598f rq-79e07d8e
def parametrize(
    structure: PolymerStructure,
    force_field: str = "gaff2",
    *,
    name: str | None = None,
    work_dir: str | os.PathLike | None = None,
    charges: str | os.PathLike | Mapping | None = None,
    params: str | os.PathLike | Mapping | None = None,
    strict: bool = True,
    dft_charges=None,
    dft_parts: str = "whole",
    charge_scale: float = 1.0,
    builtin_charges: bool = True,
) -> Topology:
    """Parametrize a built chain with a backend and assemble its F/R/E topology"""
    backend = get_backend(force_field)
    # user files first: fail before slow backend runs
    charge_path = None
    if charges is None:
        user_charges: dict[str, dict[str, float]] = {}
    elif isinstance(charges, Mapping):
        user_charges = check_charge_mapping(charges)
    else:
        charge_path = str(charges)
        user_charges = load_charge_file(charges)
    # rq-5a64d41e
    param_path = None
    if params is None:
        user_params: dict = {}
    elif isinstance(params, Mapping) and "version" in params:
        user_params = check_parameter_mapping(params)
    elif isinstance(params, Mapping):
        user_params = dict(params)  # already checked (load_parameter_file output)
    else:
        param_path = str(params)
        user_params = load_parameter_file(params)
    # rq-905bf2b0
    if dft_charges is not None and not isinstance(dft_charges, DftCharges):
        dft_charges = load_dft_charges(dft_charges)
    if dft_parts not in ("whole", "ions"):
        raise ForceFieldError(f"dft_parts must be 'whole' or 'ions', got {dft_parts!r}")
    if dft_parts == "ions" and dft_charges is None:
        raise ForceFieldError("dft_parts='ions' needs DFT charges (a folder from JOB = dft with DFT_PARTS = ions)")
    if isinstance(charge_scale, bool) or not isinstance(charge_scale, (int, float)) or not 0 < charge_scale <= 1:
        raise ForceFieldError(f"CHARGE_SCALE must be in (0, 1], got {charge_scale!r}")
    backend.check_available()
    if not structure.units:
        raise ForceFieldError("The structure has no monomer units (build it with build_polymer).")
    dft_blocks, dft_missing = _dft_coverage(structure, dft_charges)
    ion_blocks, ion_missing = _ion_coverage(structure, dft_charges) if dft_parts == "ions" else ({}, ())
    need_backend_charges = dft_charges is None or bool(dft_missing)

    contexts = _contexts(structure)
    with tempfile.TemporaryDirectory(prefix="polymer_ff_") as tmp:
        base = Path(work_dir) if work_dir is not None else Path(tmp)
        results: list[FragmentParameters] = []
        for i, ctx in enumerate(contexts):
            wd = base / f"context{i + 1}"
            wd.mkdir(parents=True, exist_ok=True)
            results.append(_call_backend(backend, ctx.structure.mol, ctx.structure.chain.net_charge, wd,
                                         need_backend_charges))
            if len(results[-1].atom_types) != ctx.structure.mol.GetNumAtoms():
                raise ForceFieldError(f"Backend '{backend.name}' returned the wrong number of atoms.")
    merged = merge_parameters([r.parameters for r in results])
    if user_params:
        full = backend.full_parameters() if hasattr(backend, "full_parameters") else None
        merged, user_keys = apply_user_parameters(merged, user_params, param_path, full)
    else:
        user_keys = ()

    chain_mol = structure.mol
    chain_res = _residue_atoms(chain_mol)
    chain_names = _residue_names(chain_mol)
    first_of: dict[str, int] = {}
    for resnum in sorted(chain_names):
        first_of.setdefault(chain_names[resnum], resnum)

    builtin = _builtin_entries(structure) if builtin_charges else {}
    builtin_used: list[str] = []
    ignored = sorted(set(user_charges) - set(first_of))
    if ignored:
        log.warning("Charge file residues not in this chain were ignored: %s", ", ".join(ignored))
    reports = []

    # residue name -> (names, types, charges)
    block_data: dict[str, tuple[tuple[str, ...], tuple[str, ...], tuple[float, ...], int]] = {}
    for ctx, res in zip(contexts, results, strict=True):
        cmol = ctx.structure.mol
        catoms = _residue_atoms(cmol)
        for res_name, cres in ctx.residues.items():
            if res_name in block_data:
                continue
            chain_idx = chain_res[first_of[res_name]]
            ctx_idx = catoms[cres]
            names = tuple(chain_mol.GetAtomWithIdx(i).GetProp("name") for i in chain_idx)
            if names != tuple(cmol.GetAtomWithIdx(i).GetProp("name") for i in ctx_idx):
                raise ForceFieldError(f"Internal error: context atoms of residue {res_name} do not match the chain.")
            target = sum(chain_mol.GetAtomWithIdx(i).GetFormalCharge() for i in chain_idx)
            raw = [res.charges[i] for i in ctx_idx]
            if res_name in dft_blocks:
                # DFT charges by template atom, shifted to the target
                dft = dft_charges.by_template[structure.templates[dft_blocks[res_name]]]
                mapped = [dft[chain_mol.GetAtomWithIdx(i).GetIntProp("template_atom")]
                          if chain_mol.GetAtomWithIdx(i).GetIntProp("template_atom") >= 0 else 0.0
                          for i in chain_idx]
                raw = list(normalize_charges(mapped, target))
            fixed: dict[str, float] = {}
            hydrogens_of = {
                chain_mol.GetAtomWithIdx(i).GetProp("name"): [
                    n.GetProp("name") for n in chain_mol.GetAtomWithIdx(i).GetNeighbors()
                    if n.GetAtomicNum() == 1 and n.GetIdx() in set(chain_idx)]
                for i in chain_idx}
            if res_name in builtin and res_name not in dft_blocks and res_name not in ion_blocks:
                # built-in RESP charges
                fixed.update(expand_entries(res_name, builtin[res_name], names, hydrogens_of, "built-in charges"))
                builtin_used.append(res_name)
            if res_name in ion_blocks:
                # rq-d4bd3142: ion-fragment charges
                ion_q = dft_charges.ions_by_template[structure.templates[ion_blocks[res_name]]]
                for i in chain_idx:
                    k = chain_mol.GetAtomWithIdx(i).GetIntProp("template_atom")
                    if k in ion_q:
                        fixed[chain_mol.GetAtomWithIdx(i).GetProp("name")] = ion_q[k]
            user_atoms = {}
            if res_name in user_charges:
                # rq-3541c593: user charges stay fixed
                user_atoms = expand_entries(res_name, user_charges[res_name], names, hydrogens_of, charge_path)
                fixed.update(user_atoms)
            if fixed:
                charges, report = apply_user_charges(res_name, names, raw, target, fixed, charge_path)
                if user_atoms:
                    reports.append(dataclasses.replace(report, user_atoms=len(user_atoms)))
            else:
                charges = normalize_charges(raw, target)
            block_data[res_name] = (names, tuple(res.atom_types[i] for i in ctx_idx), charges, target)

    missing = [r for r in first_of if r not in block_data]
    if missing:
        raise ForceFieldError(f"Internal error: no context for residue(s) {', '.join(missing)}.")
    if charge_scale != 1:
        block_data = _scale_blocks(block_data, chain_mol, chain_res, first_of, float(charge_scale))

    n_atoms = chain_mol.GetNumAtoms()
    types: list[str] = [""] * n_atoms
    charges: list[float] = [0.0] * n_atoms
    for resnum, idx in chain_res.items():
        names, btypes, bcharges, _ = block_data[chain_names[resnum]]
        if tuple(chain_mol.GetAtomWithIdx(i).GetProp("name") for i in idx) != names:
            raise ForceFieldError(f"Residue {resnum} differs from its {chain_names[resnum]} block.")
        for i, t, q in zip(idx, btypes, bcharges, strict=True):
            types[i], charges[i] = t, q

    missing_types = sorted(set(types) - set(merged.atomtypes))
    if missing_types:
        raise ForceFieldError(f"No atom-type definition for {', '.join(missing_types)}.")
    bonds, angles, propers, pairs = bonded_terms(chain_mol)
    explicit = _explicit_terms(chain_mol, merged, types, (bonds, angles, propers)) if backend.name == "oplsaa" else {}
    missing = find_missing(merged, types,
                           *([t for t in terms if t not in explicit] for terms in (bonds, angles, propers)))
    if missing and strict:
        raise MissingParameterError(missing, parameter_template(missing, merged.dihedral_funct))
    impropers = improper_terms(chain_mol, types, merged)

    blocks = tuple(_make_block(res_name, first_of[res_name], block_data[res_name], chain_mol, chain_res, impropers)
                   for res_name in first_of)
    return Topology(
        name=name or "polymer", force_field=backend.name, structure=structure, blocks=blocks,
        atom_types=tuple(types), charges=tuple(charges), parameters=merged, bonds=bonds, pairs=pairs,
        angles=angles, propers=propers, impropers=impropers, net_charge=structure.chain.net_charge,
        charge_report=tuple(sorted(reports, key=lambda r: list(first_of).index(r.residue))),
        dft_blocks=tuple(r for r in first_of if r in dft_blocks), dft_missing=dft_missing + ion_missing,
        ion_blocks=tuple(r for r in first_of if r in ion_blocks), charge_scale=float(charge_scale),
        builtin_blocks=tuple(r for r in first_of if r in builtin_used), explicit=explicit,
        user_parameters=user_keys, missing_parameters=missing,
    )


def _dft_coverage(structure: PolymerStructure, dft) -> tuple[dict[str, int], tuple[str, ...]]:
    """Residue name -> template index for residues with DFT charges"""
    if dft is None:
        return {}, ()
    covered: dict[str, int] = {}
    missing: list[str] = []
    for atom in structure.mol.GetAtoms():
        res, tid = atom.GetProp("res"), atom.GetIntProp("template")
        if res in covered or tid < 0:
            continue
        template = structure.templates[tid]
        if template in dft.by_template:
            covered[res] = tid
        else:
            label = dft.names.get(template) or _template_label(structure, tid)
            if label not in missing:
                missing.append(label)
    return covered, tuple(missing)


# rq-d4bd3142
def _ion_coverage(structure: PolymerStructure, dft: DftCharges) -> tuple[dict[str, int], tuple[str, ...]]:
    """Residues of monomer templates with complete ion-fragment charges"""
    from ..structure._build import _make_template

    unit_templates = {z.repeat_unit_smiles() for z in structure.units}
    covered: dict[str, int] = {}
    missing: list[str] = []
    for atom in structure.mol.GetAtoms():
        res, tid = atom.GetProp("res"), atom.GetIntProp("template")
        if res in covered or tid < 0 or structure.templates[tid] not in unit_templates:
            continue
        template = structure.templates[tid]
        if template in dft.ions_by_template:
            covered[res] = tid
            continue
        names = dft.ion_names.get(template) or fragment_names(
            _template_label(structure, tid), ionic_groups(_make_template(template, seed=0).mol))
        missing += [n for n in names if n not in missing]
    return covered, tuple(missing)


# rq-12060285
def _scale_blocks(block_data, chain_mol, chain_res, first_of, s: float):
    """Ionic-group atoms x s"""
    # groups stay within their charged atoms' residues
    res_of = {a.GetIdx(): a.GetIntProp("resnum") for a in chain_mol.GetAtoms()}
    ionic = set()
    for g in ionic_groups(chain_mol):
        charged = {res_of[i] for i in g.atoms if chain_mol.GetAtomWithIdx(i).GetFormalCharge() != 0}
        ionic |= {i for i in g.atoms if res_of[i] in charged}
    out = {}
    for res_name, (names, types, charges, target) in block_data.items():
        idx = chain_res[first_of[res_name]]
        scaled = [q * s if i in ionic else q for i, q in zip(idx, charges, strict=True)]
        new_target = round(target * s, 4)
        rest = [k for k, i in enumerate(idx) if i not in ionic] or list(range(len(idx)))
        shift = (new_target - sum(scaled)) / len(rest)
        for k in rest:
            scaled[k] += shift
        out[res_name] = (names, types, round_to_target(scaled, new_target, rest), new_target)
    return out


@functools.cache
def _data(name: str) -> dict:
    with resources.files("zpoly.forcefield").joinpath(f"data/{name}").open(encoding="utf-8") as fh:
        return yaml.safe_load(fh)


def _builtin_entries(structure: PolymerStructure) -> dict[str, dict]:
    """Built-in charge entries per residue name."""
    from ..zwitterions import normalize_name

    data = _data("builtin_charges.yaml")
    monomers = {normalize_name(k): v for k, v in data["monomers"].items()}
    by_template = {}
    for z in structure.units:  # by name or alias (A3361 is an alias of SBAA)
        hit = next((monomers[normalize_name(n)] for n in (z.name, *z.aliases) if normalize_name(n) in monomers), None)
        m = re.fullmatch(r"(.+)-sp(\d+)", z.name)
        if hit is None and m and int(m[2]) >= 2 and normalize_name(m[1]) in monomers:
            hit = _spacer_variant(monomers[normalize_name(m[1])], m[1], z)
        by_template[z.repeat_unit_smiles()] = hit
    out: dict[str, dict] = {}
    for atom in structure.mol.GetAtoms():
        res, tid = atom.GetProp("res"), atom.GetIntProp("template")
        if res in out or tid < 0:
            continue
        template = structure.templates[tid]
        if template in by_template:
            entry = by_template[template]
        else:
            entry = data["initiators"].get(structure.initiator.name) if structure.initiator is not None else None
        if entry is not None:
            out[res] = entry["atoms"]
    return out


# rq-c612302d: OPLS bonded type -> GAFF2 type
TO_GAFF2 = {"CT": "c3", "HC": "hc", "C_2": "c", "C": "c", "CA": "ca", "O_2": "o", "O": "o", "OS": "os", "N": "n"}
_SULFATE = Chem.MolFromSmarts("[SX4;$([SX4][OX2][#6])](=O)(=O)[O-]")
_CHEBYSHEV = ((1,), (0, 1), (-1, 0, 2), (0, -3, 0, 4), (1, 0, -8, 0, 8), (0, 5, 0, -20, 0, 16))


# rq-a67a239e
def periodic_to_rb(terms) -> tuple[float, ...]:
    """Periodic dihedral terms -> Ryckaert-Bellemans C0..C5."""
    rb = [0.0] * 6
    for phase, k, n in terms:
        if n > 5 or round(phase) % 180:
            raise ForceFieldError(f"Cannot convert the dihedral term (phase {phase}, n {n}) to Ryckaert-Bellemans.")
        sign = (1 if round(phase) % 360 == 0 else -1) * (-1) ** n
        rb[0] += k
        for power, coeff in enumerate(_CHEBYSHEV[n]):
            rb[power] += k * sign * coeff
    return tuple(round(c, 6) for c in rb)


def _explicit_terms(mol, params: ParameterSet, types, terms_by_kind) -> dict:
    """Values written on each bonded line"""
    keys = [params.atomtypes[t].bonded_key for t in types]
    gaff2 = {i: "s6" for (i, *_r) in mol.GetSubstructMatches(_SULFATE)}
    gaff2.update({n.GetIdx(): "o" if n.GetDegree() == 1 else "os" for i in list(gaff2)
                  for n in mol.GetAtomWithIdx(i).GetNeighbors()})
    g = _data("gaff2_so4.yaml")
    b = _data("oplsaa_params.yaml")
    sources = (("bonds", params.bondtypes, match_exact), ("angles", params.angletypes, match_exact),
               ("dihedrals", params.dihedraltypes, match_wildcard))
    out = {}
    for terms, (section, table, match) in zip(terms_by_kind, sources, strict=True):
        gaff = {tuple(k.split("-")): v for k, v in g[section].items()}
        own = {tuple(e["types"]): tuple(v for k, v in e.items() if k != "types") for e in b.get(section, [])}
        for term in terms:
            if any(i in gaff2 for i in term):
                value = match(gaff, tuple(gaff2.get(i) or TO_GAFF2.get(keys[i], "?") for i in term))
                if value is not None and section == "dihedrals":
                    value = periodic_to_rb(value)
            elif match(table, tuple(keys[i] for i in term)) is None:
                value = match_exact(own, tuple(keys[i] for i in term))
            else:
                value = None
            if value is not None:
                out[term] = tuple(value)
    return out


def _spacer_variant(entry: dict, parent: str, variant) -> dict | None:
    """The parent's built-in charges on a spacer variant"""
    mols = [Chem.RemoveHs(build_polymer([z], cap="H").mol) for z in (parent, variant)]
    pmol, vmol = mols
    index = {a.GetProp("name"): a.GetIdx() for a in pmol.GetAtoms()}
    listed = {index[n] for n in entry["atoms"] if n in index}
    out, seen = {}, set()
    for start in listed:
        if start in seen:
            continue
        group, todo = set(), [start]
        while todo:
            i = todo.pop()
            if i not in group:
                group.add(i)
                todo += [n.GetIdx() for n in pmol.GetAtomWithIdx(i).GetNeighbors() if n.GetIdx() in listed]
        seen |= group
        symbols = [f"[#{a.GetAtomicNum()}&H{a.GetTotalNumHs()}&D{a.GetDegree()}&{a.GetFormalCharge():+d}]"
                   for a in pmol.GetAtoms()]
        order = sorted(group)
        query = Chem.MolFromSmarts(Chem.MolFragmentToSmiles(pmol, order, atomSymbols=symbols, canonical=False))
        match = vmol.GetSubstructMatch(query)
        if len(match) != len(order):
            return None
        names = [pmol.GetAtomWithIdx(i).GetProp("name") for i in order]
        for name, j in zip(names, match, strict=True):
            out[vmol.GetAtomWithIdx(j).GetProp("name")] = entry["atoms"][name]
    return {"atoms": out}


def _template_label(structure: PolymerStructure, tid: int) -> str:
    """A readable name for a template (monomer name, or <INITIATOR>_end) for notices"""
    template = structure.templates[tid]
    for z in structure.units:
        if z.repeat_unit_smiles() == template:
            return z.name
    init = structure.initiator
    return f"{init.name}_end" if init is not None else template


def _call_backend(backend, mol, net_charge, work_dir, charges: bool) -> FragmentParameters:
    """Call the backend"""
    if charges:
        return backend.parametrize(mol, net_charge, work_dir)
    try:
        return backend.parametrize(mol, net_charge, work_dir, charges=False)
    except TypeError:
        return backend.parametrize(mol, net_charge, work_dir)


def _make_block(res_name, resnum, data, mol, chain_res, impropers) -> ResidueBlock:
    names, types, charges, target = data
    own = chain_res[resnum]
    own_set = set(own)
    prev_set = set(chain_res.get(resnum - 1, ()))
    next_set = set(chain_res.get(resnum + 1, ()))

    def label(i: int) -> str:
        nm = mol.GetAtomWithIdx(i).GetProp("name")
        if i in own_set:
            return nm
        return ("-" if i in prev_set else "+") + nm

    bonds = []
    for i in own:
        for nb in mol.GetAtomWithIdx(i).GetNeighbors():
            j = nb.GetIdx()
            if j in own_set and i < j:
                bonds.append((label(i), label(j)))
            elif j in prev_set:
                bonds.append((label(j), label(i)))
    imps = tuple(tuple(label(i) for i in imp) for imp in impropers
                 if imp[2] in own_set and all(i in own_set | prev_set | next_set for i in imp))
    return ResidueBlock(res_name, names, types, charges, tuple(bonds), imps, target)
