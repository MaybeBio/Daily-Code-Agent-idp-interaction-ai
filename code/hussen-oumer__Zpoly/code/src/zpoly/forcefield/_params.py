"""User bonded parameters"""

# rq-5a64d41e

from __future__ import annotations

import math
import os
from collections.abc import Mapping
from dataclasses import replace
from pathlib import Path

import yaml

from ._errors import ParameterFileError
from ._model import ParameterSet

SECTIONS = {  # section -> (kind, number of types, value keys)
    "bonds": ("bond", 2, ("b0", "kb")),
    "angles": ("angle", 3, ("theta0", "k")),
    "dihedrals": ("dihedral", 4, None),  # rb or terms
    "impropers": ("improper", 4, ("phase", "k", "mult")),
}
TABLES = {"bond": "bondtypes", "angle": "angletypes", "dihedral": "dihedraltypes", "improper": "impropertypes"}


def _num(value, where, path) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ParameterFileError(f"{where}: expected a number, got {value!r}", path, where.split(":")[0])
    return float(value)


def _types(value, n, where, path, key="types") -> tuple[str, ...]:
    if not isinstance(value, list) or len(value) != n or not all(isinstance(t, str) and t for t in value):
        raise ParameterFileError(f"{where}: '{key}' must be a list of {n} bond types", path, where)
    return tuple(value)


def check_parameter_mapping(data, path: str | None = None) -> dict:
    """Validate a parameter-file mapping"""
    if not isinstance(data, Mapping) or data.get("version") != 1:
        raise ParameterFileError("a parameter file needs 'version: 1'", path)
    unknown = set(data) - set(SECTIONS) - {"version"}
    if unknown:
        raise ParameterFileError(f"unknown section(s) {', '.join(sorted(map(str, unknown)))}", path)
    out: dict = {}
    for section, (kind, n, value_keys) in SECTIONS.items():
        entries = data.get(section) or []
        if not isinstance(entries, list):
            raise ParameterFileError(f"'{section}' must be a list", path)
        for k, rec in enumerate(entries, 1):
            where = f"{section}[{k}]"
            if not isinstance(rec, Mapping):
                raise ParameterFileError(f"{where}: each entry must be a mapping", path, where)
            keys = value_keys or ("rb", "terms")
            extra = set(rec) - {"types", "same_as"} - set(keys)
            if extra:
                raise ParameterFileError(f"{where}: unknown key(s) {', '.join(sorted(map(str, extra)))}", path, where)
            types = _types(rec.get("types"), n, where, path)
            given = [key for key in keys if key in rec]
            if "same_as" in rec:
                if given:
                    raise ParameterFileError(f"{where}: give either values or same_as, not both", path, where)
                entry = {"types": types, "same_as": _types(rec["same_as"], n, where, path, "same_as"), "where": where}
            elif value_keys is None:  # dihedral: rb or terms
                if len(given) != 1:
                    raise ParameterFileError(f"{where}: give 'rb' (6 coefficients) or 'terms' "
                                             "([[phase, k, mult], ...])", path, where)
                if given[0] == "rb":
                    rb = rec["rb"]
                    if not isinstance(rb, list) or len(rb) != 6:
                        raise ParameterFileError(f"{where}: 'rb' needs 6 coefficients", path, where)
                    entry = {"types": types, "style": 3,
                             "value": tuple(_num(x, where, path) for x in rb), "where": where}
                else:
                    terms = rec["terms"]
                    if not isinstance(terms, list) or not terms or not all(
                            isinstance(t, list) and len(t) == 3 for t in terms):
                        raise ParameterFileError(f"{where}: 'terms' must be [[phase, k, mult], ...]", path, where)
                    value = tuple((_num(t[0], where, path), _num(t[1], where, path), int(_num(t[2], where, path)))
                                  for t in terms)
                    entry = {"types": types, "style": 9, "value": value, "where": where}
            else:
                missing = [key for key in value_keys if key not in rec]
                if missing:
                    raise ParameterFileError(f"{where}: missing {', '.join(missing)} (or use same_as)", path, where)
                vals = [_num(rec[key], where, path) for key in value_keys]
                if kind == "improper":
                    vals[2] = int(vals[2])
                entry = {"types": types, "value": tuple(vals), "where": where}
            out.setdefault(kind, []).append(entry)
    return out


# rq-5a64d41e
def load_parameter_file(path: str | os.PathLike) -> dict:
    p = str(path)
    try:
        data = yaml.safe_load(Path(path).read_text(encoding="utf-8"))
    except OSError as exc:
        raise ParameterFileError(f"cannot read the parameter file ({exc.strerror or exc})", p) from exc
    except yaml.YAMLError as exc:
        raise ParameterFileError(f"not valid YAML ({exc})", p) from exc
    return check_parameter_mapping(data, p)


def apply_user_parameters(params: ParameterSet, user: Mapping, path: str | None = None,
                          lookup: ParameterSet | None = None):
    """Apply checked user entries"""
    from ._core import match_exact, match_wildcard

    tables = {kind: dict(getattr(params, attr)) for kind, attr in TABLES.items()}
    used = []
    for kind, entries in user.items():
        table = tables[kind]
        matcher = match_exact if kind in ("bond", "angle") else match_wildcard
        for e in entries:
            key = e["types"]
            if "same_as" in e:
                value = matcher(getattr(params, TABLES[kind]), e["same_as"])
                if value is None and lookup is not None:
                    value = matcher(getattr(lookup, TABLES[kind]), e["same_as"])
                if value is None:
                    raise ParameterFileError(f"{e['where']}: same_as {'-'.join(e['same_as'])} has no {kind} "
                                             "parameter in this force field", path, e["where"])
            else:
                if kind == "dihedral" and e["style"] != params.dihedral_funct:
                    want = "rb (6 RB coefficients)" if params.dihedral_funct == 3 else "terms ([[phase, k, mult]])"
                    raise ParameterFileError(f"{e['where']}: this force field's dihedrals need {want}",
                                             path, e["where"])
                value = e["value"]
            table.pop(tuple(reversed(key)), None)
            table[key] = value
            used.append((kind, key))
    new = replace(params, **{TABLES[k]: v for k, v in tables.items()})
    return new, tuple(used)


def parameter_template(missing: Mapping[str, list], dihedral_funct: int) -> str:
    """Ready-to-edit parameter-file text for the missing combinations"""
    lines = ["version: 1"]
    fields = {"bond": "b0: 0.0, kb: 0.0", "angle": "theta0: 0.0, k: 0.0",
              "dihedral": "rb: [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]" if dihedral_funct == 3
              else "terms: [[0.0, 0.0, 1]]",
              "improper": "phase: 180.0, k: 0.0, mult: 2"}
    for kind, section in (("bond", "bonds"), ("angle", "angles"), ("dihedral", "dihedrals")):
        if missing.get(kind):
            lines.append(f"{section}:")
            lines += [f"  - {{types: [{', '.join(t)}], {fields[kind]}}}" for t in missing[kind]]
    return "\n".join(lines) + "\n"


def _ionic_atoms(mol) -> set[int]:
    """Charged atoms plus the heavy atoms within two bonds of them"""
    ionic = {a.GetIdx() for a in mol.GetAtoms() if a.GetFormalCharge() != 0}
    frontier = set(ionic)
    for _ in range(2):
        frontier = {n.GetIdx() for i in frontier for n in mol.GetAtomWithIdx(i).GetNeighbors()
                    if n.GetAtomicNum() > 1} - ionic
        ionic |= frontier
    return ionic


def _fmt(kind: str, value, funct: int) -> str:
    if value is None:
        zeros = {"bond": "b0: 0.0, kb: 0.0", "angle": "theta0: 0.0, k: 0.0",
                 "dihedral": "rb: [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]" if funct == 3 else "terms: [[0.0, 0.0, 1]]"}
        return zeros[kind]
    if kind == "bond":
        return f"b0: {value[0]:.6f}, kb: {value[1]:.4f}"
    if kind == "angle":
        return f"theta0: {value[0]:.4f}, k: {value[1]:.6f}"
    if funct == 3:
        return "rb: [" + ", ".join(f"{x:.5f}" for x in value) + "]"
    return "terms: [" + ", ".join(f"[{p:.4f}, {k:.6f}, {m}]" for p, k, m in value) + "]"


# rq-702d4a58
def ionic_parameter_template(topology) -> str:
    """Parameter file for the ionic parts."""
    from ._core import bonded_classes, match_exact, match_wildcard

    mol = topology.structure.mol
    params = topology.parameters
    classes = bonded_classes(topology.atom_types, params)
    ionic = _ionic_atoms(mol)

    def label(i: int) -> str:
        a = mol.GetAtomWithIdx(i)
        return f"{a.GetProp('res')}:{a.GetProp('name')}"

    selections = (
        ("bond", "bonds", [t for t in topology.bonds if ionic.intersection(t)], params.bondtypes, match_exact),
        ("angle", "angles", [t for t in topology.angles if ionic.intersection(t)], params.angletypes, match_exact),
        ("dihedral", "dihedrals", [t for t in topology.propers if ionic.intersection(t)],
         params.dihedraltypes, match_wildcard),
    )
    lines = [
        "# Bonded parameters of the ionic parts (charged atoms and their bonded heavy atoms).",
        "# Values are the force field's current ones; edit them (e.g. with DFT-derived values) and pass",
        "# this file with PARAMS = <file>. Lines marked 'missing' must be filled in. GROMACS units: nm, deg, kJ/mol.",
        "version: 1",
    ]
    for kind, section, terms, table, matcher in selections:
        seen: dict[tuple, tuple[int, ...]] = {}
        for t in terms:
            key = tuple(classes[i] for i in t)
            canon = min(key, tuple(reversed(key)))
            if canon not in seen:
                seen[canon] = t if key == canon else tuple(reversed(t))
        if not seen:
            continue
        lines.append(f"{section}:")
        for key in sorted(seen):
            value = matcher(table, key)
            example = "-".join(label(i) for i in seen[key])
            note = "missing in the force field" if value is None else f"e.g. {example}"
            if value is None:
                note += f" (e.g. {example})"
            lines.append(f"  - {{types: [{', '.join(key)}], {_fmt(kind, value, params.dihedral_funct)}}}   # {note}")
    return "\n".join(lines) + "\n"
