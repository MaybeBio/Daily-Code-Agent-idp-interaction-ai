"""OPLS-AA backend"""

# rq-f6d54673

from __future__ import annotations

import functools
import logging
import os
import shutil
from collections.abc import Sequence
from dataclasses import dataclass
from importlib import resources
from pathlib import Path

import yaml
from rdkit import Chem, rdBase

from ._errors import ForceFieldError, ForceFieldUnavailableError, TypingError
from ._model import AtomType, FragmentParameters, ParameterSet

log = logging.getLogger(__name__)
OPLS_RTP_BONDEDTYPES = (1, 1, 3, 4, 1, 3, 1, 0)
LOCATE_HINT = "Set OPLSAA_FF = PATH to your oplsaa.ff folder (or set ZPOLY_OPLSAA_FF), or install GROMACS."
REQUIRED_FILES = ("forcefield.itp", "ffnonbonded.itp", "ffbonded.itp")
_PTYPES = {"A", "D", "S", "V"}


# Reading a GROMACS force-field folder
@dataclass(frozen=True)
class GromacsForceField:
    path: Path
    defaults: tuple[int, int, str, float, float]
    atomtypes: dict  # name -> (bond_type, atomic_number, mass, charge, sigma, epsilon)
    bondtypes: dict  # (bt, bt) -> (b0, kb)
    angletypes: dict  # (bt, bt, bt) -> (theta0, k)
    rb_dihedraltypes: dict  # (bt, bt, bt, bt) -> (C0..C5)
    defines: dict  # name -> (phase, k, multiplicity)
    descriptions: dict  # name -> atomtypes.atp description ({} when the file is absent)


def _lines(path: Path, defined: set[str] | None = None):
    """(line number, tokens, section) for each active data line"""
    defined = set() if defined is None else defined
    active: list[bool] = []  # one entry per open #if block: is that branch being read?
    section = None
    try:
        text = path.read_text(encoding="utf-8", errors="replace")
    except OSError as exc:
        raise ForceFieldError(f"Cannot read {path}: {exc}") from exc
    for n, raw in enumerate(text.splitlines(), 1):
        line = raw.split(";", 1)[0].strip()
        if not line:
            continue
        if line.startswith("#"):
            tok = line.split()
            word = tok[0]
            reading = all(active)
            if word in ("#ifdef", "#ifndef"):
                name = tok[1] if len(tok) > 1 else ""
                active.append((name in defined) == (word == "#ifdef"))
            elif word == "#else" and active:
                active[-1] = not active[-1]
            elif word == "#endif" and active:
                active.pop()
            elif word == "#define" and reading:
                if len(tok) > 1:
                    defined.add(tok[1])
                yield n, tok, "#define"
            continue  # #include and anything else
        if not all(active):
            continue
        if line.startswith("[") and line.endswith("]"):
            section = line[1:-1].strip()
            continue
        yield n, line.split(), section


def _parse(path: Path, handler, defined: set[str] | None = None) -> None:
    for n, tok, section in _lines(path, defined):
        try:
            handler(tok, section)
        except (ValueError, IndexError) as exc:
            raise ForceFieldError(f"{path}:{n}: cannot parse line ({exc})") from None


# rq-fb24e397 rq-9a642c3a
def read_gromacs_forcefield(path: str | os.PathLike) -> GromacsForceField:
    """Parse an oplsaa.ff folder."""
    folder = Path(path)
    for name in REQUIRED_FILES:
        if not (folder / name).is_file():
            raise ForceFieldError(f"{folder} is not a force-field folder: {name} is missing")
    defaults: list = []
    atomtypes: dict = {}
    bondtypes: dict = {}
    angletypes: dict = {}
    rb: dict = {}
    defines: dict = {}

    def on_ff(tok, section):
        if section == "defaults" and not defaults:
            defaults.extend([int(tok[0]), int(tok[1]), tok[2], float(tok[3]), float(tok[4])])

    def on_nonbonded(tok, section):
        if section != "atomtypes":
            return
        p = next(i for i, t in enumerate(tok) if t in _PTYPES and i >= 4)
        name = tok[0]
        if p >= 5:  # name bond_type at.num mass charge ptype sigma epsilon (ptype is column 6)
            bt, atnum = tok[1], int(tok[p - 3])
        else:  # name at.num mass charge ptype sigma epsilon (ptype is column 5)
            bt, atnum = name, int(tok[1])
        atomtypes[name] = (bt, atnum, float(tok[p - 2]), float(tok[p - 1]), float(tok[p + 1]), float(tok[p + 2]))

    def on_bonded(tok, section):
        if section == "#define":
            if len(tok) >= 5:
                try:
                    defines[tok[1]] = (float(tok[-3]), float(tok[-2]), int(float(tok[-1])))
                except ValueError:
                    pass  # non-numeric macro: not an improper parameter
        elif section == "bondtypes" and tok[2] == "1":
            bondtypes[(tok[0], tok[1])] = (float(tok[3]), float(tok[4]))
        elif section == "angletypes" and tok[3] == "1":
            angletypes[(tok[0], tok[1], tok[2])] = (float(tok[4]), float(tok[5]))
        elif section == "dihedraltypes":
            if len(tok) >= 11 and tok[4] == "3":
                rb[(tok[0], tok[1], tok[2], tok[3])] = tuple(float(x) for x in tok[5:11])
            elif len(tok) >= 9 and tok[2] == "3":  # old 2-type format: j k funct C0..C5
                rb[("X", tok[0], tok[1], "X")] = tuple(float(x) for x in tok[3:9])

    defined: set[str] = set()  # macros #defined so far, shared like one preprocessed file
    _parse(folder / "forcefield.itp", on_ff, defined)
    _parse(folder / "ffnonbonded.itp", on_nonbonded, defined)
    _parse(folder / "ffbonded.itp", on_bonded, defined)
    if not defaults:
        raise ForceFieldError(f"{folder / 'forcefield.itp'} has no [ defaults ] line")
    descriptions: dict = {}
    atp = folder / "atomtypes.atp"
    if atp.is_file():
        for raw in atp.read_text(encoding="utf-8", errors="replace").splitlines():
            data, _, comment = raw.partition(";")
            tok = data.split()
            if tok:
                descriptions[tok[0]] = comment.strip()
    return GromacsForceField(folder, tuple(defaults), atomtypes, bondtypes, angletypes, rb, defines,
                             descriptions)


@functools.lru_cache(maxsize=4)
def _cached_ff(path: str, stamp: tuple) -> GromacsForceField:
    return read_gromacs_forcefield(path)


def _load_ff(folder: Path) -> GromacsForceField:
    stamp = tuple((folder / n).stat().st_mtime for n in REQUIRED_FILES)
    return _cached_ff(str(folder), stamp)


# Typing rules
@dataclass(frozen=True)
class TypingRule:
    name: str
    smarts: str
    pattern: Chem.Mol
    type: str
    expect: tuple[str, ...] = ()


@dataclass(frozen=True)
class ImproperRule:
    name: str
    smarts: str
    pattern: Chem.Mol
    define: str


@dataclass(frozen=True)
class TypingRules:
    rules: tuple[TypingRule, ...]
    impropers: tuple[ImproperRule, ...]

    def __add__(self, other: TypingRules) -> TypingRules:
        return TypingRules(self.rules + other.rules, self.impropers + other.impropers)


def _compile(smarts, where: str, source: str) -> Chem.Mol:
    if not isinstance(smarts, str) or not smarts:
        raise TypingError(f"{source}: {where} needs a 'smarts' pattern", where)
    with rdBase.BlockLogs():
        patt = Chem.MolFromSmarts(smarts)
    if patt is None or patt.GetNumAtoms() == 0:
        raise TypingError(f"{source}: {where} has an invalid SMARTS {smarts!r}", where)
    return patt


# rq-779e8915 rq-87cd66f0
def load_typing_rules(path: str | os.PathLike | None = None) -> TypingRules:
    """The starter rules (path=None) or a user rule file"""
    if path is None:
        source = "starter rules"
        text = resources.files("zpoly.forcefield").joinpath("data/oplsaa_rules.yaml").read_text()
    else:
        source = str(path)
        try:
            text = Path(path).read_text(encoding="utf-8")
        except OSError as exc:
            raise TypingError(f"{source}: cannot read the rule file ({exc.strerror or exc})") from exc
    try:
        data = yaml.safe_load(text)
    except yaml.YAMLError as exc:
        raise TypingError(f"{source}: not valid YAML ({exc})") from exc
    if not isinstance(data, dict) or data.get("version") != 1:
        raise TypingError(f"{source}: a rule file needs 'version: 1'")
    unknown = set(data) - {"version", "rules", "impropers"}
    if unknown:
        raise TypingError(f"{source}: unknown key(s) {', '.join(sorted(map(str, unknown)))}")
    rules, impropers = [], []
    for kind, out, field_name in (("rules", rules, "type"), ("impropers", impropers, "define")):
        entries = data.get(kind) or []
        if not isinstance(entries, list):
            raise TypingError(f"{source}: '{kind}' must be a list")
        for k, rec in enumerate(entries, 1):
            if not isinstance(rec, dict):
                raise TypingError(f"{source}: {kind} entry {k} must be a mapping")
            where = str(rec.get("name") or f"{kind} entry {k}")
            allowed = {"name", "smarts", field_name} | ({"expect"} if kind == "rules" else set())
            extra = set(rec) - allowed
            if extra:
                raise TypingError(f"{source}: {where} has unknown key(s) {', '.join(sorted(map(str, extra)))}", where)
            value = rec.get(field_name)
            if not isinstance(value, str) or not value:
                raise TypingError(f"{source}: {where} needs a '{field_name}'", where)
            patt = _compile(rec.get("smarts"), where, source)
            if kind == "rules":
                expect = rec.get("expect", [])
                if not isinstance(expect, list) or not all(isinstance(e, str) and e for e in expect):
                    raise TypingError(f"{source}: {where}: 'expect' must be a list of keywords", where)
                out.append(TypingRule(where, rec["smarts"], patt, value, tuple(expect)))
            else:
                out.append(ImproperRule(where, rec["smarts"], patt, value))
    return TypingRules(tuple(rules), tuple(impropers))


def assign_types(mol: Chem.Mol, rules: TypingRules) -> list[str | None]:
    """First matching rule types each atom (the rule pattern's first atom)"""
    types: list[str | None] = [None] * mol.GetNumAtoms()
    for rule in rules.rules:
        for match in mol.GetSubstructMatches(rule.pattern, uniquify=False, maxMatches=100000):
            if types[match[0]] is None:
                types[match[0]] = rule.type
        if all(t is not None for t in types):
            break
    return types


def _describe(mol: Chem.Mol, i: int) -> str:
    atom = mol.GetAtomWithIdx(i)
    nbrs = ", ".join(f"{n.GetSymbol()}{n.GetIdx() + 1}" for n in atom.GetNeighbors())
    charge = f" (charge {atom.GetFormalCharge():+d})" if atom.GetFormalCharge() else ""
    return f"atom {i + 1} {atom.GetSymbol()}{charge} bonded to {nbrs or 'nothing'}"


_warned_no_atp: set[str] = set()


def _check_expectations(ff: GromacsForceField, rules: TypingRules) -> None:
    """Check rule types against atomtypes.atp."""
    checked = [r for r in rules.rules if r.expect]
    if not checked:
        return
    if not ff.descriptions:
        if str(ff.path) not in _warned_no_atp:
            _warned_no_atp.add(str(ff.path))
            log.warning("%s has no atomtypes.atp; the 'expect' checks of %d typing rules were skipped.",
                        ff.path, len(checked))
        return
    for rule in checked:
        desc = ff.descriptions.get(rule.type, "")
        if not any(k.casefold() in desc.casefold() for k in rule.expect):
            raise TypingError(
                f"Rule '{rule.name}' expects {rule.type} to be one of {list(rule.expect)}, but "
                f"{ff.path / 'atomtypes.atp'} describes it as {desc!r}. Override the rule with "
                "OPLS_RULES if your force field numbers this group differently.", rule.name)


# Backend
# rq-a47637a7 rq-040346db rq-57822d06 rq-32f43f89
class OplsaaBackend:
    """OPLS-AA via SMARTS typing on an oplsaa.ff folder."""

    name = "oplsaa"

    def __init__(self, oplsaa_ff: str | os.PathLike | None = None, rules: str | os.PathLike | None = None) -> None:
        self.oplsaa_ff = oplsaa_ff
        self.rules_path = rules

    def locate(self) -> Path:
        candidates: list[Path] = []
        if self.oplsaa_ff is not None:
            candidates = [Path(self.oplsaa_ff)]
        else:
            if os.environ.get("ZPOLY_OPLSAA_FF"):
                candidates.append(Path(os.environ["ZPOLY_OPLSAA_FF"]))
            if os.environ.get("GMXDATA"):
                candidates.append(Path(os.environ["GMXDATA"]) / "top" / "oplsaa.ff")
            for exe in ("gmx", "gmx_mpi"):
                found = shutil.which(exe)
                if found:
                    candidates.append(Path(found).resolve().parent.parent / "share" / "gromacs" / "top" / "oplsaa.ff")
        for folder in candidates:
            if all((folder / n).is_file() for n in REQUIRED_FILES):
                return folder
        where = ", ".join(str(c) for c in candidates) or "no candidate locations"
        raise ForceFieldUnavailableError(self.name, f"no oplsaa.ff folder found ({where})", LOCATE_HINT)

    def _setup(self) -> tuple[GromacsForceField, TypingRules]:
        ff = _load_ff(self.locate())
        rules = load_typing_rules()
        if self.rules_path is not None:
            rules = load_typing_rules(self.rules_path) + rules  # user rules first
        for rule in rules.rules:
            if rule.type not in ff.atomtypes:
                raise TypingError(f"Rule '{rule.name}' uses {rule.type}, which is not in "
                                  f"{ff.path / 'ffnonbonded.itp'}", rule.name)
        _check_expectations(ff, rules)
        for imp in rules.impropers:
            if imp.define not in ff.defines:
                raise TypingError(f"Improper rule '{imp.name}' uses {imp.define}, which is not defined in "
                                  f"{ff.path / 'ffbonded.itp'}", imp.name)
        return ff, rules

    def check_available(self) -> None:
        self._setup()

    def full_parameters(self) -> ParameterSet:
        """All bonded parameters of the oplsaa.ff folder."""
        ff = _load_ff(self.locate())
        return ParameterSet(ff.defaults, OPLS_RTP_BONDEDTYPES, {}, dict(ff.bondtypes), dict(ff.angletypes),
                            dict(ff.rb_dihedraltypes), {}, dihedral_funct=3)

    def parametrize(self, mol: Chem.Mol, net_charge: int, work_dir: Path, *, charges: bool = True
                    ) -> FragmentParameters:
        ff, rules = self._setup()
        types = assign_types(mol, rules)
        untyped = [i for i, t in enumerate(types) if t is None]
        if untyped:
            atoms = [_describe(mol, i) for i in untyped]
            smiles = Chem.MolToSmiles(Chem.RemoveHs(mol, sanitize=False))
            raise TypingError(
                f"No OPLS-AA typing rule matches {len(untyped)} atom(s) of fragment {smiles}:\n  "
                + "\n  ".join(atoms[:20]) + ("\n  ..." if len(atoms) > 20 else "")
                + "\nAdd rules for them in a rule file (OPLS_RULES = FILE).", None, atoms)
        typed: Sequence[str] = [t for t in types if t is not None]
        charges = tuple(ff.atomtypes[t][3] for t in typed)

        present = sorted(set(typed))
        atomtypes = {}
        for t in present:
            bt, atnum, mass, _q, sigma, eps = ff.atomtypes[t]
            atomtypes[t] = AtomType(t, atnum, mass, sigma, eps, bt)
        classes = {at.bonded_key for at in atomtypes.values()}

        def keep(key) -> bool:
            return all(k in classes or k == "X" for k in key)

        impropers: dict = {}
        for imp in rules.impropers:
            for match in mol.GetSubstructMatches(imp.pattern, uniquify=False):
                key = ("X", "X", atomtypes[typed[match[0]]].bonded_key, "X")
                impropers.setdefault(key, ff.defines[imp.define])
        params = ParameterSet(
            ff.defaults, OPLS_RTP_BONDEDTYPES, atomtypes,
            {k: v for k, v in ff.bondtypes.items() if keep(k)},
            {k: v for k, v in ff.angletypes.items() if keep(k)},
            {k: v for k, v in ff.rb_dihedraltypes.items() if keep(k)},
            impropers, dihedral_funct=3,
        )
        return FragmentParameters(tuple(typed), charges, params)
