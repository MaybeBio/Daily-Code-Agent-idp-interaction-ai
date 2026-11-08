"""Loading, validating and querying the zwitterion library"""

from __future__ import annotations

import difflib
import functools
import os
import re
from collections.abc import Iterator, Mapping
from importlib import resources
from pathlib import Path
from types import MappingProxyType
from typing import Any

import yaml

from . import _chem
from ._errors import (
    DuplicateEntryError,
    InvalidStructureError,
    LibraryFileError,
    NoPolymerizableGroupError,
    UnknownBlockError,
    UnknownEntryError,
    UnknownInitiatorError,
    ZwitterionError,
)
from ._model import Block, BlockCategory, ForceField, Initiator, Zwitterion

BUILTIN_SOURCE = "<built-in>"
SCHEMA_VERSION = 1
STATUSES = frozenset({"common", "specialty", "combinatorial", "computed"})

_TOP_KEYS = frozenset({"version", "blocks", "monomers", "comonomers", "small_molecules", "initiators"})
_INITIATOR_FIELDS = frozenset(
    {"name", "aliases", "full_name", "head", "tail", "counterion", "source", "notes", "parameters"}
)
_BLOCK_FIELDS = frozenset({"name", "smiles", "smarts", "notes"})
_ENTRY_FIELDS = frozenset(
    {"name", "aliases", "full_name", "family", "status", "blocks", "smiles", "source", "notes",
     "properties", "parameters", "residue"}
)
_COMONOMER_FIELDS = _ENTRY_FIELDS | {"counterion"}
_SMARTS_REQUIRED = frozenset({BlockCategory.POLYMERIZABLE, BlockCategory.CATION, BlockCategory.ANION})
_SPACER_SHORTHAND = re.compile(r"C([1-9][0-9]*)")
_RESIDUE_NAME = re.compile(r"[A-Z0-9]{2}")


# rq-beab5e79
def normalize_name(name: str) -> str:
    """Case-insensitive key ignoring '-', '_', spaces and '/'"""
    return re.sub(r"[-_ /]", "", name).casefold()


# rq-226ae0c8
class Library:
    """An immutable, validated collection of building blocks and zwitterion entries"""

    def __init__(
        self,
        blocks: Mapping[str, Block],
        entries: Mapping[str, Zwitterion],
        initiators: Mapping[str, Initiator] | None = None,
    ) -> None:
        self._blocks = MappingProxyType(dict(blocks))
        inits = dict(initiators or {})
        self._initiator_index = MappingProxyType(inits)
        self._initiators = tuple(sorted({id(i): i for i in inits.values()}.values(), key=lambda i: i.name))
        # each entry once
        unique = {id(z): z for z in entries.values()}
        self._entries = tuple(sorted(unique.values(), key=lambda z: z.name))
        self._index = MappingProxyType(dict(entries))

    def get(self, name: str) -> Zwitterion:
        try:
            return self._index[normalize_name(name)]
        except KeyError:
            pool = {normalize_name(z.name): z.name for z in self._entries}
            close = difflib.get_close_matches(normalize_name(name), list(pool), n=3, cutoff=0.6)
            raise UnknownEntryError(name, [pool[c] for c in close]) from None

    def block(self, key: str) -> Block:
        if key in self._blocks:
            return self._blocks[key]
        match = _SPACER_SHORTHAND.fullmatch(key)
        if match:
            n = int(match.group(1))
            return Block(key, BlockCategory.SPACER, f"(CH2)n, n={n}",
                         "[*:1]" + "C" * n + "[*:2]", None, 0, False)
        raise UnknownBlockError(key)

    def blocks(self, category: BlockCategory | str | None = None) -> list[Block]:
        wanted = None if category is None else BlockCategory(category)
        return [b for b in self._blocks.values() if wanted is None or b.category is wanted]

    def entries(
        self,
        *,
        family: str | None = None,
        status: str | None = None,
        polymerizable: bool | None = None,
        cation: str | None = None,
        anion: str | None = None,
        force_field: ForceField | str | None = None,
    ) -> list[Zwitterion]:
        ff = None if force_field is None else ForceField(force_field)
        out = []
        for z in self._entries:
            if family is not None and family.casefold() not in z.family.casefold():
                continue
            if status is not None and z.status != status:
                continue
            if polymerizable is not None and z.polymerizable != polymerizable:
                continue
            if cation is not None and cation not in z.cation_classes:
                continue
            if anion is not None and anion not in z.anion_classes:
                continue
            if ff is not None and ff not in z.parameters:
                continue
            out.append(z)
        return out

    def __contains__(self, name: object) -> bool:
        return isinstance(name, str) and normalize_name(name) in self._index

    def __len__(self) -> int:
        return len(self._entries)

    def __iter__(self) -> Iterator[Zwitterion]:
        return iter(self._entries)

    # rq-c7899aba
    def initiator(self, name: str) -> Initiator:
        try:
            return self._initiator_index[normalize_name(name)]
        except KeyError:
            pool = {normalize_name(i.name): i.name for i in self._initiators}
            close = difflib.get_close_matches(normalize_name(name), list(pool), n=3, cutoff=0.6)
            raise UnknownInitiatorError(name, [pool[c] for c in close]) from None

    def initiators(self) -> list[Initiator]:
        return list(self._initiators)


# Structure validation
def _build(
    blocks: Mapping[str, Block],
    lookup_block,
    *,
    name: str,
    smiles: str | None,
    block_keys: list[str] | None,
    require_polymerizable: bool,
    vinyl_atoms: tuple[int, int] | None = None,
    zwitterionic: bool = True,
    counterion: str | None = None,
    **fields: Any,
) -> Zwitterion:
    if block_keys is not None:
        mol = _chem.assemble([lookup_block(k) for k in block_keys])
        if smiles is not None:
            given = _chem.parse_smiles(smiles)
            if _chem.canonical(given)[1] != _chem.canonical(mol)[1]:
                raise ZwitterionError(
                    f"'{name}': SMILES {smiles} differs from the structure assembled from blocks "
                    f"{block_keys}"
                )
    else:
        mol = _chem.parse_smiles(smiles)
    source_smiles = smiles or ".".join(block_keys or [])
    if zwitterionic:
        _chem.check_zwitterion(mol, source_smiles)
    else:
        _chem.check_comonomer(mol, source_smiles, counterion)

    poly_blocks = [b for b in blocks.values() if b.category is BlockCategory.POLYMERIZABLE]
    group, vinyl = _chem.find_polymerizable_group(mol, source_smiles, poly_blocks, vinyl_atoms)
    if require_polymerizable and group is None:
        raise NoPolymerizableGroupError(f"'{name}' has no recognised polymerizable group.")

    mol_c, smiles_c = _chem.canonical(mol)
    if vinyl is not None:
        order = list(map(int, mol.GetProp("_smilesAtomOutputOrder").strip("[],").split(",")))
        vinyl = (order.index(vinyl[0]), order.index(vinyl[1]))
    formula, weight = _chem.formula_and_weight(mol_c)
    return Zwitterion(
        name=name,
        smiles=smiles_c,
        blocks=tuple(block_keys) if block_keys is not None else None,
        polymerizable_group=group,
        vinyl_atoms=vinyl,
        cation_classes=_chem.matching_keys(mol_c, [b for b in blocks.values() if b.category is BlockCategory.CATION]),
        anion_classes=_chem.matching_keys(mol_c, [b for b in blocks.values() if b.category is BlockCategory.ANION]),
        formula=formula,
        molecular_weight=weight,
        zwitterionic=zwitterionic,
        counterion=counterion,
        **fields,
    )


# File loading
class _Builder:
    def __init__(self) -> None:
        self.blocks: dict[str, Block] = {}
        self.block_sources: dict[str, str] = {}
        self.entries: dict[str, Zwitterion] = {}
        self.entry_sources: dict[str, str] = {}
        self.initiators: dict[str, Initiator] = {}
        self.initiator_sources: dict[str, str] = {}

    def lookup_block(self, key: str) -> Block:
        return Library(self.blocks, {}).block(key)

    def add_file(self, path: Path, label: str) -> None:
        try:
            text = path.read_text(encoding="utf-8")
        except OSError as exc:
            raise LibraryFileError(path, f"cannot read file ({exc.strerror or exc})") from exc
        try:
            data = yaml.safe_load(text)
        except yaml.YAMLError as exc:
            raise LibraryFileError(path, f"not valid YAML ({exc})") from exc
        if not isinstance(data, dict):
            raise LibraryFileError(path, "top level must be a mapping")
        unknown = set(data) - _TOP_KEYS
        if unknown:
            raise LibraryFileError(path, f"unknown top-level key(s): {', '.join(sorted(map(str, unknown)))}")
        if data.get("version") != SCHEMA_VERSION:
            raise LibraryFileError(path, f"unsupported version {data.get('version')!r}; expected {SCHEMA_VERSION}")

        self._add_blocks(path, label, data.get("blocks") or {})
        for section, needs_poly in (("monomers", True), ("comonomers", True), ("small_molecules", False)):
            records = data.get(section) or []
            if not isinstance(records, list):
                raise LibraryFileError(path, f"'{section}' must be a list")
            for record in records:
                self._add_entry(path, label, record, needs_poly, zwitterionic=section != "comonomers")
        records = data.get("initiators") or []
        if not isinstance(records, list):
            raise LibraryFileError(path, "'initiators' must be a list")
        for record in records:
            self._add_initiator(path, label, record)

    def _add_initiator(self, path: Path, label: str, rec: Any) -> None:
        if not isinstance(rec, dict) or not isinstance(rec.get("name"), str) or not rec["name"]:
            raise LibraryFileError(path, f"initiator {rec!r} must be a mapping with a 'name'")
        name = rec["name"]
        where = f"initiator '{name}'"
        unknown = set(rec) - _INITIATOR_FIELDS
        if unknown:
            raise LibraryFileError(path, f"{where}: unknown field(s) {', '.join(sorted(map(str, unknown)))}")
        if not isinstance(rec.get("full_name"), str) or not rec["full_name"]:
            raise LibraryFileError(path, f"{where}: missing 'full_name'")
        aliases = rec.get("aliases") or []
        if not isinstance(aliases, list) or not all(isinstance(a, str) and a for a in aliases):
            raise LibraryFileError(path, f"{where}: 'aliases' must be a list of strings")
        raw_params = rec.get("parameters") or {}
        if not isinstance(raw_params, dict):
            raise LibraryFileError(path, f"{where}: 'parameters' must be a mapping")
        keys = [normalize_name(n) for n in [name, *aliases]]
        for raw, key in zip([name, *aliases], keys, strict=True):
            if key in self.initiators or keys.count(key) > 1:
                raise DuplicateEntryError(raw, [self.initiator_sources.get(key, label), label])
        try:
            parameters = {ForceField(k): str(v) for k, v in raw_params.items()}
            initiator = make_initiator(
                name, head=rec.get("head"), tail=rec.get("tail"), counterion=rec.get("counterion"),
                full_name=rec["full_name"], aliases=tuple(aliases), source=rec.get("source", "") or "",
                notes=rec.get("notes", "") or "", parameters=parameters,
            )
        except (ZwitterionError, ValueError) as exc:
            raise LibraryFileError(path, f"{where}: {exc}") from exc
        for key in keys:
            self.initiators[key] = initiator
            self.initiator_sources[key] = label

    def _add_blocks(self, path: Path, label: str, blocks: Any) -> None:
        if not isinstance(blocks, dict):
            raise LibraryFileError(path, "'blocks' must be a mapping of category to blocks")
        for cat_name, items in blocks.items():
            try:
                category = BlockCategory(cat_name)
            except ValueError:
                raise LibraryFileError(path, f"unknown block category '{cat_name}'") from None
            if not isinstance(items, dict):
                raise LibraryFileError(path, f"blocks of category '{cat_name}' must be a mapping")
            for key, rec in items.items():
                self._add_block(path, label, category, str(key), rec)

    def _add_block(self, path: Path, label: str, category: BlockCategory, key: str, rec: Any) -> None:
        where = f"block '{key}'"
        if not isinstance(rec, dict):
            raise LibraryFileError(path, f"{where} must be a mapping")
        unknown = set(rec) - _BLOCK_FIELDS
        if unknown:
            raise LibraryFileError(path, f"{where}: unknown field(s) {', '.join(sorted(unknown))}")
        if key in self.blocks:
            raise DuplicateEntryError(key, [self.block_sources[key], label])
        for req in ("name", "smiles") + (("smarts",) if category in _SMARTS_REQUIRED else ()):
            if not isinstance(rec.get(req), str) or not rec[req]:
                raise LibraryFileError(path, f"{where}: missing '{req}'")
        try:
            mol = _chem.parse_smiles(rec["smiles"])
            if rec.get("smarts"):
                _chem.parse_smarts(rec["smarts"])
        except ZwitterionError as exc:
            raise LibraryFileError(path, f"{where}: {exc}") from exc
        labels = _chem.block_attachment_labels(mol)
        expected = [2] if category is BlockCategory.POLYMERIZABLE else None
        if expected is not None and labels != expected:
            raise LibraryFileError(path, f"{where}: a polymerizable block needs exactly one [*:2]")
        if expected is None and labels not in ([1], [1, 2]):
            raise LibraryFileError(path, f"{where}: needs [*:1] and optionally [*:2] attachment points")
        self.blocks[key] = Block(
            key=key, category=category, name=rec["name"], smiles=rec["smiles"],
            smarts=rec.get("smarts"), formal_charge=_chem.net_charge(mol),
            terminal=labels == [1], notes=rec.get("notes", "") or "",
        )
        self.block_sources[key] = label

    def _add_entry(self, path: Path, label: str, rec: Any, needs_poly: bool, zwitterionic: bool = True) -> None:
        if not isinstance(rec, dict) or not isinstance(rec.get("name"), str) or not rec["name"]:
            raise LibraryFileError(path, f"entry {rec!r} must be a mapping with a 'name'")
        name = rec["name"]
        where = f"entry '{name}'"
        unknown = set(rec) - (_ENTRY_FIELDS if zwitterionic else _COMONOMER_FIELDS)
        if unknown:
            raise LibraryFileError(path, f"{where}: unknown field(s) {', '.join(sorted(map(str, unknown)))}")
        for req in ("full_name", "family", "status"):
            if not isinstance(rec.get(req), str) or not rec[req]:
                raise LibraryFileError(path, f"{where}: missing '{req}'")
        if rec["status"] not in STATUSES:
            raise LibraryFileError(path, f"{where}: status must be one of {sorted(STATUSES)}")
        if "blocks" not in rec and "smiles" not in rec:
            raise LibraryFileError(path, f"{where}: needs 'blocks' or 'smiles'")
        aliases = rec.get("aliases") or []
        if not isinstance(aliases, list) or not all(isinstance(a, str) and a for a in aliases):
            raise LibraryFileError(path, f"{where}: 'aliases' must be a list of strings")
        properties = rec.get("properties") or {}
        if not isinstance(properties, dict):
            raise LibraryFileError(path, f"{where}: 'properties' must be a mapping")
        raw_params = rec.get("parameters") or {}
        if not isinstance(raw_params, dict):
            raise LibraryFileError(path, f"{where}: 'parameters' must be a mapping")
        try:
            parameters = {ForceField(k): str(v) for k, v in raw_params.items()}
        except ValueError as exc:
            raise LibraryFileError(path, f"{where}: {exc}") from exc
        residue = rec.get("residue")
        if residue is not None and (not isinstance(residue, str) or not _RESIDUE_NAME.fullmatch(residue)):
            raise LibraryFileError(
                path, f"{where}: 'residue' must be a two-character prefix of uppercase letters or digits")
        block_keys = rec.get("blocks")
        if block_keys is not None and (not isinstance(block_keys, list) or not block_keys):
            raise LibraryFileError(path, f"{where}: 'blocks' must be a non-empty list")

        # duplicates first (cheap)
        keys = [normalize_name(n) for n in [name, *aliases]]
        for raw, key in zip([name, *aliases], keys, strict=True):
            if key in self.entries or keys.count(key) > 1:
                other = self.entry_sources.get(key, label)
                raise DuplicateEntryError(raw, [other, label])
        try:
            entry = _build(
                self.blocks, self.lookup_block, name=name, smiles=rec.get("smiles"),
                block_keys=[str(k) for k in block_keys] if block_keys is not None else None,
                require_polymerizable=needs_poly,
                zwitterionic=zwitterionic, counterion=rec.get("counterion"),
                aliases=tuple(aliases), full_name=rec["full_name"], family=rec["family"],
                status=rec["status"], source=rec.get("source", "") or "",
                notes=rec.get("notes", "") or "",
                residue=residue,
                properties=MappingProxyType(dict(properties)),
                parameters=MappingProxyType(parameters),
            )
        except ZwitterionError as exc:
            raise LibraryFileError(path, f"{where}: {exc}") from exc
        if not needs_poly and entry.polymerizable:
            raise LibraryFileError(path, f"{where}: small molecules must not be polymerizable")
        for key in keys:
            self.entries[key] = entry
            self.entry_sources[key] = label


def _builtin_path():
    return resources.files("zpoly.zwitterions").joinpath("data/zwitterion_library.yaml")


# rq-54e8d4f2
def load_library(*paths: str | os.PathLike[str], include_builtin: bool = True) -> Library:
    """Load the built-in library (optionally) merged with user YAML files, in order"""
    builder = _Builder()
    if include_builtin:
        with resources.as_file(_builtin_path()) as builtin:
            builder.add_file(Path(builtin), BUILTIN_SOURCE)
    for p in paths:
        builder.add_file(Path(p), os.fspath(p))
    return Library(builder.blocks, builder.entries, builder.initiators)


# rq-26d9d524
@functools.lru_cache(maxsize=1)
def default_library() -> Library:
    """The built-in library, loaded and validated once per process"""
    return load_library()


# rq-3c5af765 rq-458386a5
def compose(*block_keys: str, name: str | None = None, library: Library | None = None) -> Zwitterion:
    """Assemble a zwitterion from building blocks, in backbone-to-chain-end order"""
    lib = library if library is not None else default_library()
    keys = [str(k) for k in block_keys]
    return _build(
        {b.key: b for b in lib.blocks()}, lib.block, name=name or "-".join(keys), smiles=None,
        block_keys=keys, require_polymerizable=True, full_name=name or "-".join(keys),
    )


# rq-d2f618ff
def from_smiles(
    smiles: str,
    name: str,
    *,
    require_polymerizable: bool = True,
    vinyl_atoms: tuple[int, int] | None = None,
    library: Library | None = None,
) -> Zwitterion:
    """Validate an arbitrary zwitterion SMILES and identify its polymerizable group"""
    lib = library if library is not None else default_library()
    return _build(
        {b.key: b for b in lib.blocks()}, lib.block, name=name, smiles=smiles, block_keys=None,
        require_polymerizable=require_polymerizable, vinyl_atoms=vinyl_atoms, full_name=name,
    )


def make_initiator(
    name: str,
    *,
    head: str | None = None,
    tail: str | None = None,
    counterion: str | None = None,
    full_name: str | None = None,
    aliases: tuple[str, ...] = (),
    source: str = "",
    notes: str = "",
    parameters: Mapping[ForceField, str] | None = None,
) -> Initiator:
    """Validate end-group fragments and build an Initiator"""
    if head is None and tail is None:
        raise InvalidStructureError("", f"initiator '{name}' needs a head and/or tail fragment")
    charges = {}
    for side, frag in (("head", head), ("tail", tail)):
        if frag is None:
            charges[side] = 0
            continue
        if not isinstance(frag, str):
            raise InvalidStructureError(str(frag), f"{side} fragment must be a SMILES string")
        mol = _chem.parse_smiles(frag)
        if _chem.block_attachment_labels(mol) != [1]:
            raise InvalidStructureError(frag, f"{side} fragment needs exactly one [*:1] attachment point")
        charges[side] = _chem.net_charge(mol)
    if counterion is not None and not isinstance(counterion, str):
        raise InvalidStructureError(str(counterion), "counterion must be text such as 'K+'")
    return Initiator(
        name=name, head_smiles=head, tail_smiles=tail, full_name=full_name or name,
        aliases=tuple(aliases), head_charge=charges["head"], tail_charge=charges["tail"],
        counterion=counterion, source=source, notes=notes,
        parameters=MappingProxyType(dict(parameters or {})),
    )


# rq-44802ad6
def custom_initiator(
    name: str,
    *,
    head: str | None = None,
    tail: str | None = None,
    counterion: str | None = None,
    full_name: str | None = None,
) -> Initiator:
    """Define an initiator in code, validated like a library-file initiator"""
    return make_initiator(name, head=head, tail=tail, counterion=counterion, full_name=full_name)
