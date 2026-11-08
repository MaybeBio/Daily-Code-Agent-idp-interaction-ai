"""Variants of library monomers"""

# rq-6516039f

from __future__ import annotations

import re
from types import MappingProxyType

from ._errors import UnknownBlockError, VariantError
from ._library import Library, _build, default_library
from ._model import Block, BlockCategory, Zwitterion

_LINEAR = re.compile(r"C([1-9][0-9]*)")
# ammonium block -> (alkyl length, terminal?)
_ALKYL_AMMONIUM = {"NMe2": (1, False), "NEt2": (2, False), "NMe3": (1, True)}


def _alkyl_cation_key(n: int, terminal: bool) -> str:
    return f"NR3:{n}" if terminal else f"NR2:{n}"


def _alkyl_cation_smiles(n: int, terminal: bool) -> str:
    r = "C" * n
    return f"[*:1][N+]({r})({r}){r}" if terminal else f"[*:1][N+]({r})({r})[*:2]"


# rq-c9621362
def make_variant(
    entry: Zwitterion | str,
    *,
    spacer: int | None = None,
    alkyl: int | None = None,
    library: Library | None = None,
) -> Zwitterion:
    """Derive a monomer with a different cation-anion spacer and/or N-alkyl length"""
    lib = library if library is not None else default_library()
    z = lib.get(entry) if isinstance(entry, str) else entry
    for option, value in (("spacer", spacer), ("alkyl", alkyl)):
        if value is not None and (not isinstance(value, int) or isinstance(value, bool) or value < 1):
            raise ValueError(f"{option} must be an integer >= 1, got {value!r}")
    if spacer is None and alkyl is None:
        return z
    option = "spacer" if spacer is not None else "alkyl"
    if not z.zwitterionic:
        raise VariantError(z.name, option, "it is a comonomer, not a zwitterion")
    if z.blocks is None:
        raise VariantError(z.name, option, "the entry has no building-block decomposition")

    keys = list(z.blocks)
    categories = []
    for k in keys:
        try:
            categories.append(lib.block(k).category)
        except UnknownBlockError:
            categories.append(None)
    charged = [i for i, c in enumerate(categories) if c in (BlockCategory.CATION, BlockCategory.ANION)]
    suffix = []
    extra_blocks = {}

    if spacer is not None:
        if len(charged) < 2:
            raise VariantError(z.name, "spacer", "it has no separate cation and anion blocks")
        between = [i for i in range(charged[0] + 1, charged[1]) if categories[i] is BlockCategory.SPACER]
        if len(between) != 1:
            raise VariantError(z.name, "spacer", "there is no single spacer block between the charges")
        idx = between[0]
        match = _LINEAR.fullmatch(keys[idx])
        if not match:
            raise VariantError(z.name, "spacer", f"its spacer '{keys[idx]}' is not a plain (CH2)n chain")
        if int(match.group(1)) != spacer:
            keys[idx] = f"C{spacer}"
            suffix.append(f"sp{spacer}")

    if alkyl is not None:
        idx = next((i for i in charged if keys[i] in _ALKYL_AMMONIUM), None)
        if idx is None:
            raise VariantError(z.name, "alkyl", "its cation is not a methyl/ethyl ammonium (NMe2, NEt2, NMe3)")
        current, terminal = _ALKYL_AMMONIUM[keys[idx]]
        if current != alkyl:
            key = _alkyl_cation_key(alkyl, terminal)
            extra_blocks[key] = _alkyl_cation_smiles(alkyl, terminal)
            keys[idx] = key
            suffix.append(f"alk{alkyl}")

    if not suffix:
        return z

    def lookup(key: str):
        if key in extra_blocks:
            smi = extra_blocks[key]
            return Block(key, BlockCategory.CATION, f"N-alkyl ammonium (n={key.split(':')[1]})", smi,
                         None, 1, "[*:2]" not in smi)
        return lib.block(key)

    name = "-".join([z.name, *suffix])
    return _build(
        {b.key: b for b in lib.blocks()}, lookup, name=name, smiles=None, block_keys=keys,
        require_polymerizable=z.polymerizable, full_name=f"{z.full_name} ({', '.join(suffix)})",
        family=z.family, status="user", notes=f"Variant of {z.name}: {', '.join(suffix)}.",
        properties=MappingProxyType({}), parameters=MappingProxyType({}),
    )
