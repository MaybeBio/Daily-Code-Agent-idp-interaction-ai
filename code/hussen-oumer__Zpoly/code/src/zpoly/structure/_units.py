"""Monomer or sequence -> repeat units."""

# rq-35c1237e

from __future__ import annotations

from ..zwitterions import Library, Zwitterion, make_variant
from ._sequence import homopolymer


def variant(unit, spacer: int | None, alkyl: int | None, lib: Library) -> Zwitterion:
    if spacer is None and alkyl is None:
        return lib.get(unit) if isinstance(unit, str) else unit
    return make_variant(unit, spacer=spacer, alkyl=alkyl, library=lib)


def resolve_units(monomer, n: int | None, spacer: int | None, alkyl: int | None,
                  lib: Library) -> tuple[list[Zwitterion], dict[str, str]]:
    """Returns (units in chain order, {name as given"""
    if isinstance(monomer, (str, Zwitterion)):
        if n is None:
            raise ValueError("Give the number of units n for a homopolymer.")
        name = monomer if isinstance(monomer, str) else monomer.name
        sequence = [monomer] * len(homopolymer(name, n))  # homopolymer() checks n
    else:
        sequence = list(monomer)
        if n is not None and n != len(sequence):
            raise ValueError(f"n={n} does not match the sequence length {len(sequence)}.")
    cache: dict[str, Zwitterion] = {}
    for item in sequence:
        key = item if isinstance(item, str) else item.name
        if key not in cache:
            cache[key] = variant(item, spacer, alkyl, lib)
    units = [cache[item if isinstance(item, str) else item.name] for item in sequence]
    return units, {key: z.name for key, z in cache.items()}
