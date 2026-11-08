"""Sequence helpers"""

# rq-b08ea041

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence

import numpy as np

from ._errors import SequenceError


def _check_n(n: int) -> None:
    if not isinstance(n, int) or isinstance(n, bool) or n < 1:
        raise SequenceError(f"The number of units must be an integer >= 1, got {n!r}.")


# rq-d34a145b
def homopolymer(monomer: str, n: int) -> list[str]:
    """n copies of one monomer"""
    _check_n(n)
    return [monomer] * n


# rq-9c67798c
def alternating(a: str, b: str, n: int) -> list[str]:
    """n units alternating a, b, a, ... (an odd n has one more a)."""
    _check_n(n)
    return [a if i % 2 == 0 else b for i in range(n)]


# rq-88a61c9d
def block(blocks: Sequence[tuple[str, int]]) -> list[str]:
    """Blocks in order, e.g. [("DMAPS", 10), ("CBMA", 20)]."""
    if not blocks:
        raise SequenceError("A block copolymer needs at least one block.")
    out: list[str] = []
    for item in blocks:
        try:
            monomer, count = item
        except (TypeError, ValueError):
            raise SequenceError(f"Each block must be (monomer, count), got {item!r}.") from None
        if not isinstance(count, int) or isinstance(count, bool) or count < 1:
            raise SequenceError(f"Block '{monomer}' needs a count >= 1, got {count!r}.")
        out.extend([monomer] * count)
    return out


# rq-ff223b41
def random_copolymer(composition: Mapping[str, float], n: int, seed: int = 0) -> list[str]:
    """Statistical copolymer"""
    _check_n(n)
    if not composition:
        raise SequenceError("A random copolymer needs at least one monomer.")
    names = list(composition)
    fractions = np.array([composition[k] for k in names], dtype=float)
    if not np.all(np.isfinite(fractions)) or np.any(fractions <= 0):
        raise SequenceError(f"All fractions must be positive, got {dict(composition)}.")
    if not math.isclose(fractions.sum(), 1.0, abs_tol=1e-6):
        raise SequenceError(f"Fractions must sum to 1, got {fractions.sum():.6g}.")
    exact = fractions * n
    counts = np.floor(exact).astype(int)
    remainder = n - counts.sum()
    # largest remainder
    order = sorted(range(len(names)), key=lambda i: (-(exact[i] - counts[i]), i))
    for i in order[:remainder]:
        counts[i] += 1
    units = [name for name, c in zip(names, counts, strict=True) for _ in range(c)]
    np.random.default_rng(seed).shuffle(units)
    return units


# rq-6652f857
def pattern(pattern: str, mapping: Mapping[str, str], n: int) -> list[str]:
    """Repeat a symbol pattern such as "AAB" and cut it at n units"""
    _check_n(n)
    if not pattern:
        raise SequenceError("The pattern is empty.")
    missing = sorted(set(pattern) - set(mapping))
    if missing:
        raise SequenceError(f"Pattern symbol(s) {', '.join(missing)} have no monomer in the mapping.")
    return [mapping[pattern[i % len(pattern)]] for i in range(n)]
