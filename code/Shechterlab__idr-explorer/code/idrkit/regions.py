"""Intervals: merging, and the disordered regions long enough to call an IDR.

Positions here are 1-based inclusive [start, end], as MobiDB and UniProt give them.
"""
from __future__ import annotations


def merge_intervals(spans) -> list[list[int]]:
    """Overlapping 1-based inclusive intervals merged (PTM_Pipeline scripts/36c merge())."""
    out: list[list[int]] = []
    for a, b in sorted(spans):
        if out and a <= out[-1][1]:
            out[-1][1] = max(out[-1][1], b)
        else:
            out.append([a, b])
    return out


def idr_intervals(spans, length: int, min_len: int = 20) -> list[tuple[int, int]]:
    """IDRs: disordered intervals clipped to the sequence, merged where they overlap, of at least min_len residues."""
    clipped = [(max(1, int(a)), min(length, int(b))) for a, b in spans if int(b) >= 1 and int(a) <= length]
    return [(a, b) for a, b in merge_intervals(clipped) if b - a + 1 >= min_len]
