"""Sequence grammars of the two poles: acidic stretches and runs, glutamate in acidic runs, the CK2 consensus, RG/GR
grammar and RG/RGG clusters.

acidic_stretches, residue_runs, longest_run and rg_clusters are glu_proteomics' glu.sequence functions, moved here
unchanged; e_tract_count, ck2_motif_count and longest_acidic_run are MicroDomains' acidic_dtidr_score definitions and
rg_motif_count its methyl_microdomain_score one. Positions are 0-based half-open [start, end).
"""
from __future__ import annotations

import re

import numpy as np

from .charge import runs

E_TRACT_FORWARD = re.compile(r"(?=E[DE])")
E_TRACT_REVERSE = re.compile(r"(?=[DE]E)")
CK2_CONSENSUS = re.compile(r"(?=[ST]..[DE])")
RG_FORWARD = re.compile(r"(?=RG)")
RG_REVERSE = re.compile(r"(?=GR)")


def acidic_stretches(seq: str, window: int = 10, min_acidic: int = 6, max_basic: int = 1,
                     merge_gap: int = 3, min_length: int = 8) -> list[tuple[int, int]]:
    """Acidic stretches as 0-based [start, end) intervals.

    A window is acidic if it has >= min_acidic D/E and <= max_basic K/R. The union of acidic windows is merged across
    gaps <= merge_gap, trimmed so both ends are D/E, and kept if >= min_length.
    """
    n = len(seq)
    if n == 0:
        return []
    w = min(window, n)
    ac = np.fromiter((c in "DE" for c in seq), dtype=int, count=n)
    ba = np.fromiter((c in "KR" for c in seq), dtype=int, count=n)
    cs_a = np.concatenate([[0], np.cumsum(ac)])
    cs_b = np.concatenate([[0], np.cumsum(ba)])
    na = cs_a[w:] - cs_a[:-w]
    nb = cs_b[w:] - cs_b[:-w]
    ok = np.flatnonzero((na >= min(min_acidic, w)) & (nb <= max_basic))
    diff = np.zeros(n + 1, dtype=int)
    np.add.at(diff, ok, 1)
    np.add.at(diff, ok + w, -1)
    mark = np.cumsum(diff[:-1]) > 0
    merged: list[list[int]] = []
    for s, e in runs(mark):
        if merged and s - merged[-1][1] <= merge_gap:
            merged[-1][1] = e
        else:
            merged.append([s, e])
    out = []
    for s, e in merged:
        while s < e and seq[s] not in "DE":
            s += 1
        while e > s and seq[e - 1] not in "DE":
            e -= 1
        if e - s >= min_length:
            out.append((int(s), int(e)))
    return out


def longest_run(sub: str, chars: str) -> int:
    best = cur = 0
    for c in sub:
        cur = cur + 1 if c in chars else 0
        best = max(best, cur)
    return best


def residue_runs(seq: str, chars: str, min_run: int) -> list[tuple[int, int]]:
    """Uninterrupted runs of `chars` of at least min_run residues."""
    return [(m.start(), m.end()) for m in re.finditer(f"[{chars}]{{{min_run},}}", seq)]


def acidic_runs(seq: str, min_run: int = 6) -> list[tuple[int, int]]:
    """Uninterrupted D/E runs of at least min_run (6: the MicroDomains prospective list and PTM_Pipeline
    scripts/52_preview_acidic_tracts.py)."""
    return residue_runs(seq, "DE", min_run)


def longest_acidic_run(seq) -> int:
    """Longest uninterrupted D/E stretch."""
    return longest_run(seq, "DE") if isinstance(seq, str) else 0


def e_tract_count(seq) -> int:
    """Glutamates in an acidic run, counted in both directions (E[DE] and [DE]E, overlapping), as the acidic dtIDR
    score counts them."""
    if not isinstance(seq, str) or not seq:
        return 0
    return len(E_TRACT_FORWARD.findall(seq)) + len(E_TRACT_REVERSE.findall(seq))


def ck2_motif_count(seq) -> int:
    """S or T with an acidic residue at +3 ([ST]xx[DE], overlapping)."""
    return len(CK2_CONSENSUS.findall(seq)) if isinstance(seq, str) and seq else 0


def ck2_motif_sites(seq: str) -> list[int]:
    """0-based positions of the S/T of every [ST]xx[DE]."""
    return [m.start() for m in CK2_CONSENSUS.finditer(seq)]


def rg_motif_count(seq) -> int:
    """Overlapping RG and GR dipeptides, so RGRGR scores 4 dipeptides the way the dtIDR score counts them."""
    if not isinstance(seq, str) or not seq:
        return 0
    return len(RG_FORWARD.findall(seq)) + len(RG_REVERSE.findall(seq))


def rg_clusters(seq: str, window: int = 25, min_rg: int = 3) -> list[tuple[int, int]]:
    """Regions with >= min_rg RG dipeptides within `window` aa (RGG/RG boxes)."""
    pos = [m.start() for m in re.finditer(r"(?=RG)", seq)]
    out: list[list[int]] = []
    for i in range(len(pos)):
        j = i
        while j + 1 < len(pos) and pos[j + 1] - pos[i] <= window - 2:
            j += 1
        if j - i + 1 >= min_rg:
            s, e = pos[i], pos[j] + 2
            if out and s <= out[-1][1]:
                out[-1][1] = max(out[-1][1], e)
            else:
                out.append([s, e])
    return [(a, b) for a, b in out]
