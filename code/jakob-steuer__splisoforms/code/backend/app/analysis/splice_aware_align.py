"""
Splice-aware pairwise alignment for canonical-vs-novel isoform sequences.

Port of the frontend's spliceAwareAlign() (frontend/src/utils/alignment.ts).
Anchors on exact k-mer matches between canonical and novel; chains them by
longest-weight non-crossing subsequence; the resulting alignment looks like:

  * anchor regions          → exact match  (identity)
  * canonical-only gaps     → deletions    (skipped exons)
  * isoform-only insertions → novel seq    (novel exons)
  * short symmetric mismatches (≤ SNP_WIDTH) → substitution block

This is the right tool for AS isoforms because global Needleman-Wunsch
force-aligns biologically unrelated exon sequences as mismatches, producing
random PTM-position mappings.  The k-mer anchor approach lets entirely-novel
exons stay novel (no spurious alignment), so canonical PTMs in those regions
are correctly classified as `lost`.

The Python implementation is line-for-line equivalent to the TS source so
backend PTM mapping and frontend MSA display always agree.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

K = 5                    # k-mer seed size (mirrors frontend default)
MAX_ISO_OVERLAP = 5      # soft isoform-overlap tolerance during chaining
SNP_WIDTH = 3            # short bilateral gap → treat as substitution


@dataclass
class _Anchor:
    c0: int
    i0: int
    len: int


def _needleman_wunsch_fallback(seq1: str, seq2: str) -> Tuple[str, str]:
    """Semi-global NW with free end-gaps on seq1.  Used only when sequences
    are shorter than K (rare).  Mirrors the JS fallback path."""
    match, mismatch, gap = 2, -1, -3
    n, m = len(seq1), len(seq2)
    if n == 0 or m == 0:
        return seq1, seq2

    score = [[0] * (m + 1) for _ in range(n + 1)]
    for j in range(m + 1):
        score[0][j] = j * gap
    # No leading-gap penalty on seq1 (semi-global):
    # row 0 left as zeros for column 0; init column 0 of subsequent rows with gap
    for i in range(1, n + 1):
        score[i][0] = 0  # free leading gaps on canonical
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            s = match if seq1[i - 1] == seq2[j - 1] else mismatch
            score[i][j] = max(
                score[i - 1][j - 1] + s,
                score[i - 1][j] + gap,
                score[i][j - 1] + gap,
            )

    # Traceback from the best cell in the last column
    i = max(range(n + 1), key=lambda x: score[x][m])
    j = m
    a1: list[str] = []
    a2: list[str] = []
    # Trailing canonical residues if any
    while i < n:
        a1.append(seq1[i])
        a2.append("-")
        i += 1
    i = max(range(n + 1), key=lambda x: score[x][m])
    while i > 0 and j > 0:
        s = match if seq1[i - 1] == seq2[j - 1] else mismatch
        if score[i][j] == score[i - 1][j - 1] + s:
            a1.append(seq1[i - 1]); a2.append(seq2[j - 1]); i -= 1; j -= 1
        elif score[i][j] == score[i - 1][j] + gap:
            a1.append(seq1[i - 1]); a2.append("-"); i -= 1
        else:
            a1.append("-"); a2.append(seq2[j - 1]); j -= 1
    while j > 0:
        a1.append("-"); a2.append(seq2[j - 1]); j -= 1
    while i > 0:
        a1.append(seq1[i - 1]); a2.append("-"); i -= 1
    return "".join(reversed(a1)), "".join(reversed(a2))


def splice_aware_align(seq1: str, seq2: str) -> Tuple[str, str]:
    """Return (canon_aln, query_aln) — aligned strings of equal length."""
    seq1 = seq1.strip().replace(" ", "")
    seq2 = seq2.strip().replace(" ", "")
    n, m = len(seq1), len(seq2)

    if n == 0:
        return "", seq2
    if m == 0:
        return seq1, "-" * n
    if n < K or m < K:
        return _needleman_wunsch_fallback(seq1, seq2)

    # ── 1. k-mer index on canonical ──────────────────────────────────────────
    canon_index: Dict[str, List[int]] = {}
    for i in range(n - K + 1):
        kmer = seq1[i:i + K]
        canon_index.setdefault(kmer, []).append(i)

    # ── 2. Collect seeds ─────────────────────────────────────────────────────
    raw_seeds: List[_Anchor] = []
    for j in range(m - K + 1):
        kmer = seq2[j:j + K]
        for c in canon_index.get(kmer, ()):
            raw_seeds.append(_Anchor(c0=c, i0=j, len=K))

    if not raw_seeds:
        # No shared k-mers at all — isoform comes from an entirely different
        # exon set.  Render the whole canonical as deletion, the whole isoform
        # as insertion after it.
        return seq1 + "-" * m, "-" * n + seq2

    # ── 3. Extend seeds to maximal exact-match blocks ────────────────────────
    extended: List[_Anchor] = []
    for s in raw_seeds:
        c0, i0 = s.c0, s.i0
        # Extend right
        length = K
        while c0 + length < n and i0 + length < m and seq1[c0 + length] == seq2[i0 + length]:
            length += 1
        # Extend left
        left = 0
        while c0 - left > 0 and i0 - left > 0 and seq1[c0 - left - 1] == seq2[i0 - left - 1]:
            left += 1
        extended.append(_Anchor(c0=c0 - left, i0=i0 - left, len=length + left))

    # ── 4. Dedup by (c0, i0); sort by canonical position then isoform ────────
    seen: set[Tuple[int, int]] = set()
    unique: List[_Anchor] = []
    for a in extended:
        key = (a.c0, a.i0)
        if key not in seen:
            seen.add(key)
            unique.append(a)
    unique.sort(key=lambda a: (a.c0, a.i0))

    # ── 5. Longest-weight non-crossing chain (dp) ───────────────────────────
    if not unique:
        return seq1 + "-" * m, "-" * n + seq2

    dp = [a.len for a in unique]
    prev = [-1] * len(unique)
    for i in range(1, len(unique)):
        for j in range(i):
            # Canonical strict non-overlap
            if unique[j].c0 + unique[j].len > unique[i].c0:
                continue
            # Isoform allows MAX_ISO_OVERLAP soft overlap
            if unique[j].i0 + unique[j].len > unique[i].i0 + MAX_ISO_OVERLAP:
                continue
            candidate = dp[j] + unique[i].len
            if candidate > dp[i]:
                dp[i] = candidate
                prev[i] = j

    best = max(range(len(dp)), key=lambda x: dp[x])

    # Reconstruct chain
    chain: List[_Anchor] = []
    cur = best
    while cur != -1:
        chain.insert(0, unique[cur])
        cur = prev[cur]

    # ── 6. Build alignment from chain ────────────────────────────────────────
    canon_aln: list[str] = []
    query_aln: list[str] = []
    c_pos = 0
    i_pos = 0

    for anchor in chain:
        iso_overlap   = max(0, i_pos - anchor.i0)
        canon_overlap = max(0, c_pos - anchor.c0)
        skip = max(iso_overlap, canon_overlap)

        eff_c0  = anchor.c0 + skip
        eff_i0  = anchor.i0 + skip
        eff_len = anchor.len - skip

        if eff_len <= 0:
            c_pos = max(c_pos, anchor.c0 + anchor.len)
            i_pos = max(i_pos, anchor.i0 + anchor.len)
            continue

        c_gap = eff_c0 - c_pos
        i_gap = eff_i0 - i_pos
        shared_mismatch = (min(c_gap, i_gap)
                           if c_gap > 0 and i_gap > 0 and min(c_gap, i_gap) <= SNP_WIDTH
                           else 0)

        if shared_mismatch > 0:
            canon_aln.append(seq1[c_pos:c_pos + shared_mismatch])
            query_aln.append(seq2[i_pos:i_pos + shared_mismatch])

        if c_gap > shared_mismatch:
            d = seq1[c_pos + shared_mismatch: eff_c0]
            canon_aln.append(d)
            query_aln.append("-" * len(d))

        if i_gap > shared_mismatch:
            ins = seq2[i_pos + shared_mismatch: eff_i0]
            canon_aln.append("-" * len(ins))
            query_aln.append(ins)

        match_seq = seq1[eff_c0: eff_c0 + eff_len]
        canon_aln.append(match_seq)
        query_aln.append(match_seq)

        c_pos = eff_c0 + eff_len
        i_pos = eff_i0 + eff_len

    if c_pos < n:
        canon_aln.append(seq1[c_pos:])
        query_aln.append("-" * (n - c_pos))
    if i_pos < m:
        canon_aln.append("-" * (m - i_pos))
        query_aln.append(seq2[i_pos:])

    return "".join(canon_aln), "".join(query_aln)


def map_canonical_to_novel(canonical_seq: str, novel_seq: str) -> Dict[int, Optional[int]]:
    """
    Build a {1-based canonical pos → 1-based novel pos or None} map from a
    splice-aware alignment.  None means the canonical residue is deleted in
    the novel (skipped exon).  Used by ptm_accessibility.py to project known
    canonical PTM sites onto the novel isoform's coordinate space.
    """
    canon_aln, query_aln = splice_aware_align(canonical_seq, novel_seq)
    position_map: Dict[int, Optional[int]] = {}
    c_pos = 0
    n_pos = 0
    for c_char, n_char in zip(canon_aln, query_aln):
        c_is_gap = (c_char == "-")
        n_is_gap = (n_char == "-")
        if not c_is_gap and not n_is_gap:
            position_map[c_pos + 1] = n_pos + 1
            c_pos += 1
            n_pos += 1
        elif not c_is_gap and n_is_gap:
            position_map[c_pos + 1] = None
            c_pos += 1
        elif c_is_gap and not n_is_gap:
            n_pos += 1
    return position_map
