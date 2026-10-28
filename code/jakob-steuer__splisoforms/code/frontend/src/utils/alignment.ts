/**
 * Centralized alignment utilities for Splisoforms.
 * Ensures 3D structure viewers align exactly with the Sequence Alignment (MSA) panel.
 */

export interface AlignedPair {
  canonAln: string;
  queryAln: string;
}

// ─────────────────────────────────────────────────────────────────────────────
// Splice-aware alignment  (primary — used for MSA panel & structure superposition)
// ─────────────────────────────────────────────────────────────────────────────

/**
 * Splice-aware alignment: aligns a novel isoform to the canonical protein sequence
 * using exact k-mer anchors (shared exon segments).
 *
 * Motivation
 * ----------
 * Global/semi-global Needleman-Wunsch treats every sequence difference as a
 * substitution.  For splice isoforms this produces misleading "mutations" (amber)
 * wherever a novel exon or alternative C-terminus is force-aligned against an
 * unrelated canonical region.
 *
 * Instead this function:
 *  1. Builds a k-mer index on seq1 (canonical).
 *  2. Finds all exact 5-mer seed hits in seq2 (isoform).
 *  3. Extends each seed to its maximal exact-match block.
 *  4. Deduplicates and chains the blocks in canonical order
 *     (longest-weight increasing subsequence).
 *  5. Builds the final alignment:
 *       • anchor regions          → exact match (shown as identity in MSA)
 *       • canonical-only gaps     → deletions   (skipped exons, shown as pink dashes)
 *       • isoform-only insertions → novel seq   (novel exons, shown in green)
 *
 * This means only genuine shared exon sequence is ever aligned as "identity",
 * and novel/alternative exon sequence appears as green insertions rather than
 * amber mutations.
 *
 * Falls back to semi-global NW when sequences are shorter than k (rare edge case).
 */
export function spliceAwareAlign(seq1: string, seq2: string): AlignedPair {
  seq1 = seq1.trim().replace(/\s/g, "");
  seq2 = seq2.trim().replace(/\s/g, "");

  const n = seq1.length;
  const m = seq2.length;

  if (n === 0) return { canonAln: "", queryAln: seq2 };
  if (m === 0) return { canonAln: seq1, queryAln: "-".repeat(n) };

  const K = 5; // k-mer seed size

  // Fall back to NW for very short sequences
  if (n < K || m < K) return needlemanWunsch(seq1, seq2);

  // ─── 1. K-mer index on canonical ─────────────────────────────────────────
  const canonIndex = new Map<string, number[]>();
  for (let i = 0; i <= n - K; i++) {
    const kmer = seq1.substring(i, i + K);
    const list = canonIndex.get(kmer);
    if (list) list.push(i);
    else canonIndex.set(kmer, [i]);
  }

  // ─── 2. Collect seeds ─────────────────────────────────────────────────────
  type Anchor = { c0: number; i0: number; len: number };
  const rawSeeds: Anchor[] = [];

  for (let j = 0; j <= m - K; j++) {
    const kmer = seq2.substring(j, j + K);
    const hits = canonIndex.get(kmer);
    if (!hits) continue;
    for (const c of hits) {
      rawSeeds.push({ c0: c, i0: j, len: K });
    }
  }

  // No shared k-mers → isoform is entirely novel (different exon set).
  // Show canonical as "skipped" (pink dashes in isoform row) then the isoform
  // as a novel green insertion after the canonical ends.
  if (rawSeeds.length === 0) {
    return {
      canonAln: seq1 + "-".repeat(m),
      queryAln: "-".repeat(n) + seq2,
    };
  }

  // ─── 3. Extend seeds to maximal exact-match blocks ────────────────────────
  // All seeds within the same exact-match block will extend to the same
  // (c0, i0, len), so dedup-by-(c0,i0) collapses them to one entry per block.
  const extended: Anchor[] = rawSeeds.map(({ c0, i0 }) => {
    // Extend right
    let len = K;
    while (c0 + len < n && i0 + len < m && seq1[c0 + len] === seq2[i0 + len]) len++;
    // Extend left
    let leftExt = 0;
    while (
      c0 - leftExt > 0 &&
      i0 - leftExt > 0 &&
      seq1[c0 - leftExt - 1] === seq2[i0 - leftExt - 1]
    ) leftExt++;
    return { c0: c0 - leftExt, i0: i0 - leftExt, len: len + leftExt };
  });

  // ─── 4. Deduplicate by (c0, i0), then sort by canonical position ──────────
  const seen = new Set<string>();
  const unique: Anchor[] = [];
  for (const a of extended) {
    const key = `${a.c0},${a.i0}`;
    if (!seen.has(key)) {
      seen.add(key);
      unique.push(a);
    }
  }
  unique.sort((a, b) => a.c0 - b.c0 || a.i0 - b.i0);

  // ─── 5. Chain: longest-weight non-crossing subsequence ───────────────────
  // Weight = total aligned residues (longer anchors are more reliable).
  //
  // Soft isoform-overlap tolerance (MAX_ISO_OVERLAP):
  //   Greedy seed extension can extend an anchor a few residues past an exon
  //   boundary when the junction sequence happens to be shared by two exons
  //   (e.g. "QA" at the tail of one exon AND the head of the next).  The
  //   canonical-space positions never cross, but the isoform positions can
  //   overlap by a few residues.  We allow up to MAX_ISO_OVERLAP residues of
  //   soft isoform overlap so both anchors can still be chained; the overlap
  //   is resolved during alignment building by advancing the effective start
  //   of the later anchor past the collision.
  const MAX_ISO_OVERLAP = 5;

  const dp = unique.map(a => a.len);   // best total weight ending at anchor i
  const prev = unique.map(() => -1);

  for (let i = 1; i < unique.length; i++) {
    for (let j = 0; j < i; j++) {
      // Canonical space: strictly non-overlapping (no tolerance needed)
      if (unique[j].c0 + unique[j].len > unique[i].c0) continue;
      // Isoform space: allow small overlaps caused by junction ambiguity
      if (unique[j].i0 + unique[j].len > unique[i].i0 + MAX_ISO_OVERLAP) continue;

      const candidate = dp[j] + unique[i].len;
      if (candidate > dp[i]) {
        dp[i] = candidate;
        prev[i] = j;
      }
    }
  }

  // Best chain end
  let bestIdx = 0;
  for (let i = 1; i < dp.length; i++) {
    if (dp[i] > dp[bestIdx]) bestIdx = i;
  }

  // Reconstruct chain (in canonical order)
  const chain: Anchor[] = [];
  let cur = bestIdx;
  while (cur !== -1) {
    chain.unshift(unique[cur]);
    cur = prev[cur];
  }

  // ─── 6. Build alignment from chain ───────────────────────────────────────
  let canonAln = "";
  let queryAln = "";
  let cPos = 0;
  let iPos = 0;

  for (const anchor of chain) {
    // Resolve any soft overlap from the previous anchor by advancing this
    // anchor's effective start past the already-consumed positions.
    // isoOverlap > 0 when the greedy extension of the previous anchor
    // "consumed" a few residues that belong to this anchor's exon.
    const isoOverlap  = Math.max(0, iPos - anchor.i0);
    const canonOverlap = Math.max(0, cPos - anchor.c0);
    const skip = Math.max(isoOverlap, canonOverlap);

    const effC0  = anchor.c0 + skip;
    const effI0  = anchor.i0 + skip;
    const effLen = anchor.len - skip;

    if (effLen <= 0) {
      // Entire anchor already consumed — skip it
      cPos = Math.max(cPos, anchor.c0 + anchor.len);
      iPos = Math.max(iPos, anchor.i0 + anchor.len);
      continue;
    }

    // Gap between previous position and this anchor.
    //
    // When BOTH sequences have residues to fill the gap, we need to decide
    // whether this represents:
    //   (a) A point mutation / short substitution (SNP) → align as amber mismatches
    //   (b) Two different exons sitting in the gap    → canonical gap = pink (skipped),
    //                                                   isoform gap   = green (novel)
    //
    // Heuristic: if both sides of the gap are ≤ SNP_WIDTH residues, treat as
    // a substitution (amber). Longer gaps almost certainly represent different
    // exon sequences and must NOT be force-aligned as mismatches — that would
    // produce long amber columns for unrelated sequences and misclassify the
    // biological event (splice, not point mutation).
    const SNP_WIDTH = 3;
    const cGap = effC0 - cPos;
    const iGap = effI0 - iPos;
    const sharedMismatch =
      cGap > 0 && iGap > 0 && Math.min(cGap, iGap) <= SNP_WIDTH
        ? Math.min(cGap, iGap)
        : 0;

    if (sharedMismatch > 0) {
      // Short substitution block (likely SNP) — shown amber by the MSA renderer
      canonAln += seq1.substring(cPos, cPos + sharedMismatch);
      queryAln  += seq2.substring(iPos,  iPos  + sharedMismatch);
    }

    // Remaining canonical-only residues → deletions (skipped exon, pink dashes)
    if (cGap > sharedMismatch) {
      const del = seq1.substring(cPos + sharedMismatch, effC0);
      canonAln += del;
      queryAln += "-".repeat(del.length);
    }

    // Remaining isoform-only residues → novel insertions (novel exon, green)
    if (iGap > sharedMismatch) {
      const ins = seq2.substring(iPos + sharedMismatch, effI0);
      canonAln += "-".repeat(ins.length);
      queryAln += ins;
    }

    // Anchor (exact match, trimmed for any overlap)
    const matchSeq = seq1.substring(effC0, effC0 + effLen);
    canonAln += matchSeq;
    queryAln += matchSeq;

    cPos = effC0 + effLen;
    iPos = effI0 + effLen;
  }

  // Trailing canonical residues (isoform ends before canonical)
  if (cPos < n) {
    canonAln += seq1.substring(cPos);
    queryAln += "-".repeat(n - cPos);
  }

  // Trailing isoform residues (novel C-terminus / alternative last exon)
  if (iPos < m) {
    canonAln += "-".repeat(m - iPos);
    queryAln += seq2.substring(iPos);
  }

  return { canonAln, queryAln };
}

// ─────────────────────────────────────────────────────────────────────────────
// Needleman-Wunsch  (kept for pairwise panels: DomainsPanel, PtmAccessibilityPanel)
// ─────────────────────────────────────────────────────────────────────────────

/**
 * Semi-global Needleman-Wunsch alignment: free end-gaps on seq1 (canonical).
 * Used for pairwise comparisons that need continuous coordinate mapping,
 * e.g. the compare page.  DomainsPanel and PtmAccessibilityPanel now use
 * spliceAwareAlign instead to avoid fragmented coverage bins.
 */
export function needlemanWunsch(seq1: string, seq2: string): AlignedPair {
  // Defensive cleaning: ensure no hidden whitespace shifts indices
  seq1 = seq1.trim().replace(/\s/g, "");
  seq2 = seq2.trim().replace(/\s/g, "");

  const match = 2;
  const mismatch = -1;
  const gap = -3;

  const n = seq1.length;
  const m = seq2.length;

  if (n === 0 || m === 0) return { canonAln: seq1, queryAln: seq2 };

  const score: number[][] = Array.from({ length: n + 1 }, () => new Array(m + 1).fill(0));

  // Initialize: seq2 must be fully consumed, but seq1 (canonical) can have free leading gaps
  for (let j = 0; j <= m; j++) score[0][j] = j * gap;

  for (let i = 1; i <= n; i++) {
    for (let j = 1; j <= m; j++) {
      const s = seq1[i - 1] === seq2[j - 1] ? match : mismatch;
      score[i][j] = Math.max(
        score[i - 1][j - 1] + s,
        score[i - 1][j] + gap,
        score[i][j - 1] + gap
      );
    }
  }

  // Find best score in last column (free trailing gaps on seq1)
  let bestI = n;
  let maxScore = score[n][m];
  for (let i = 0; i <= n; i++) {
    if (score[i][m] > maxScore) {
      maxScore = score[i][m];
      bestI = i;
    }
  }

  // Traceback from (bestI, m)
  let res1 = "", res2 = "";
  let i = bestI, j = m;

  // Trailing canonical residues
  for (let k = n; k > bestI; k--) {
    res1 = seq1[k - 1] + res1;
    res2 = "-" + res2;
  }

  while (i > 0 && j > 0) {
    const s = seq1[i - 1] === seq2[j - 1] ? match : mismatch;
    if (score[i][j] === score[i - 1][j - 1] + s) {
      res1 = seq1[i - 1] + res1;
      res2 = seq2[j - 1] + res2;
      i--; j--;
    } else if (score[i][j] === score[i - 1][j] + gap) {
      res1 = seq1[i - 1] + res1;
      res2 = "-" + res2;
      i--;
    } else {
      res1 = "-" + res1;
      res2 = seq2[j - 1] + res2;
      j--;
    }
  }

  // Remaining seq2 residues
  while (j > 0) {
    res1 = "-" + res1;
    res2 = seq2[j - 1] + res2;
    j--;
  }

  // Leading canonical residues
  while (i > 0) {
    res1 = seq1[i - 1] + res1;
    res2 = "-" + res2;
    i--;
  }

  return { canonAln: res1, queryAln: res2 };
}

// ─────────────────────────────────────────────────────────────────────────────
// Structure superposition indices
// ─────────────────────────────────────────────────────────────────────────────

/**
 * Returns 1-based residue indices for Molstar superposition.
 * Uses splice-aware alignment so only genuinely shared exon residues
 * (exact k-mer anchored matches) are used as superposition anchors.
 * Novel exon sequences are never included, preventing bogus superpositions.
 */
export function getAlignedIndices(
  canonicalSeq: string,
  querySeq: string
): { idxCanon: number[]; idxQuery: number[] } {
  const { canonAln, queryAln } = spliceAwareAlign(canonicalSeq, querySeq);

  const idxCanon: number[] = [];
  const idxQuery: number[] = [];

  let cPos = 1; // 1-based
  let qPos = 1;

  for (let i = 0; i < canonAln.length; i++) {
    const c = canonAln[i];
    const q = queryAln[i];

    if (c !== "-" && q !== "-") {
      if (c === q) {
        // Exact match in a shared-exon anchor → use as superposition anchor
        idxCanon.push(cPos);
        idxQuery.push(qPos);
      }
      cPos++;
      qPos++;
    } else if (c !== "-") {
      cPos++;
    } else if (q !== "-") {
      qPos++;
    }
  }

  return { idxCanon, idxQuery };
}
