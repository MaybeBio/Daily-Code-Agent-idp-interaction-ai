"""Statistics added after the first run (DESIGN §7): bootstrap intervals over entries and an
alignment-free near-duplicate check on the sequences, used for collapsed statistics and for the
grouping of cross-validation folds.
"""
from __future__ import annotations

from collections.abc import Callable

import numpy as np

KMER = 5
# near duplicates: Jaccard similarity of the 5-mer sets. Two chains of equal length at 95% identity
# share about 0.95^5 = 0.77 of their 5-mers, a Jaccard of 0.63; 0.5 (about 93% identity) leaves room
# for a tag or a short truncation
NEAR_DUPLICATE = 0.5
# cross-validation groups: single linkage on containment (shared 5-mers over the smaller set), so a
# domain construct and the full chain, or two truncations of one protein, land in one group
GROUP_SIMILARITY = 0.3


def bootstrap_ci(values: np.ndarray, stat: Callable[[np.ndarray], float] = np.median, n_boot: int = 2000,
                 seed: int = 0, alpha: float = 0.05) -> tuple[float, float]:
    """Percentile bootstrap interval of `stat` over entries (resampling with replacement)."""
    v = np.asarray(values, float)
    v = v[np.isfinite(v)]
    if len(v) < 2:
        return (float("nan"), float("nan"))
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, len(v), size=(n_boot, len(v)))
    s = np.array([stat(v[i]) for i in idx])
    return (float(np.quantile(s, alpha / 2)), float(np.quantile(s, 1 - alpha / 2)))


def share(x: np.ndarray) -> float:
    return float(np.mean(x))


def kmer_set(seq: str, k: int = KMER) -> frozenset[str]:
    seq = seq.upper()
    return frozenset(seq[i:i + k] for i in range(max(0, len(seq) - k + 1)))


def containment(a: frozenset[str], b: frozenset[str]) -> float:
    """Shared k-mers over the smaller set: 1 for identical sequences and for a construct that is a
    sub-sequence of the other; robust to tags and truncations where Jaccard is not."""
    if not a or not b:
        return 0.0
    return len(a & b) / min(len(a), len(b))


def jaccard(a: frozenset[str], b: frozenset[str]) -> float:
    return len(a & b) / len(a | b) if (a or b) else 0.0


def cluster_sequences(seqs: list[str], threshold: float, k: int = KMER,
                      measure: Callable[[frozenset, frozenset], float] = containment) -> np.ndarray:
    """Single-linkage clusters of sequences whose k-mer similarity (`measure`, containment by
    default) is at least `threshold`. Returns an integer cluster id per sequence, numbered in order
    of first appearance."""
    sets = [kmer_set(s, k) for s in seqs]
    n = len(sets)
    parent = list(range(n))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for i in range(n):
        for j in range(i + 1, n):
            if measure(sets[i], sets[j]) >= threshold:
                ri, rj = find(i), find(j)
                if ri != rj:
                    parent[max(ri, rj)] = min(ri, rj)
    roots = [find(i) for i in range(n)]
    ids, out = {}, np.empty(n, int)
    for i, r in enumerate(roots):
        out[i] = ids.setdefault(r, len(ids))
    return out


def near_duplicate_pairs(labels: list[str], seqs: list[str], threshold: float = NEAR_DUPLICATE,
                         k: int = KMER) -> list[dict]:
    """Every pair whose Jaccard similarity is at least the threshold, with both measures."""
    sets = [kmer_set(s, k) for s in seqs]
    out = []
    for i in range(len(sets)):
        for j in range(i + 1, len(sets)):
            jc = jaccard(sets[i], sets[j])
            if jc >= threshold:
                out.append({"a": labels[i], "b": labels[j], "jaccard": round(jc, 3),
                            "containment": round(containment(sets[i], sets[j]), 3),
                            "identical": seqs[i].upper() == seqs[j].upper()})
    return out
