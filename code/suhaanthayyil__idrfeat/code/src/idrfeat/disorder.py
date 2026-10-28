"""Per-residue disorder annotation and IDR segment calling.

metapredict is the primary predictor and AIUPred the secondary one; AIUPred also supplies the
fold-on-binding propensity. Both are established published tools. metapredict is required, since
IDR segments are defined by a real disorder predictor rather than a local stand-in. AIUPred
models, when used, load once into a process-level singleton.
"""

from __future__ import annotations

import numpy as np

from .constants import AA

_predictor = None

_STANDARD = set(AA)
_TRANS = {ord(k): v for k, v in {"U": "C", "O": "K", "B": "D", "Z": "E", "J": "L", "X": "A"}.items()}


def standardize_sequence(seq: str) -> str:
    s = seq.upper().translate(_TRANS)
    if all(c in _STANDARD for c in s):
        return s
    return "".join(c if c in _STANDARD else "A" for c in s)


def segments_from_scores(
    scores: np.ndarray, threshold: float, min_len: int
) -> list[tuple[int, int]]:
    """Contiguous runs of scores >= threshold at least min_len long, 1-based inclusive."""
    above = np.asarray(scores) >= threshold
    segments = []
    start = None
    for i, hit in enumerate(above):
        if hit and start is None:
            start = i
        elif not hit and start is not None:
            if i - start >= min_len:
                segments.append((start + 1, i))
            start = None
    if start is not None and len(above) - start >= min_len:
        segments.append((start + 1, len(above)))
    return segments


def disordered(scores: np.ndarray, threshold: float) -> np.ndarray:
    """Boolean mask of residues with score >= threshold."""
    return np.asarray(scores) >= threshold


def metapredict_scores(seq: str) -> np.ndarray:
    """Per-residue metapredict disorder scores in [0, 1] for one sequence."""
    import metapredict as meta

    return np.asarray(meta.predict_disorder(standardize_sequence(seq), return_numpy=True), dtype=float)


def _get_predictor():
    global _predictor
    if _predictor is None:
        from aiupred import AIUPred

        _predictor = AIUPred(force_cpu=True)
    return _predictor


def aiupred_disorder_scores(seq: str) -> np.ndarray:
    """AIUPred per-residue disorder propensity in [0, 1]."""
    return np.asarray(_get_predictor().predict_disorder(standardize_sequence(seq)), dtype=float)


def aiupred_binding_scores(seq: str) -> np.ndarray:
    """AIUPred per-residue binding (fold-on-binding) propensity in [0, 1]."""
    return np.asarray(_get_predictor().predict_binding(standardize_sequence(seq)), dtype=float)


def metapredict_available() -> bool:
    import importlib.util

    return importlib.util.find_spec("metapredict") is not None


def aiupred_available() -> bool:
    import importlib.util

    return importlib.util.find_spec("aiupred") is not None


def default_backend() -> str:
    """Primary disorder backend that defines IDR segments."""
    return "metapredict"


def primary_disorder(seq: str, backend: str | None = None) -> tuple[np.ndarray, str]:
    """Per-residue disorder scores from the chosen backend, with the backend name."""
    backend = backend or default_backend()
    if backend == "metapredict":
        return metapredict_scores(seq), backend
    if backend == "aiupred":
        return aiupred_disorder_scores(seq), backend
    raise ValueError(f"unknown disorder backend: {backend}")
