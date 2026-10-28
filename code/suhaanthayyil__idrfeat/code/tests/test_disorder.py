"""Tests for disorder segment calling."""

from __future__ import annotations

import numpy as np
import pytest

from idrfeat.disorder import (
    default_backend,
    disordered,
    primary_disorder,
    segments_from_scores,
    standardize_sequence,
)


def test_segments_are_one_based_inclusive() -> None:
    scores = np.array([0.0, 0.6, 0.6, 0.6, 0.0])
    assert segments_from_scores(scores, threshold=0.5, min_len=3) == [(2, 4)]
    assert segments_from_scores(scores, threshold=0.5, min_len=4) == []


def test_segments_at_sequence_end() -> None:
    scores = np.array([0.0, 0.9, 0.9, 0.9])
    assert segments_from_scores(scores, threshold=0.5, min_len=2) == [(2, 4)]


def test_disordered_mask() -> None:
    scores = np.array([0.2, 0.8, 0.5])
    assert list(disordered(scores, threshold=0.5)) == [False, True, True]


def test_default_backend_is_metapredict() -> None:
    assert default_backend() == "metapredict"


def test_unknown_backend_raises() -> None:
    with pytest.raises(ValueError):
        primary_disorder("SEKPGSEKPG", backend="heuristic")


def test_standardize_preserves_length_and_maps_nonstandard() -> None:
    seq = "MUOXBZJ"
    out = standardize_sequence(seq)
    assert len(out) == len(seq)
    assert set(out) <= set("ACDEFGHIKLMNPQRSTVWY")
    assert out == "MCKADEL"


def test_metapredict_handles_selenocysteine() -> None:
    from idrfeat.disorder import metapredict_available, metapredict_scores

    if not metapredict_available():
        return
    scores = metapredict_scores("MKTAYIAKQRUQISFVKSHFSRQLEERLGLIEVQ")
    assert len(scores) == 34
