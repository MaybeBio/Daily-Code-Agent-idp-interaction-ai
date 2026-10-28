from __future__ import annotations

import math

import numpy as np

from idrfeat.esm import (
    attention_position_features,
    segment_attention_features,
    segment_embedding_mean,
)


def _uniform_att(layers: int = 2, heads: int = 2, length: int = 3) -> np.ndarray:
    return np.full((layers, heads, length, length), 1.0 / length)


def test_attention_features_uniform_rows() -> None:
    att = _uniform_att()
    feats = attention_position_features(att)
    assert feats["received"].shape == (3,)
    assert feats["entropy"].shape == (3,)
    assert np.allclose(feats["received"], 1.0)
    assert np.allclose(feats["entropy"], math.log(3))


def test_segment_attention_features_keys() -> None:
    out = segment_attention_features(_uniform_att(), 1, 3)
    assert math.isclose(out["attn_received_mean"], 1.0)
    assert math.isclose(out["attn_entropy_mean"], math.log(3))


def test_segment_embedding_mean() -> None:
    emb = np.array([[0.0, 0.0], [2.0, 4.0], [4.0, 8.0]])
    out = segment_embedding_mean(emb, 1, 3)
    assert np.allclose(out, [2.0, 4.0])
