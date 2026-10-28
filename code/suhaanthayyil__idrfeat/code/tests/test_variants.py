from __future__ import annotations

import pytest

from idrfeat.variants import apply_mutation, segment_sequence_features, variant_delta


def test_apply_mutation_replaces_one_residue() -> None:
    assert apply_mutation("AAAA", 2, "D") == "ADAA"


def test_apply_mutation_out_of_range() -> None:
    with pytest.raises(ValueError):
        apply_mutation("AAAA", 5, "D")


def test_segment_features_keys_present() -> None:
    f = segment_sequence_features("ASKGDESKGDESKGDE")
    assert "net_charge" in f and "omega" in f and "hydropathy_mean" in f


def test_variant_delta_charge_shift() -> None:
    seq = "ASKGDESKGDASKGDE"
    d = variant_delta(seq, pos=1, alt="D", start=1, end=len(seq))
    assert d["d_net_charge"] == -1.0
    assert d["d_frac_negative"] > 0.0


def test_variant_delta_position_must_be_in_segment() -> None:
    with pytest.raises(ValueError):
        variant_delta("AAAAAA", pos=1, alt="D", start=3, end=6)
