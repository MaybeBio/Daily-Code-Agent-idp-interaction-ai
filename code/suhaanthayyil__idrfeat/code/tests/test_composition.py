"""Tests for composition and hydropathy features, all sourced from localCIDER."""

from __future__ import annotations

import math

from idrfeat.composition import composition_features, mean_hydropathy


def test_aa_fractions_sum_to_one() -> None:
    feats = composition_features("ACDEFGHIKLMNPQRSTVWY")
    total = sum(v for k, v in feats.items() if k.startswith("aa_frac_"))
    assert math.isclose(total, 1.0, abs_tol=1e-9)


def test_each_residue_fraction_is_correct() -> None:
    feats = composition_features("AAAACDEF")  # 4 A out of 8
    assert math.isclose(feats["aa_frac_A"], 0.5)
    assert math.isclose(feats["aa_frac_C"], 0.125)
    assert math.isclose(feats["aa_frac_Y"], 0.0)


def test_grouped_fractions() -> None:
    # D E K R are charged; F W Y aromatic; P proline; G glycine.
    feats = composition_features("DEKRFWYPG")  # length 9
    assert math.isclose(feats["frac_charged"], 4 / 9)
    assert math.isclose(feats["frac_aromatic"], 3 / 9)
    assert math.isclose(feats["frac_proline"], 1 / 9)
    assert math.isclose(feats["frac_glycine"], 1 / 9)


def test_disorder_promoting_fraction() -> None:
    assert math.isclose(composition_features("ARGQSEKP")["frac_disorder_promoting"], 1.0)
    assert math.isclose(composition_features("ILVFWYCM")["frac_disorder_promoting"], 0.0)


def test_mean_hydropathy_localcider_scale() -> None:
    # localCIDER normalized Kyte-Doolittle scale, 0 to 9. Isoleucine is most hydrophobic (9),
    # arginine most hydrophilic (0).
    assert math.isclose(mean_hydropathy("I"), 9.0, abs_tol=1e-9)
    assert math.isclose(mean_hydropathy("R"), 0.0, abs_tol=1e-9)
    assert 0.0 < mean_hydropathy("IR") < 9.0


def test_matches_localcider() -> None:
    import pytest

    pytest.importorskip("localcider")
    from localcider.sequenceParameters import SequenceParameters

    from idrfeat.disorder import standardize_sequence

    for s in ["DEKRFWYPG", "ARGQSEKP", "MSKGEEDNMAIIKEFMRFKVHMEGSVNGHEF"]:
        sp = SequenceParameters(standardize_sequence(s))
        feats = composition_features(s)
        assert math.isclose(feats["frac_charged"], sp.get_FCR(), abs_tol=1e-9)
        assert math.isclose(
            feats["frac_disorder_promoting"], sp.get_fraction_disorder_promoting(), abs_tol=1e-9
        )
        assert math.isclose(mean_hydropathy(s), sp.get_mean_hydropathy(), abs_tol=1e-9)


def test_empty_sequence_is_safe() -> None:
    feats = composition_features("")
    assert feats["aa_frac_A"] == 0.0
    assert feats["frac_charged"] == 0.0
    assert mean_hydropathy("") == 0.0
