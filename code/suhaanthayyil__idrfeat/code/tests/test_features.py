"""Tests for assembling the per-segment feature table."""

from __future__ import annotations

import math

import numpy as np
import pandas as pd
import pytest

from idrfeat.disorder import metapredict_available
from idrfeat.features import TABLE_COLUMNS, Annotations, build_feature_table
from idrfeat.io import load_config

pytestmark = pytest.mark.skipif(
    not metapredict_available(), reason="metapredict is required to call IDR segments"
)

# A protein with a clear charged/polar IDR flanked by hydrophobic caps.
SEQS = {
    "EX1": "ILVFILVFIL" + "SEKPGQSEKPGQSEKPGQSEKPGQSEKPGQSEKPGQ" + "ILVFILVFIL",
}


def _cfg() -> dict:
    return load_config()


def test_table_has_documented_columns_in_order() -> None:
    df = build_feature_table(SEQS, _cfg())
    assert list(df.columns) == TABLE_COLUMNS


def test_at_least_one_segment_called() -> None:
    df = build_feature_table(SEQS, _cfg())
    assert len(df) >= 1
    row = df.iloc[0]
    assert row["accession"] == "EX1"
    assert row["disorder_backend"] == "metapredict"
    assert row["idr_start"] >= 1 and row["idr_end"] <= len(SEQS["EX1"])
    assert row["seg_length"] == row["idr_end"] - row["idr_start"] + 1


def test_sequence_features_present_and_sane() -> None:
    df = build_feature_table(SEQS, _cfg())
    row = df.iloc[0]
    aa_total = sum(row[c] for c in df.columns if c.startswith("aa_frac_"))
    assert math.isclose(aa_total, 1.0, abs_tol=1e-9)
    assert np.isfinite(row["hydropathy_mean"])
    assert np.isfinite(row["disorder_mean"])


def test_annotations_absent_are_missing_not_zero_confusion() -> None:
    df = build_feature_table(SEQS, _cfg(), annotations=None)
    row = df.iloc[0]
    # With no databases or SEG masks, those features are NaN, not fabricated.
    assert pd.isna(row["elm_motif_residues"])
    assert pd.isna(row["ptm_count"])
    assert pd.isna(row["phasepro_overlap"])
    assert pd.isna(row["lowcomplexity_frac"])
    assert pd.isna(row["plddt_mean"])


def test_annotations_present_are_computed() -> None:
    df0 = build_feature_table(SEQS, _cfg())
    start = int(df0.iloc[0]["idr_start"])
    # Put an ELM motif and PTM sites inside the called IDR.
    anns = Annotations(
        elm={"EX1": [(start, start + 4, "LIG", "LIG_TEST_1")]},
        ptm={"EX1": [(start + 1, "phospho"), (start + 2, "acetyl")]},
        phasepro={"EX1": [(start, start + 3)]},
        drllps={"EX1"},
    )
    df = build_feature_table(SEQS, _cfg(), annotations=anns)
    row = df.iloc[0]
    assert row["elm_motif_residues"] == 5
    assert row["overlaps_binding_motif"] == True  # noqa: E712
    assert row["ptm_count"] == 2
    assert row["ptm_phospho_count"] == 1
    assert row["phasepro_overlap"] == True  # noqa: E712
    assert row["drllps_member"] == True  # noqa: E712


def test_deterministic() -> None:
    a = build_feature_table(SEQS, _cfg())
    b = build_feature_table(SEQS, _cfg())
    pd.testing.assert_frame_equal(a, b)
