from __future__ import annotations

import math

import numpy as np
import pandas as pd

from idrfeat.labels import (
    aggregate_segment_labels,
    check_coordinate_concordance,
    idr_segments_from_residues,
    load_residue_table,
)


def _residues() -> pd.DataFrame:
    seq = "MASKDEFG"
    idr = [0, 0, 1, 1, 1, 1, 0, 0]
    peptide = [0, 0, 2, 2, 3, 1, 0, 0]
    depleted = [0, 0, 1, 1, 0, 0, 0, 0]
    lfc = [np.nan, np.nan, -2.0, -1.5, -0.2, 0.1, np.nan, np.nan]
    rows = []
    for i in range(8):
        rows.append(
            {
                "accession": "P1",
                "aa_loc": i + 1,
                "aa": seq[i],
                "idr": idr[i],
                "depleted": depleted[i],
                "lfc": lfc[i],
                "pvalue": 0.01,
                "peptide_count": peptide[i],
            }
        )
    return pd.DataFrame(rows)


def test_load_residue_table_renames_schema(tmp_path) -> None:
    raw = pd.DataFrame(
        {
            "Uniprot_ID": ["P1", "P1"],
            "AA_loc": [1, 2],
            "AA": ["M", "A"],
            "IDR": [0, 1],
            "Fitness_depleted": [0, 1],
            "LFC": [np.nan, -1.0],
            "P_value": [np.nan, 0.01],
            "Peptide_count": [0, 3],
        }
    )
    path = tmp_path / "sample.csv"
    raw.to_csv(path, index=False)
    df = load_residue_table(path)
    assert list(df.columns) == [
        "accession",
        "aa_loc",
        "aa",
        "idr",
        "depleted",
        "lfc",
        "pvalue",
        "peptide_count",
    ]
    assert df["peptide_count"].iloc[0] == 0


def test_coordinate_concordance_match_and_mismatch() -> None:
    df = _residues()
    assert check_coordinate_concordance(df, {"P1": "MASKDEFG"})["match_frac"] == 1.0
    bad = check_coordinate_concordance(df, {"P1": "MASKDEFX"})
    assert bad["match_frac"] < 1.0
    assert "P1" in bad["mismatched_accessions"]


def test_idr_segments_from_residue_flags() -> None:
    seg = idr_segments_from_residues(_residues())
    assert list(seg[["accession", "idr_start", "idr_end", "seg_length"]].iloc[0]) == ["P1", 3, 6, 4]
    assert len(seg) == 1


def test_aggregate_segment_labels() -> None:
    df = _residues()
    seg = idr_segments_from_residues(df)
    out = aggregate_segment_labels(df, seg, min_covered=3, frac_cutoff=0.5).iloc[0]
    assert out["n_covered"] == 4
    assert math.isclose(out["frac_depleted"], 0.5)
    assert out["fitness_label"] == 1.0


def test_aggregate_below_min_coverage_is_nan() -> None:
    df = _residues()
    seg = idr_segments_from_residues(df)
    out = aggregate_segment_labels(df, seg, min_covered=10).iloc[0]
    assert math.isnan(out["fitness_label"])


def test_load_residue_table_tab_delimited_txt(tmp_path) -> None:
    raw = pd.DataFrame(
        {
            "Uniprot_ID": ["P1", "P1"],
            "AA_loc": [1, 2],
            "AA": ["M", "A"],
            "IDR": [0, 1],
            "Fitness_depleted": [0, 1],
            "LFC": [np.nan, -1.0],
            "P_value": [np.nan, 0.01],
            "Peptide_count": [0, 3],
        }
    )
    path = tmp_path / "screen.txt"
    raw.to_csv(path, sep="\t", index=False)
    df = load_residue_table(path)
    assert len(df) == 2
    assert df["accession"].tolist() == ["P1", "P1"]
    assert df["idr"].tolist() == [0, 1]
