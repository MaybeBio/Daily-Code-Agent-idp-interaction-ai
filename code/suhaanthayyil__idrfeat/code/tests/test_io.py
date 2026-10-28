"""Tests for FASTA reading, config defaults, table writing, and the bundled example."""

from __future__ import annotations

import pandas as pd

from idrfeat.io import (
    example_fasta_path,
    load_config,
    read_fasta,
    write_table,
)


def test_read_fasta_uniprot_and_plain_headers(tmp_path) -> None:
    path = tmp_path / "p.fasta"
    path.write_text(">sp|P12345|NAME_HUMAN desc\nMKAEILV\nFASD\n>plainid foo\nGGGG\n")
    seqs = read_fasta(path)
    assert seqs["P12345"] == "MKAEILVFASD"
    assert seqs["plainid"] == "GGGG"


def test_load_config_defaults_and_override(tmp_path) -> None:
    cfg = load_config()
    assert cfg["disorder"]["score_threshold"] == 0.5
    override = tmp_path / "c.yaml"
    override.write_text("disorder:\n  score_threshold: 0.7\n")
    cfg2 = load_config(override)
    assert cfg2["disorder"]["score_threshold"] == 0.7
    # Unspecified defaults survive the merge.
    assert cfg2["disorder"]["min_segment_length"] == 10


def test_write_table_writes_parquet_and_csv(tmp_path) -> None:
    df = pd.DataFrame({"a": [1, 2], "b": ["x", "y"]})
    out = tmp_path / "feats.parquet"
    parquet_path, csv_path = write_table(df, out)
    assert parquet_path.exists() and csv_path.exists()
    assert csv_path.suffix == ".csv"
    pd.testing.assert_frame_equal(pd.read_parquet(parquet_path), df)


def test_example_fasta_is_bundled_and_parses() -> None:
    path = example_fasta_path()
    assert path.exists()
    seqs = read_fasta(path)
    assert len(seqs) >= 1
    assert all(seq for seq in seqs.values())
