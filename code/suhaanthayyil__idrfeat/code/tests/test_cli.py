"""End-to-end test of the `idrfeat features` command on the bundled example."""

from __future__ import annotations

import pandas as pd
import pytest

from idrfeat.cli import main
from idrfeat.disorder import metapredict_available
from idrfeat.features import TABLE_COLUMNS
from idrfeat.io import example_fasta_path

pytestmark = pytest.mark.skipif(
    not metapredict_available(), reason="metapredict is required to call IDR segments"
)


def test_features_command_writes_table(tmp_path) -> None:
    out = tmp_path / "features.parquet"
    rc = main(
        [
            "features",
            "--fasta",
            str(example_fasta_path()),
            "--out",
            str(out),
        ]
    )
    assert rc == 0
    assert out.exists()
    assert out.with_suffix(".csv").exists()
    df = pd.read_parquet(out)
    assert list(df.columns) == TABLE_COLUMNS
    assert len(df) >= 1
    assert (df["disorder_backend"] == "metapredict").all()
