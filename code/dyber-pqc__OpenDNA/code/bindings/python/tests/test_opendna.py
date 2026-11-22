"""Tests for the Python bindings.

Split by what they need. The refusals need only the binary, so they run
anywhere it can be built. The predictions need trained heads, which are not in
the repository, so they skip with a message rather than failing on every
machine that has not downloaded them -- the same trade the Rust exit-code test
makes, and for the same reason: a test that fails for everybody gets disabled,
and then it checks nothing.
"""

from __future__ import annotations

import os
import shutil
from pathlib import Path

import pytest

import opendna

UBIQUITIN = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"


def models_dir() -> Path | None:
    """Where the heads are, if they are anywhere."""
    candidate = Path(os.environ.get("OPENDNA_MODELS", "models"))
    for path in (candidate, Path(__file__).resolve().parents[3] / "models"):
        if (path / "disorder.safetensors").is_file():
            return path
    return None


needs_heads = pytest.mark.skipif(models_dir() is None, reason="no trained heads installed")


def test_the_binary_is_found_or_the_error_says_how_to_get_one():
    try:
        assert Path(opendna.binary()).exists()
    except opendna.ToolMissing as error:
        assert "releases" in str(error) or "OPENDNA_BIN" in str(error)


def test_an_unknown_task_is_refused_before_anything_runs():
    with pytest.raises(opendna.InputRejected) as caught:
        opendna.predict("folding", {"u": UBIQUITIN})
    assert "folding" in str(caught.value)


def test_stability_without_mutations_says_what_it_needs():
    with pytest.raises(opendna.InputRejected) as caught:
        opendna.predict("stability", {"u": UBIQUITIN})
    # The message has to name the argument and show its shape, because this is
    # the mistake every first-time caller makes.
    assert "mutations" in str(caught.value)
    assert "L8P" in str(caught.value)


def test_mutations_are_refused_on_a_task_that_cannot_use_them():
    with pytest.raises(opendna.InputRejected):
        opendna.predict("disorder", {"u": UBIQUITIN}, mutations=[("u", "L8P")])


@pytest.mark.parametrize(
    "sequences",
    [
        {},
        {"u": ""},
        {"has\nnewline": UBIQUITIN},
        {"has>bracket": UBIQUITIN},
    ],
)
def test_sequences_that_would_corrupt_the_fasta_are_refused(sequences):
    with pytest.raises(opendna.InputRejected):
        opendna.predict("disorder", sequences)


def test_a_missing_head_raises_rather_than_estimating(tmp_path):
    """The binding must inherit the tool's refusal, not soften it."""
    empty = tmp_path / "no-models"
    empty.mkdir()
    with pytest.raises(opendna.HeadMissing) as caught:
        opendna.predict("disorder", {"u": UBIQUITIN}, models=empty)
    assert "not installed" in str(caught.value)


def test_models_reports_an_empty_directory_without_failing(tmp_path):
    empty = tmp_path / "no-models"
    empty.mkdir()
    report = opendna.models(models_dir=empty)
    assert "not installed" in report


@needs_heads
def test_disorder_returns_a_row_per_residue_with_its_provenance():
    run = opendna.predict("disorder", {"ubiquitin": UBIQUITIN}, models=models_dir())
    assert len(run) == len(UBIQUITIN)
    assert {"name", "position", "residue", "probability", "disordered"} <= set(run.rows[0])
    # The point of the binding: the numbers do not arrive without their source.
    assert run.build, "a run with no build id cannot be traced to a model"
    assert run.limitations, "the card's limitations must travel with the numbers"
    assert run.manifest["input"]["sha256"], "the manifest must record what went in"
    assert run.head_version in run.cite()


@needs_heads
def test_stability_scores_the_mutations_it_was_given():
    run = opendna.predict(
        "stability",
        {"ubiquitin": UBIQUITIN},
        mutations=[("ubiquitin", "L8P"), ("ubiquitin", "I44D")],
        models=models_dir(),
    )
    assert [row["mutation"] for row in run.rows] == ["L8P", "I44D"]
    # Every regressor reports an interval beside the point, and a caller who
    # reads only the point is the reason it is there.
    for row in run.rows:
        assert row["interval_low"] < row["ddg_kcal_per_mol"] < row["interval_high"]


@needs_heads
def test_a_mutation_the_sequence_contradicts_is_refused():
    with pytest.raises(opendna.InputRejected) as caught:
        opendna.predict(
            "stability",
            {"ubiquitin": UBIQUITIN},
            mutations=[("ubiquitin", "Q8P")],
            models=models_dir(),
        )
    assert "position 8" in str(caught.value)


@needs_heads
@pytest.mark.skipif(shutil.which("python") is None, reason="needs an interpreter")
def test_the_frame_has_a_row_per_result_when_pandas_is_installed():
    pytest.importorskip("pandas")
    run = opendna.predict("burial", {"ubiquitin": UBIQUITIN}, models=models_dir())
    frame = run.frame
    assert len(frame) == len(UBIQUITIN)
    assert "accessibility" in frame.columns
