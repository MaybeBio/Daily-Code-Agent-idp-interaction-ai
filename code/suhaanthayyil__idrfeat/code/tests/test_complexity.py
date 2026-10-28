"""Tests for SEG low complexity and the longest residue run."""

from __future__ import annotations

import pytest

from idrfeat.complexity import longest_single_run, seg_lowcomplexity_frac, segmasker_available


def test_longest_single_run() -> None:
    assert longest_single_run("AABBBBA") == 4
    assert longest_single_run("ABCDE") == 1
    assert longest_single_run("") == 0


def test_seg_lowcomplexity_high_for_homopolymer() -> None:
    if not segmasker_available():
        pytest.skip("segmasker (NCBI SEG) not installed")
    assert seg_lowcomplexity_frac("A" * 40) > 0.9
    assert seg_lowcomplexity_frac("ACDEFGHIKLMNPQRSTVWY" * 3) < 0.5


def test_seg_empty_is_zero() -> None:
    assert seg_lowcomplexity_frac("") == 0.0
