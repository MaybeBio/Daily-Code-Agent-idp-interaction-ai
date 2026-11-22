"""Checks that reported coordinates actually locate the predicted domains.

Run with:

    python -m pytest tests/
    python tests/test_coordinates.py      # also works without pytest

These cover the contract that start and end positions are 0-based and
inclusive, so that sequence[start:end + 1] returns exactly the reported
subsequence. An earlier version returned the start positions of the
below-threshold stretches as the fourth value, which did not line up with the
activation domains whenever a sequence began below threshold.
"""

from pathlib import Path

from madrigal_predictor import analyze_fasta_file, read_sequence
from madrigal_predictor.package_script import sanitize

HERE = Path(__file__).resolve().parent
FASTA = HERE.parent / "demo" / "demo_seqs.fasta"


def _records(path):
    """Yield (title, sequence) from a FASTA file."""
    title, seq = None, ""
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if title is not None:
                yield title, sanitize(seq)
            title, seq = line[1:], ""
        else:
            seq += line
    if title is not None:
        yield title, sanitize(seq)


def test_coordinates_locate_the_domains():
    """start/end must slice the reported subsequence back out."""
    titles, above_runs, _below, starts, ends = analyze_fasta_file(FASTA)
    sequences = [s for _t, s in _records(FASTA)]

    assert len(titles) == len(sequences)
    assert len(above_runs) == len(starts) == len(ends) == len(sequences)

    for seq, runs, ss, ee in zip(sequences, above_runs, starts, ends):
        assert len(runs) == len(ss) == len(ee)
        for run_seq, start, end in zip(runs, ss, ee):
            assert 0 <= start <= end < len(seq)
            assert end - start + 1 == len(run_seq)
            assert seq[start:end + 1] == run_seq


def test_domains_are_ordered_and_disjoint():
    _titles, above_runs, _below, starts, ends = analyze_fasta_file(FASTA)
    for runs, ss, ee in zip(above_runs, starts, ends):
        for j in range(1, len(runs)):
            assert ss[j] > ee[j - 1], "domains must not overlap or go backwards"


def test_single_sequence_matches_fasta_result():
    """read_sequence and analyze_fasta_file must agree on the same sequence."""
    titles, above_runs, _below, starts, ends = analyze_fasta_file(FASTA)
    sequences = [s for _t, s in _records(FASTA)]

    for seq, runs, ss, ee in zip(sequences, above_runs, starts, ends):
        r2, _b2, s2, e2 = read_sequence(seq)
        assert r2 == runs
        assert s2 == ss
        assert e2 == ee


def test_threshold_is_monotonic():
    """A stricter threshold cannot call more residues than a looser one."""
    seq = next(s for _t, s in _records(FASTA))
    loose, _b1, _s1, _e1 = read_sequence(seq, threshold=0.4)
    strict, _b2, _s2, _e2 = read_sequence(seq, threshold=0.8)
    assert sum(len(r) for r in strict) <= sum(len(r) for r in loose)


if __name__ == "__main__":
    failures = 0
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            try:
                fn()
                print(f"[PASS] {name}")
            except AssertionError as exc:
                failures += 1
                print(f"[FAIL] {name}: {exc}")
    print("\n" + ("all tests passed" if not failures
                  else f"{failures} test(s) failed"))
    raise SystemExit(1 if failures else 0)
