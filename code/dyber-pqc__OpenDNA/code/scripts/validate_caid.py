"""Validate a .caid submission, independently of the Rust that produced it.

    python scripts/validate_caid.py predictions.caid reference.fasta

This reimplements the rules in `vectorized_metrics/parsers.py` from
BioComputingUP/CAID rather than executing that repository, and recomputes
AUROC in numpy. Two implementations of one measurement is how an error in
either gets noticed: on CAID3 Disorder-PDB this reports 0.9409 against the
0.9410 from `opendna-train bench`, a difference in tie handling in the
trapezoid rule rather than in the model.

It is deliberately Python and deliberately not sharing code with the CLI. A
checker written in the same language, against the same types, by the same
author, mostly checks that the code equals itself. The point here is that the
file on disk is readable by something else.

Checks, each of which corresponds to a way a submission fails silently:

  - blocks group on lines starting with '>' (so a stray comment desynchronises
    every later target)
  - the accession is header[1:].strip(), matched against the reference (an
    unmatched target is skipped without complaint, scoring nothing)
  - body lines split on tabs into 3 columns or 4
  - positions run 1..n contiguously (their scorer joins on position)
  - the residue column equals the reference sequence
  - the binary column equals the published score thresholded at 0.5
"""

import sys
from itertools import groupby

import numpy as np

THRESHOLD = 0.5


def read_reference(path):
    """Reference records: header, sequence, states. '-' is not assessed."""
    pattern = {"0": 0.0, "1": 1.0, "-": np.nan}
    lines = [l.rstrip("\n") for l in open(path) if l.strip()]
    if len(lines) % 3:
        sys.exit(f"{path}: {len(lines)} non-empty lines is not a multiple of 3")
    ref = {}
    for i in range(0, len(lines), 3):
        acc = lines[i][1:].strip().split()[0]
        ref[acc] = (lines[i + 1], np.array([pattern.get(c, np.nan) for c in lines[i + 2]]))
    return ref


def read_predictions(path):
    """Prediction blocks, grouped the way their reader groups them."""
    parsed = {}
    with open(path) as fh:
        faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
        for acc_group in faiter:
            acc = next(acc_group).strip()[1:]
            rows = [l.strip("\n").rsplit("\t") for l in next(faiter)]
            widths = {len(r) for r in rows}
            if not widths <= {3, 4}:
                sys.exit(f"{acc}: column counts {sorted(widths)}; their parser accepts 3 or 4")
            positions = [int(r[0]) for r in rows]
            if positions != list(range(1, len(rows) + 1)):
                sys.exit(f"{acc}: positions are not contiguous from 1; the scorer joins on them")
            parsed[acc] = {
                "seq": "".join(r[1] for r in rows),
                "scores": np.array([float(r[2]) for r in rows]),
                "states": np.array([float(r[3]) for r in rows]) if len(rows[0]) == 4 else None,
            }
    return parsed


def auroc(scores, truth):
    """Rank-based AUROC, written out rather than imported, so this does not
    agree with the Rust by both calling the same library."""
    order = np.argsort(-scores, kind="mergesort")
    t = truth[order]
    positives, negatives = t.sum(), (1 - t).sum()
    if positives == 0 or negatives == 0:
        sys.exit("one class is empty; AUROC is undefined")
    tpr = np.concatenate([[0.0], np.cumsum(t) / positives])
    fpr = np.concatenate([[0.0], np.cumsum(1 - t) / negatives])
    return float(np.trapezoid(tpr, fpr))


def main(argv):
    if len(argv) != 3:
        sys.exit(__doc__)
    pred_path, ref_path = argv[1], argv[2]
    ref = read_reference(ref_path)
    pred = read_predictions(pred_path)

    print(f"targets in file:     {len(pred)}")
    matched = [a for a in pred if a in ref]
    print(f"matched to reference: {len(matched)} / {len(pred)}")
    if not matched:
        sys.exit("no accession matched the reference: this submission would score nothing")
    if len(matched) < len(ref):
        print(f"note: {len(ref) - len(matched)} reference targets have no prediction")

    for acc in matched:
        ref_seq, _ = ref[acc]
        got = pred[acc]
        if len(got["scores"]) != len(ref_seq):
            sys.exit(f"{acc}: {len(got['scores'])} scores for {len(ref_seq)} reference residues")
        if got["seq"] != ref_seq:
            sys.exit(f"{acc}: the residue column disagrees with the reference sequence")
        states = got["states"]
        if states is None:
            continue
        if not set(np.unique(states)) <= {0.0, 1.0}:
            sys.exit(f"{acc}: the binary column is not 0/1")
        if not np.array_equal(states, (got["scores"] >= THRESHOLD).astype(float)):
            sys.exit(
                f"{acc}: the binary column is not the published score thresholded at "
                f"{THRESHOLD}. The file contradicts itself, which reads as malformed even "
                f"though their parser prefers the state column."
            )

    scores, truth = [], []
    for acc in matched:
        _, states = ref[acc]
        keep = ~np.isnan(states)
        scores.append(pred[acc]["scores"][keep])
        truth.append(states[keep])
    scores, truth = np.concatenate(scores), np.concatenate(truth)

    print(f"residues assessed:    {len(truth)}")
    print(f"positive rate:        {truth.mean():.4f}")
    print(f"AUROC (numpy):        {auroc(scores, truth):.4f}")
    print("every check passed")


if __name__ == "__main__":
    main(sys.argv)
