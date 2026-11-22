# Submitting the disorder head to CAID

[CAID](https://caid.idpcentral.org) is the Critical Assessment of Intrinsic
Disorder: the organisers run participants' software on blind targets and score
it themselves. That is the point. Our own CAID3 number is a reimplementation
of their protocol against their published reference file — honest, but we held
the answers while computing it. A CAID ranking is the same measurement made
by someone who does not.

The disorder head is the one worth entering. On CAID3 Disorder-PDB it scores
AUROC 0.941 by our reckoning, against published figures of 0.925 for AIUPred
v3 and 0.954 for MUSE-IDR.

## What CAID needs, and what we hand over

CAID wants runnable software, not a file of predictions, which is why the CLI
is the submission. Given a FASTA of targets it writes their format directly:

```bash
opendna --models models predict --task disorder \
    -i targets.fasta --format caid -o OpenDNA.caid
```

The format, verified against `vectorized_metrics/parsers.py` in
[BioComputingUP/CAID](https://github.com/BioComputingUP/CAID) and against
`AIUPred.caid`, a real CAID3 entry, rather than inferred from prose:

```text
>DP02732
1	M	0.727	1
2	E	0.774	1
```

One `>accession` block per target; tab-separated 1-based position, residue,
score, and the binary call. Their parser accepts three columns and derives the
call by thresholding at 0.5 — which is the threshold this head was evaluated
at — but we emit four so the call is ours on the record.

### Two things that will silently ruin a submission

**The accession must be the first token of the FASTA header.** Their reference
parser does `acc, *desc = header[1:].split()`, and `parse_prediction` skips
any target whose accession is not in the reference *without complaining*. Emit
the whole header and the submission scores zero targets while looking like it
merely failed. The CLI takes the first token.

**This format carries no provenance comments.** Every other output the CLI
writes has a `#` block with the model version and the card's limitations. Here
it cannot: their reader groups lines on whether they start with `>`, so a
leading comment becomes a block of its own, is read as an accession, and
shifts every later target's body onto the wrong protein. It would parse
without error and score nonsense. Provenance goes to a sidecar
`<output>.provenance.txt` instead.

## Before submitting

1. **Check the deadline and the current categories** on
   <https://caid.idpcentral.org>. CAID4 has launched, with disorder, binding
   and linker categories, but the dates are on the site and are not reproduced
   here because a stale deadline in a repository is worse than none.

2. **Enter only the disorder category.** There is no binding or linker head,
   and entering a category you cannot predict wastes the organisers' compute.

3. **Confirm the targets parse.** Real CAID targets contain residues the
   strict FASTA reader refuses — `DP03745` has a selenocysteine at position
   48. The CLI reads sequences as written for exactly this reason, and ESM-2
   has a distinct token for `U`, so the residue is used rather than collapsed
   to `X`. If a future target carries something stranger, the run should still
   complete; check the output covers every target.

4. **Expect windowing on long targets.** 34 of the 319 CAID3 Disorder-PDB
   proteins exceed ESM-2's 1024-residue limit and are covered by overlapping
   windows. No published CAID entrant does this, so those targets are not
   scored the way another method would score them. It is in the benchmark
   report's caveats and should be mentioned in the submission notes.

5. **State what the head is.** A small network on frozen ESM-2 150M
   embeddings — not a fine-tuned language model. Trained on DisProt 2025_12
   with CAID3 proteins held out of training entirely, which was verified
   rather than assumed: 318 of the 319 appear only in our test split, and
   none in train or validation.

## Validating the output without running their code

`parse_prediction` is short enough to reimplement, and doing so avoids
executing a cloned repository. The checks worth making are: every accession
matches the reference; the residue column equals the reference sequence;
positions are contiguous from 1; the binary column equals the scores
thresholded at 0.5; and the recomputed AUROC over assessed residues matches
what `opendna-train bench` reported. If the last two agree, the file and the
benchmark describe the same model.
