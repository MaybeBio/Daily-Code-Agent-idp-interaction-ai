# Demo

A small, self-contained demo that runs in about ten seconds on a normal desktop
computer, with no GPU. Everything it needs is in this folder.

It scores sixteen *Saccharomyces cerevisiae* transcription factors (12,981 total
residues) and reports the activation domains the model predicts in each. The
model was never trained on these sequences.

## Running the demo

Install the package first — see
[Installation guide](https://github.com/staller-lab/consensus_model_package/blob/HEAD/README.md#2-installation-guide) in the main README.
Then, from the repository root:

```
python demo/run_demo.py
```

## Expected runtime

**About 10 seconds** on a Windows 11 desktop (Intel i7, 16 GB RAM, no GPU).
Scoring all sixteen proteins takes **0.2 seconds**; the rest is importing
PyTorch and loading the model, which is only 13 KB.

## Expected output

Predicted domains are printed to the terminal and written to
`demo/output/predicted_ADs.csv`, one row per domain with its protein, start and
end coordinates (1-based, inclusive), length and sequence.

A correct run reports **17 predicted activation domains across 12 of the 16
proteins**, scored at the package's default threshold of 0.6.

That figure is not the whole of the model's output. On these sixteen proteins
the package returns **52** runs above the threshold, and `run_demo.py` lists
only the 17 that are at least 10 aa long. The other 35 are too short to be real
domains: 19 are single residues, 6 are 2 aa, 6 are 3 aa, and the rest are 6–9 aa.
The shortest domain listed is 11 aa and the longest is 62 aa.

This filter belongs to the demo script, not to the package.
`madrigal-predictor` has no minimum-length setting and returns every run above
the threshold. Set `MIN_LENGTH = 1` at the top of `run_demo.py` to list all 52.

The output begins:

```
ARG81     1 predicted activation domain
          residues  209- 240  ( 32 aa)  VPPTPNLLDYDWNNLNITGYEWISSELRDDAL
CUP2      1 predicted activation domain
          residues  171- 187  ( 17 aa)  LEQSPSISMEPVSINET
```

### Checking the result

Four things can be checked without any external data:

- **Totals.** 17 domains across 12 of 16 proteins. No domain is predicted in
  RDS2, ASG1, RFX1 or FAP1.
- **RFX1 is a repressor.** RFX1 (CRT1) is the major repressor of
  DNA-damage-inducible genes in *S. cerevisiae*, and MADRIGAL predicts no
  activation domain in it. This is the clearest negative control in the set.
- **ECM22 and UPC2 are paralogues**, and the model recovers the same conserved
  motif in both without being told they are related: `FDESIADLYPVEID` at
  residues 685–698 of ECM22, and `FDESIADLYPVGLDSP` within the domain at
  757–798 of UPC2.
- **Composition.** The predicted domains are acidic and hydrophobic and sit
  mostly toward the C-terminus, as expected of yeast activation domains.
  Several end in a characteristic acidic run, for example `DFFDYDLLFQND` in
  PIP2 and `FGGLDLFDYDFLFGNDFA` in OAF1.

If those hold, the installation is working correctly.

## Parameters

Both are set at the top of `run_demo.py`:

| Setting | Value | Meaning |
|---|---|---|
| `THRESHOLD` | 0.6 | Passed to `madrigal-predictor`. A residue is called part of an activation domain when the predicted consensus — the fraction of the six source predictors expected to call it — is at or above this. 0.6 is the package default. |
| `MIN_LENGTH` | 10 | Used by this script only, after the model has run. Domains shorter than this are left out of the printed listing and the CSV. The package itself applies no length filter. Set to 1 to list everything. |

## Input file

`demo_seqs.fasta` — the sixteen yeast transcription factors of the held-out
test set.
