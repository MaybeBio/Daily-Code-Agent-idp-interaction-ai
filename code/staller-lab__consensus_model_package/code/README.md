# madrigal-predictor

**MADRIGAL** — Multi-model Activation Domain Recurrent-neural-network Inferring
Grouped Activity Labels.

MADRIGAL predicts transcriptional **activation domains** in protein sequences.
It is a bidirectional LSTM, trained with [PARROT](https://github.com/idptools/parrot),
that takes a protein sequence of any length and predicts, for each residue, the
fraction of six published activation-domain predictors that would call that
residue part of an activation domain. Residues scoring above a threshold are
reported as predicted activation domains.

The six predictors it summarises are PADDLE, TADA, ADpred, ADhunterLite,
ADhunterEnsemble and the Kotha composition model. Running all six over a
proteome is slow; MADRIGAL reproduces their consensus in a single fast pass.

- Web interface: <https://www.madrigal-predictor.net>
- Package: <https://pypi.org/project/madrigal-predictor/>
- Model development and evaluation: [staller-lab/Madrigal-development](https://github.com/staller-lab/Madrigal-development)
- Experimental analysis: [staller-lab/designed-ads](https://github.com/staller-lab/designed-ads)

---

## 1. System requirements

### Software dependencies

| Package | Version |
|---|---|
| python | >= 3.9, < 3.12 |
| idptools-parrot | >= 1.7.5 |
| torch | >= 2.0 |
| numpy | >= 2.0 |

Installing the package pulls these in automatically. A clean install on
2026-11-19 resolved `idptools-parrot` 1.7.5, `torch` 2.15.0 and `numpy` 2.4.6.

### Operating systems tested

Tested on Windows 11 Home (build 26200), where both the installation and the
demo were verified. The package is pure Python and depends only on
cross-platform libraries, so it should run on any system meeting the
requirements above, macOS and Linux included.

### Required non-standard hardware

**None.** MADRIGAL runs on a standard desktop computer using the CPU only. The
trained model is 13 KB, and scoring sixteen full-length transcription factors
(12,981 residues) takes under a second.

A GPU was used to *train* the model and to run one of the six source
predictors, but neither is needed to use this package. That work is documented
in [Madrigal-development](https://github.com/staller-lab/Madrigal-development).

---

## 2. Installation guide

### Instructions

MADRIGAL requires Python 3.9, 3.10 or 3.11. PyTorch does not yet support 3.12
for this dependency set, so create a dedicated environment:

```
conda create -n madrigal python=3.11
conda activate madrigal
pip install madrigal-predictor
```

Or, with an existing Python 3.9–3.11 installation:

```
python -m venv madrigal
source madrigal/bin/activate        # Windows: madrigal\Scripts\activate
pip install madrigal-predictor
```

Nothing else needs to be installed by hand. PyTorch, NumPy and PARROT are
pulled in automatically, and the trained model ships inside the package. To
check that the installation worked, run the demo in section 3.

### Typical install time

**About 8 minutes** on a standard desktop computer with a broadband connection
— roughly 1 minute to create the environment and 7 minutes to download and
install the package with its PyTorch dependency. Measured on Windows 11
(Intel i7, 16 GB RAM) with no packages previously cached. Subsequent installs
from a warm pip cache take under a minute.

---

## 3. Demo

The `demo/` folder contains everything needed. No data has to be downloaded.

### Instructions to run

```
python demo/run_demo.py
```

The demo scores `demo/demo_seqs.fasta` — sixteen *Saccharomyces cerevisiae*
transcription factors (12,981 residues) that form the held-out test set used to
evaluate MADRIGAL in the paper. The model was never trained on them.

### Expected output

Predicted activation domains are printed to the terminal and written to
`demo/output/predicted_ADs.csv`, one row per domain with its protein, start and
end coordinates (1-based, inclusive), length and sequence.

A correct run reports **17 predicted activation domains across 12 of the 16
proteins**, scored at the package's default threshold of 0.6. The demo script
leaves runs shorter than 10 aa out of its listing — a choice made by
`run_demo.py` when printing, not by the package, which applies no length filter.
The output begins:

```
ARG81     1 predicted activation domain
          residues  209- 240  ( 32 aa)  VPPTPNLLDYDWNNLNITGYEWISSELRDDAL
CUP2      1 predicted activation domain
          residues  171- 187  ( 17 aa)  LEQSPSISMEPVSINET
```

Four results can be checked without any external data:

- **Totals.** 17 domains across 12 of 16 proteins. No domain is predicted in
  RDS2, ASG1, RFX1 or FAP1.
- **RFX1 is a repressor.** RFX1 (CRT1) is the major repressor of
  DNA-damage-inducible genes in *S. cerevisiae*, and MADRIGAL predicts no
  activation domain in it.
- **ECM22 and UPC2 are paralogues**, and the same conserved motif is recovered
  in both: `FDESIADLYPVEID` at residues 685–698 of ECM22, and
  `FDESIADLYPVGLDSP` within the domain at 757–798 of UPC2.
- **Composition.** The predicted domains are acidic and hydrophobic and sit
  mostly toward the C-terminus, as expected of yeast activation domains.

### Expected run time

**About 10 seconds** on a standard desktop computer (Windows 11, Intel i7,
16 GB RAM, no GPU). Scoring all sixteen proteins takes **0.2 seconds**; the
remainder is importing PyTorch and loading the model.

---

## 4. Instructions for use

### Running on your own data

The package exposes two functions.

**A FASTA file:**

```python
from madrigal_predictor import analyze_fasta_file

titles, above_runs, below_runs, above_starts, above_ends = analyze_fasta_file(
    "my_sequences.fasta", threshold=0.6
)
```

**A single sequence:**

```python
from madrigal_predictor import read_sequence

above_runs, below_runs, above_starts, above_ends = read_sequence(
    "MSEYQPSLFALNPMGFSPLDGSKSTNENVSAST...", threshold=0.6
)
```

Returned values, one list per input sequence for `analyze_fasta_file`:

| Value | Contents |
|---|---|
| `titles` | FASTA headers, without the leading `>` |
| `above_runs` | Subsequences predicted to be activation domains |
| `below_runs` | The remaining subsequences |
| `above_starts` | Start position of each predicted domain |
| `above_ends` | End position of each predicted domain, inclusive |

**Coordinates are 0-based and inclusive**, so a domain is
`sequence[start:end + 1]`. To report 1-based residue numbers, as used in the
literature, add 1 to each start and end — this is what `demo/run_demo.py`
prints.

### Choosing a threshold

The threshold is the fraction of the six source predictors expected to call a
residue part of an activation domain. The package default is **0.6**.

For a result with high confidence, we recommend a threshold of 0.8. For a
confident result, 0.6 is appropriate. For a more inclusive result, set your
threshold to 0.4.

| Threshold | Result |
|---|---|
| 0.8 | High confidence |
| 0.6 | Confident — the package default |
| 0.4 | More inclusive |

Predicted domains are more likely to be genuine when they fall in intrinsically
disordered regions, so it is worth checking disorder — with
[metapredict](https://github.com/idptools/metapredict), for example — alongside
these predictions.

### Practical notes

- MADRIGAL takes a protein sequence of arbitrary length and returns a
  per-residue prediction, so there is no minimum or maximum input length. The
  sequences it was trained on were between 53 and 1,442 residues.
- Input is case-insensitive and whitespace is stripped.

### Reproducing the published results

Model training, evaluation and all figures are in
[Madrigal-development](https://github.com/staller-lab/Madrigal-development)
and [designed-ads](https://github.com/staller-lab/designed-ads).

---

## Citation

> Zintel, M. A. *et al.* Machine learning guided design of strong
> transcriptional activation domains. *bioRxiv* 2026.10.30.763198 (2026).
> doi:[10.64898/2026.10.30.763198](https://doi.org/10.64898/2026.10.30.763198)

Please also cite PARROT, which MADRIGAL is built on, and the six source
predictors where appropriate.

## Credits

Developed by the [Staller Lab](https://www.stallerlab.com/): Claire Josephine
LeBlanc, Marissa Ann Zintel, Sanjana R. Kotha and Rishaan Jayce Bansal, using
the PARROT package from the Holehouse Lab
(<https://github.com/idptools/parrot>).

## License

MIT — see [`LICENSE`](https://github.com/staller-lab/consensus_model_package/blob/HEAD/LICENSE).
