---
license: other
license_name: cc-by-sa-2.1-jp
license_link: https://creativecommons.org/licenses/by-sa/2.1/jp/deed.en
library_name: opendna
tags:
  - protein
  - biology
  - protein-language-model
  - solubility
  - expression
  - regression
datasets:
  - dyber-pqc/opendna-solubility
base_model: facebook/esm2_t30_150M_UR50D
---

# OpenDNA solubility head

Predicts protein solubility from sequence.

This is a small head trained on **frozen** [facebook/esm2_t30_150M_UR50D](https://huggingface.co/facebook/esm2_t30_150M_UR50D) embeddings. It is not a fine-tuned language model: the embedding model is unchanged, and only the head below was trained.

## Using it

These are not transformers weights and no `AutoModel` will load them. They are read by [OpenDNA](https://github.com/dyber-pqc/OpenDNA), which is a single binary with no runtime to install — take it from the [releases page](https://github.com/dyber-pqc/OpenDNA/releases) or build it with `cargo build --release -p opendna-cli`.

```bash
hf download dyber-pqc/opendna-solubility --local-dir models
export OPENDNA_MODELS=$PWD/models

# confirm the download is intact: the weights are hashed against the
# digest recorded in their metadata
opendna models

opendna predict --task solubility -i candidates.fasta -o solubility.csv
```

The output carries the head version, the embedding model, the device and the limitations below as `#` comments, and a `solubility.manifest.json` beside it records the digests of what went in and what came out. The first run downloads the embedding model; after that `opendna stage` makes a directory that works with no network at all.

`opendna models` checks that the file is the one it claims to be. `opendna selftest` goes further and re-runs the heads over a known sequence to check they reproduce recorded values — worth doing, and it needs all four heads installed.

## Inputs and outputs

**Input.** The mean-pooled ESM-2 embedding of the sequence.

**Output.** A solubility percentage on the source assay's scale. Values above 100 occur in the source data and are not clipped.

## Results

| Metric | Validation | Test | Baseline (test) |
|---|---|---|---|
| MAE | 18.313 | 18.823 | 28.416 |
| Pearson r | 0.740 | 0.691 | — |
| RMSE | 23.015 | 23.388 | 32.256 |
| Spearman ρ | 0.737 | 0.701 | — |

A — in the baseline column means the metric is undefined for that baseline, not unmeasured. This baseline predicts the same value for every input, so it has no correlation with the truth.

Baseline: Predicting the training-set mean for every protein.

Validation chose the epoch that was kept (epoch 35 of 45 run). The test split was evaluated once, afterwards. Splits come from the dataset builder, which clusters by sequence identity so near-identical proteins cannot appear on both sides.

### Against published benchmarks

None. The external benchmarks for this task come from sources whose terms do not permit redistribution, or from services that have been retired, so there is no set this head can be scored on alongside other published methods. The figures above are therefore not evidence of where it stands relative to the field — only of how it does on data drawn the way its training data was drawn.

## Training

| | |
|---|---|
| Embedding model | `facebook/esm2_t30_150M_UR50D` (640 dimensions, frozen) |
| Optimiser | AdamW, learning rate 0.001 |
| Batch size | 256 |
| Seed | 20260101 |
| Dataset | [`dyber-pqc/opendna-solubility`](https://huggingface.co/datasets/dyber-pqc/opendna-solubility) |
| Data licence | CC BY-SA 2.1 JP (eSOL, Niwa et al. 2009); the UniProtKB sequences it is joined against are CC BY 4.0 |

### Licence of these weights

The eSOL measurements are licensed CC BY-SA 2.1 JP, which is share-alike. Whether trained weights are a derivative work of their training data is an open legal question, and this project is not going to answer it in its own favour, so these weights are released under the same share-alike terms rather than the Apache-2.0 used for the other heads.

### Reproducing this run

Download the dataset from [`dyber-pqc/opendna-solubility`](https://huggingface.co/datasets/dyber-pqc/opendna-solubility) into `datasets/out/solubility/`, or rebuild it from its upstream sources with the `datasets/` workspace of the [OpenDNA repository](https://github.com/dyber-pqc/OpenDNA) — the dataset card lists every source file with its SHA-256. Then, from `training/`:

```bash
cargo run --release -p opendna-train -- embed --task solubility
cargo run --release -p opendna-train -- train --task solubility \
    --batch-size 256 --learning-rate 0.001 --seed 20260101
```

The hyperparameters are spelled out above even where they match the current defaults, so the command reproduces *this* run rather than whatever the defaults happen to be later.

Both commands take `--datasets <dir>` if the dataset is not at the default `../datasets/out`. Embedding is cached, so the second run of either command is cheap.

## Limitations

- Trained on E. coli cell-free expression measurements. Solubility in other hosts, or under different expression conditions, is a different quantity.
- Sequences are UniProtKB canonical sequences, not the exact expressed constructs; tags and truncations in the real experiments are not represented.
- A few thousand proteins is a small training set. Treat the output as a ranking signal for prioritising constructs, not as a predicted assay value.

## Intended use

Ranking and triage during protein engineering: comparing candidate mutations or constructs against each other. It is not a substitute for measurement, and it should not be used to make claims about a specific protein without experimental confirmation.

## Citation

The training data comes from:

> Niwa, T. et al. Bimodal protein solubility distribution revealed by an aggregation analysis of the entire ensemble of Escherichia coli proteins. PNAS 106, 4201-4206 (2009).

Please cite that work as well as this model.

## How wrong the number is likely to be

A point estimate invites being read as a measurement. These intervals say how wide the range around it has to be, fitted by split conformal prediction on the validation split: the half-width is whatever the residuals require, with no assumption about the shape of the error.

| Confidence | Interval | Coverage achieved on held-out data |
|---|---|---|
| 50% | ±15.675 percentage points | 50.2% |
| 80% | ±29.515 percentage points | 77.0% |
| 90% | ±38.945 percentage points | 91.1% |

Mean absolute error on the same split is 18.823 percentage points.

The achieved column is there because it differs from the nominal one, and it can differ either way. A conformal guarantee assumes the calibration and held-out sets are exchangeable, and this project's splits are deliberately not: the dataset builder clusters by sequence identity so that near-identical proteins cannot straddle a split, which makes the held-out set systematically more distant than a random sample. That is the right way to build the splits, and the cost is that the guarantee is approximate. The measured column is the one to believe, in whichever direction it points. On this head it falls short at some levels and not others; read the column.

Each protein contributes exactly one residual, so unlike the per-residue heads these observations are independent. There are simply few of them: the split the quantile is fitted on holds a few hundred proteins, which is why the interval is wide and why it should be read as the honest width rather than a tight one.

The half-width is constant, which is the method's main limitation: the interval is wider than it needs to be where the head is reliable and narrower than it should be where it is not.
