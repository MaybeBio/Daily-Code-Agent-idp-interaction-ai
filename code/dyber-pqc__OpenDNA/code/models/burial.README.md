---
license: apache-2.0
library_name: opendna
tags:
  - protein
  - biology
  - protein-language-model
  - solvent-accessibility
  - burial
  - rsa
  - per-residue
  - regression
datasets:
  - dyber-pqc/opendna-burial
base_model: facebook/esm2_t30_150M_UR50D
---

# OpenDNA burial head

Predicts, for each residue, how solvent-exposed it is.

This is a small head trained on **frozen** [facebook/esm2_t30_150M_UR50D](https://huggingface.co/facebook/esm2_t30_150M_UR50D) embeddings. It is not a fine-tuned language model: the embedding model is unchanged, and only the head below was trained.

## Using it

These are not transformers weights and no `AutoModel` will load them. They are read by [OpenDNA](https://github.com/dyber-pqc/OpenDNA), which is a single binary with no runtime to install — take it from the [releases page](https://github.com/dyber-pqc/OpenDNA/releases) or build it with `cargo build --release -p opendna-cli`.

```bash
hf download dyber-pqc/opendna-burial --local-dir models
export OPENDNA_MODELS=$PWD/models

# confirm the download is intact: the weights are hashed against the
# digest recorded in their metadata
opendna models

opendna predict --task burial -i candidates.fasta -o burial.csv
```

The output carries the head version, the embedding model, the device and the limitations below as `#` comments, and a `burial.manifest.json` beside it records the digests of what went in and what came out. The first run downloads the embedding model; after that `opendna stage` makes a directory that works with no network at all.

`opendna models` checks that the file is the one it claims to be. `opendna selftest` goes further and re-runs the heads over a known sequence to check they reproduce recorded values — worth doing, and it needs all four heads installed.

## Inputs and outputs

**Input.** Per-residue ESM-2 embeddings. Sequences longer than the embedding model's limit are covered by overlapping windows.

**Output.** Relative solvent accessibility per residue: the residue's accessible area over the theoretical maximum for that amino acid, so roughly 0 in the core and roughly 1 fully exposed. The output is squashed into [0, 1.25], because a negative area is not a thing and the theoretical maxima are occasionally exceeded.

## Results

| Metric | Validation | Test | Baseline (test) |
|---|---|---|---|
| MAE | 0.102 | 0.115 | 0.190 |
| Pearson r | 0.767 | 0.721 | — |
| RMSE | 0.141 | 0.154 | 0.222 |
| Spearman ρ | 0.782 | 0.737 | — |

A — in the baseline column means the metric is undefined for that baseline, not unmeasured. This baseline predicts the same value for every input, so it has no correlation with the truth.

Baseline: Predicting the training-set mean accessibility for every residue.

Validation chose the epoch that was kept (epoch 4 of 14 run). The test split was evaluated once, afterwards. Splits come from the dataset builder, which clusters by sequence identity so near-identical proteins cannot appear on both sides.

### Against published benchmarks

None. The external benchmarks for this task come from sources whose terms do not permit redistribution, or from services that have been retired, so there is no set this head can be scored on alongside other published methods. The figures above are therefore not evidence of where it stands relative to the field — only of how it does on data drawn the way its training data was drawn.

## Training

| | |
|---|---|
| Embedding model | `facebook/esm2_t30_150M_UR50D` (640 dimensions, frozen) |
| Optimiser | AdamW, learning rate 0.001 |
| Batch size | 256 |
| Seed | 20260101 |
| Dataset | [`dyber-pqc/opendna-burial`](https://huggingface.co/datasets/dyber-pqc/opendna-burial) |
| Data licence | CC BY 4.0 (CATH v4.4.0); the coordinates it is computed from are CC0 1.0 (wwPDB) |

### Licence of these weights

CATH is CC BY 4.0 and wwPDB coordinates are CC0, neither of which places a condition on the licence of derived work, so the weights are Apache-2.0. Cite the dataset.

### Reproducing this run

Download the dataset from [`dyber-pqc/opendna-burial`](https://huggingface.co/datasets/dyber-pqc/opendna-burial) into `datasets/out/burial/`, or rebuild it from its upstream sources with the `datasets/` workspace of the [OpenDNA repository](https://github.com/dyber-pqc/OpenDNA) — the dataset card lists every source file with its SHA-256. Then, from `training/`:

```bash
cargo run --release -p opendna-train -- embed --task burial
cargo run --release -p opendna-train -- train --task burial \
    --batch-size 256 --learning-rate 0.001 --seed 20260101
```

The hyperparameters are spelled out above even where they match the current defaults, so the command reproduces *this* run rather than whatever the defaults happen to be later.

Both commands take `--datasets <dir>` if the dataset is not at the default `../datasets/out`. Embedding is cached, so the second run of either command is cheap.

## Limitations

- Labels come from a single chain of a crystal or NMR model, so a residue that forms an inter-chain interface in an obligate oligomer is labelled exposed. The head learns that convention and will do the same.
- Trained on single-segment CATH domains whose sequence could be located in current coordinates, which excludes discontinuous domains — about one in eight of CATH — and entries substantially revised since CATH v4.4.0. The training distribution is biased toward simple, stable structures.
- Accessibility is a property of a fold. For a sequence with no stable structure the number is still produced and still means nothing; read it alongside the disorder head rather than alone.
- The metrics are pooled over residues, and residues within a domain are strongly correlated because neighbours share a local environment. The effective sample size is nearer the number of domains than the number of residues, so small differences between splits are noise.
- Trained on frozen embeddings, so accuracy is bounded by what ESM-2 already encodes about the position.

## Intended use

Ranking and triage during protein engineering: comparing candidate mutations or constructs against each other. It is not a substitute for measurement, and it should not be used to make claims about a specific protein without experimental confirmation.

## Citation

The training data comes from:

> Sillitoe, I. et al. CATH: increased structural coverage of functional space. Nucleic Acids Research 49, D266-D273 (2021).
>
> The labels are computed from wwPDB coordinates and normalised by the theoretical maxima of: Tien, M. Z. et al. Maximum allowed solvent accessibilites of residues in proteins. PLOS ONE 8, e80635 (2013).

Please cite that work as well as this model.

## How wrong the number is likely to be

A point estimate invites being read as a measurement. These intervals say how wide the range around it has to be, fitted by split conformal prediction on the validation split: the half-width is whatever the residuals require, with no assumption about the shape of the error.

| Confidence | Interval | Coverage achieved on held-out data |
|---|---|---|
| 50% | ±0.072 RSA | 43.3% |
| 80% | ±0.170 RSA | 76.0% |
| 90% | ±0.237 RSA | 87.5% |

Mean absolute error on the same split is 0.115 RSA.

The achieved column is there because it differs from the nominal one, and it can differ either way. A conformal guarantee assumes the calibration and held-out sets are exchangeable, and this project's splits are deliberately not: the dataset builder clusters by sequence identity so that near-identical proteins cannot straddle a split, which makes the held-out set systematically more distant than a random sample. That is the right way to build the splits, and the cost is that the guarantee is approximate. The measured column is the one to believe, in whichever direction it points. On this head it points the same way at every level: the intervals cover less than they say, so read each one as the lower bound of a range that ought to be wider.

It rests on fewer observations than the residual count suggests. Residuals from one domain are correlated — neighbouring residues share a local environment and whatever the head has got wrong about that fold — so the effective number of independent observations behind each quantile is nearer the number of domains in the split than the number of residues. The dataset card records both.

The half-width is constant, which is the method's main limitation: the interval is wider than it needs to be where the head is reliable and narrower than it should be where it is not.
