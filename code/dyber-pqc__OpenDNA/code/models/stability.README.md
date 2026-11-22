---
license: apache-2.0
library_name: opendna
tags:
  - protein
  - biology
  - protein-language-model
  - stability
  - ddg
  - mutation-effect
  - variant-effect
  - regression
datasets:
  - dyber-pqc/opendna-stability
base_model: facebook/esm2_t30_150M_UR50D
---

# OpenDNA stability head

Predicts the change in folding free energy (ΔΔG) caused by a single amino-acid substitution.

This is a small head trained on **frozen** [facebook/esm2_t30_150M_UR50D](https://huggingface.co/facebook/esm2_t30_150M_UR50D) embeddings. It is not a fine-tuned language model: the embedding model is unchanged, and only the head below was trained.

## Using it

These are not transformers weights and no `AutoModel` will load them. They are read by [OpenDNA](https://github.com/dyber-pqc/OpenDNA), which is a single binary with no runtime to install — take it from the [releases page](https://github.com/dyber-pqc/OpenDNA/releases) or build it with `cargo build --release -p opendna-cli`.

```bash
hf download dyber-pqc/opendna-stability --local-dir models
export OPENDNA_MODELS=$PWD/models

# confirm the download is intact: the weights are hashed against the
# digest recorded in their metadata
opendna models

# a CSV of `name,mutation` rows, for example `ubiquitin,L8P`
opendna predict --task stability -i candidates.fasta \
    --mutations muts.csv -o stability.csv
```

The output carries the head version, the embedding model, the device and the limitations below as `#` comments, and a `stability.manifest.json` beside it records the digests of what went in and what came out. The first run downloads the embedding model; after that `opendna stage` makes a directory that works with no network at all.

`opendna models` checks that the file is the one it claims to be. `opendna selftest` goes further and re-runs the heads over a known sequence to check they reproduce recorded values — worth doing, and it needs all four heads installed.

## Inputs and outputs

**Input.** The wild-type sequence's ESM-2 embedding at the mutated position, one-hot encodings of the wild-type and mutant residues, and the masked log-odds `log P(mutant) - log P(wild type)` at that position.

**Output.** A single ΔΔG value in kcal/mol, following the dataset's sign convention: negative is destabilising.

## Results

| Metric | Validation | Test | Baseline (test) |
|---|---|---|---|
| MAE | 0.605 | 0.532 | — |
| Pearson r | 0.678 | 0.717 | 0.497 |
| RMSE | 0.802 | 0.699 | — |
| Spearman ρ | 0.674 | 0.674 | 0.445 |

A — in the baseline column means the metric is undefined for that baseline, not unmeasured. This baseline is a score rather than a prediction in the same units as the measurements, so an absolute or squared error against them is not defined for it — only how well it orders them.

Baseline: The language model's masked log-odds used directly, with no trained head.

Validation chose the epoch that was kept (epoch 1 of 11 run). The test split was evaluated once, afterwards. Splits come from the dataset builder, which clusters by sequence identity so near-identical proteins cannot appear on both sides.

### Against published benchmarks

**proteingym-stability-out-of-domain** — 13747 items scored.

| Metric | Value |
|---|---:|
| spearman::ESTA_BACSU_Nutschel_2020 | 0.3163 |
| spearman::NUD15_HUMAN_Suiter_2020 | 0.5351 |
| spearman::PTEN_HUMAN_Matreyek_2021 | 0.4253 |
| spearman::TPMT_HUMAN_Matreyek_2018 | 0.4911 |
| spearman::mean_over_assays | 0.4419 |

- These assays measure variant abundance or thermostability in whole proteins, not ΔΔG by calorimetry. They are correlates of folding stability, so a positive Spearman means the head orders destabilising substitutions correctly; the magnitude is not comparable to a ΔΔG-against-ΔΔG correlation.
- Every ProteinGym assay derived from Tsuboyama 2023 is excluded. Those 64 assays are the head's own training data, and scoring them would measure memorisation.
- ProteinGym scores are higher-is-fitter, and the head predicts higher-is-more-stable, so the expected sign is positive. A negative correlation here would mean the head is ordering mutations backwards on this protein, not that the convention is flipped.
- 4 assays, 13747 mutants scored. This is a handful of proteins, so the mean across assays has no useful confidence interval; read the per-assay numbers.
- ESTA_BACSU_Nutschel_2020: thermostability in a 212-residue protein, 2172 single substitutions.
- NUD15_HUMAN_Suiter_2020: VAMP-seq, drug sensitivity in a 164-residue protein, 2844 single substitutions.
- PTEN_HUMAN_Matreyek_2021: Protein stability in a 403-residue protein, 5083 single substitutions.
- TPMT_HUMAN_Matreyek_2018: Protein stability in a 245-residue protein, 3648 single substitutions.

## Training

| | |
|---|---|
| Embedding model | `facebook/esm2_t30_150M_UR50D` (640 dimensions, frozen) |
| Optimiser | AdamW, learning rate 0.001 |
| Batch size | 512 |
| Seed | 20260101 |
| Dataset | [`dyber-pqc/opendna-stability`](https://huggingface.co/datasets/dyber-pqc/opendna-stability) |
| Data licence | CC BY 4.0 (Tsuboyama et al. 2023) |

### Licence of these weights

The training data is CC BY 4.0, which asks for attribution and places no condition on the licence of derived work, so the weights are Apache-2.0. Cite the dataset.

### Reproducing this run

Download the dataset from [`dyber-pqc/opendna-stability`](https://huggingface.co/datasets/dyber-pqc/opendna-stability) into `datasets/out/stability/`, or rebuild it from its upstream sources with the `datasets/` workspace of the [OpenDNA repository](https://github.com/dyber-pqc/OpenDNA) — the dataset card lists every source file with its SHA-256. Then, from `training/`:

```bash
cargo run --release -p opendna-train -- embed --task stability
cargo run --release -p opendna-train -- train --task stability \
    --batch-size 512 --learning-rate 0.001 --seed 20260101
```

The hyperparameters are spelled out above even where they match the current defaults, so the command reproduces *this* run rather than whatever the defaults happen to be later.

Both commands take `--datasets <dir>` if the dataset is not at the default `../datasets/out`. Embedding is cached, so the second run of either command is cheap.

## Limitations

- Trained almost entirely on small, single-domain proteins of 32-72 residues. Behaviour on large multi-domain proteins, membrane proteins and complexes is untested.
- The source assay infers ΔΔG from protease resistance over a range of roughly -1 to 5 kcal/mol for ΔG. Mutations outside that range are absent from training, so strongly destabilising substitutions are under-represented.
- Only single substitutions. The model says nothing about combined mutations, insertions or deletions.
- Trained on frozen embeddings, so accuracy is bounded by what ESM-2 already encodes about the position.
- Trained with a squared-error loss, so predictions are compressed toward the training mean and span a much narrower range than the measurements do. The rank order is the usable signal, which is why Spearman is the headline metric; treat an individual predicted value as a weak estimate.

## Intended use

Ranking and triage during protein engineering: comparing candidate mutations or constructs against each other. It is not a substitute for measurement, and it should not be used to make claims about a specific protein without experimental confirmation.

## Citation

The training data comes from:

> Tsuboyama, K. et al. Mega-scale experimental analysis of protein folding stability in biology and design. Nature 620, 434-444 (2023).

Please cite that work as well as this model.

## How wrong the number is likely to be

A point estimate invites being read as a measurement. These intervals say how wide the range around it has to be, fitted by split conformal prediction on the validation split: the half-width is whatever the residuals require, with no assumption about the shape of the error.

| Confidence | Interval | Coverage achieved on held-out data |
|---|---|---|
| 50% | ±0.465 kcal/mol | 54.6% |
| 80% | ±0.970 kcal/mol | 85.2% |
| 90% | ±1.327 kcal/mol | 93.7% |

Mean absolute error on the same split is 0.532 kcal/mol.

The achieved column is there because it differs from the nominal one, and it can differ either way. A conformal guarantee assumes the calibration and held-out sets are exchangeable, and this project's splits are deliberately not: the dataset builder clusters by sequence identity so that near-identical proteins cannot straddle a split, which makes the held-out set systematically more distant than a random sample. That is the right way to build the splits, and the cost is that the guarantee is approximate. The measured column is the one to believe, in whichever direction it points. On this head it points the other way at every level: the intervals cover more than they say, so they are conservative rather than optimistic.

It rests on fewer observations than the residual count suggests. Residuals from one wild-type protein are correlated — substitutions in the same fold share whatever the head has got wrong about it — so the effective number of independent observations behind each quantile is nearer the number of distinct wild types in the split than the number of substitutions. The dataset card records both.

The half-width is constant, which is the method's main limitation: the interval is wider than it needs to be where the head is reliable and narrower than it should be where it is not.
