---
license: apache-2.0
library_name: opendna
pipeline_tag: token-classification
tags:
  - protein
  - biology
  - protein-language-model
  - intrinsic-disorder
  - idp
  - per-residue
  - token-classification
datasets:
  - dyber-pqc/opendna-disorder
base_model: facebook/esm2_t30_150M_UR50D
---

# OpenDNA disorder head

Predicts, for each residue, whether it lies in an intrinsically disordered region.

This is a small head trained on **frozen** [facebook/esm2_t30_150M_UR50D](https://huggingface.co/facebook/esm2_t30_150M_UR50D) embeddings. It is not a fine-tuned language model: the embedding model is unchanged, and only the head below was trained.

## Using it

These are not transformers weights and no `AutoModel` will load them. They are read by [OpenDNA](https://github.com/dyber-pqc/OpenDNA), which is a single binary with no runtime to install — take it from the [releases page](https://github.com/dyber-pqc/OpenDNA/releases) or build it with `cargo build --release -p opendna-cli`.

```bash
hf download dyber-pqc/opendna-disorder --local-dir models
export OPENDNA_MODELS=$PWD/models

# confirm the download is intact: the weights are hashed against the
# digest recorded in their metadata
opendna models

opendna predict --task disorder -i candidates.fasta -o disorder.csv
```

The output carries the head version, the embedding model, the device and the limitations below as `#` comments, and a `disorder.manifest.json` beside it records the digests of what went in and what came out. The first run downloads the embedding model; after that `opendna stage` makes a directory that works with no network at all.

`opendna models` checks that the file is the one it claims to be. `opendna selftest` goes further and re-runs the heads over a known sequence to check they reproduce recorded values — worth doing, and it needs all four heads installed.

For a CAID submission, add `--format caid`.

## Inputs and outputs

**Input.** Per-residue ESM-2 embeddings. Sequences longer than the embedding model's limit are covered by overlapping windows.

**Output.** One calibrated probability per residue. The threshold for a disorder call is recorded in the head's metadata and stated under Calibration below, because it is a property of the trained head rather than a convention.

## Results

| Metric | Validation | Test | Baseline (test) |
|---|---|---|---|
| AUPRC | 0.540 | 0.611 | 0.179 |
| AUROC | 0.840 | 0.883 | 0.500 |
| Balanced accuracy | 0.772 | 0.812 | 0.500 |
| F1 | 0.535 | 0.573 | 0.304 |
| MCC | 0.434 | 0.493 | — |

A — in the baseline column means the metric is undefined for that baseline, not unmeasured. This baseline is a score rather than a prediction in the same units as the measurements, so an absolute or squared error against them is not defined for it — only how well it orders them.

Baseline: Labelling every residue disordered.

F1, balanced accuracy and MCC above are computed at a probability of 0.5 on the head's raw output, which is how the training report records them. The head's own decision threshold is lower and is stated under Calibration, along with the Matthews correlation it achieves — a different and better number, measured a different way. AUROC and AUPRC are rank statistics and do not depend on either threshold.

Validation chose the epoch that was kept (epoch 5 of 15 run). The test split was evaluated once, afterwards. Splits come from the dataset builder, which clusters by sequence identity so near-identical proteins cannot appear on both sides.

### Against published benchmarks

**caid3-disorder-pdb** — 99239 items scored.

| Metric | Value |
|---|---:|
| auroc | 0.9410 |
| auprc | 0.8996 |
| mcc | 0.7096 |
| f1 | 0.7788 |
| balanced_accuracy | 0.8250 |
| auroc (mean over proteins) | 0.9252 |

- Residues the reference did not assess ('-') are excluded, which is the CAID protocol. 61183 of 160422 residues were excluded.
- 34 of 319 proteins are longer than the embedding model's 1024-residue limit and were covered by overlapping windows. No published CAID3 entrant uses this windowing, so those proteins are not scored the way another method would score them.
- The threshold for MCC, F1 and balanced accuracy is 0.34, taken from the head rather than tuned on this set: tuning it here would make the number incomparable to a published one and would be fitting to the test set. It therefore comes from a different distribution than this one -- a calibrated head's threshold is fitted on its own validation split, whose class balance and definition of disorder are not this reference set's -- so a threshold chosen here would score better. The rank metrics above are unaffected by any of this, which is why AUROC and AUPRC are the figures to compare.
- CAID3 rankings are computed by the organisers' own code on a server-side ground truth. This is our reimplementation of the protocol against the published reference file, so a small difference from an official number is expected and this figure should not be presented as an official CAID3 result.

## Training

| | |
|---|---|
| Embedding model | `facebook/esm2_t30_150M_UR50D` (640 dimensions, frozen) |
| Optimiser | AdamW, learning rate 0.001 |
| Batch size | 256 |
| Seed | 20260101 |
| Dataset | [`dyber-pqc/opendna-disorder`](https://huggingface.co/datasets/dyber-pqc/opendna-disorder) |
| Data licence | CC BY 4.0 (DisProt 2025_12); the CAID3 reference set used only to hold proteins out is CC BY 3.0 |

### Licence of these weights

Both sources ask for attribution and place no condition on the licence of derived work, so the weights are Apache-2.0. Cite the dataset.

### Reproducing this run

Download the dataset from [`dyber-pqc/opendna-disorder`](https://huggingface.co/datasets/dyber-pqc/opendna-disorder) into `datasets/out/disorder/`, or rebuild it from its upstream sources with the `datasets/` workspace of the [OpenDNA repository](https://github.com/dyber-pqc/OpenDNA) — the dataset card lists every source file with its SHA-256. Then, from `training/`:

```bash
cargo run --release -p opendna-train -- embed --task disorder
cargo run --release -p opendna-train -- train --task disorder \
    --batch-size 256 --learning-rate 0.001 --seed 20260101
```

The hyperparameters are spelled out above even where they match the current defaults, so the command reproduces *this* run rather than whatever the defaults happen to be later.

Both commands take `--datasets <dir>` if the dataset is not at the default `../datasets/out`. Embedding is cached, so the second run of either command is cheap.

## Limitations

- DisProt records disorder that has been experimentally demonstrated. Residues that are simply unannotated are treated as ordered during training, so the negative class contains an unknown number of undiscovered disordered residues.
- Annotation is biased toward proteins people chose to study, and toward longer disordered regions. Short linkers are under-represented.
- The CAID3 reference proteins are held out of training, and their whole identity clusters with them, which is what makes the CAID3 figure above a measurement rather than a leak. The results table beside it is this dataset's own test split and is not comparable to anybody's published number: our negative class contains undiscovered disorder and CAID masks what it did not assess, so the two are scoring different things.
- The metrics are pooled over residues, but residues within one protein are strongly correlated — neighbouring positions almost always share a label. The effective sample size is therefore closer to the few hundred proteins in each split than to the hundred thousand or more residues the counts show, so small differences between splits should not be read as real. On this run the test split scored higher than validation on every metric, which is what that noise looks like rather than evidence the model generalises unusually well.

## Intended use

Ranking and triage during protein engineering: comparing candidate mutations or constructs against each other. It is not a substitute for measurement, and it should not be used to make claims about a specific protein without experimental confirmation.

## Citation

The training data comes from:

> Aspromonte, M. C. et al. DisProt in 2024: improving function annotation of intrinsically disordered proteins. Nucleic Acids Research 52, D434-D441 (2024).

Please cite that work as well as this model.

## Calibration

The probabilities are corrected before they are reported. This head is trained with weighted cross-entropy so that it does not predict "ordered" everywhere, which helps the ranking and inflates the probabilities; a Platt map fitted on the validation split removes the inflation.

| | |
|---|---|
| Map | `p = sigmoid(0.7237 * logit + -1.5665)` |
| Expected calibration error, held-out split | 0.2167 before, 0.0321 after |
| Threshold for a disorder call | 0.34 |
| Matthews correlation at that threshold, held-out split | 0.5213 |

The map is monotonic, so it changes no ranking and no rank statistic: AUROC is identical before and after. What it changes is the meaning of the number, and therefore where the decision threshold belongs -- which is why the threshold is recorded here rather than left at the conventional 0.5. The threshold was chosen on the validation split, so the correlation above is an honest estimate rather than the best figure findable with the held-out labels in hand.

Calibration is against the labels. DisProt records disorder that has been demonstrated experimentally and unannotated residues are treated as ordered, so these probabilities are calibrated against the annotation rather than against nature. If the annotation under-reports disorder, and it does by an unknown amount, the calibrated probabilities inherit that.

The threshold is a property of a distribution, not of the head. It was chosen on the validation split of this dataset, and it does not suit every set: on the CAID3 reference above, MCC, F1 and balanced accuracy are lower at this threshold than they were at the conventional 0.5, because 0.34 calibrated is about 0.78 raw and that set prefers a lower cut. The rank metrics are identical either way. If you threshold these probabilities for a decision of your own, the right cut depends on your class balance and on what a false positive costs you; the calibration is what makes it possible to reason about that at all.

## Confidence, and whether to abstain

How far the calibrated probability sits from the decision threshold orders residues by how likely the call is to be wrong, so a caller can trade coverage for accuracy. Measured on the held-out split:

| Minimum confidence | Residues kept | Error among them |
|---|---|---|
| 0.00 | 100% | 15.7% |
| 0.05 | 89% | 12.0% |
| 0.10 | 80% | 9.3% |
| 0.15 | 73% | 7.8% |
| 0.20 | 66% | 6.5% |
| 0.25 | 57% | 5.3% |
| 0.30 | 31% | 4.8% |

Nothing abstains by default. The confidence is reported for every residue and the cut-off is the caller's, because the right one depends on what a wrong call costs them. `opendna predict --min-confidence` marks the rows that do not clear a cut-off rather than removing them: a table quietly missing rows looks exactly like a complete one.

This signal captures 55% of the error reduction an oracle shows is available -- an oracle being a rule that discards exactly the wrong calls. Half of what is there, which is worth having and is not a guarantee.

Measured and rejected: k-mer novelty. Ordering residues by this was no better than discarding them at random, so no rule is built on it. It is recorded because the obvious next idea is to try it.
