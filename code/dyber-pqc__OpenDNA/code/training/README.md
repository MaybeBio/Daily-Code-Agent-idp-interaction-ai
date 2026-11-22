# Training

Builds the embedding caches and trains OpenDNA's prediction heads. A separate
workspace, so the library crates and the desktop app do not pull in Arrow,
Parquet or the CUDA toolchain.

The heads themselves are defined in [`opendna-ml`](https://github.com/dyber-pqc/OpenDNA/blob/HEAD/crates/opendna-ml), not
here, so training and inference share one definition and cannot drift apart.

## What gets trained

Each head is a small network on top of **frozen** ESM-2 embeddings. The
language model is never fine-tuned, which is what keeps a run to minutes
rather than days and makes it reproducible from a seed.

| Head | Input | Output |
|---|---|---|
| `stability` | wild-type embedding at the mutated position, one-hot wild-type and mutant residues, and the model's masked log-odds for the substitution | ΔΔG |
| `solubility` | mean-pooled sequence embedding | solubility percentage |
| `disorder` | per-residue embeddings, then one convolution along the chain | per-residue logit |

## Prerequisites

1. Build the datasets first — see [`datasets/`](https://github.com/dyber-pqc/OpenDNA/blob/HEAD/datasets). The default
   `--datasets ../datasets/out` expects them there.
2. For GPU training, a CUDA toolkit. On Windows use
   [`scripts/cargo-cuda.ps1`](https://github.com/dyber-pqc/OpenDNA/blob/HEAD/scripts/cargo-cuda.ps1) in place of `cargo`;
   it picks a toolkit the driver supports and puts MSVC and the CUDA runtime
   on PATH. Without the `cuda` feature everything runs on the CPU, which is
   fine for `stability` and `solubility` but slow for `disorder`.

## Running

Embedding is the slow step and is cached, so it is done once per dataset and
model. Re-running skips sequences that are already cached, which also means
an interrupted run can simply be restarted.

```bash
cargo run --release -p opendna-train --features cuda -- embed --task stability
cargo run --release -p opendna-train --features cuda -- train --task stability
```

On Windows:

```powershell
..\scripts\cargo-cuda.ps1 run --release -p opendna-train --features cuda '--' embed --task stability
```

Useful flags: `--embedder` chooses the ESM-2 checkpoint (default
`esm2-large`, 150M parameters), `--epochs`, `--batch-size`,
`--learning-rate`, `--patience` and `--seed` control the run, and `--out`
sets where weights go (default `../models`).

## Reproducibility

A run is reproducible from its seed: re-running `train` with the same
`--seed`, `--batch-size` and `--learning-rate` reproduces the weights byte
for byte. That was checked by retraining the stability head and comparing
SHA-256 against the committed file.

Two caveats worth knowing:

- The embedding caches must match too. They are keyed by sequence and by
  embedding model, so a cache built with a different checkpoint is rejected
  rather than silently reused.
- On a CUDA device, embedding is bit-reproducible only because this pipeline
  runs inference on one thread. Several threads sharing one model on a GPU
  perturb results at the level of f32 rounding; see the note on
  [`Esm2`](https://github.com/dyber-pqc/OpenDNA/blob/HEAD/crates/opendna-ml/src/esm2/mod.rs).

## Outputs

Per task, in `--out`:

- `<task>.safetensors` — the trained weights
- `<task>.metadata.json` — which embedding model and width the head expects,
  so inference cannot pair a head with the wrong embeddings
- `<task>.metrics.json` — validation and test metrics, and the baseline

## How the numbers are produced

- The **validation** split chooses when to stop and which epoch to keep. The
  **test** split is read once, at the end. No tuning against test.
- Every report includes a **baseline** that a useful model has to beat:
  the language model's own log-odds for stability, predicting the training
  mean for solubility, and always answering "disordered" for disorder.
- Splits come from the dataset builder, which clusters by sequence identity
  so that near-identical proteins cannot land on both sides.
- Metrics are implemented in [`metrics.rs`](https://github.com/dyber-pqc/OpenDNA/blob/HEAD/training/opendna-train/src/metrics.rs) and
  tested against worked examples rather than trusted.
