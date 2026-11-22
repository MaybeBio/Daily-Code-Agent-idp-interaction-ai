# Models

Model weights are never committed to this repository; they are downloaded or
built on demand.

This directory holds OpenDNA's own prediction heads, trained by the
[training workspace](https://github.com/dyber-pqc/OpenDNA/blob/HEAD/training) and described below.

It used to hold a second thing as well: `manifest.yaml`, a list of third-party
models — ESMFold, ProteinMPNN and others — that the Python backend downloaded
into `~/.opendna/models/`. That backend is gone, nothing reads the manifest,
and none of those models can be loaded by the Rust application, so the
manifest has been removed rather than left to imply they are available. The
only third-party model in use is ESM-2, fetched by `opendna-ml` into the
standard Hugging Face cache.

## OpenDNA prediction heads

Each is a small network on top of **frozen** ESM-2 embeddings — the language
model itself is never fine-tuned.

## What a trained head consists of

| File | Purpose |
|---|---|
| `<task>.safetensors` | the weights |
| `<task>.metadata.json` | which embedding model and width the head expects |
| `<task>.metrics.json` | validation and test metrics, and the baseline |
| `<task>.README.md` | the model card, generated from the metrics |

The metadata is checked at load time, so a head can never be paired with
embeddings it was not trained on — the application refuses rather than
returning numbers computed from the wrong input.

## Building them

```bash
cd training
cargo run --release -p opendna-train -- embed --task stability
cargo run --release -p opendna-train -- train --task stability
```

See [training/README.md](https://github.com/dyber-pqc/OpenDNA/blob/HEAD/training/README.md) for GPU setup and options.

## Published repositories

| Task | Model | Dataset |
|---|---|---|
| Stability (ΔΔG) | [`dyber-pqc/opendna-stability`](https://huggingface.co/dyber-pqc/opendna-stability) | [`datasets/dyber-pqc/opendna-stability`](https://huggingface.co/datasets/dyber-pqc/opendna-stability) |
| Solubility | [`dyber-pqc/opendna-solubility`](https://huggingface.co/dyber-pqc/opendna-solubility) | [`datasets/dyber-pqc/opendna-solubility`](https://huggingface.co/datasets/dyber-pqc/opendna-solubility) |
| Disorder | [`dyber-pqc/opendna-disorder`](https://huggingface.co/dyber-pqc/opendna-disorder) | [`datasets/dyber-pqc/opendna-disorder`](https://huggingface.co/datasets/dyber-pqc/opendna-disorder) |
| — (no head yet) | — | [`datasets/dyber-pqc/opendna-structures`](https://huggingface.co/datasets/dyber-pqc/opendna-structures) |

A model repository holds `<task>.safetensors`, `<task>.metadata.json`,
`<task>.metrics.json` and the card as `README.md`. The weights keep their
task-prefixed names so that all three repositories can be downloaded into one
directory and used as-is:

```bash
hf download dyber-pqc/opendna-stability --local-dir models
```

Then point the application at it with `OPENDNA_MODELS=models`.

### Licensing

The heads do not all carry the same licence, because their training data does
not. Stability and disorder are Apache-2.0: their sources are CC BY, which
asks for attribution and sets no condition on derived work. Solubility is
CC BY-SA 2.1 JP, following eSOL — whether weights are a derivative work of
their training data is unsettled, and a share-alike condition is not something
to assume away in our own favour. Each card states its own licence and why.

## Publishing

Publishing is a manual step, done by a maintainer. Nothing in this repository
uploads automatically, and no token is stored here.

`scripts/publish_hf.py` does steps 3 and 4 for all nine repositories, after
checking each head's weights against the SHA-256 in its own metadata and each
card's front matter. It does nothing without a flag:

```bash
python scripts/publish_hf.py             # check, and print what it would send
python scripts/publish_hf.py --upload    # send the cards, weights and manifests
python scripts/publish_hf.py --public    # make the repositories public
```

Upload before making anything public. A repository holds whatever was uploaded
last, so opening one whose card predates the current one publishes a superseded
document. The steps below are what the script does and why; read them before
running it, and step 5 is still yours.

1. Log in with your own account:

   ```bash
   hf auth login
   ```

2. Review the generated card. The metrics table, hyperparameters and kept
   epoch come from the run, but the intended-use and limitations sections are
   written by hand and should be re-read before each release — particularly
   if the training data changed. Re-render a card without retraining with:

   ```bash
   cargo run --release -p opendna-train -- card --task stability --out ../models
   ```

3. Publish the dataset that produced the model first, or at least alongside
   it, so the card's dataset link resolves. A dataset directory from
   `datasets/out/` uploads whole:

   ```bash
   hf repos create dyber-pqc/opendna-stability --repo-type dataset
   hf upload dyber-pqc/opendna-stability datasets/out/stability . --repo-type dataset
   ```

4. Create the model repository and upload the four files, with the card as
   `README.md`:

   ```bash
   hf repos create dyber-pqc/opendna-stability --repo-type model
   hf upload dyber-pqc/opendna-stability models/stability.safetensors
   hf upload dyber-pqc/opendna-stability models/stability.metadata.json
   hf upload dyber-pqc/opendna-stability models/stability.metrics.json
   hf upload dyber-pqc/opendna-stability models/stability.README.md README.md
   ```

5. Check what landed rather than trusting the upload summary. Download into an
   empty directory, verify each Parquet file against the `sha256` in the
   dataset's own `manifest.json`, and run one prediction from the downloaded
   weights:

   ```bash
   cargo run --release -p opendna-ml --example predict_stability -- <dir> <sequence> 8
   ```

   The numbers must match a local run exactly. They are the same weights; any
   difference means the upload is not what was trained.

## Honesty requirements

These apply to anything published from here:

- **Report the baseline.** Every card states what the model is compared
  against — the language model's zero-shot signal, the training mean, or the
  majority class. A model that does not beat its baseline is not published.
- **Test once.** The validation split chooses the epoch; the test split is
  evaluated afterwards, once. Nothing is tuned against test numbers.
- **Say what the data does not cover.** The limitation sections exist because
  each of these datasets is narrower than the task name suggests: stability
  is small single domains, solubility is one organism's expression assay, and
  disorder is annotated rather than exhaustive.
- **No silent fallbacks.** If a head is missing, the application says the
  prediction is unavailable. It never substitutes a heuristic and presents
  the result as a prediction.
