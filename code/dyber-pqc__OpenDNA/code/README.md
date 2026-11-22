# OpenDNA

Protein property prediction that runs on your own machine. No server, no
account, no network: sequences stay where you opened them.

Four trained models score candidate sequences for per-residue disorder,
per-residue solvent accessibility (`burial`), solubility, and the stability
effect of a point mutation. Each is a small head
on frozen [ESM-2](https://huggingface.co/facebook/esm2_t30_150M_UR50D)
embeddings. When a model is missing the answer is an error saying so — there
is no heuristic fallback anywhere, and no number in the interface comes from
anything but a trained model.

## Who this is for

**If you have candidate sequences you cannot upload, and you want to triage
them before committing bench time**, this is built for you. Of the commercial
protein-design platforms, essentially all are cloud-only; the reason to use
OpenDNA is usually not that it predicts better but that it predicts at all
inside an environment where sending unpublished constructs to a vendor is not
allowed.

Concretely, it is useful today if you are:

- **Working on intrinsic disorder.** The disorder head is competitive with the
  CAID3 field (below), runs over a whole proteome from one command, and costs
  nothing.
- **Screening constructs under IP or data-governance constraints.** Nothing
  leaves the machine, so there is no review to pass.
- **Somewhere that has to show its work.** Every prediction carries the model
  version, the embedding model, and the limitations from that model's card;
  every dataset file is checksummed; the published models and this repository
  provably describe the same build.

It is **not** yet useful if you need docking, molecular dynamics, sequence
design, a trained antibody model, a signed installer, or a support contract.
Those are either absent or unfinished, and the section below says which.
Structures it will predict, but only by driving an ESMFold or ColabFold you
installed yourself. For antibodies it will tell you how often each 9-residue
window of a heavy chain occurs in human repertoires, which is a count over a
published table rather than a model — see below.

## Quick start

Take `opendna` for your platform from the
[releases page](https://github.com/dyber-pqc/OpenDNA/releases), check it
against the `SHA256SUMS` published beside it, and put it on your path. There
is nothing to install and no runtime to bring along. To build it instead, see
[Build](#build).

Predictions need the trained heads, which are not bundled:

```bash
hf download dyber-pqc/opendna-disorder --local-dir models
export OPENDNA_MODELS=$PWD/models
```

Check that this machine reproduces the numbers the heads were published with:

```bash
opendna selftest
```

Then score a whole FASTA file:

```bash
opendna predict --task disorder --input candidates.fasta --output disorder.csv
```

`--task` takes `disorder`, `burial`, `solubility` or `stability`; `--format`
takes `csv`, `tsv`, `json`, or `caid` for a CAID submission. A run that will
take hours can be made interruptible with `--checkpoint`. See
[crates/opendna-cli](https://github.com/dyber-pqc/OpenDNA/blob/HEAD/crates/opendna-cli) for the command, and
[CAID.md](https://github.com/dyber-pqc/OpenDNA/blob/HEAD/crates/opendna-cli/CAID.md) for the assessment format.

The desktop application is the same models with a structure viewer; see
[Build](#build).

From a notebook instead, with the provenance attached to the frame rather than
left behind in a CSV:

```python
import opendna
run = opendna.predict("disorder", {"ubiquitin": sequence})
run.frame.head()
run.cite()
```

See [bindings/python](https://github.com/dyber-pqc/OpenDNA/blob/HEAD/bindings/python). What changed between releases, and what
each one does not cover, is in [docs/CHANGELOG.md](https://github.com/dyber-pqc/OpenDNA/blob/HEAD/docs/CHANGELOG.md).

## How good the models are

Two kinds of number, and the distinction matters more than the values.
**In-domain** is our own held-out test split: how a head does on data drawn
the way its training data was drawn. **External** is a reference set that
other published methods also report on, which is the only kind that says where
a head stands relative to the field.

| Head | In-domain test | External |
|---|---|---|
| **Disorder** | AUROC 0.883 | **AUROC 0.941** on CAID3 Disorder-PDB |
| Accessibility | Spearman 0.737 | not benchmarked externally yet |
| Stability (ΔΔG) | Spearman 0.674 | Spearman 0.442 over 4 ProteinGym assays |
| Solubility | Spearman 0.701 | none exists under usable terms |

**Disorder is the strong one.** Published CAID3 Disorder-PDB figures include
AIUPred v3 at 0.925 and MUSE-IDR at 0.954, so 0.941 sits inside that field.
Two caveats, both ours: this is a reimplementation of the CAID protocol
against their published reference file, not an official CAID entry; and 34 of
the 319 targets exceed ESM-2's length limit and were covered by overlapping
windows, which no published entrant does. The number was independently
recomputed in numpy at 0.9409 by `scripts/validate_caid.py`, so it is not an
artefact of one implementation.

It scores *higher* externally than in-domain, which is not an improvement but
a different target: our split labels unannotated residues ordered, because
DisProt records only demonstrated disorder, so the negative class contains
disordered residues nobody has shown yet. CAID masks what it did not assess.

**Stability is a ranking signal, not a measurement.** It falls from 0.674 to
0.442 out of domain because it trains on single domains of 32–74 residues and
the external assays are whole proteins of 164–403. The in-domain figure
overstates generalisation by roughly 0.23 Spearman. All four assays correlate
positively, so the ordering is right; treat an individual ΔΔG as weak.

**Solubility has no external number** because the benchmarks for it come from
sources whose terms forbid redistribution, or from services that have been
retired. Saying so beats quoting the in-domain figure as though it were
comparable.

Both suites are reproducible with `opendna-train bench`; the reports and the
reasoning are in [training/benchmarks](https://github.com/dyber-pqc/OpenDNA/blob/HEAD/training/benchmarks). A bigger
embedding model was tried and rejected: ESM-2 650M made the disorder head
*worse* on every external metric, for 2.61 GB of weights against 0.60 GB.

## What it does not do

Earlier versions of this README advertised all of these. None of them is
implemented here:

- docking, molecular dynamics, or free-energy calculation
- sequence or backbone design
- any trained antibody model, light-chain or paired-chain scoring, or
  developability beyond humanness
- any HTTP API or hosted service

Two of those have changed shape rather than arriving. OpenDNA still folds
nothing and designs nothing itself; what it has is an adapter that finds,
identifies and pins external programs, and one command that uses it:

- **`opendna fold`** predicts structures by driving an ESMFold or ColabFold
  that you installed. With neither installed it exits 3 naming what it looked
  for. It checks that the chain each engine returns is the sequence it was
  sent, and records which engine produced each structure. See
  [the CLI reference](https://github.com/dyber-pqc/OpenDNA/blob/HEAD/crates/opendna-cli/README.md#folding-a-sequence).
- **`opendna engines`** says which external programs are present, which
  version each reports, and pins them in an `engines.lock`. Docking and
  molecular dynamics are declared there and not yet driven.

The heads call a notebook through [Python bindings](https://github.com/dyber-pqc/OpenDNA/blob/HEAD/bindings/python), which is
the only way OpenDNA is reachable from outside a terminal. There is still no
server. [The roadmap](https://github.com/dyber-pqc/OpenDNA/blob/HEAD/docs/ROADMAP.md) says what is next and in what order.

## Datasets

Built by [datasets/](https://github.com/dyber-pqc/OpenDNA/blob/HEAD/datasets) from pinned upstream sources, every file
checksummed, with a licence review recorded per source — including the sources
**excluded** because their terms forbid redistribution.

| Dataset | Rows | Source |
|---|---:|---|
| stability | 389,068 | Tsuboyama et al. 2023 |
| solubility | 3,127 | eSOL, UniProtKB |
| disorder | 3,186 | DisProt 2025_12 |
| structures | 106,910 | CATH v4.4.0, wwPDB |

Splits are by 30% sequence-identity cluster, so no close homolog spans two
splits. Validation chooses the epoch; the test split is scored once,
afterwards; external benchmarks are never tuned against at all. The published
datasets and heads are listed in [models/README.md](https://github.com/dyber-pqc/OpenDNA/blob/HEAD/models/README.md).

## Build

Rust 1.95 or newer. `unsafe` is forbidden across the workspace.

```bash
cargo test --workspace
cargo run --release -p opendna-cli -- models
```

The desktop application is a separate workspace:

```bash
cd app
rustup target add wasm32-unknown-unknown
cargo install trunk
cargo run -p xtask -- fetch-molstar
cd ui && trunk build --release && cd ..
cargo run -p opendna-desktop
```

GPU builds need a `cuda` feature and a toolkit candle supports; on Windows use
`scripts/cargo-cuda.ps1`, which loads nvcc and MSVC and documents the version
bounds and the PowerShell argument traps.

## Layout

| Directory | Workspace | Contents |
|---|---|---|
| [crates/](https://github.com/dyber-pqc/OpenDNA/blob/HEAD/crates) | root | core types, I/O, analysis, storage, ESM-2, the heads, the CLI |
| [app/](https://github.com/dyber-pqc/OpenDNA/blob/HEAD/app) | `app` | Tauri desktop shell, Leptos UI, IPC protocol |
| [datasets/](https://github.com/dyber-pqc/OpenDNA/blob/HEAD/datasets) | `datasets` | the dataset builder |
| [training/](https://github.com/dyber-pqc/OpenDNA/blob/HEAD/training) | `training` | embedding cache, training, model cards, benchmarks |
| [pipelines/](https://github.com/dyber-pqc/OpenDNA/blob/HEAD/pipelines) | — | Nextflow and Snakemake modules for running parts in parallel |
| [bindings/](https://github.com/dyber-pqc/OpenDNA/blob/HEAD/bindings/python) | — | Python package, for calling the heads from a notebook |
| [docs/](https://github.com/dyber-pqc/OpenDNA/blob/HEAD/docs) | — | the documentation site, including [VALIDATION.md](https://github.com/dyber-pqc/OpenDNA/blob/HEAD/docs/VALIDATION.md) |
| [models/](https://github.com/dyber-pqc/OpenDNA/blob/HEAD/models) | — | trained heads and their cards; weights are not committed |

Each workspace is built, linted with `clippy -D warnings`, and tested in CI.
[DEVELOPMENT.md](https://github.com/dyber-pqc/OpenDNA/blob/HEAD/DEVELOPMENT.md) has the ground rules a change is held to.

## Licence

The code is [Apache-2.0](https://github.com/dyber-pqc/OpenDNA/blob/HEAD/LICENSE), which is what the `LICENSE` file and every
crate's metadata say. An earlier README claimed "Apache 2.0 + Commons Clause",
a different and more restrictive grant that appears nowhere in the repository;
it was removed rather than acted on.

**Model weights are not all under the same licence, because their training
data is not.** The stability and disorder heads are Apache-2.0: their sources
are CC BY, which asks for attribution and sets no condition on derived work.
The solubility head is **CC BY-SA 2.1 JP**, following eSOL — whether weights
are a derivative work of their training data is unsettled, and a share-alike
condition is not something to assume away in our own favour. If you are
embedding OpenDNA commercially, that distinction is the one to take to
counsel. Each model card states its own licence and the reasoning.
