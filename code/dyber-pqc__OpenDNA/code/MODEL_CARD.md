# Model cards

Every model OpenDNA runs has a card. This page says where each one is and what
it is for; it does not restate their numbers.

That is deliberate. Each head's card is **generated from its own measured
training report** by `opendna-train`, and lives beside its weights, so the card
and the model cannot drift apart. A second hand-written copy of the metrics
here would be a copy that goes stale, and the version of this file before
`v0.6.0-rc1` is the proof: it described ESMFold, ESM-IF1 and ProteinMPNN, none
of which OpenDNA has run since the rewrite.

## OpenDNA's own heads

Four, each a small network on **frozen** ESM-2 embeddings. The language model
is never fine-tuned.

| Head | Predicts | Card |
|---|---|---|
| `disorder` | per-residue intrinsic disorder | [models/disorder.README.md](https://github.com/dyber-pqc/OpenDNA/blob/HEAD/models/disorder.README.md) |
| `burial` | per-residue relative solvent accessibility | [models/burial.README.md](https://github.com/dyber-pqc/OpenDNA/blob/HEAD/models/burial.README.md) |
| `stability` | ΔΔG of a point mutation, kcal/mol | [models/stability.README.md](https://github.com/dyber-pqc/OpenDNA/blob/HEAD/models/stability.README.md) |
| `solubility` | percent soluble expression | [models/solubility.README.md](https://github.com/dyber-pqc/OpenDNA/blob/HEAD/models/solubility.README.md) |

Each card states the head's architecture, the dataset it was trained on and
how that dataset was split, its validation and test metrics, the **baseline it
is compared against**, its external benchmark figure or the fact that it has
none, its calibration or lack of it, its licence, and its limitations. The same
cards are published with the weights on Hugging Face, and
[models/README.md](https://github.com/dyber-pqc/OpenDNA/blob/HEAD/models/README.md) describes what a trained head consists of.

Three things hold for all four:

- **A missing head is an error**, naming the files it looked for. Nothing falls
  back to a heuristic and reports the result as a prediction.
- **Every prediction carries its card's limitations** into the output file as
  comments, read from the card rather than written twice.
- **Every head has a build id** derived from its weights, its training
  configuration and its data, recorded in each run's manifest, so a number can
  be traced to the exact model that produced it.

[docs/VALIDATION.md](https://github.com/dyber-pqc/OpenDNA/blob/HEAD/docs/VALIDATION.md) is the page for somebody who has to
sign off on using these, and includes what has not been checked.

## The one third-party model

**ESM-2**, [`facebook/esm2_t30_150M_UR50D`](https://huggingface.co/facebook/esm2_t30_150M_UR50D)
— 150M parameters, 640-dimensional embeddings, used frozen as the input
representation for all four heads.

- **Developer**: Meta AI (FAIR)
- **Reference**: Lin, Z. et al. (2023). Evolutionary-scale prediction of
  atomic-level protein structure. *Science* 379, 1123–1130.
  [doi:10.1126/science.ade2574](https://doi.org/10.1126/science.ade2574)
- **Licence**: MIT
- **Limits that reach the user**: a 1,024-residue window, being its 1,026
  position embeddings less the two special tokens. Disorder and
  accessibility cover longer sequences with overlapping windows, taking each
  residue's embedding from the window it sits most centrally in. Solubility
  pools over a single pass and has no windowed equivalent, so it refuses a
  longer sequence rather than reporting a number computed differently from
  training.

ESM-2 650M was tried in place of it and rejected: it made the disorder head
worse on every external metric, for four times the weight size. The comparison
is in [training/benchmarks](https://github.com/dyber-pqc/OpenDNA/blob/HEAD/training/benchmarks).

## Models OpenDNA no longer runs

Earlier versions of this file carried cards for ESMFold, ESM-IF1 and
ProteinMPNN, which the removed Python backend downloaded. None of them is
loaded by the Rust application, and none of them is a model OpenDNA ships,
trains, versions or can make a card for.

ESMFold is nonetheless reachable again: `opendna fold` runs an installation of
it that you provide, as an external program. That is a different relationship
from the one this file describes, and deliberately so — OpenDNA does not have
its weights, cannot say what they are beyond the version string the program
prints, and makes no claim about their accuracy. What it records instead is
which program ran: its version, the digest of the file invoked, and how much
that digest identifies. ESM-IF1 and ProteinMPNN are not driven at all; see
[docs/ROADMAP.md](https://github.com/dyber-pqc/OpenDNA/blob/HEAD/docs/ROADMAP.md).

They are named here so somebody who remembers those cards gets an answer
rather than silence.

## Framework

The per-head cards follow
[Model Cards for Model Reporting](https://arxiv.org/abs/1810.03993) in the
sections they carry: intended use, training data, evaluation data, metrics
with a baseline, and limitations. What they add beyond that framework is the
build id and the digests, because a card that cannot be tied to a specific file
of weights describes a model rather than the one that produced your numbers.
