# External benchmarks

Each file here is the report one `opendna-train bench` run produced, written by
the run rather than by hand, and named `<suite>.<embedder>.json`.

These exist because the numbers in the model cards were unsituated. A test
Spearman of 0.674 or an AUROC of 0.883 says how a head does on data drawn the
way its training data was drawn, and nothing at all about how it compares to
published methods — because no other method reports a number on our splits.

## Does predicted burial help the stability head?

**No.** It makes it worse, in domain and out, and the point of asking was that
it might not.

The stability head's weakness is that it has no structural context: an
embedding at the mutated position, two one-hots, and a log-odds. Structure-
aware methods beat sequence-only ones precisely there, so the obvious move is
to train a solvent-accessibility head from deposited coordinates and feed its
prediction in at the mutated position.

| | plain | with burial | Δ |
|---|---:|---:|---:|
| **Out of domain** (ProteinGym, mean over 4 assays) | **0.4419** | 0.4122 | **−0.0297** |
| ESTA_BACSU_Nutschel_2020 | 0.3163 | 0.2897 | −0.0266 |
| NUD15_HUMAN_Suiter_2020 | 0.5351 | 0.4570 | −0.0782 |
| PTEN_HUMAN_Matreyek_2021 | 0.4253 | 0.3965 | −0.0287 |
| TPMT_HUMAN_Matreyek_2018 | 0.4911 | 0.5057 | +0.0146 |
| In domain, Spearman | 0.6742 | 0.6697 | −0.0044 |
| In domain, Pearson | 0.7166 | 0.6990 | −0.0176 |
| In domain, RMSE | 0.6994 | 0.7507 | +0.0514 |

Worse on three of four external assays and on every in-domain metric. Every
hyperparameter was held at the published run's value — batch 512, learning
rate 1e-3, seed 20260101 — so the feature is the only difference.

The reason was predictable and was written down before the experiment ran: a
burial value predicted from an ESM-2 embedding is a *deterministic function*
of input the head already receives. It cannot add information. The hope was
that distilling structural supervision into one scalar would help a small MLP
that has to learn from 389,068 mutations what the burial head learned from
deposited coordinates. Instead it added a parameter and a noisy feature to
overfit, on domains of 32–74 residues where predicted burial is least
reliable.

Burial is not irrelevant to stability — predicted burial alone correlates
0.394 Spearman with ΔΔG on the in-domain test split, which is most of the way
to the zero-shot log-odds baseline at 0.445. The signal is real. It is simply
already present in the embedding, which is the whole point.

**What would plausibly work instead is true structure**, not predicted burial:
accessibility computed from a real backbone for the protein being scored, the
way ThermoMPNN uses one. That needs a structure for the query, which the
application does not have and cannot predict. It is a real limitation rather
than a thing to work around with a proxy.

The published stability head is unchanged. The burial variant is kept here as
evidence, and the regression check below confirms the published head still
scores 0.4419 exactly after all the plumbing this experiment required.

### The published head is byte-for-byte unaffected

Supporting an optional feature meant touching the head's constructor, its
metadata and the benchmark path, any of which could have changed the model
that is actually released. Re-running the out-of-domain suite against
`models/` after all of it gives 0.3163, 0.5351, 0.4253, 0.4911 and a mean of
0.4419 — identical to the figures above it. `with_burial` defaults to false
in the metadata, so a head written before the field existed loads as the
plain variant, which is what the published one is.

## Is a bigger embedding model worth it?

Asked and answered: **mostly no**, which is not what we expected and is the
reason these suites were built before the experiment rather than after.

All three heads were retrained on ESM-2 650M with every hyperparameter held at
the value the 150M run used, so the embedder is the only variable. Reports are
in this directory as `*.esm2-650m.json`, with the in-domain metrics under
`in-domain/`.

| Suite | Metric | 150M | 650M | Δ |
|---|---|---:|---:|---:|
| CAID3 Disorder-PDB | AUROC | 0.9410 | 0.9384 | **−0.0025** |
| CAID3 Disorder-PDB | MCC | 0.7292 | 0.7018 | **−0.0274** |
| CAID3 Disorder-PDB | AUPRC | 0.8996 | 0.8916 | **−0.0080** |
| ProteinGym stability | Spearman, mean | 0.4419 | 0.4607 | +0.0187 |
| In domain: stability | Spearman | 0.6742 | 0.6964 | +0.0222 |
| In domain: solubility | Spearman | 0.7009 | 0.7105 | +0.0096 |
| In domain: disorder | AUROC | 0.8833 | 0.8836 | +0.0003 |

The disorder head gets **worse** on every external metric with an embedder
four times the size, and gains nothing in domain. The stability head improves
consistently, in domain and out. Solubility gains 0.01 Spearman in domain with
a worse MAE and has no external check at all.

The 650M checkpoint is 2.61 GB against 0.60 GB. Promoting only the head that
benefits would mean the application downloading and running two embedding
models — 3.2 GB of weights and two forward passes per prediction — to buy
0.019 external Spearman on one head whose honest external figure is 0.46
either way. For a local-first desktop application that is a bad trade, so
**all three published heads stay on 150M**.

The useful part of this is the negative result. "Bigger embeddings are a
near-free win" is the kind of claim that is obviously true until measured, and
measuring it cost one afternoon because the suites already existed. Had the
retrain come first, all three heads would have been promoted on the strength
of the in-domain numbers, and the disorder head would have quietly got worse.

## Results

Heads trained on frozen ESM-2 150M embeddings, which is what is published.

| Suite | Metric | Ours | In-domain equivalent |
|---|---|---:|---:|
| CAID3 Disorder-PDB | AUROC, pooled | 0.941 | 0.883 |
| CAID3 Disorder-PDB | AUROC, mean over proteins | 0.925 | — |
| ProteinGym stability (out of domain) | Spearman, mean over 4 assays | 0.442 | 0.674 |

Both differences are informative, in opposite directions.

**Disorder scores higher externally, and this is not an improvement.** Our
split labels every unannotated residue ordered, because DisProt records
disorder that has been demonstrated; the negative class therefore contains an
unknown number of disordered residues nobody has shown yet. CAID masks the
residues its assessment did not cover and scores only the rest. The two
numbers measure different things, and the CAID one is the comparable one.

For context, published CAID3 Disorder-PDB figures include AIUPred v3 at 0.825
on Disorder-NOX and 0.925 on Disorder-PDB, and MUSE-IDR at 0.954 on
Disorder-PDB. Those are self-reported by their authors against the official
assessment. Ours is a reimplementation of the protocol against the published
reference file, so it is not an official CAID3 result and a small difference
from one is expected.

**Stability falls from 0.674 to 0.442, and that gap is the real finding.** The
head is trained on single domains of 32–74 residues. The external assays are
whole proteins of 164–403 residues, measuring variant abundance and
thermostability. The in-domain figure overstates generalisation by roughly
0.23 Spearman. All four assays correlate positively, so the head orders
substitutions correctly even where the strength is modest, and the sign
agreement across four independent proteins is itself a check that the feature
pipeline is wired up correctly end to end.

**Solubility has no entry.** The external benchmarks for it come from sources
whose terms forbid redistribution, or from services that have been retired
(see the exclusions in `datasets/README.md`). An absent row is more honest
than quoting the in-domain number as though it were comparable.

## Running them

The data is **not** in this repository and is never built into a dataset. Some
of it carries no licence at all — `protddg-bench` has no `LICENSE` file and
GitHub reports `"license": null` — which permits evaluation but not
republication. That is the same reasoning that kept FireProtDB out of the
stability dataset.

CAID3 Disorder-PDB is fetched by the `datasets` workspace, which records its
licence and checksum, so it is already in that cache:

```bash
cargo run --release -p opendna-train -- bench \
    --suite caid3-disorder-pdb \
    --reference ../datasets/cache/caid3-disorder-pdb/disorder_pdb.fasta \
    --models ../models \
    --out benchmarks/caid3-disorder-pdb.esm2-150m.json
```

ProteinGym is MIT-licensed. Fetch the assay index and the per-assay score
files, then:

```bash
curl -o DMS_substitutions.csv \
    https://raw.githubusercontent.com/OATML-Markslab/ProteinGym/main/reference_files/DMS_substitutions.csv

cargo run --release -p opendna-train -- bench \
    --suite proteingym-stability \
    --reference DMS_substitutions.csv \
    --assays <dir of per-assay csv files> \
    --models ../models \
    --out benchmarks/proteingym-stability.esm2-150m.json
```

## What the suites refuse to do

- **No threshold tuning.** The disorder threshold is 0.5, the value training
  used. Choosing it against a benchmark would be fitting to the benchmark and
  would make the result incomparable to a published one.
- **No scoring of masked residues.** CAID ships a mask for a reason, and
  scoring the residues it excluded would inflate every metric.
- **No training data.** Every ProteinGym assay using `cDNA display proteolysis`
  is excluded — 64 of its 217 assays — because that is Tsuboyama 2023, which is
  the stability head's training set. Scoring them would measure memorisation.
- **No silent renumbering.** Every mutant label is checked against the
  wild-type residue the sequence actually has at that position, and a
  disagreement is an error rather than a warning. A benchmark run on a
  numbering convention that is off by one produces a plausible number and a
  worthless one.
- **No windowing where training did not window.** An assay longer than the
  embedding model's limit is skipped and named, because the stability head was
  trained on single-pass embeddings.

## Provenance

Each report records the embedding model, the head version, how many items were
scored, how many were excluded and why, and the caveats that apply. A report
whose `head_version` does not match the head you are holding describes a
different model.
