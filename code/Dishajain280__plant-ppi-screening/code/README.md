# plant-ppi-screening

[![tests](https://github.com/Dishajain280/plant-ppi-screening/actions/workflows/ci.yml/badge.svg)](https://github.com/Dishajain280/plant-ppi-screening/actions/workflows/ci.yml)
[![python](https://img.shields.io/badge/python-3.10%2B-blue)](https://www.python.org/)
[![license](https://img.shields.io/badge/licence-MIT-green)](LICENSE)

A reproducible harness for screening plant protein–protein interactions with
structure-prediction models. It takes a curated candidate list, resolves the
proteins, builds AlphaFold Server batch files within the free quota, parses the
returned confidence metrics, and scores the ranking against known positives and
honest negatives.

**Status: the harness is complete and tested (36 unit tests). The predictions have
not been run yet, so this repository contains no biological result.** The figure
below is synthetic demo output, generated to prove the pipeline works. What remains,
and how long it takes, is in [ROADMAP.md](https://github.com/Dishajain280/plant-ppi-screening/blob/HEAD/ROADMAP.md).

![Ranked screening scores, demo data](https://github.com/Dishajain280/plant-ppi-screening/raw/HEAD/demo/report/figure.svg)

Method, metric definitions and the honest limitations are in
[docs/method.md](https://github.com/Dishajain280/plant-ppi-screening/blob/HEAD/docs/method.md). Uploading jobs and downloading results step by
step is in [docs/running-predictions.md](https://github.com/Dishajain280/plant-ppi-screening/blob/HEAD/docs/running-predictions.md). CI runs the
test suite and this demo on every push, across Python 3.10–3.12.

## Why this exists

Interactions found by yeast two-hybrid or BiFC cost weeks per pair, so nobody can
screen a transcription factor family that way. Structure prediction reverses the
order of work: rank hundreds of pairs cheaply, then spend bench time on the top
few percent. The catch is that naive use of these models produces confident
nonsense on plant proteins — few plant structures are solved, so alignments are
thin, and transcription factors carry disordered regions that depress every
confidence score. This harness exists to make the ranking defensible rather than
decorative.

## Quickstart

```bash
python cli.py resolve            # gene symbols -> UniProt accessions + FASTA
python cli.py check-labels       # validate the label file before spending quota
python cli.py batch              # candidate pairs  -> batches/batch_001.json
# upload batches/*.json at https://alphafoldserver.com, download the outputs
python cli.py ingest --results results/
python cli.py evaluate --results results/scores.csv
python cli.py demo               # offline end-to-end check on synthetic scores
python -m unittest discover -s tests
```

No third-party packages: Python 3.10+ standard library only.

## The honesty guard

`data/candidates.csv` carries a `verified` column. `evaluate` **refuses to compute
metrics** while any included pair is marked `verified=no`, and tells you which
ones. A curated "positive" that nobody checked produces an impressive AUC and a
worthless claim, so the pipeline blocks that by default:

```
$ python cli.py evaluate
refusing to score: these pairs are marked unverified in candidates.csv:
  - VIP1_BZIP29
check each label against its source, set verified=yes, or pass --allow-unverified for a dry run
```

## Server constraints this accounts for

| Constraint | Value | Handling |
| --- | --- | --- |
| Jobs per day | 30 | `QuotaTracker` records submissions in `.state/quota.json` |
| Tokens per job | 5,000 | pairs above the cap are skipped with a reason, not silently dropped |
| Jobs per batch file | 100 | jobs chunked into `batch_001.json`, `batch_002.json`, … |
| Output layout | one folder per job: summary JSON, mmCIF, per-token data | `ingest` walks the tree for `*summary_confidences*.json` |

Full-length plant proteins routinely exceed the token cap, which is why
`range_a` / `range_b` in the candidate file accept domain boundaries (`1-150`).

## Reading the scores

- **ipTM** — confidence in the relative placement of the two chains. Above ~0.8 is
  a strong interface, below ~0.6 is usually a failed model, 0.6–0.8 is the grey zone.
- **min chain-pair PAE** — smallest predicted alignment error across the interface.
  This is the metric to trust when a protein has disordered regions dragging its
  global confidence down.
- **pLDDT** — per-residue local confidence; necessary but not evidence of binding.

A full-length pair can score far below its isolated domains (published plant
example: BRI1–BAK1 at ipTM ≈ 0.38 full-length versus ≈ 0.90 for the interacting
domains, with a true negative near 0.2). `BRI1_BAK1_domains` in the candidate file
exists to reproduce that gap, which is the calibration step before trusting
anything else.

## Data schema

`data/candidates.csv`:

| Column | Meaning |
| --- | --- |
| `pair_id` | unique name; also the job name on the server |
| `gene_a`, `gene_b` | Arabidopsis gene symbols (resolved via UniProt) |
| `acc_a`, `acc_b` | optional explicit accessions, required when a symbol is ambiguous |
| `range_a`, `range_b` | optional 1-based inclusive domain boundaries |
| `label` | 1 for a known interactor, 0 for a negative |
| `label_type` | `positive_curated`, `negative_curated`, or `decoy_localization` |
| `source` | where the label comes from — fill this in or the label is not usable |
| `verified` | `yes` only after someone checked the source |

Two naming traps this repository documents, both found by running it:

- **VIP1** matches two Arabidopsis genes (a transcription factor, Q9MA75, and an
  inositol kinase, F4J8C6), so the accession is pinned explicitly.
- **AGB1** is stored in UniProt under the gene name `GB1` (P49177, At4g34460),
  which is why the resolver falls back from exact-name to synonym search.

## Repository layout

```
cli.py                   resolve | batch | ingest | evaluate | demo
ppi/uniprot.py           symbol resolution, FASTA cache, sequence slicing
ppi/batch.py             token budgeting, chunking, daily quota tracking
ppi/results.py           parse summary confidences into a scores table
ppi/metrics.py           AUC, recall in the top slice, precision, bootstrap CIs
ppi/curate.py            label validation: duplicates, placeholders, decoy share
ppi/report.py            markdown table and dependency-free SVG figure
data/candidates.csv      curated pairs with sources and verification flags
docs/method.md           protocol, metric definitions, limitations
docs/curating-labels.md  seven-step protocol for building a verified label set
docs/running-predictions.md  upload, quota, download and scoring runbook
tests/                   36 offline tests, no network required
.github/workflows/ci.yml tests plus the demo, on every push
CITATION.cff             citation metadata (GitHub "Cite this repository")
```

## Limits

- **Manual upload.** The AlphaFold Server has no public API, so the pipeline stops
  at a JSON file you submit by hand. That is a terms-of-use constraint, not a gap
  in the code.
- **Labels are the hard part.** Predicting pairs is easy; knowing which pairs truly
  interact is not. Curated positives from the literature and decoys built from
  subcellular localisation are honest; "random pairs assumed not to interact" are not.
- **No negative experimental validation.** A high ipTM is a hypothesis. The intended
  loop is prediction → Y2H or BiFC in the lab → updated labels.

## Next steps

1. Curate 30–60 labelled pairs for one protein family, with a source per row.
2. Calibrate on the BRI1–BAK1 full-length versus domain gap.
3. Run at 30 pairs/day and record every batch submitted.
4. Report AUC with a bootstrap interval, plus recall at a fixed precision.

## Contact

Disha Pirodiya · dishapirodiya@gmail.com · [github.com/Dishajain280](https://github.com/Dishajain280)
