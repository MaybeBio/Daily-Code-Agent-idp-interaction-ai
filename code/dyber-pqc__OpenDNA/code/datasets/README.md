# OpenDNA datasets

This directory builds the training and evaluation datasets for the OpenDNA
prediction heads (stability, solubility, per-residue disorder), a
structure/fold dataset, antibody evaluation cohorts, and the peptide humanness
prevalence table. Its Cargo workspace is separate from the application
workspace, so the application never depends on arrow, parquet or an HTTP client.

| Dataset | Rows | Row is | Sources | License of the output | Split unit |
|---|---:|---|---|---|---|
| `structures` | 106,910 | a CATH domain | CATH v4.4.0 S100, wwPDB entries index | CC BY 4.0 | CATH topology merged with 30% identity clusters |
| `stability` | 389,068 | a single substitution | Tsuboyama et al. 2023 | CC BY 4.0 | 30% identity cluster of the wild type |
| `solubility` | 3,127 | an E. coli protein | eSOL (NBDC archive), UniProtKB | CC BY-SA 2.1 JP | 30% identity cluster |
| `disorder` | 3,186 | a protein, with per-residue labels | DisProt 2025_12, CAID3 reference (hold-out only) | CC BY 4.0 | 30% identity cluster; CAID3 clusters fixed to test |
| `burial` | 939 | a CATH domain, with per-residue labels | CATH v4.4.0 S100, wwPDB entries index | CC BY 4.0 | inherited from `structures`, so no topology or cluster spans two splits |
| `antibody` | 246,543 | a heavy-chain V region | Observed Antibody Space (54 data units) | CC BY 4.0 | none: four cohorts, not a split. See below |
| `humanness` | 1,096,088 peptides | a 9-mer and the number of the 141 human donors holding it | Observed Antibody Space (150 donor units) | CC BY 4.0 | none: a reference table, nothing is trained |

Row counts are from the committed manifests in `cards/`, which describe the
builds that produced them. The first five are published on the Hub and trained
the heads; `antibody` and `humanness` are built and their cards committed, but
not yet published.

The last two have no split, and that is a finding rather than an omission.
`antibody` is four evaluation cohorts: all 18 of its mouse donors share
identical rearrangements and form one indivisible group, because inbred
laboratory strains are not independent donors, so no grouped split can put
mouse rows on both sides. `humanness` is a single reference table of peptide
prevalence; a score computed from it is an exact count, not a prediction, so
there is nothing to hold out. Both cards explain this where a reader meets the
numbers, and the antibody module records the measurements behind it.

## Published builds

| Dataset | On Hugging Face |
|---|---|
| `structures` | [`dyber-pqc/opendna-structures`](https://huggingface.co/datasets/dyber-pqc/opendna-structures) |
| `stability` | [`dyber-pqc/opendna-stability`](https://huggingface.co/datasets/dyber-pqc/opendna-stability) |
| `solubility` | [`dyber-pqc/opendna-solubility`](https://huggingface.co/datasets/dyber-pqc/opendna-solubility) |
| `disorder` | [`dyber-pqc/opendna-disorder`](https://huggingface.co/datasets/dyber-pqc/opendna-disorder) |
| `burial` | [`dyber-pqc/opendna-burial`](https://huggingface.co/datasets/dyber-pqc/opendna-burial) |

`cards/<dataset>/manifest.json` records the SHA-256 of every Parquet file in
the published build, so what is on the Hub can be checked against what is in
git without rebuilding anything.

Each built dataset gets a Hugging Face dataset card (`README.md`) with sources,
licenses, construction steps, schema, split sizes, cleaning statistics, known
limitations and citations. The card is generated from the same tables that are
written to Parquet.

## Requirements

- Rust 1.95 or newer, the same minimum as the rest of the project.
- [MMseqs2](https://github.com/soedinglab/MMseqs2) on `PATH`, or pointed to with
  `MMSEQS=<path>` or `--mmseqs <path>`, for the five clustered datasets. The
  pipeline stops with an error if MMseqs2 is missing rather than falling back to
  a weaker method. `antibody` and `humanness` do not cluster and do not need it.
- About 1.2 GB of disk for the download cache of the first five datasets (most
  of it the Tsuboyama archive), plus 1.4 GB for the antibody data units and
  0.4 GB for the humanness donor units.

## Usage

Run from this directory:

```bash
cargo run --release -- sources            # list sources, licenses, exclusions
cargo run --release -- download all       # fetch and checksum all upstream files
cargo run --release -- verify             # re-hash the cache against sources.lock
cargo run --release -- build stability    # or: structures, solubility, disorder,
                                          #     burial, antibody, humanness, all
```

`humanness` has one more step, which nothing else does:

```bash
cargo run --release -- evaluate humanness   # needs `build antibody` first
```

It scores the four antibody cohorts against the table and writes
`out/humanness/EVALUATION.md`. Not a validation score -- the measure is a
count, so there is nothing to validate -- but whether the table separates
cohorts it has never seen is a checkable claim, and an uncheckable one would
not be worth making.

Output goes to `out/<dataset>/`:

```text
out/stability/
  README.md            dataset card
  manifest.json        parameters, source checksums, file checksums, split sizes
  data/train.parquet
  data/validation.parquet
  data/test.parquet
```

`cards/<dataset>/` holds a committed copy of the `README.md` and `manifest.json` from the
last local build, so the cards and the checksums of the Parquet files can be reviewed in
git. The Parquet files themselves are not committed.

`cache/`, `out/` and `target/` are git-ignored. Nothing in this pipeline uploads
anything; publishing to Hugging Face is a separate, manual step.

## Reproducibility

- **Pinned inputs.** `sources.lock` records the URL, size and SHA-256 of every
  upstream file. A cached file is reused only if its hash matches. A changed
  upstream file is an error until it is accepted explicitly with
  `--update-lock`. URLs point at versioned artefacts where the upstream offers
  them (Zenodo record, DisProt release, CATH release, git commit). UniProt REST
  and the wwPDB index are not versioned, so for those the lockfile is what
  detects drift.
- **Deterministic splits.** Splits use SplitMix64 with a fixed seed (default
  `20261114`, override with `--seed`) over groups sorted by id. The generator is
  implemented in `src/rng.rs` and tested against reference outputs, so a
  dependency upgrade cannot change the splits.
- **Deterministic output.** Rows are sorted before writing. Running the same build
  twice produces byte-identical Parquet files; `manifest.json` lists their hashes.
- **Offline rebuilds.** `--offline` refuses to touch the network.

## License review

Every source's license was checked on the maintainer's own site, not on mirrors,
on 2026-11-14.

| Source | License | Used |
|---|---|---|
| Tsuboyama et al. 2023 (Zenodo 7992926) | CC BY 4.0 | yes |
| DisProt | CC BY 4.0 | yes |
| CAID3 reference (BioComputingUP/CAID) | CC BY 3.0 | protein ids and sequences only, to choose the test split; nothing is redistributed |
| eSOL (NBDC LSDB Archive) | CC BY-SA 2.1 JP | yes; the share-alike term makes the solubility dataset CC BY-SA 2.1 JP |
| UniProtKB | CC BY 4.0 | yes (sequences for eSOL) |
| CATH | CC BY 4.0 | yes |
| wwPDB | CC0 1.0 | yes |
| FireProtDB 2.0 | custom terms: non-transferable license, no sublicensing of output, no derivative works | **excluded** |
| NESG solubility scores (Price et al. 2011) | no public data release | **excluded** |
| PSI:Biology / TargetTrack | service retired, no licensed download | **excluded** |

The full reasons are in `src/sources.rs` (`EXCLUDED_SOURCES`) and are repeated
in the affected dataset cards.

## Development

```bash
cargo fmt --all -- --check
cargo clippy --all-targets -- -D warnings
cargo test
```

Tests use small fixtures in `opendna-data/tests/fixtures/` and never touch the
network. They cover parsing of every upstream format, sign and unit
normalisation, duplicate merging, label construction, the lockfile rules, split
assignment, and Parquet output.

Code layout (`opendna-data/src/`):

| Module | Role |
|---|---|
| `sources.rs` | catalogue of upstream files, licenses and exclusions |
| `cache.rs` | cached downloads with retries and the SHA-256 lockfile |
| `cluster.rs` | MMseqs2 invocation and output parsing |
| `split.rs`, `rng.rs`, `unionfind.rs` | group-aware deterministic splitting |
| `parquet_out.rs` | table type and Parquet writer (column descriptions stored as field metadata) |
| `card.rs`, `context.rs` | dataset cards, manifests, output layout |
| `datasets/*.rs` | one module per dataset: parsing, cleaning, build |
