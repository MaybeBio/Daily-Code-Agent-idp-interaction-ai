# idr_embed_divergence

Protein language model analyses of functional divergence among intrinsically disordered regions (IDRs).

## Overview

This repository contains code to reproduce selected analyses from **Using Protein Language Model Embeddings to Identify Intrinsically Disordered Regions Driving Protein Family Diversification** ([manuscript in preparation]). The project uses protein language model embeddings to investigate functional divergence among IDRs within protein families, including cyclin-dependent kinases (CDKs).

The code supports generating sequence embeddings, calculating neighbour distances and Family Neighbour Divergence (FND), inferring functional associations through Gene Ontology (GO) enrichment, and training IDR-LM.


## Contents

- [Repository structure](#repository-structure)
- [Installation](#installation)
- [Input data](#input-data)
- [Analysis workflows](#analysis-workflows)
- [Outputs](#outputs)
- [Reproducibility and known limitations](#reproducibility-and-known-limitations)
- [Contributions and acknowledgments](#contributions-and-acknowledgments)
- [Citation](#citation)
- [License](#license)

## Repository structure

| Location | Contents |
| --- | --- |
| `src/idr_diverge/` | Core Python modules for embedding generation and downstream analyses. |
| `src/idr_diverge/IDR_LM_train/` | IDR-LM training code and model configuration files. |
| `scripts/` | Entry-point scripts for running workflows. |
| `configs/` | YAML configuration files specifying inputs, outputs, and analysis settings. |
| `envs/` | Conda environments and dependencies for different workflows |
| `data/` | Input datasets; see availability information below. |
| `res/` or `results/` | Output locations specified in configuration files. |
| `docs/` | Detailed workflow instructions. |

Detailed workflow instructions should be placed in `docs/`; proposed documentation filenames are listed below.

## Installation

The analyses were originally run on Linux. [Add the tested distribution/version and hardware requirements.]

### 1. Clone the repository

```bash
git clone <repo-link>
cd idr_embed_divergence
```

### 2. Create the analysis environment

```bash
conda env create -f envs/env.yml
conda activate diverge_env
pip install -e .
```

Run these commands from the repository root. The editable installation makes the `idr_diverge` package available to the workflow scripts.

### 3. Create workflow-specific environments as needed

For neighbour distance and FND analyses:

```bash
conda env create -f envs/fnd_env.yml
conda activate fnd_env
pip install -e .
```

For IDR-LM training:

```bash
conda env create -f envs/train_env.yml
conda activate lm_train
pip install -e .
```

## Input data

Workflows require some or all of the following:

- Protein or protein-segment sequences in FASTA format.
- Sequence embeddings in CSV or HDF5 (`.h5`) format, as supported by each script.
- Ortholog embeddings and, when filtering orthologs, an Ensembl homolog annotation file.
- A JSON mapping of protein-family names to UniProt IDs.
- A GO ontology file (`.obo`) and GO annotation table.
- A trained model checkpoint when generating embeddings with `IDR_LM`.

Most scripts expect sequence identifiers in the following format:

```text
{SPECIES_OR_ENSEMBL_ID}_{UNIPROT_ID}_{SEGMENT_TYPE}_{START}_{END}
```

Examples:

```text
HUMAN_Q14004_IDR_1_694
ENSNLEP00000006574_P12345_PFAM_45_120
```

Use consistent identifiers across related files. GO-enrichment query IDs must exactly match those in the embedding file.

Detailed specifications for identifiers, coordinate conventions, FASTA headers, and embedding formats belong in `docs/input-data.md` [to be added].


## Data availability

- **Toy datasets:** Small example inputs for testing the workflows are
  provided in `data/toy-data/`.
- **Processed examples:** Example processed datasets are provided in
  `data/examples/`.
- **Raw data:** Available from [Zenodo](https://github.com/denisele99/IDR-embed-diverge/blob/HEAD/YOUR_ZENODO_LINK). Includes annotation data, plm embeddings (human IDRs, human PFAM domains, CDK IDR and PFAM domain orthologs), IDR-LM training data, IDR-LM model weights and configs.


## Analysis workflows

Commands below assume that the entry-point scripts are located in `scripts/` and are run from the repository root. [Verify these paths against the repository.]

| Workflow | Purpose | Instructions |
|---|---|---|
| Data preparation | Prepare sequences, segment annotations, and protein families | [Data preparation](https://github.com/denisele99/IDR-embed-diverge/blob/HEAD/docs/data-preparation.md) |
| Embedding generation | Generate representations of IDRs and domains | [Embeddings](https://github.com/denisele99/IDR-embed-diverge/blob/HEAD/docs/embeddings.md) |
| Neighbour distance/divergence | Identify nearby segments in embedding space | [Nearest neighbours](https://github.com/denisele99/IDR-embed-diverge/blob/HEAD/docs/nearest-neighbours.md) |
| GO enrichment and similarity | Infer functional associations and compare enrichment profiles | [GO analysis](https://github.com/denisele99/IDR-embed-diverge/blob/HEAD/docs/go-analysis.md) |
| Family divergence | Calculate FND for protein families | [Family divergence](https://github.com/denisele99/IDR-embed-diverge/blob/HEAD/docs/family-divergence.md) |
| Visualization | Generate embedding visualizations and figures | [Visualization](https://github.com/denisele99/IDR-embed-diverge/blob/HEAD/docs/visualization.md) |

Each workflow document describes its required inputs, commands, parameters, and outputs.


### 1. Generate protein language model embeddings

Generate embeddings for sequences in a single FASTA file or a directory of FASTA files. Supported embedding types are:

- `esm`: ESM-1b.
- `IDR_LM`: trained IDR-LM; requires a model checkpoint.
- `IDR_LM_random`: randomly initialized IDR-LM.

```bash
python scripts/get_embeddings.py --config configs/embed_config.yaml
```

Set the embedding type, input path, output directory, and model/script paths in the configuration file. Paths in the embedding configuration are resolved relative to that configuration file.

**Output:** One HDF5 embedding file per input FASTA file.

**Detailed documentation:** `docs/embeddings.md`.

### 2. Calculate neighbour distances and FND

Compare neighbour distances between paralogous segments with an ortholog-based background. FND measures the difference between the average paralog neighbour distance and the expected ortholog/background neighbour distance.

**Required inputs:** Background embeddings, ortholog embeddings, and a family-to-UniProt mapping. An ortholog annotation file is also required when ortholog filtering is enabled. The family mapping is not required when running only the random-background calculation.

```bash
python scripts/run_ndivergence.py configs/ndist_config.yml
```

The configuration controls three analysis components:

1. `neighbour_divergence`: calculate neighbour distances between paralogous segments (`between`) or among orthologs of the same segment (`within`).
2. `random_neighbour_divergence`: estimate background divergence using orthologs of randomly sampled segments.
3. `calculate_FND`: combine family and background distance results to calculate FND.

For a staged run, generate the family and random-background distances first. Then set their output paths under `calculate_FND`, enable that component, and run the pipeline again. The configuration also supports recomputing background distances when no background distance file is supplied.

**Outputs:** Pickled neighbour-distance dictionaries and a CSV of FND results. Higher positive FND values indicate greater paralog divergence relative to the ortholog/background baseline.

**Detailed documentation:** `docs/divergence.md`.

### 3. Perform GO enrichment

Identify GO terms enriched among the neighbours of each query segment. Neighbours can be selected using:

- `knn`: nearest neighbours in embedding space (default).
- `BLAST`: top BLAST hits.
- `RANDOM`: randomly selected neighbours as a control.

The default neighbourhood size is 50.

```bash
python scripts/run_go_enrichment.py --config configs/go_config.yaml
```

**Required inputs:** An embedding file containing the query and background regions, GO ontology, and GO annotation table. The BLAST option additionally requires a BLAST output file.

**Output:** A CSV containing query identifiers, enriched GO terms with adjusted p-values, and term counts. Example columns are `query_id`, `term_adj_p`, and `term_counts`.

The related script `scripts/run_enrich_eval.py` evaluates or summarizes enrichment results. =

**Detailed documentation:** `docs/go-analysis.md`.

### 4. Train IDR-LM (optional)

Training code and configuration files are located in `src/idr_diverge/IDR_LM_train/`.

- `bert_config.json`: model architecture configuration.
- `idrlm_pretrain.yaml`: training arguments and file paths.
- `idrlm_train.py`: training script.
- `idrlm_embed.py`: embedding generation from trained or randomly initialized weights.

Activate `lm_train`, then configure the training FASTA, model configuration, optional resume checkpoint, and checkpoint output directories in `idrlm_pretrain.yaml`. Training configuration paths are resolved relative to the YAML file.

**Training data:** MobiDB-derived IDR sequences. [Add the exact dataset citation, release, and Zenodo download link.]

**Logging:** Training is configured to log to Weights & Biases. Configure the account, project, and run name, or disable that integration before running.

**Outputs:** Model checkpoints and training logs. If `load_best_model_at_end=True`, the final saved model uses the best checkpoint according to the selected evaluation metric.

**Detailed documentation:** `docs/training.md`.

## Outputs

| Analysis | Output |
| --- | --- |
| Embedding generation | HDF5 files containing sequence identifiers and embeddings. |
| Neighbour distance | `.pkl` dictionaries; structure depends on `dist_type` and `return_mode`. |
| FND | CSV containing group identifiers, FND values, and statistical results. |
| GO enrichment | CSV containing enriched GO terms, adjusted p-values, and term counts per query. |
| IDR-LM training | Model checkpoints and training logs. |

Output paths are specified in the workflow configurations. Full output schemas and examples belong in the corresponding workflow documentation.

