# IDPInterp: Towards Interpretable Intrinsic Disorder in Protein Language Models

We train sparse autoencoders (SAEs) on PLMs with contrasting training objectives, evolutionary covariation (ESM-2) and biophysical dynamics (SeqDance), and interpret the resulting disorder features with region-aware evaluation and held-out-validated automated descriptions. Layer-wise linear probes, shuffled-region nulls, sequence-composition baseline, and paired region bootstraps are used to discover differences in model representations irrespective of unsupervised features.


## Repository layout

```
idpinterp/   SAE architectures, trainer, region-aware evaluator, embedders, data utils
scripts/     Numbered pipeline stages (01–19) plus make_figures.py, make_probe_figure.py
conf/        Hydra configs (model/, sae/, dataset/, paths/)
tests/       Unit tests
outputs/     Curated result snapshot: feature atlases, probe results, figures
```

## Installation

```bash
python -m venv .venv && source .venv/bin/activate
pip install -e .          # or: pip install -r requirements.txt
```

Requires Python ≥ 3.9. Training uses a CUDA GPU; analysis and figure scripts run on CPU.

## Feature atlas

The interpreted atlases live under `outputs/atlas/`, one directory per model
(`esm2_8m_topk_32x_k64_layer4/`, `seqdance_topk_32x_k64_layer10/`). Each
`atlas_table.csv` has one row per evaluated feature:

| column | meaning |
|---|---|
| `feature_idx` | SAE feature index |
| `triplet` | best-matching disorder axis (MobiDB/DisProt annotation) |
| `f1`, `precision`, `region_recall` | region-aware scores on the held-out test split |
| `best_threshold` | activation threshold selected on validation |
| `average_precision` | threshold-free precision |
| `description` | language-model description of the driving sequence property |
| `pearson_r` | held-out validation of the description |
| `family` | decoder-clustering family assignment |

`feature_descriptions.csv` is a human-readable description table and `atlas_summary.json`
holds per-model aggregate statistics.

## Pipeline

```bash
python scripts/01_fetch_mobidb.py                                # disorder corpus
python scripts/02_generate_embeddings.py model=esm2_8m           # per-residue embeddings
python scripts/03_train_sae.py      model=esm2_8m sae=topk_32x
python scripts/04_evaluate_atlas.py model=esm2_8m sae=topk_32x   # region-aware F1
python scripts/05_describe_features.py                           # descriptions + validation
python scripts/11_compare_models.py                              # activation correlation, CKA
python scripts/12b_bootstrap_cv.py                               # cross-validated paired gap
```

SAE hyperparameters, layers, and dataset are set through the Hydra configs in `conf/`. Swap to SeqDance with `model=seqdance`.

### Probes and controls

```bash
python scripts/16_layerwise_probe.py model=esm2_8m '+probe.layers=[1,2,3,4,5,6]' \
    '+probe.axes=[disprot-molten_globule]' +probe.dump_regions=true +probe.null_control=true
python scripts/18_rebuild_paired.py        # main paired tests
python scripts/19_paired_by_layer.py       # paired test at every matched layer
python scripts/make_probe_figure.py        # layer-wise probe figure
```

## Data sources

Residue- and region-level disorder annotations come from **DisProt** and **MobiDB**;
sequences from **UniProt/Swiss-Prot**. SAEs are trained without disorder labels.

SeqDance weights and code come from the public SeqDance repository (to local, `export SEQDANCE_REPO=/path/to/SeqDance`). Automated interpretation uses `ANTHROPIC_API_KEY` or `GEMINI_API_KEY`.
