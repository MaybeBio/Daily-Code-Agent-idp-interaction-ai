# PriorSLiM

SLiM prediction with biological priors and sparse pretrained features.

Data and model resources: [Zenodo, DOI: 10.5281/zenodo.24190614](https://doi.org/10.5281/zenodo.24190614). The repository supplies the code, default prediction model and small example; the Zenodo record supplies experiment inputs, model archives, predictions and figure source data.

## Install

```bash
conda env create -f environment.yml
conda activate priorslim
```

Run commands from this directory. `src/` contains feature generation, training and prediction; `reproduce/` contains experiments, plotting code and small input tables.

## Prediction model

Use the released `model.pt` for prediction. It contains the trained predictor, normalisation parameters, selected MotifAE feature IDs, positional weights and a validation-derived classification threshold. No training or model selection is required. ESM2, MotifAE and GussIDR are separate resources used to prepare inputs.

The repository includes `data/models/model.pt`, the prediction program's default checkpoint. No separate download is needed. To restore a missing copy from Zenodo, run:

```bash
python src/download.py
```

The downloader uses the published Zenodo URL and verifies SHA-256 against `src/configs/resources.json`; no URL argument is needed. If the model is already installed, continue below. Models used for individual paper experiments are bundled in [models.tar.gz](https://zenodo.org/records/24190614/files/models.tar.gz?download=1); see [the reproduction guide](https://github.com/Zamprini/PriorSLiM/blob/HEAD/reproduce/README.md#experiment-models).

## Quick start: run the included example

This checks that the included model reproduces the expected scores for 32 synthetic windows. It uses prepared features, so ESM2, MotifAE and GussIDR do not need to run:

```bash
python src/predict.py --example
```

The command exits successfully if all scores match the supplied reference within an absolute tolerance of `1e-5`. This example checks model loading and numerical prediction; it is not a biological benchmark.

## Rebuild the paper's experiment features

Download [feature_inputs.tar.gz](https://zenodo.org/records/24190614/files/feature_inputs.tar.gz?download=1) with the command below. It contains the experiment PSSMs, scoring windows, actual sequences (including variants), region annotations, frozen GussIDR predictions and ESM2 crop mappings. Large explicit tensors and MotifAE activation caches are **regenerated locally**, not distributed.

```bash
mkdir -p DATA
python src/download.py --asset feature_inputs --output-dir data/archives
tar -xzf data/archives/feature_inputs.tar.gz -C DATA
python src/download_esm2.py
python src/download_motifae.py --output-dir data/motifae
export MOTIFAE_CHECKPOINT="$(pwd)/data/motifae/MotifAE_k3_step_80000.pt"

# Start with the smallest external cohort. Use cuda, mps or cpu.
python src/rebuild_features.py --data-root DATA \
  --experiment localization_viral --device mps
```

Use `--experiment benchmark` for the training/test windows, or `--experiment all` for all experiments. `--list` lists available stores; `--check-inputs` verifies inputs without computation; `--explicit-only` builds the five explicit arrays without ESM2 or MotifAE. Rerunning the command reuses verified completed stores. Temporary ESM2 embeddings are removed after sparse encoding unless `--keep-esm` is set.

The supplied GussIDR tables and training-derived accessibility reference avoid rerunning GussIDR or fitting new references. ESM2 and MotifAE inference is still required, and full reconstruction is computationally substantial. For archive extraction, experiment names and training setup, follow the [reproduction guide](https://github.com/Zamprini/PriorSLiM/blob/HEAD/reproduce/README.md#rebuild-experiment-features).

## Predict sites in your own proteins

PriorSLiM scores a 29-aa window centred on each residue in an annotated IDR. Scanning an IDR produces a score profile that can be used to rank candidate sites. The workflow is:

```text
Protein sequences + IDR/ordered-region annotations + PSSMs
    -> windows.tsv
    -> explicit features and sparse MotifAE activations
    -> one score per window centre
```

The prediction command consumes features, not FASTA directly. Follow steps 1-5 to prepare them. If you already have the explicit feature files and sparse MotifAE store, go directly to step 5 and substitute their paths. Training data are not needed for this workflow.

### Obtain ESM2 weights

ESM2 weights are downloaded separately from [Meta's official distribution](https://github.com/facebookresearch/esm#available-models-and-datasets), not included in the Zenodo archive. Before generating features, run:

```bash
python src/download_esm2.py
```

This downloads the [2.60 GB checkpoint](https://dl.fbaipublicfiles.com/fair-esm/models/esm2_t33_650M_UR50D.pt) and its [small contact-regression companion](https://dl.fbaipublicfiles.com/fair-esm/regression/esm2_t33_650M_UR50D-contact-regression.pt), which the ESM loader also requires. It uses `curl`, resumes interrupted downloads on rerun, and verifies sizes and SHA-256 checksums before installing the files. Existing valid files are reused.

The default cache is `~/.cache/torch/hub/checkpoints/`; `TORCH_HOME` and `XDG_CACHE_HOME` are respected. To use another disk, set `TORCH_HOME` before both downloading and generating features. Direct browser downloads can also be placed in this cache; rerun the command to verify them. Official URLs and reference checksums are recorded in `src/configs/resources.json`. Once features have been generated, prediction and training do not reload ESM2.

### 1. Prepare sequences and regions

Place your inputs under `inputs/`:

| Input | Contents |
|---|---|
| `inputs/proteins.fasta` | Full protein sequences, with identifiers such as `>P12345`. |
| `inputs/regions.xlsx` | One row per protein, using the columns shown below. |
| `inputs/pssm/` | A PSI-BLAST ASCII PSSM for each protein, named by its identifier, e.g. `P12345.pssm`. PSSMs are not generated by this repository. |

Example region-table format (replace these illustrative coordinates with your annotations):

| UniProt_ID | Disorder_Regions | Ordered_Regions |
|---|---|---|
| P12345 | 31-120;181-240 | 1-30;121-180;241-300 |

Coordinates are **1-based and inclusive**. Use the same identifiers in FASTA, the region table and PSSM filenames. Separate intervals with semicolons; leave `Ordered_Regions` blank if no ordered regions are annotated.

Generate the scoring windows and extract the training-derived feature reference from the downloaded checkpoint:

```bash
export DEVICE=cpu  # alternatively: cuda or mps

python - <<'PY'
import json
from pathlib import Path
from src.features.chemical_window_idppi_fast import parse_fasta
from src.features.build_contextual_esm2 import load_idr_intervals, make_idr_scan_windows
import torch

Path("work").mkdir(exist_ok=True)
windows = make_idr_scan_windows(
    parse_fasta("inputs/proteins.fasta"),
    load_idr_intervals(Path("inputs/regions.xlsx")),
)
windows.to_csv("work/windows.tsv", sep="\t", index=False)
bundle = torch.load("data/models/model.pt", map_location="cpu", weights_only=False)
Path("work/training_feature_manifest.json").write_text(
    json.dumps(bundle["feature_manifest"], indent=2) + "\n"
)
print(f"Prepared {len(windows)} windows")
PY
```

`work/windows.tsv` contains `window_id`, `uniprot`, `window_start`, `window_end`, `window_center`, `sequence`, `idr_start` and `idr_end`. Each IDR residue is used as a centre, with a stride of one. Windows may extend beyond the IDR; only positions outside the protein are padded with `X`. No labels are required.

### 2. Prepare the IDR accessibility cache

The explicit-feature builder needs `work/cache/gussidr_ensemble_features.tsv`, containing GussIDR predictions for the input IDRs. If this cache has already been prepared for these proteins and regions, place it there and continue to step 3. Otherwise, use the cache-generation commands below.

<details>
<summary>Generate the cache from GussIDR resources</summary>

Set `GUSSIDR_DIR` to the directory containing the `As/`, `ree/`, `rg/`, `v/` and `spr/` model folders. Each folder must include its prediction module, hyperparameter JSON and weights. These MIT-licensed resources are included under `GussIDR_zsq/` in [models.tar.gz](https://zenodo.org/records/24190614/files/models.tar.gz?download=1).

The shared cache builder also prepares regional ESM2 and NARDINI caches. Its NARDINI step uses a separate environment; these additional caches are not additional inputs to the PriorSLiM predictor.

```bash
conda env create -f src/configs/environment-nardini.yml
export NARDINI_PYTHON="$(conda info --base)/envs/nardini/bin/python"
python src/download.py --asset models --output-dir data/archives
mkdir -p DATA
tar -xzf data/archives/models.tar.gz -C DATA
export GUSSIDR_DIR="$(pwd)/DATA/GussIDR_zsq"

python -m src.features.build_all_features \
  --windows work/windows.tsv \
  --pssm-dir inputs/pssm --fragment inputs/regions.xlsx \
  --cache-dir work/cache --gussidr-dir "$GUSSIDR_DIR" \
  --nardini-python "$NARDINI_PYTHON" \
  --esm-device "$DEVICE" --gussidr-device "$DEVICE" \
  --prepare-cache-only
```

This step prepares IDR-level resources only; it does not train PriorSLiM. Continue in the `priorslim` environment for the commands below.

</details>

### 3. Build explicit biological features

```bash
python -m src.features.build_priorslim_features \
  --windows work/windows.tsv \
  --pssm-dir inputs/pssm --regions inputs/regions.xlsx \
  --cache-dir work/cache \
  --training-manifest work/training_feature_manifest.json \
  --output-prefix work/features
```

This creates `work/features_tensor.npz`, `work/features_metadata.tsv` and `work/features_feature_manifest.json`. They contain the conservation, chemistry, motif and accessibility inputs in model order. The reference extracted in step 1 supplies the training-derived accessibility normalisation.

### 4. Build sparse pretrained features

Download `MotifAE_k3_step_80000.pt` (419.6 MB) directly from the [authors' versioned Zenodo record](https://doi.org/10.5281/zenodo.20547228). It is not included in our `models.tar.gz`. The command below resumes interrupted downloads, checks SHA-256 and reuses an existing verified file. Use the **k3** checkpoint, not k4, k5 or the PriorSLiM `model.pt`.

```bash
python src/download_motifae.py --output-dir data/motifae
export MOTIFAE_CHECKPOINT="$(pwd)/data/motifae/MotifAE_k3_step_80000.pt"
```

For a copy downloaded earlier, set `MOTIFAE_CHECKPOINT` to its absolute path; `python src/download_motifae.py --output-dir /path/to/checkpoint_directory --check-only` verifies it without downloading. ESM2 uses `esm2_t33_650M_UR50D`, layer 33, from the cache prepared above. Both encoders remain frozen. Neither upstream checkpoint is needed when using prepared features.

```bash

python -m src.features.build_contextual_esm2 \
  --windows work/windows.tsv \
  --pssm-dir inputs/pssm --regions inputs/regions.xlsx \
  --crop-mode idr_coverage --context-margin 128 \
  --output-dir work/contextual_esm2 \
  --metadata-output work/contextual_windows.tsv \
  --chunk-metadata-output work/contextual_chunks.tsv \
  --device "$DEVICE"

python -m src.features.build_motifae_sparse_windows \
  --context-metadata work/contextual_windows.tsv \
  --chunk-metadata work/contextual_chunks.tsv \
  --checkpoint "$MOTIFAE_CHECKPOINT" \
  --output-dir work/motifae_sparse --device "$DEVICE"
```

The second command writes the sparse activations and their window mapping to `work/motifae_sparse/`. All feature-generation commands above use the same `work/windows.tsv`.

### 5. Predict and read the scores

```bash
python src/predict.py \
  --tensor work/features_tensor.npz \
  --metadata work/features_metadata.tsv \
  --feature-manifest work/features_feature_manifest.json \
  --store-dir work/motifae_sparse \
  --output-prefix results/my_proteins --device cpu
```

The command loads `data/models/model.pt` automatically. Change `cpu` to `cuda` or `mps` to use a GPU.

Open `results/my_proteins.tsv` for the predictions:

| Column | Meaning |
|---|---|
| `window_id`, `uniprot`, `window_start`, `window_end` | The scored window and its protein coordinates. Its centre is `(window_start + window_end) / 2`. |
| `probability` | The sigmoid-transformed candidate score; higher values indicate stronger predicted SLiM evidence. |
| `logit` | The score before the sigmoid transformation. |
| `explicit_logit`, `motifae_logit` | The explicit and pretrained contributions; their sum is `logit`. |
| `prediction_at_validation_threshold` | Whether the score is at least the checkpoint's validation-derived window-classification threshold. |

Plot or rank `probability` against window centre within each IDR. These are candidate scores, not binding affinities or probabilities for a specified partner. The supplied threshold classifies individual windows; this command does not merge windows into peaks or apply an IDR-level peak threshold. The accompanying `.npz` stores numerical arrays and `.json` records the model and input provenance. Use a new output prefix for another run.

## Reproduce and test

```bash
python reproduce/plot.py --figure all --output results/figures
python -m unittest discover -s tests -v
python reproduce/test_data.py
python -m unittest discover -s reproduce -p 'test_experiments.py'
```

See [reproduce/README.md](https://github.com/Zamprini/PriorSLiM/blob/HEAD/reproduce/README.md) for experiment models, retraining commands and figure-input mappings. The prediction regression test runs automatically when `data/models/model.pt` is installed. Plotting uses Arial; install it to avoid font substitution.

## Citation and license

Citation metadata: [CITATION.cff](https://github.com/Zamprini/PriorSLiM/blob/HEAD/CITATION.cff).

Data and model record: [10.5281/zenodo.24190614](https://doi.org/10.5281/zenodo.24190614).

Code is released under the [MIT License](https://github.com/Zamprini/PriorSLiM/blob/HEAD/LICENSE). Third-party code and resources retain their own terms.
