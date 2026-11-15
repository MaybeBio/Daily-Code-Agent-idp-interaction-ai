# Reproduction

The versioned data and model release is [Zenodo 24190614](https://doi.org/10.5281/zenodo.24190614). Download only the resources needed for your analysis:

| File | Purpose |
|---|---|
| [benchmark.tar.gz](https://zenodo.org/records/24190614/files/benchmark.tar.gz?download=1) | Benchmark sequences, labels and fixed splits |
| [models.tar.gz](https://zenodo.org/records/24190614/files/models.tar.gz?download=1) | Experiment checkpoints, training assets and GussIDR |
| [feature_inputs.tar.gz](https://zenodo.org/records/24190614/files/feature_inputs.tar.gz?download=1) | PSSMs and sequence-level inputs for feature reconstruction |
| [external_datasets.tar.gz](https://zenodo.org/records/24190614/files/external_datasets.tar.gz?download=1) | External evaluation inputs and annotations |
| [predictions.tar.gz](https://zenodo.org/records/24190614/files/predictions.tar.gz?download=1) | Saved predictions and contribution-analysis inputs |
| [figure_source_data.tar.gz](https://zenodo.org/records/24190614/files/figure_source_data.tar.gz?download=1) | Numerical figure tables and statistical inputs |

Inspecting supplied results does not require feature reconstruction or training. ESM2 and MotifAE weights are obtained separately from their authors, as described below.

## Experiment models

The application checkpoint `data/models/model.pt` has identical contents to the fold-8, K619 checkpoint in `models.tar.gz`, including weights and fold/feature metadata. The included synthetic example uses this checkpoint.

Paper reproduction uses models from ten cross-validation folds. Each fold has two validation-selected operating points:

- **Performance:** the converged positive-penalty Stage C solution with the highest validation average precision (AP).
- **Efficiency:** the smallest feature set within 0.005 AP of that solution.

Both are subsequently refined in Stage D. `K` is the number of retained MotifAE features, not the window length. Folds 1, 3 and 8 select the same model for both operating points, giving 17 distinct experiment checkpoints. The application model is an alias of one of them, not an additional trained model or an ensemble.

| Fold | Performance K | Efficiency K |
|---|---|---|
| 1 | 929 | 929 |
| 2 | 695 | 501 |
| 3 | 490 | 490 |
| 4 | 902 | 610 |
| 5 | 825 | 573 |
| 6 | 1109 | 781 |
| 7 | 1039 | 722 |
| 8 | 619 | 619 |
| 9 | 500 | 265 |
| 10 | 1393 | 1033 |

Experiment checkpoints are distributed in `models.tar.gz`, not as separate downloads. After extraction into `DATA`, each checkpoint is at `DATA/train/outputs_jwarm_cv10/seed42/fold_XX/KNNN/jwarm_bundle.pt`, where `XX` is the two-digit fold and `NNN` is the selected K. The `experiment_models` entries in `src/configs/resources.json` record their archive paths and checksums.

For example, download the archive and use the fold-6 efficiency checkpoint:

```bash
python src/download.py --version 2026-11-13 --asset models \
  --output-dir data/archives
mkdir -p DATA
tar -xzf data/archives/models.tar.gz -C DATA

python src/predict.py \
  --checkpoint DATA/train/outputs_jwarm_cv10/seed42/fold_06/K781/jwarm_bundle.pt \
  --tensor INPUT.npz --metadata INPUT.tsv --feature-manifest MANIFEST.json \
  --store-dir SPARSE_FEATURE_DIR --output-prefix results/fold06_predictions
```

Skip the download and extraction if the archive is already available. Ordinary application users only need `model.pt` as their PriorSLiM checkpoint; upstream feature-generation resources are described in the main README.

## Rebuild experiment features

The release supplies `feature_inputs.tar.gz` instead of precomputed feature caches. Its PSSMs are accompanied by the exact scoring windows, actual protein/variant sequences, IDR and ordered-region annotations, frozen GussIDR predictions and contextual crop mappings. These inputs preserve the experiment cohorts and sequence context; do not rescan the proteins into a different set of windows.

1. Download the required archives from the links above to one directory. The resource downloader verifies SHA-256 automatically for its supported assets. Zenodo also displays an MD5 checksum for each published file. Extract the resources into a new `DATA` directory:

```bash
export RELEASE=/path/to/release_downloads
export PRIORSLIM_DATA_ROOT=/absolute/path/to/DATA
mkdir -p "$PRIORSLIM_DATA_ROOT"
tar -xzf "$RELEASE/benchmark.tar.gz" --strip-components=1 -C "$PRIORSLIM_DATA_ROOT"
for name in models feature_inputs external_datasets predictions figure_source_data; do
  tar -xzf "$RELEASE/$name.tar.gz" -C "$PRIORSLIM_DATA_ROOT"
done
python "$PRIORSLIM_DATA_ROOT/restore_splits.py"
python "$PRIORSLIM_DATA_ROOT/restore_analysis_paths.py"
```

2. From the GitHub repository, download the two frozen encoders and check the reconstruction inputs:

```bash
python src/download_esm2.py
python src/download_motifae.py --output-dir data/motifae
export MOTIFAE_CHECKPOINT="$(pwd)/data/motifae/MotifAE_k3_step_80000.pt"
python src/rebuild_features.py --data-root "$PRIORSLIM_DATA_ROOT" --check-inputs
```

3. Rebuild the required experiments. For example:

```bash
python src/rebuild_features.py --data-root "$PRIORSLIM_DATA_ROOT" \
  --experiment benchmark --device mps

# All experiment cohorts, including calibration and external evaluations:
python src/rebuild_features.py --data-root "$PRIORSLIM_DATA_ROOT" \
  --experiment all --device mps
```

| Experiment name | Reconstructed inputs |
|---|---|
| `benchmark` | Training, validation and fixed-test windows; Figures 2-3 |
| `calibration` | Full-IDR validation scans used for peak-threshold calibration |
| `localization_ASHI`, `localization_MOTIF`, `localization_phospho`, `localization_viral` | The four external localization cohorts; Figure 4 |
| `mutation` | Wild-type and variant entities for the reported mutation evaluation; Figure 5 |
| `functional_annotations`, `functional_background` | Annotated IDRs and matched-background scans; Figure 6 |
| `structural` | IDR interaction structural analysis; Figure 6 |

Use `cuda`, `mps` or `cpu` as supported by your machine. `--list` shows store IDs and window counts; `--job ID` rebuilds one store. `--explicit-only` needs no encoder weights. Completed outputs are verified and reused on rerun, and ESM2 crop files are removed after MotifAE encoding by default. A full run covers approximately 4.6 million windows and requires substantial compute and disk space; run only the cohorts needed for your analysis. ESM2/MotifAE inference can vary slightly across hardware and software versions, so newly computed activations are not promised to be byte-identical to the original caches. The published predictions and figure tables remain available for checking reported results without recomputation.

4. For analyses using the released fold-specific projections, build them after `benchmark` finishes:

```bash
python "$PRIORSLIM_DATA_ROOT/rebuild_compact.py" --data-root "$PRIORSLIM_DATA_ROOT"
```

This applies the released training-derived feature IDs and P95 scales; it does not refit them on external data. A fresh training run uses its own output directory and prepares its own fold-local assets.

## Retrain the paper models

Complete the `benchmark` reconstruction above. `PRIORSLIM_DATA_ROOT` must point to the resulting directory containing `new_data/`, `train/` and `feature_cache/`.

```bash
export PRIORSLIM_DATA_ROOT=/path/to/extracted_training_resources

# Check the supplied inputs and all ten fold splits; do not train.
python src/train.py check --device cpu

# Train stages A-C for fold 8 and select its validation operating points.
python src/train.py frontier --fold 8 --device cpu --output-root results/training

# Refine those selected models in stage D.
python src/train.py refine --fold 8 --device cpu --output-root results/training --resume
```

Replace `8` with the fold to reproduce; replace `cpu` with `cuda` or `mps` as appropriate. Results are written under `results/training/seed42/fold_08/`, including `automatic_k_selection.tsv` and the refined `K*/jwarm_bundle.pt` checkpoints. The training entry point selects models using validation data and does not evaluate the fixed test set.

## Regenerate figure panels

Run from the repository root in the project environment:

```bash
python reproduce/test_data.py
python reproduce/plot.py --figure all --output results/figures
```

`--figure` accepts `2`, `3`, `4`, `5`, `6`, or `all`. The output directory receives 27 SVG panels and matching 300-dpi PNGs. Manually drawn panels (Figure 1, 4a and 6e–f) are not generated by this command. The scripts use NumPy, pandas, Matplotlib and scikit-learn. Install Arial to match the plotting font; Matplotlib otherwise substitutes an installed font.

Numerical panels are regenerated from the included tables; final PowerPoint typography and spacing require separate layout postprocessing.

## Inputs and commands

Every row uses `python reproduce/plot.py --figure N --output results/figures`, with `N` equal to its figure number. Input paths below are relative to `data/fig0N/`.

| Panels | Calculation or display | Included input tables |
|---|---|---|
| 2a | Test discrimination curves | `aligned_test_predictions.tsv`, `a_metrics.tsv` |
| 2b | Negative-type false-positive rates and recall | `b_negative_types.tsv`, `aligned_test_predictions.tsv` |
| 2c | Corrected/still-wrong feature effects | `e_feature_effects.tsv` |
| 2d–e | Discovery feature selection and positional evidence | `discovery_f_candidates.tsv`, `locked_candidates.json`, `f_cases.tsv`, `f_position_profiles.tsv` |
| 2f–g | Dependence, activation coverage, retained support | `g_all_candidate_dependence.tsv`, `g_broad_candidate_evidence.tsv`, `all_candidate_fold_metrics.tsv` |
| 3a | Joint deletion and matched controls | `joint_deletion_endpoint.tsv`, `joint_deletion_curves.tsv` |
| 3b | Retained target differences at matched support | `retention_near600_display.tsv` |
| 3c | Contribution/composition partial correlations | `c_signed_contribution_composition.tsv` |
| 3d | Chemistry-control contributions and deletion effects | `d_original_absolute_values.tsv` |
| 3e | Training-strategy validation AP | `e_training_strategy_seed_values.tsv` |
| 3f | Five-seed support-budget curves | `f_budget_envelope.tsv` |
| 4b | Pooled Top-1 through Top-10 recovery | `pooled_summary.tsv`, `pooled_per_fold.tsv` |
| 4c | Thresholded localization and Top-K recovery | `peak_summary.tsv`, `topk_summary.tsv` |
| 4d | Same-pattern target/background ranking | `pattern_summary.tsv`, `pattern_per_fold.tsv` |
| 4e | Center-aligned profiles and locus availability | `combined_profiles.tsv` |
| 5a–c | Mutation direction and component changes | `../experiments/mutation.tsv` |
| 5d–f | NKX2-5 A61V, CAPRIN1 F372A, AMOTL2 R71I traces | `../experiments/mutation.tsv`, `curve_plot_data.tsv` |
| 6a–b | Functional-annotation enrichment/depletion | `main_results.tsv`, `plotted_values.tsv` |
| 6c–d | Domain-partner density ratios and adjusted odds ratios | `primary360_tests.tsv` |

## Recompute statistics

Run from the repository root; each output directory must be new.

```bash
python reproduce/experiments.py cv10 --output-dir results/cv10
python reproduce/experiments.py localization --output-dir results/localization
python reproduce/experiments.py mutation --output-dir results/mutation
python reproduce/experiments.py functional --output-dir results/functional
python reproduce/experiments.py structural --output-dir results/structural
python -m unittest discover -s reproduce -p 'test_experiments.py'
```

| Command | Included input under `reproduce/data/experiments/` | Output |
|---|---|---|
| `cv10` | `cv10.tsv.gz` | `cv10_test_per_fold.tsv` |
| `localization` | `localization.npz` | `pooled_summary.tsv`, `pooled_per_fold.tsv` |
| `mutation` | `mutation.tsv` | `mutation_summary.tsv` |
| `functional` | `functional/{pairs.tsv.gz,hits.npz,planned_tests.json}` | `main_results.tsv`, `per_fold_summary.tsv` |
| `structural` | `structural.tsv.gz` | `primary360_tests.tsv` |

Use `--input-dir` to select another input directory. These commands consume sample-level predictions, hits or peak counts. Model inference and training use `src/`.

## PowerPoint spacing

```bash
python reproduce/patch_spacing.py INPUT.pptx --check
python reproduce/patch_spacing.py INPUT.pptx OUTPUT.pptx --arial-font-dir ARIAL_DIR --inkscape INKSCAPE
```

Adjusts spacing in Figures 2 and 4. Requires `lxml`, `fonttools`, Inkscape and installed Arial fonts. The supported input is `Figures_1-6_Arial18-20_排版优化版.pptx`, SHA-256 `3cc4169cf0cec73dd34799d3c2dea03e0b7d17a69a1831d75ba699430fc31764`. Other input versions are rejected; the output must be a new file. Font rendering may differ across platforms.

## Contribution and ablation runs

`DATA_ROOT` is the extracted training-resource directory containing `new_data/`, `train/`, `feature_cache/` and `prior_compression_experiments/`. Outputs must be outside that directory.

Decoder-based controls also require the complete official MotifAE k3 checkpoint, downloaded separately from the [authors' versioned release](https://doi.org/10.5281/zenodo.20547228):

```bash
python src/download_motifae.py --output-dir data/motifae
export MOTIFAE_CHECKPOINT="$(pwd)/data/motifae/MotifAE_k3_step_80000.pt"
```

An existing copy can be used by setting `MOTIFAE_CHECKPOINT` to its absolute path. Alternatively, pass `--decoder-checkpoint /path/to/MotifAE_k3_step_80000.pt` to `run_ablation.py`. The full checkpoint contains both the encoder and decoder; an encoder-only export is not suitable for these controls. Contribution analysis using prepared activations does not require this download.

```bash
python reproduce/run_contributions.py check --data-root DATA_ROOT --output-dir results/contributions
python reproduce/run_contributions.py discovery --data-root DATA_ROOT --output-dir results/contributions
python reproduce/run_contributions.py followup --data-root DATA_ROOT --output-dir results/contributions
python reproduce/run_contributions.py candidates --data-root DATA_ROOT --output-dir results/contributions
python reproduce/run_ablation.py check --data-root DATA_ROOT --output-dir results/ablation
python reproduce/run_ablation.py smoke --data-root DATA_ROOT --output-dir results/ablation --arm C --seed 2026
```

| Entry | Phases | Inputs |
|---|---|---|
| `run_contributions.py` | `discovery`, `followup`, `candidates` | CV10 splits, fold assets/profiles, locked model checkpoints and validation predictions, MotifAE-only comparison outputs |
| `run_ablation.py` | `prepare`, `train`, `compress`, `refine` | Combined features, sparse MotifAE stores, fold-6 assets; `--arm`, `--seed`, `--device`, optional `--resume` |
| `run_ablation.py` | `decoder`, `profile`, `controls`, `refit` | Trained outputs in the same output directory and original feature resources |
| `run_ablation.py` | `analyze`, `analysis-check` | Completed ablation outputs; `analysis-check` runs numerical unit checks only |

Use `--help` for available arms and phases.

```bash
python reproduce/run_descriptors.py effects --data-root DATA_ROOT --output-dir results/descriptor_effects
python reproduce/run_descriptors.py composition --data-root DATA_ROOT --output-dir results/composition
```

`effects` reads CV10 test predictions and combined features. `composition` reads development-window contributions and descriptor resources. `composition-extract` first reconstructs contributions from the fold checkpoints and sparse stores. Each output directory must be new.
