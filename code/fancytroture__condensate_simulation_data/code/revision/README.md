# Revised manuscript analysis and figures — 2026-10-18

This release accompanies *Biomolecular Condensation Buffers Global Chain Structure
while Enhancing Conformational State Remodeling*, revised manuscript version
2026-10-17 v15. The original analysis directory layout is retained so the archived
scripts and their dependencies remain identifiable. Dates in directory names
identify the calculation, not the current manuscript figure number.

## Reproduce the final figures

From the repository root, using the environment in `environment.yml`:

```bash
python revision/reproduce_figures.py
```

This command uses the included compact numerical inputs and writes PDF/PNG files
to `outputs/revision_20261018/`. It also writes the response-figure aliases R1–R7.
Individual groups can be selected with `--only main`, `phase`, `states`, `models`,
`sensitivity`, `mobility`, `tj`, or `events`.

| Current manuscript figure | Response alias | Source within this workspace |
|---|---|---|
| Figure 3 | — | `make_figure3.py`; full-series builder: `revision_2026_10_08/JCP_revision/update_dynamics_figures.py` |
| S2: coexistence phase diagrams | — | `revision_2026_10_08/JCP_revision/phase_diagram_20261013/plot_phase_diagram.py` |
| S8: common-state transitions and populations | R2 | `revision_2026_10_08/R2C2_state_remodeling/response_package/build_package.py` |
| S9: model populations and assignment overlap | R3 | `revision_2026_10_08/R2C3_model_consistency/analyze.py` |
| S11: crossing-threshold ratios | R4 | `revision_2026_10_08/JCP_revision/response_revision_20261016/plot_figures.py` |
| S12: finite-window integral ratios | R5 | Same as S11 |
| S13: selected/all-chain lateral mobility | R6 | `revision_2026_10_08/JCP_revision/remove_R6B_20261016/plot_mobility.py` |
| S14: matched 100 ns residence comparison | R7 | Same as S11 |
| S15: partner loss/gain at a 2 ns lag | — | `reproduce_figures.py --only events`; full-series calculation in `update_dynamics_figures.py` |
| S16: TJ dense/dilute regional densities | — | `revision_2026_10_08/JCP_revision/tj_notation_20261017/make_figures.py` |
| S17: four TJ crossing-time ratios | R1 | Same as S16 |

Figures 1–2 and inherited SI illustrations retain their original builders. The
table above gives the current numbering for the revised analyses. The plotting
entry point preserves the final TJ superscripts, four-observable ratio panels,
inside-panel legends/titles, and single-panel selected/all-chain mobility figure.
Older helper modules retain their original plotting functions for provenance;
the entry point calls only the final selected panels.

## Analysis entry points

All paths below are relative to this `revision/` workspace.

| Analysis | Entry points |
|---|---|
| Original selected-window intermolecular series, correlations, synchronized block bootstrap | `revision_2026_09_30/interchain_recompute/compute.py`, `analyze.py`, `sensitivity.py` |
| Paired TJ observables and paired-snapshot bootstrap | `revision_2026_10_08/paired_slab_temperature_jump/extract_full20ns.py`, `analyze_full20ns.py` |
| Extended TJ partner-retention lag | `revision_2026_10_08/JCP_revision/response_revision_20261015/extend_partner_lags.py` |
| Coexistence densities at added temperatures | `revision_2026_10_08/critical_temperature_slab_simulations/check_density_3to7us.py` |
| Three-temperature critical fit from included densities | `revision_2026_10_08/JCP_revision/refit_critical_temperature.py` |
| Common-state matching and trajectory-sampling uncertainty | `revision_2026_10_08/R2C2_state_remodeling/no_population_filter/analyze.py`, `followup/analyze_followup.py`, `response_package/build_package.py` |
| Model populations, assignment agreement, and all-entry comparisons | `revision_2026_10_08/R2C3_model_consistency/analyze.py`, `all_edges_recalculation_2026_10_09/recalculate.py`, `add_sampling.py` |
| Four-observable threshold/integral uncertainty | `revision_2026_10_08/JCP_revision/response_revision_20261015_v2/rebuild_intramolecular.py`, `sensitivity_statistics.py` |
| Selected/all-chain mobility and strict residence | `revision_2026_10_08/R2C5_chain_selection/extract.py`, `central60.py`, `analyze.py` |
| Add Rg/Nintra to the matched residence analysis | `revision_2026_10_08/JCP_revision/response_revision_20261015_v2/complete_residence.py`, `refine_residence_intra.py` |

The small-data critical-temperature fit can be run directly:

```bash
python revision/revision_2026_10_08/JCP_revision/refit_critical_temperature.py
```

Full reanalysis requires the archived trajectories, model outputs, selection
metadata, and full bootstrap arrays. Those large inputs are not included here.
The scripts retain workspace-relative paths under `fus/`, `A1/`, `TDP43/`,
`revision_2026_09_30/`, and `revision_2026_10_08/`; restore the corresponding archive
with that directory layout before invoking raw-data stages. Simulation stages
also require their configuration, serialized systems, starting snapshots, and
completion metadata. The plotting-only command does not require those raw inputs.
Use a working copy when rerunning analysis stages, because they write calculated
tables beside their scripts.

## Numerical inputs and provenance

`MANIFEST.json` records SHA-256 hashes of the selected original scripts and
numerical inputs, along with any path-only release adaptations. The release adds
small entry-point wrappers but does not change the scientific estimators.

`figure3_inputs/` stores the displayed curves, conditional interval endpoints,
and Rg temperature-response intervals. It is derived from the existing original
intramolecular caches and the physical-time intermolecular results. Its export
script keeps the original 5,000-resample Rg bootstrap and the 5,000 synchronized
250 ns origin-block intermolecular bootstrap:

```bash
python revision/export_figure3_inputs.py --workspace /path/to/archived_workspace
```

The exporter is for use with the full archive. It computes the compact plotting
inputs; it does not launch simulations. The full-series analysis and figure
scripts are also retained so the interval calculation is not represented only by
its plotted endpoints.

Additional numerical tables include all five populations for each of the twenty
models per condition, complete transition matrices, common-state assignments,
sampling/model-initialization comparisons, and the nine sampled coexistence
density pairs with fit parameters. These tables are computation inputs or outputs;
their archived names do not change the three-table organization of the current SI.

## Interpretation and definitions

- TJ simulations branch from paired starting configurations; their regional
  densities evolve. They are transient controls, not equilibrium fixed-density
  simulations.
- Threshold times are first crossings of the stated levels. Noncrossing samples
  remain censored, and partner retention permits departure and return.
- The primary edge-selected distance, all-entry distance, and common-state
  distance have distinct definitions. Trajectory-sampling and model-initialization
  intervals remain separate.
- The matched-residence comparison includes Rg, Nintra, Ninter, and partner
  retention, with the documented cadence for each observable.

## Verification

The release validation records figure comparisons and numerical checks in
`validation.json`. Existing tests can be run without launching MD:

```bash
python -m unittest discover -s revision/revision_2026_09_30/interchain_recompute -p 'test_*.py' -v
python -m unittest discover -s revision/revision_2026_10_08/R2C5_chain_selection -p 'test_analysis.py' -v
```

These tests require the scientific dependencies listed in the root environment,
including MDAnalysis. Rendering the figures checks the supplied plotting inputs;
it does not repeat full trajectory extraction, VAMPNet training, or MD production.
