# Beyond Stickers and Spacers

Code for "Beyond Stickers and Spacers: Residue-level Contact Networks Predict IDR Organization and Phosphorylation Rewiring"
(T. Rahman and V. N. Uversky, submitted to RECOMB 2027).

## What is here

| Folder | What it holds |
|---|---|
| `original/` | The scripts exactly as they were run, with `SHA256SUMS`. Check with `cd original && shasum -a 256 -c SHA256SUMS`. They contain the absolute paths of the machines used and are a record, not meant to be run as they stand. |
| `runnable/` | The same scripts with only the path prefixes changed, by `tools/build_runnable.py`. Every change is listed in `audit/runnable_changes.csv`. |
| `helpers/` | The few functions the phosphorylation analyses need from a separate unpublished project, copied verbatim by `tools/extract_helpers.py`. `tools/test_helpers.py` confirms that they give results identical to the source. The rest of that project is not released, and those five source files are the only scripts missing from `original/`. |
| `audit/` | How the contents were chosen. `manifest.csv` lists every script, why it is included, and its hash. `runnable_changes.csv` lists every path edit. |
| `tools/` | The scripts that built all of the above. |

Simulation inputs and run scripts are in [Sabian98/BENDER](https://github.com/Sabian98/BENDER). Trajectories are on Hugging Face:
[CALVADOS-2](https://huggingface.co/datasets/taseef/BENDER), [Mpipi-Recharged](https://huggingface.co/datasets/taseef/MPIPI-Recharged), and
[Mpipi-GG](https://huggingface.co/datasets/taseef/BENDER-MPIPI). Labels, trained models, and analysis results are in the [v1.0 release](https://github.com/Sabian98/beyond-stickers-spacers/releases/tag/v1.0).

## Running

```bash
cd runnable
for f in labels models_and_tables results_and_inputs; do
  curl -LO https://github.com/Sabian98/beyond-stickers-spacers/releases/download/v1.0/$f.tar.gz
  tar -xzf $f.tar.gz
done
export BSS_DATA=$PWD/data
source env.sh
cd redo_evidence/recharged_last16 && python compare.py
```

The archives unpack into `runnable/`. Some steps also need public files that are not in the deposit, placed under `$BSS_DATA` or `runnable/`:

| File | Source | Place at | Needed by |
|---|---|---|---|
| `BENDER_BIO.csv` | `data/bender_bio/` of [Sabian98/BENDER](https://github.com/Sabian98/BENDER) | `$BSS_DATA/BENDER-BIO/` | Many analyses and figures |
| CALVADOS-2 trajectories at 300 K | [taseef/BENDER](https://huggingface.co/datasets/taseef/BENDER) | `$BSS_DATA/BENDER-MT/BENDER_300K/` | Phosphorylation reweighting |
| Mpipi-GG trajectories at 300 K | [taseef/BENDER-MPIPI](https://huggingface.co/datasets/taseef/BENDER-MPIPI) | `$BSS_DATA/BENDER-BIO/MPIPI/mpipi_sim/300k_runs/` | Fig. 3 |
| `missense_variants.csv.gz` | The ClinVar table of the IDRome study | `runnable/external_data/idrome_clinvar/` | Fig. S15 and the ClinVar analyses |
| `variant_summary_2026-10-05.txt.gz` | The ClinVar release of 5 October 2026, from NCBI | `runnable/external_data/idrome_clinvar/` | Fig. S15 and the ClinVar analyses |

`runnable/env.sh` sets `BSS_REPO` and `PYTHONPATH`. The scripts need numpy, scipy, pandas, scikit-learn, networkx, matplotlib, and mdtraj, and the simulation scripts also need OpenMM and OpenMpipi.
The trained models load with Python 3.9.6, scikit-learn 1.6.1, and numpy 2.0.2.
We ran the figure scripts in two environments: Python 3.9 with scikit-learn for most, and Python 3.12 with pandas for Figs. 5, S22, S23, S24, and S25. From a clean copy with only the deposit and the files above, all 20 figure scripts run. That checks that they run, not that each output matches the paper pixel by pixel.
Commands printed for the computing cluster keep the original cluster paths, so change those to suit your own system.

## Figures

| Figure | Script |
|---|---|
| Fig. 1 | `methods_paper/figures/fig1_overview.py` |
| Fig. 2 | `methods_paper/figures/fig2_readouts.py` |
| Fig. 3 | `methods_paper/figures/fig3_arrangement.py` |
| Fig. 4 | `methods_paper/figures/fig4_validation.py` |
| Fig. 5 | `methods_paper/figures/fig6_phospho.py` |
| Fig. S1 | `methods_paper/figures/fig1_split.py` |
| Fig. S2 | `methods_paper/figures/figS_extra.py` |
| Fig. S3 | `methods_paper/figures/figS_null.py` |
| Fig. S4 | `methods_paper/figures/figS_context.py` |
| Fig. S5 | `methods_paper/figures/figS_predictor.py` |
| Fig. S6 | `methods_paper/figures/figS_extra.py` |
| Fig. S7 | `methods_paper/figures/figS_taxa.py` |
| Fig. S8 | `methods_paper/figures/figS_context.py` |
| Fig. S9 | `methods_paper/figures/figS_context.py` |
| Fig. S10 | No script survives (see Known limits) |
| Fig. S11 | No script survives (see Known limits) |
| Fig. S12 | No script survives (see Known limits) |
| Fig. S13 | `methods_paper/figures/figS_extra.py` |
| Fig. S14 | No script survives (see Known limits) |
| Fig. S15 | `methods_paper/figures/figS_clinvar.py` |
| Fig. S16 | `methods_paper/figures/figS_design.py` |
| Fig. S17 | `methods_paper/figures/fig5_case_fus.py` |
| Fig. S18 | `methods_paper/figures/figS_case_fus.py` |
| Fig. S19 | `methods_paper/figures/fig_atlas.py` |
| Fig. S20 | `methods_paper/figures/fig_atlas.py` |
| Fig. S21 | `methods_paper/figures/figS_extra.py` |
| Fig. S22 | `methods_paper/figures/figS_degree.py` |
| Fig. S23 | `methods_paper/figures/figS_phospho_recharged.py` |
| Fig. S24 | `methods_paper/figures/figS_reweighting.py` |
| Fig. S25 | `methods_paper/figures/figS_nmr_maps.py` |

## Analyses by section

| Part of the paper | Main scripts, in `runnable/` |
|---|---|
| Main text, Sections 3.1 and 3.2 | `redo_evidence/sep2_vs_sep3.py`, `redo_evidence/pitfall_theory.py`, `redo_evidence/frames_reliability.py` |
| S1 Simulation data and labels | `redo_evidence/hub_frames.py`, `redo_evidence/gaivi/hub_labels_gaivi.py`, `redo_evidence/gaivi/submit_hub_labels.sh`, `redo_evidence/recompute_v3.sh`, `redo_evidence/label_recharged_v3.sh`, `redo_evidence/update_recharged_v3_all.sh`, `redo_evidence/recharged_last16/compare.py` |
| S2 Analysis choices | `redo_evidence/sep2_vs_sep3.py`, `redo_evidence/pitfall_theory.py`, `redo_evidence/frames_reliability.py` |
| S3 Residue-neutral null | `redo_evidence/sim_null.py`, `redo_evidence/sim_null_v3.py`, `redo_evidence/sim_null_gg_pairs.py`, `redo_evidence/null_multi.py`, `redo_evidence/positional_null_check.py` |
| S4 Arrangement and interaction parameters | `redo_evidence/arrangement_test.py`, `redo_evidence/arrangement_test_amend_h4.py`, `redo_evidence/arrangement_test_degree.py`, `redo_evidence/arrangement_ff_table.py`, `redo_evidence/crossff_aa_degree.py`, `redo_evidence/label_shuffles_v3.sh`, `redo_evidence/label_yf_v3.sh`, `redo_evidence/agreement_ceiling.py` |
| S5, S6 Predictor and taxa | `redo_evidence/hub_predictor` |
| S7 Sticker and spacer context | `redo_evidence/hub_predictor`, `redo_evidence/context_stickers.py` |
| S8 Ensembles from experiment | `redo_evidence/fresh_ped_test.py`, `redo_evidence/ped_blind11_report.py`, `redo_evidence/ped_control.py`, `redo_evidence/ped_residue_classes.py`, `redo_evidence/ped_restraint_information.py`, `redo_evidence/fetch_ped.sh` |
| S9 Second test set, all-atom force fields | `redo_evidence/ped_v2_build.py`, `redo_evidence/ped_v2_labels.py`, `redo_evidence/ped_v2_score.py`, `redo_evidence/ped_extra_test.py`, `redo_evidence/helicity_allatom.py` |
| S10 PRE reweighting | `redo_evidence/pre_reweight.py`, `redo_evidence/pre_forward_check.py`, `redo_evidence/pre_direct_check.py` |
| S11 ClinVar variants | `redo_evidence/clinvar_hub.py`, `redo_evidence/clinvar_redo_predict.py`, `redo_evidence/clinvar_redo_test.py`, `redo_evidence/clinvar_hub_conservation.py` |
| S12 Design batch | `redo_evidence/gaivi/design_sequences.py`, `redo_evidence/gaivi/design_glycine.py`, `redo_evidence/gaivi/analyze_glycine.py`, `redo_evidence/analyze_design.py`, `redo_evidence/design_vs_experiment.py`, `redo_evidence/structure_vs_agreement.py`, `redo_evidence/dssp_numpy.py` |
| S13 FUS case study | `methods_paper/figures/case_fus_prep.py`, `methods_paper/figures/case_fus_raster.py` |
| S14 Predicted hubs across IDRs | `redo_evidence/hub_atlas.py`, `redo_evidence/consensus_hubs.py` |
| S15 Phase separation | `redo_evidence/llps_profiles.py`, `redo_evidence/llps_csat.py`, `redo_evidence/llps_unevenness.py`, `redo_evidence/llps_mutations.py`, `redo_evidence/llps_sticker_hub.py`, `redo_evidence/llps_csat_degree.py`, `redo_evidence/llps_unevenness_degree.py`, `redo_evidence/llps_mutations_degree.py`, `redo_evidence/llps_degree_family.py` |
| S16 Degree beside hub occupancy | `redo_evidence/hub_vs_degree_sim.py`, `redo_evidence/hub_vs_degree_models.py`, `redo_evidence/hub_vs_degree_ped.py`, `redo_evidence/hub_vs_degree_ext.py`, `redo_evidence/hub_vs_degree_e49.py`, `redo_evidence/degree_wetlab.py` |
| S17 Phosphorylation, Fig. 5 | `redo_evidence/phospho_network.py`, `redo_evidence/phospho_network_test.py`, `redo_evidence/phospho_direct_inputs.py`, `redo_evidence/phospho_direct_labels.py`, `redo_evidence/phospho_direct_test.py`, `redo_evidence/phospho_direct_recharged_test.py`, `redo_evidence/phospho_charge_test.py`, `redo_evidence/phospho_charge_posthoc_scd.py`, `redo_evidence/phospho_charge_posthoc_kappa.py`, `redo_evidence/phospho_function_test.py` |
| S17 NMR case series | `redo_evidence/phospho_nmr_inputs.py`, `redo_evidence/phospho_nmr_test.py`, `redo_evidence/phospho_nmr_recharged_test.py`, `redo_evidence/phospho_nmr_amend2/phospho_nmr_amend2_test.py`, `redo_evidence/phospho_nmr_window/phospho_nmr_window.py`, `redo_evidence/phospho_nmr_ceiling/phospho_nmr_ceiling.py`, `redo_evidence/phospho_rauh/phospho_rauh.py`, `redo_evidence/phospho_rauh/phospho_rauh_direct.py` |

This table lists the main entry points. `audit/manifest.csv` has the full list, including what each script imports.

## Known limits

- Nine panels of Figs. S10, S11, S12, and S14 come from images whose scripts no longer exist, so they cannot be regenerated. Their labels were later edited to match the paper's terms.
- The phosphorylation reweighting reads per-frame energies from a separate project. The fields it needs will be deposited with the labels.
- The Mpipi-Recharged labels of the paper cover 10,266 of the 10,282 runs, as Supplementary Section S1 explains. `redo_evidence/recharged_last16/` holds the check with all 10,282.

## License

MIT, see `LICENSE`.
