# October 2026 revision inputs

The current revision uses the additional workspace documented in
[`revision/README.md`](https://github.com/fancytroture/condensate_simulation_data/blob/HEAD/revision/README.md). Its tracked CSV/NPZ inputs are sufficient
to regenerate revised Figure 3 and the ten added or updated SI figures with
`python revision/reproduce_figures.py`.

Raw trajectories, trained-model caches, and full bootstrap ensembles remain
external inputs for a complete reanalysis. The revision guide identifies these
dependencies; the compact plotting inputs do not substitute for them.

# Archival data layout

Set `JPCL_DATA_ROOT` to the directory containing the extracted archival data. The
figure scripts use paths relative to that root so the analysis is independent of the
authors' workstation layout.

The archive retains the following top-level organization:

```text
fus/fus_protein_single_simulation/
A1/A1_protein_single_simulation/
TDP43/tdp43_protein_single_simulation/
project_figure_set_Tl_Th/source_data/intrachain_relaxation/
project_figure/source_data/interchain_ratio_bootstrap/
ppt_5_11_figures/results/
ppt_5_11_figures/paper_story_figures/
storyline_figures_for_ppt_corrected_v3/
toc_assets/
```

## Figure dependencies

- **Figure 1:** per-condition `R_g`/`R_ee` caches, state-population tables, and
  state-resolved `avg_contacts.npy` arrays.
- **Figure 2:** aligned transition-network tables, state contact maps, and the
  20-model mapped-transition ensembles.
- **Figure 3:** `R_g`, intrachain-contact, interchain-contact, and partner-survival
  correlation summaries and bootstrap tables.
- **Figure S1:** slab density profiles.
- **Figure S2:** stratified intrachain-contact samples and trajectory-level means.
- **Figures S3-S6:** aligned state contact maps, TICA projections, implied-timescale
  ensembles, and Chapman-Kolmogorov validation arrays.
- **Figure S7:** `R_ee` correlation caches.
- **Figure S8:** partner-turnover, partner-number-change, and compensation summaries.
- **TOC:** FUS condensate renderings in `toc_assets/`, FUS structural/contact inputs,
  and the measured intrachain/interchain ACF arrays.

The small, publication-facing CSV/NPZ summaries that are practical to version are
also mirrored under `source_data/` in this repository. Raw trajectories and trained
network checkpoints are intentionally distributed through the archival dataset.

If an expected file is absent, the builders raise `FileNotFoundError`; verify that
`JPCL_DATA_ROOT` points to the extracted archive root rather than to an individual
system directory.
