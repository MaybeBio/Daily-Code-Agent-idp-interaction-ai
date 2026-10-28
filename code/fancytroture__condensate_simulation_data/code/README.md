# Temperature-dependent structure and dynamics of biomolecular condensates

This repository contains the simulation, preprocessing, VAMPNet/Markov-state-model,
analysis, and figure-generation code accompanying the manuscript by Xuezhen Cai and
Xiakun Chu.

The study compares FUS-LCD, A1-LCD, and TDP-43-LCD as free chains and in
200-chain slab condensates at paired low and high temperatures.

## October 2026 revision

The analysis and plotting code for the revised manuscript (2026-10-17 v15) is in
[`revision/`](https://github.com/fancytroture/condensate_simulation_data/blob/HEAD/revision/README.md), released on 2026-10-18. It includes paired
temperature-jump controls, three-temperature coexistence fits, common-state
comparisons, model-initialization checks, threshold/integral sensitivity, and
matched-residence and mobility analyses.

To reproduce the revised Figure 3 and the ten added or updated SI figures from
the included compact numerical inputs:

```bash
python revision/reproduce_figures.py
```

Outputs are written to `outputs/revision_20261018/`. Current figure numbers,
analysis entry points, dependencies, and the distinction between plotting and
reanalysis from trajectories are listed in the [revision guide](https://github.com/fancytroture/condensate_simulation_data/blob/HEAD/revision/README.md).
The original-release commands below retain their historical figure numbering;
use the revision entry point for the current Figure 3 and revised SI figures.

## Repository map

```text
simulation/       OpenABC/OpenMM condensate setup and production workflows
preprocessing/    dense-phase chain extraction and contact-feature generation
msm/              VAMPNet implementation and training code
analysis/         bootstrap and partner/contact-number analyses
figures/          final main, supporting, and TOC figure builders
source_data/      lightweight processed data and model summaries
toc_assets/       rendered condensate snapshots used by the TOC builder
revision/         October 2026 analyses, final figure builders, and compact inputs
outputs/          generated locally; ignored by Git
```

## Installation

Create the analysis environment with:

```bash
conda env create -f environment.yml
conda activate jpcl-condensate-temperature
```

OpenABC is required for running new simulations and is installed by the pip section
of `environment.yml`. VMD is needed only to regenerate the repository-provided
condensate snapshots used by the TOC graphic; it is not required to assemble the
graphic itself.

## Data location

Lightweight source tables and the rendered TOC inputs are versioned in this
repository. Large trajectories, feature arrays, representative structures, and
trained-model ensembles are not versioned on GitHub because of their size.

After obtaining the accompanying data package, point the code to its extracted root:

```bash
export JPCL_DATA_ROOT=/absolute/path/to/archival_data
```

The required directory layout and the inputs used by each figure are documented in
[`DATA_LAYOUT.md`](https://github.com/fancytroture/condensate_simulation_data/blob/HEAD/DATA_LAYOUT.md). Scripts never require the authors' absolute local
paths.

## Original-release figure builders

From the repository root:

```bash
python figures/make_figure1_thermodynamic_conformational_ensemble.py
python figures/make_figure3_msm_story.py
python figures/make_figure2_dynamics_story.py
python figures/make_s2_s3_structural_controls.py
python figures/make_supplementary_figures.py
python figures/make_toc_graphic.py
```

The manuscript ordering is:

| Manuscript item | Builder |
|---|---|
| Figure 1 | `make_figure1_thermodynamic_conformational_ensemble.py` |
| Figure 2 | `make_figure3_msm_story.py` |
| Figure 3 | `make_figure2_dynamics_story.py` |
| Figures S1, S3-S6, S8 | `make_supplementary_figures.py` |
| Figures S2 and S7 | `make_s2_s3_structural_controls.py` |
| TOC graphic | `make_toc_graphic.py` |

Generated files are written below `outputs/`. To select only one standard SI builder,
use, for example:

```bash
python figures/make_supplementary_figures.py --only 1
```

## Analysis workflow

```text
simulation trajectories
  -> reconstruction across periodic boundaries
  -> dense-slab chain selection
  -> intrachain/interchain contact features
  -> correlation and partner-survival analyses
  -> VAMPNet training and state alignment
  -> bootstrap/model-ensemble summaries
  -> manuscript figures
```

The principal analysis choices used in the manuscript are encoded in the scripts,
including the 0.66 nm residue-contact cutoff, exclusion of near-sequence neighbors,
the 2 ns Markov lag, and trajectory/model-level resampling.

## Systems

| System | Low temperature | High temperature | Condensate chains |
|---|---:|---:|---:|
| FUS-LCD | 280 K | 325 K | 200 |
| A1-LCD | 230 K | 260 K | 200 |
| TDP-43-LCD | 260 K | 300 K | 200 |

The A1-LCD simulation sequence contains the documented PY-to-GS substitution in the
M9 motif. Sequence definitions and model parameters are included directly in the
simulation inputs and described in the manuscript Supporting Information.

## Citation

Please cite the associated article and this repository. Machine-readable citation
metadata are provided in [`CITATION.cff`](https://github.com/fancytroture/condensate_simulation_data/blob/HEAD/CITATION.cff).

## License

Code is released under the MIT License; see [`LICENSE`](https://github.com/fancytroture/condensate_simulation_data/blob/HEAD/LICENSE). Repository-provided
data may be reused with attribution to the associated article.
