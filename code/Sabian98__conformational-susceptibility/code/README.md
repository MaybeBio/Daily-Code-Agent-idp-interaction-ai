# Conformational susceptibility of disordered proteins

Code for "The Rules of Conformational Susceptibility of Intrinsically Disordered Proteins across the Tree of Life: Predicting
the Unpredictable" (manuscript in preparation).

The paper computes how strongly the ensembles of 10,282 highly disordered proteins respond to temperature, salt, and
modification, from single simulated ensembles with fluctuation-response relations, in CALVADOS-2, Mpipi-GG, and
Mpipi-Recharged.

## What is here

| Folder | What it holds |
|---|---|
| `original/` | The scripts exactly as they were run, with `SHA256SUMS`. Check with `cd original && shasum -a 256 -c SHA256SUMS`. They contain the absolute paths of the machines used and are a record, not meant to be run as they stand. |
| `runnable/` | The same scripts with only the path prefixes changed, by `tools/build_runnable.py`, to `$CS_REPO` (this repository) and `$CS_DATA` (the data folder). Every change is listed in `audit/runnable_changes.csv`. |
| `audit/` | How the contents were chosen. `manifest.csv` lists every script of the project, whether it is included, why, and its hash. |
| `tools/` | The scripts that built all of the above. |

Within `original/` and `runnable/`, `code/` holds the analysis and figure scripts, `code/gaivi_templates/` the run scripts that
add the hydration forcings to CALVADOS-2 and Mpipi-Recharged simulations, and `bender_reanalysis/` two shared helper modules of
the BENDER analysis code (`bender_common.py`, `graph_features.py`). The published equivalents of these helpers are in
[Sabian98/BENDER](https://github.com/Sabian98/BENDER) under `analysis/`.

Left out are slide-deck and infographic scripts, exploratory analyses that the paper does not report, file housekeeping for the
simulation runs, and planning scripts for experiments. `audit/manifest.csv` gives the reason for each.

## Figures

| Figure | Script (`code/`) |
|---|---|
| Fig. 1 | `fig_overview.py` |
| Fig. 2 | `fig_manuscript2.py` |
| Fig. 3 | `fig_forcing.py` |
| Fig. 4 | `fig_extra.py` |
| Fig. 5 | `fig_heldout.py` |
| Fig. 6 | `fig_apps.py` |
| Fig. 7 | `fig_page4.py` |
| Supplementary figures | `fig_supplementary.py`, `fig_manuscript2.py`, `fig_extra.py`, `fig_forcing.py`, `fig_supp_redraw.py`, `fig_acetyl.py`, `fig_extrapolation.py`, `fig_lineage.py`, `fig_ogt_lineup.py`, `fig_phosphosite.py` |

The figure scripts read tables written by the analysis scripts. The header of each script states what it reads and writes.

## Data

Simulation inputs and run scripts of the proteome simulations are in [Sabian98/BENDER](https://github.com/Sabian98/BENDER).
Trajectories are on Hugging Face: [CALVADOS-2](https://huggingface.co/datasets/taseef/BENDER),
[Mpipi-Recharged](https://huggingface.co/datasets/taseef/MPIPI-Recharged), and [Mpipi-GG](https://huggingface.co/datasets/taseef/BENDER-MPIPI).
The tables that the figure scripts read are in `figure_data/`, in the project's folder layout, with a manifest (see `figure_data/README.md`). Per-frame caches are not included.

## Running

Python 3 with numpy, pandas, scipy, matplotlib, scikit-learn, statsmodels, networkx, and mdtraj. Scripts that recompute energies
or set up simulations also need OpenMM, OpenMpipi, and CALVADOS. The versions used were numpy 2.5, pandas 3.0, scipy 1.18,
matplotlib 3.11, scikit-learn 1.9, statsmodels 0.15, networkx 3.7, mdtraj 1.11, and OpenMM 8.6.

```bash
export CS_DATA=/path/to/data      # the folder that holds the data, in the project's folder layout
source runnable/env.sh
python runnable/code/fig_page4.py
```

Scripts that set up simulation jobs on a cluster keep the output folders of the original runs (under `/data` and `/general`)
as they were written, so adjust them before running.

## License

MIT, see `LICENSE`.
