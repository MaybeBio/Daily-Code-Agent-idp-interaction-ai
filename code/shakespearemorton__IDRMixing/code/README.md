# Quantifying (de)Mixing of Disordered Proteins in Molecular Dynamics Simulations

[![Paper](https://img.shields.io/badge/Paper-Link-blue)](https://doi.org/PLACEHOLDER)
[![Data](https://img.shields.io/badge/Data-Zenodo-blue)](https://doi.org/10.5281/zenodo.23689084)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This repository contains the code and data for reproducing the results in:

> **[Quantifying (de)Mixing of Disordered Proteins in Molecular Dynamics Simulations](https://doi.org/PLACEHOLDER)**  
> William Morton & Robert Vácha  
> *Under review*, 2026

## Overview

Biomolecular condensates organize cellular biochemistry, and their composition is often set by weak, multivalent interactions between intrinsically disordered regions (IDRs). This work introduces a **domain decomposition method (DDM)** that measures dense- and dilute-phase concentrations in direct-coexistence simulations without locating the condensate interface, and a continuous **demixing index** (*D*) that quantifies how two IDRs are organized inside a condensate.

**Key contributions:**
- A voxel-based method that recovers dense- and dilute-phase concentrations within 1% of the Gibbs dividing surface method across 399 single-IDR simulations
- A continuous demixing index (*D*) that distinguishes homogeneous mixing, surface enrichment, core–shell architectures and complete segregation
- Simulations of all 496 pairs of 32 IDRs, 2,068 binary mixtures in total, at several concentration ratios
- Comparisons with three experimental datasets: FUS–hnRNPA1 phase diagrams (Farag et al.), MED1/NPM1 partitioning (Lyons et al.) and an IDR colocalization network (Gilat et al.)

## The domain decomposition method

`DDM.py` analyses one slab (direct-coexistence) simulation of one or two IDR species:

1. Every frame is centred on the largest contiguous density along the long axis of the box.
2. The box is divided into voxels of side *S* (20 Å by default), and the residue count of each voxel is averaged over the trajectory.
3. The dense phase, *c*<sub>den</sub>, is the tallest peak (above the lowest one) of a kernel density estimate (KDE) of the voxel concentrations.
4. The dilute phase, *c*<sub>dil</sub>, is the lowest peak of a KDE of ln *c* over the occupied voxels. Peaks produced by voxels visited only a few times are skipped: fewer than 10 residue-visits over the run *and* under 1% of the occupied voxels.
5. Uncertainties come from block averaging the per-frame concentration over the voxels in each peak's half-maximum window.
6. ΔG = *k*<sub>B</sub>*T* ln(*c*<sub>dil</sub>/*c*<sub>den</sub>), clamped at −10 *k*<sub>B</sub>*T*. For two species, the demixing index *D* is also reported.

```
python DDM.py top.pdb traj.dcd [traj_2.dcd ...] --name NAME --start 1000 --voxel-size 20 --out NAME.npz
```

Two species are detected automatically from the topology. The output is one `.npz` per simulation holding the time-averaged voxel maps, the phase concentrations with their errors, ΔG and, for mixtures, the per-species concentrations and *D*.

### Checking the result

`DDM.py` always returns a *c*<sub>dil</sub>, and returns a *c*<sub>den</sub> and ΔG whenever the KDE of *c* finds two or more peaks, which can also happen in a single-phase system (for example from transient clusters). When only one peak is found, *c*<sub>den</sub> and ΔG are NaN and the run is flagged `no_dense_peak`. Before treating a system as two-phase, check two quantities saved in the `.npz` and printed in the report:

| Field | Check | What a failure suggests |
|-------|-------|-------------------------|
| `frac_vox_ge_cden` | at least 5% of voxels at or above *c*<sub>den</sub> | no real dense phase, only transient clusters |
| `dil_visits` | at least 10 residue-visits at the dilute peak (*c*<sub>dil</sub> × frames analysed) | the dilute peak is a counting artifact of voxels entered only a few times, not a resolved dilute phase |

If either check fails, the system may be a single phase. Report `cdil_onephase`, the mean concentration over all voxels, as its concentration, with no dense phase and no ΔG. `DDM.py` flags both cases in its `warning` and `warning_reasons` fields (`one_phase`, `dil_few_visits`), along with other conditions worth a look, and prints `cdil_onephase` when one of them trips. A flag is a prompt to inspect the system, for example its concentration profile, not an automatic verdict.

In the heterotypic analysis, a map is treated as two-phase only when at least 5% of its voxels sit at or above *c*<sub>den</sub>; maps that fail this gate are excluded. A pair is considered well mixed when *D* < 0.13 and phase separating when ΔG<sub>*ij*</sub> < −3 *k*<sub>B</sub>*T*. These criteria are applied in the figure scripts, not in `DDM.py`.

## Repository structure

```
├── DDM.py               # Domain decomposition method (single simulation)
├── Blender_Renders/     # Blender scenes and renders used for the snapshots
├── Create_Subset/       # Selection of the 32 IDRs and the mixture compositions
├── DATASETS/            # Sequences, simulation set-ups and experimental data
├── DDM_Homotypic/       # DDM over the single-IDR simulations, all voxel sizes
├── DDM_Heterotypic/     # DDM over the binary-mixture simulations
├── Figure1F/            # Finite-size check (Figure 1F)
├── Figures/             # Figure scripts and the figures they produce
└── Run_MD/              # CALVADOS2 simulation scripts
```

### `Create_Subset/`
- **`create.py`**: selects IDRs from the von Bülow et al. training set (`df_training.csv`) across GIN sequence-grammar clusters and ΔG ranges, giving `gin_samples.csv`
- **`compositions.py`**: builds the mixture set-ups for every pair at 50:50 and further ratios (about 90:10, 75:25, 25:75 and 10:90), keeping about 36,000 residues per box, giving `gin_prep.csv`

### `DDM_Homotypic/`
- **`run_ddm.py`**: runs `DDM.py` on each single-IDR simulation of von Bülow et al. at every voxel size *S* = 10–60 Å, writing `results_<S>/` and collecting `slab_results_<S>.csv`
- **`frac_two_peak.py`**: fraction of frames whose instantaneous voxel concentrations resolve two peaks (Figure 1E), written to `peaks_20/` and merged into `slab_results_20.csv`
- **`sub.sh`**: PBS array submission

### `DDM_Heterotypic/`
- **`run_ddm_het.py`**: runs `DDM.py` on every binary mixture in four sets and collects one CSV per set:

  | Set | Prep CSV | Output | Used in |
  |-----|----------|--------|---------|
  | main | `gin_prep.csv` | `heterotypic_main.csv` | Figures 2–3, Supplementary Figures 2–4 |
  | known | `known_gin_prep.csv` | `heterotypic_known.csv` | Figure 4A (FUS–hnRNPA1) |
  | sabari | `sabari_gin_prep.csv` | `heterotypic_sabari.csv` | Figure 4C (MED1/NPM1) |
  | exp | `exp_gin_prep.csv` | `heterotypic_exp.csv` | Figure 4D (Gilat et al. network) |

  The exp trajectories are already centred and trimmed (frames 1000 onward), so that set is analysed from frame 0.
- **`sub.sh`**: PBS array submission

### `Figure1F/`
- **`finite.py`**: runs `DDM.py` on three IDRs simulated in 17×17×170 nm and 25×25×250 nm boxes at the same residue density (Figure 1F)

### `DATASETS/`

| File | Description |
|------|-------------|
| `df_training.csv` | Reference IDRs and their simulated ΔG from von Bülow et al. |
| `gin_samples.csv` | The 32 IDRs of this study with their sequence features |
| `gin_prep.csv` | Set-ups for the 2,068 binary-mixture simulations |
| `residues.csv` | CALVADOS2 amino-acid parameters |
| `IDRome_DB_full.csv` | IDRome database, for sequence look-ups |
| `known_*.csv` | FUS–hnRNPA1 (WT and p12D) comparison with Farag et al. |
| `sabari_*.csv` | MED1/NPM1 comparison with Lyons et al., including the pellet/supernatant data (`sabari_all_pellet.csv`) |
| `exp_*.csv` | The 18 IDRs of the Gilat et al. colocalization network |
| `*_demixing.csv`, `demixing.csv`, `homotypic_conc.csv` | Earlier analysis results, kept for reference |

### `Figures/`
Each script writes its figure to `Figures/` and the plotted data to a `figure*_data/` or `supplementary*_data/` directory.

| Script | Content |
|--------|---------|
| `Figure1.py` | DDM validation against the Gibbs dividing surface: voxel-size dependence, dense and dilute concentrations, two-peak frequency, finite-size check and the hnRNPA1 phase diagram |
| `Figure2AC.py` | Species-resolved KDEs and voxel concentration space for mixed and demixed pairs |
| `Figure2GH.py` | Distribution of the demixing index and its dependence on composition |
| `Figure3.py` | Demixing by chemical interaction category, composition dependence of ΔG<sub>*ij*</sub>, feature differences between mixed and demixed pairs, and stability of mixed condensates |
| `Figure4.py` | Comparison with experiments: FUS–hnRNPA1 tie lines, MED1/NPM1 partitioning and the Gilat et al. interaction network |
| `Supplementary1.py` | Sequence features of the 32 selected IDRs against the source dataset |
| `Supplementary2.py` | Demixing index against the local mixing order parameter of Chen & Jacobs |
| `Supplementary3.py` | Homotypic ΔG difference for mixed and demixed pairs, by interaction category |
| `Supplementary4.py` | Heterotypic ΔG<sub>*ij*</sub> against the mean homotypic ΔG, by interaction category |

### `Run_MD/`
Scripts for running the CALVADOS2 direct-coexistence simulations in OpenMM (`sim_calvados.py`, `calvados_param.py`, `run_list.py`).

## Data availability

The `.npz` files written by `DDM.py` are too large for this repository and are available on Zenodo: [10.5281/zenodo.23689084](https://doi.org/10.5281/zenodo.23689084). Each holds one simulation's time-averaged voxel maps, phase concentrations, errors, ΔG and, for mixtures, the per-species concentrations and *D*. They cover:

- `DDM_Homotypic/results_<S>/`: single-IDR simulations at every voxel size
- `DDM_Heterotypic/results/<set>/`: binary mixtures (main, known, sabari, exp)
- `Figure1F/results_20/`: finite-size check
- `A1_Phase/results/`: hnRNPA1 phase diagram

Place them at these paths to rebuild the CSVs (`--collect`) or the figures without rerunning the trajectory analysis.

## Requirements

```
numpy>=1.26
pandas>=2.1
scipy>=1.13
scikit-learn>=1.5
matplotlib>=3.9
seaborn
colormaps
pycirclize          # Figure 4D chord diagram
MDAnalysis>=2.9
OpenMM>=8.2         # for running simulations
```

Block-averaged errors use [BLOCKING](https://github.com/fpesceKU/BLOCKING). Set its path at the top of `DDM.py`.

## Citation

If you use this code or data, please cite:

```bibtex
@article{morton2026mixing,
  title={Quantifying (de)Mixing of Disordered Proteins in Molecular Dynamics Simulations},
  author={Morton, William and V{\'a}cha, Robert},
  journal={},
  year={2026},
  doi={PLACEHOLDER}
}
```

## Related Resources

- [CALVADOS2 Force Field](https://github.com/KULL-Centre/CALVADOS)
- [von Bülow et al. Dataset](https://github.com/KULL-Centre/_2024_buelow_PSpred)
- [IDRome Database](https://github.com/KULL-Centre/_2023_Tesei_IDRome)

## License

This project is licensed under the MIT License - see the [LICENSE](https://github.com/shakespearemorton/IDRMixing/blob/HEAD/LICENSE) file for details.

## Contact

William Morton — william.morton@ceitec.muni.cz  
CEITEC, Masaryk University, Brno, Czech Republic
