# TGIF-1-HD-IDR-interaction

Scripts, simulation input files, and representative structures used to investigate interactions between the TGIF-1 homeodomain (HD) and intrinsically disordered region (IDR) using mD-VcMD simulations.

Two systems were analyzed throughout this study:

- **pIDR-HD**: phosphorylated IDR model
- **non-pIDR-HD**: non-phosphorylated IDR model

---

## Repository Structure

### [01_initial_structures](https://github.com/nakatani-m/TGIF-1-HD-IDR-interaction/blob/main/01_initial_structures/README.md)/

Files used to construct the initial simulation systems.

#### pIDR-HD/
#### non-pIDR-HD/

Each directory contains representative files used for preparing the simulation systems.

##### input_structures/

Input molecular structures used during system construction.

Examples:

- `6fqp_after_tleap.pdb`
- `IDP_P_Min.pdb` or `IDP_nonP_Min.pdb`
- `spc.pdb`

##### parameter_files/

Topology and parameter files required for system generation.

Files include:

- `complex_*_edit_01.top`
- `ions.itp`
- `spce.itp`
- `water_and_ions_nonbonded.itp`
- `ions.mdp`

##### scripts/

Scripts used for system preparation.

- `01_HD_IDR_pdb.sh`
  - Preparation of HD and IDR structures.
  - Generation of `6fqp_after_tleap.pdb` and minimized IDR structures.

- `02_packmol_complex.sh`
  - Generation of initial HD-IDR configurations using Packmol.

- `03_ion.bash`
  - Addition of ions and preparation of MD-ready systems using GROMACS.


---

### [02_equilibration](https://github.com/nakatani-m/TGIF-1-HD-IDR-interaction/blob/main/02_equilibration/README.md)/

Files used for energy minimization and equilibration before mD-VcMD simulations.

#### pIDR-HD/
#### non-pIDR-HD/

Each directory contains:

##### mdp/

GROMACS parameter files.

- `minim.mdp`
  - Energy minimization.

- `posre_Berendsen.mdp`
  - Position-restrained equilibration using the Berendsen barostat.

- `posre_Parrinello-Rahman.mdp`
  - Position-restrained equilibration using the Parrinello-Rahman barostat.

##### scripts/

Execution scripts for the equilibration workflow.

- `min.bash`
  - Energy minimization.

- `berendesen_posre.bash`
  - Berendsen equilibration.

- `parrinello_rahman_posre.bash`
  - Parrinello-Rahman equilibration.

---

### [03_md_vcmd](https://github.com/nakatani-m/TGIF-1-HD-IDR-interaction/blob/main/03_md_vcmd/README.md)/

Input files and representative examples used for mD-VcMD simulations.

#### common/

Files shared among all simulations.

- `params.dat`
  - Definition of reaction-coordinate space used in mD-VcMD simulations.

#### pIDR-HD/
#### non-pIDR-HD/

Each system contains the following components.

##### parameters/

Simulation topology and restraint files.

Files include:

- `md.top`
- `ions.itp`
- `spce.itp`
- `water_and_ions_nonbonded.itp`
- `posre_HD.itp`
- `posre_IDP_P.itp`
- `posre_IDP_nonP.itp`

##### iteration1/

Representative files for the first mD-VcMD iteration.

Files include:

- `run.bash`
- `run.mdp`
- `plumed.dat`
- `state.dat`
- `posre.gro`
- `run.gro`

In this stage, simulations start from equilibrated structures obtained after MD equilibration.

##### iteration2/

Representative files for the second mD-VcMD iteration.

Files include:

- `run.bash`
- `run.mdp`
- `plumed.dat`
- `state.dat`
- `next_conf/run.gro`

In this stage, simulations are continued using structures generated in the previous mD-VcMD iteration.

##### iteration3plus/

Representative files for the third and subsequent mD-VcMD iterations.

Files include:

- `run.bash`
- `run.mdp`
- `plumed.dat`
- `state.dat`

In this stage, initial structures are reassigned from the entire ensemble according to the reconstructed canonical distribution.

##### structure_selection/

Scripts used to generate initial structures for the next mD-VcMD iteration.

Files include:

- `link_conf.bash`
  - Execution script for structure reassignment.

- `link_conf_dim3.py`
  - Generates assignments of structures among reaction-coordinate states.

---

### [04_analysis](https://github.com/nakatani-m/TGIF-1-HD-IDR-interaction/blob/main/04_analysis/README.md)/

Scripts used for trajectory analysis.

#### PCA/

Principal component analysis (PCA) based on HD-IDR distance coordinates.

Files include:

- `mda_pca_traj_selectatoms.py`
- `pca.bash`

Part of the workflow was adapted from:

- [Omegagene Toolkit](https://github.com/kotakasahara/omegagene/tree/master/toolkit)

#### Qcano/

Canonical distribution reconstruction and parameter updates for subsequent mD-VcMD iterations.

Files include:

- `merge_vcmd_qraw.py`
  - Merges distributions from all replicas.

- `power_params.py`
  - Updates VcMD parameters for the next iteration.

- `prep_next_qraw.bash`
  - Prepares files required for Qcano calculations.

This workflow was performed using tools from:

- [Omegagene Toolkit](https://github.com/kotakasahara/omegagene/tree/master/toolkit)
- [Omegagene Qcano (qcano_szfit)](https://github.com/kotakasahara/omegagene/tree/master/vcmd/qcano_szfit)

#### DBSCAN/

Density-based clustering of PCA coordinates using the DBSCAN algorithm.

##### pIDR-HD/

DBSCAN analysis of PCA coordinates obtained from the phosphorylated model.

Files include:

- `dbscan.py`

##### non-pIDR-HD/

DBSCAN analysis of PCA coordinates obtained from the non-phosphorylated model.

Files include:

- `dbscan.py`

For both systems, PCA coordinates were standardized using z-score normalization prior to clustering.

Cluster labels were reassigned according to cluster population size:

- C1: largest cluster
- C2: second largest cluster
- ...

---

### [05_clustered_structure](https://github.com/nakatani-m/TGIF-1-HD-IDR-interaction/blob/main/05_clustered_structure/README.md)/

Representative structures used for structural analysis and discussion.

For each cluster, the 10 structures closest to the cluster centroid are provided.

#### pIDR-HD/

Included clusters:

- C1
- C2

#### non-pIDR-HD/

Included clusters:

- C1

Only clusters discussed in this study are included.

---

## Workflow

```text
Initial structure preparation
            ↓
Packmol system generation
            ↓
Ion placement
            ↓
Energy minimization
            ↓
Equilibration
            ↓
mD-VcMD Iteration 1
            ↓
Qcano reconstruction
            ↓
Structure reassignment
            ↓
mD-VcMD Iteration 2
            ↓
Qcano reconstruction
            ↓
Structure reassignment
            ↓
mD-VcMD Iteration 3+
            ↓
PCA
            ↓
DBSCAN clustering
            ↓
Representative structure extraction
```

---

## Notes

This repository contains representative input files, scripts, and structures required to reproduce the computational workflow.

The actual simulations were performed using **1000 independent replicas** for each system and iteration.

Large files such as trajectory files and complete ensembles of replicated structures are not included because of their file sizes.

Representative files are provided for each stage of the workflow to illustrate the directory structure and simulation protocol.


The complete dataset associated with this study, including:

- Packmol-generated initial structures used in the simulations,
- PCA coordinate files used as inputs for DBSCAN clustering analyses,

is available on Zenodo:

[Supporting Dataset for TGIF-1 HD-IDR Interaction Simulations](https://doi.org/10.5281/zenodo.23982610)
