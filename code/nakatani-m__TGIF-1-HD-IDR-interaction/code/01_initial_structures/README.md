# 01_initial_structures

This directory contains files used to generate the initial simulation systems for the TGIF-1 HD-IDR interaction simulations.

Two systems were prepared:

- `pIDR-HD`: phosphorylated IDR model
- `non-pIDR-HD`: non-phosphorylated IDR model

The system preparation workflow consisted of:

```text
HD / IDR structure preparation
             ↓
Packmol structure generation
             ↓
Ion addition
             ↓
Initial structures for MD simulations
```

---

## Directory Structure

Each system directory contains:

```text
pIDR-HD/
├── input_structures/
├── parameter_files/
├── scripts/
├── posre.gro
└── state.dat

non-pIDR-HD/
├── input_structures/
├── parameter_files/
├── scripts/
├── posre.gro
└── state.dat
```

---

## input_structures/

Input structures used for system generation.

### 6fqp_after_tleap.pdb

TGIF-1 homeodomain structure after processing with Modeller and AmberTools.

### IDP_P_Min.pdb / IDP_nonP_Min.pdb

Energy-minimized IDR structures.

### spc.pdb

Water molecule model used for Packmol system construction.

---

## parameter_files/

Topology and parameter files required during system preparation.

Files include:

- `complex_P_edit_01.top`
- `complex_nonP_edit_01.top`
- `ions.itp`
- `spce.itp`
- `water_and_ions_nonbonded.itp`
- `ions.mdp`

These files were used during topology generation and ion addition.

---

## scripts/

### 01_HD_IDR_pdb.sh

Preparation of HD and IDR structures.

This script:

- processes the HD structure
- generates missing residues using Modeller
- prepares Amber-compatible structures
- constructs phosphorylated and non-phosphorylated IDR models
- performs IDR energy minimization

Outputs:

- `6fqp_after_tleap.pdb`
- `IDP_P_Min.pdb`
- `IDP_nonP_Min.pdb`

---

### 02_packmol_complex.sh

Generation of initial HD-IDR configurations using Packmol.

This script:

- places the HD and IDR structures in a simulation box
- generates multiple independent initial configurations

---

### 03_ion.bash

System neutralization and ion addition using GROMACS.

This script:

- prepares coordinate files
- generates topology files
- adds ions to achieve the target ionic concentration
- generates simulation-ready systems

---

## Notes

The actual simulations were performed using 1000 independently generated initial configurations for each system.

Only representative files are included in this repository.

The complete set of Packmol-generated initial structures is available in the Zenodo dataset described in the repository root README.
