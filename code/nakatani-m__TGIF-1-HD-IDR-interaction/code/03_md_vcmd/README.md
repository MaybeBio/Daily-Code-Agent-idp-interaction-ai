# 03_md_vcmd

This directory contains representative input files and execution workflows used for mD-VcMD simulations of TGIF-1 HD-IDR interactions.

Two systems were analyzed:

- `pIDR-HD`: phosphorylated IDR model
- `non-pIDR-HD`: non-phosphorylated IDR model

The mD-VcMD workflow consisted of repeated cycles of simulation, canonical distribution reconstruction, and initial structure reassignment.

```text
Equilibrated structure
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
...
```

---

## Directory Structure

```text
03_md_vcmd/
├── common/
├── pIDR-HD/
└── non-pIDR-HD/
```

---

## common/

Files shared among all simulations.

### params.dat

Definition of the reaction-coordinate space used in mD-VcMD simulations.

---

## pIDR-HD/
## non-pIDR-HD/

The directory structure is identical for both systems.

```text
system/
├── parameters/
├── iteration1/
├── iteration2/
├── iteration3plus/
└── structure_selection/
```

---

## parameters/

Topology and restraint files required for mD-VcMD calculations.

Files include:

- `md.top`
- `ions.itp`
- `spce.itp`
- `water_and_ions_nonbonded.itp`
- `posre_HD.itp`
- `posre_IDP_P.itp`
- `posre_IDP_nonP.itp`

These files were used throughout all mD-VcMD iterations.

---

## iteration1/

Representative files for the first mD-VcMD iteration.

Files include:

- `run.bash`
- `run.mdp`
- `plumed.dat`
- `state.dat`
- `posre.gro`
- `run.gro`

### Description

This stage starts from the equilibrated structures generated in `02_equilibration`.

The simulation is initiated from:

```text
posre.gro
```

which corresponds to an equilibrated structure obtained after energy minimization and NPT equilibration.

---

## iteration2/

Representative files for the second mD-VcMD iteration.

Files include:

- `run.bash`
- `run.mdp`
- `plumed.dat`
- `state.dat`
- `next_conf/run.gro`

### Description

The second iteration starts from structures generated during the previous mD-VcMD iteration.

The simulation is initiated from:

```text
run.gro
```

generated in Iteration 1.

---

## iteration3plus/

Representative files for the third and subsequent mD-VcMD iterations.

Files include:

- `run.bash`
- `run.mdp`
- `plumed.dat`
- `state.dat`

### Description

For the third and subsequent iterations, initial structures are reassigned according to the canonical distribution reconstructed by Qcano.

The reassigned structures are generated during the structure selection procedure and stored in:

```text
iteration2/next_conf/
```

The simulations are then initiated from these reassigned structures.

This reassignment procedure is repeated after each Qcano reconstruction, providing updated initial configurations for the next mD-VcMD iteration.

---

## structure_selection/

Scripts used to assign initial structures for the next mD-VcMD iteration.

### link_conf.bash

Wrapper script used to generate assignments for the next iteration.

This script:

1. executes the structure selection procedure,
2. generates temporary assignment scripts,
3. creates the initial configurations for the next iteration.

### link_conf_dim3.py

Python script used for state-based structure reassignment in three-dimensional reaction-coordinate space.

This script assigns structures to each VcMD state according to the reconstructed canonical distribution.

---

## Relationship with Qcano Analysis

After each mD-VcMD iteration:

```text
mD-VcMD simulation
          ↓
Qcano reconstruction
          ↓
Parameter update
          ↓
Structure reassignment
          ↓
Next mD-VcMD iteration
```

The Qcano calculations used for canonical distribution reconstruction are provided in:

```text
04_analysis/Qcano
```

---

## Notes

The actual simulations were performed using 1000 independent replicas for each system.

Only representative files are included in this repository to illustrate the workflow and directory structure.

Representative examples of a single replica are provided for each iteration.

Qcano outputs and structure reassignment procedures were used to iteratively improve sampling and generate the initial structures for subsequent mD-VcMD iterations.
