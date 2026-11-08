# 02_equilibration

This directory contains files used for energy minimization and equilibration before mD-VcMD simulations.

Two systems were prepared:

- `pIDR-HD`: phosphorylated IDR model
- `non-pIDR-HD`: non-phosphorylated IDR model

The equilibration workflow consisted of:

```text
Initial structure
        ↓
Energy minimization
        ↓
NPT equilibration (Berendsen)
        ↓
NPT equilibration (Parrinello-Rahman)
        ↓
Equilibrated structures for mD-VcMD
```

---

## Directory Structure

Each system directory contains:

```text
pIDR-HD/
├── mdp/
└── scripts/

non-pIDR-HD/
├── mdp/
└── scripts/
```

---

## mdp/

GROMACS parameter files used during equilibration.

### minim.mdp

Parameter file used for energy minimization.

Purpose:

- remove steric clashes
- relax unfavorable contacts
- prepare structures for equilibration

---

### posre_Berendsen.mdp

Parameter file used for position-restrained equilibration with the Berendsen barostat.

Purpose:

- stabilize temperature and pressure
- gradually relax the solvent environment
- maintain the overall protein structure during early equilibration

---

### posre_Parrinello-Rahman.mdp

Parameter file used for position-restrained equilibration with the Parrinello-Rahman barostat.

Purpose:

- produce a physically meaningful NPT ensemble
- equilibrate system density
- generate structures used as starting points for mD-VcMD simulations

---

## scripts/

Execution scripts used for the equilibration procedure.

### min.bash

Runs GROMACS energy minimization.

Input:

- initial structures generated in `01_initial_structures`

Output:

- minimized structure

---

### berendesen_posre.bash

Runs position-restrained equilibration using the Berendsen barostat.

Input:

- minimized structure

Output:

- Berendsen-equilibrated structure

---

### parrinello_rahman_posre.bash

Runs position-restrained equilibration using the Parrinello-Rahman barostat.

Input:

- Berendsen-equilibrated structure

Output:

- final equilibrated structure used for mD-VcMD simulations

---

## Notes

The actual simulations were performed using 1000 independently generated initial structures for each system.

The equilibration protocol described here was applied independently to every generated system.

Only representative parameter files and execution scripts are provided in this repository.

Representative equilibrated structures used as inputs for subsequent mD-VcMD simulations are included in:

- `03_md_vcmd/pIDR-HD/iteration1`
- `03_md_vcmd/non-pIDR-HD/iteration1`
