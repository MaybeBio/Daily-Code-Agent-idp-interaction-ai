## Preparing and Running Simulations

### 1. System Preparation

**Behavior:**

* **With ligand:**

  * Extracts ligand from the PDB file.
  * Runs `antechamber` to generate ligand `.mol2`, `.prepin`, and `.frcmod` files.
  * Uses `build_complex.leap` to combine protein and ligand, solvate, add ions, and generate `step3_input.{parm7,rst7,pdb}`.

* **Protein-only:**

  * Extracts protein without hydrogens.
  * Uses `build_protein.leap` to solvate and add ions.
  * Generates `step3_input.{parm7,rst7,pdb}`.

---

### 2. Equilibration

Two-step equilibration is implemented:

1. **step4.1_equilibration.inp** – Minimization + restrained MD

   * Steps: 5000 minimization + 125000 restrained MD
   * Time-step: 1 fs
   * Positional restraints: backbone (400 kJ/mol/nm²), side-chain (40 kJ/mol/nm²)
   * Temperature: 303.15 K, no pressure coupling

2. **step4.2_equilibration.inp** – Unrestrained MD

   * Steps: 1,000,000
   * Timestep: 2 fs
   * Backbone restraint: 40 kJ/mol/nm², side-chain: 0
   * Pressure coupling enabled (1 bar, isotropic Monte Carlo barostat)

> The equilibration runs are automatically managed in `utils/run.sh`.

---

### 3. Production Simulation

- **step5_production.inp**
   
   * Steps: 15,000,000 (30 ns)
   * Timestep: 2 fs
   * Backbone/side-chain restraints: none
   * NPT ensemble: temperature 303.15 K, pressure 1 bar, isotropic
   * Output frequency: nstout = 25000, nstdcd = 250000

- **How it works**
   
   * Managed through `utils/run.sh`.
   * Generates both **checkpoint (`.chk`)** and **restart (`.rst`)** files:
       - `.chk` files store full OpenMM internal state, including GPU-specific randomization, and allow a **true continuation** of the simulation on the same hardware.
     - `.rst` files contain coordinates and velocities only; they can be used to **start simulations on different GPUs** but do not fully preserve the GPU-specific state like `.chk`.

---
