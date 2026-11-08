# Amber Simulations of Soluble Proteins

This repository provides an **automated pipeline** to prepare and run **Amber molecular dynamics simulations** of soluble proteins, with or without a ligand, using **Python’s OpenMM library**. The workflow prepares all necessary input files, solvation, ligand parameters, and simulation configurations for equilibration and production runs.

Two force-field models are available with the same workflow; only the force field, water, ions and equilibration steps change.

---

## Force-field models

| `MODEL`                    | Protein                 | Water                          | Ions (default salt)                     | Equilibration | Typical use |
|----------------------------|-------------------------|--------------------------------|-----------------------------------------|---------------|-------------|
| `ff19SB_tip3p` (default)   | ff19SB                  | TIP3P                          | Joung-Cheatham, neutralize only (Na+/Cl−) | 2 steps       | folded proteins, protein–ligand |
| `a99SBdisp_tip4pd`         | a99SB-disp (+ extra backbone dihedrals) | TIP4P-D (a99SB-disp variant)    | CHARMM22, neutralize + 0.15 M KCl        | 3 steps       | disordered proteins, unfolding; protocol of Sarthak et al. 2023 (JCTC 19:3721) |

The a99SB-disp/TIP4P-D files come from `pipeline_amber_tip4pd` (T. Lemmin, M. McCully; MIT licence, see `utils/forcefields/a99SBdisp_tip4pd/LICENSE_pipeline_amber_tip4pd`).

---

## Build options

All options are environment variables in front of `build.sh`; anything not set uses the default.

| Variable | Default | Options | Meaning |
|----------|---------|---------|---------|
| `MODEL`  | `ff19SB_tip3p` | `ff19SB_tip3p`, `a99SBdisp_tip4pd` | force field + water + ions (table above) |
| `WATER`  | `authors` | `authors`, `amber` | `a99SBdisp_tip4pd` only. `authors`: water of pipeline_amber_tip4pd (H +0.59 / EP −1.18). `amber`: AmberTools `leaprc.water.tip4pd-a99SBdisp` (H +0.58 / EP −1.16). Same LJ. |
| `BOX`    | `rect`  | `rect`, `cubic` | box shape |
| `PAD`    | `15.0`  | number (Å) | minimum solvent padding around the solute |
| `SALT`   | `0` (ff19SB_tip3p) / `0.15` (a99SBdisp_tip4pd) | number (M) | salt added on top of neutralization (0 = neutralize only) |
| `CATION` | `Na+` (ff19SB_tip3p) / `K+` (a99SBdisp_tip4pd) | `Na+`, `K+` | cation for neutralization and salt (anion is always Cl−) |

Examples:

```bash
bash utils/build.sh input.pdb                                         # ff19SB / TIP3P, as before
MODEL=a99SBdisp_tip4pd bash utils/build.sh input.pdb                  # a99SB-disp / TIP4P-D, 0.15 M KCl
MODEL=a99SBdisp_tip4pd WATER=amber BOX=cubic bash utils/build.sh input.pdb
SALT=0.15 bash utils/build.sh input.pdb LIG -1 7.4                    # ff19SB + ligand + 0.15 M NaCl
```

---

## Usage

Copy the repository once per system, then:

**Protein-only system:**

```bash
conda activate ambertools23
bash utils/build.sh <input_pdb_file>
conda activate openmm
bash utils/run.sh <gpu_ids>
```

**Protein-ligand system:**

```bash
conda activate ambertools23
bash utils/build.sh <input_pdb_file> <ligand_resname> <ligand_charge> <ph>
conda activate openmm
bash utils/run.sh <gpu_ids>
```

* `<input_pdb_file>` – Path to the PDB file (protein or protein-ligand complex).
* `<ligand_resname>` – Three-letter residue name of the ligand in the PDB file.
* `<ligand_charge>` – Total charge of the ligand (integer or decimal).
* `<ph>` – pH value of the system (normally 7.4).
* `<gpu_ids>` – GPU index/indices for OpenMM (default `0`).

`build.sh` copies the OpenMM inputs of the chosen model from `inputs/<MODEL>/` to the main folder.
Edit those copies (`step4.*_equilibration.inp`, `step5_production.inp`) for a specific system; `inputs/` holds the templates.

---

## Repository Structure

* `inputs/<MODEL>/` – OpenMM input templates per model.
  * `step4.1_equilibration.inp` – minimization + restrained MD (NVT).
  * `step4.2_equilibration.inp` (`step4.3` for `a99SBdisp_tip4pd`) – further equilibration.
  * `step5_production.inp` – production MD (30 ns, unrestrained, NPT; identical for both models).
* `utils/` – Utility scripts for system preparation, equilibration, production runs, and server submission.
  * `build.sh` – Builds the system (ligand parameters, solvation, ions, restraints, checks).
  * `build_protein.leap`, `build_complex.leap` – tleap templates for **protein-only** and **protein-ligand** systems.
  * `forcefields/<MODEL>/` – force-field files and leaprc of each model.
  * `add8.py` – adds the a99SB-disp extra backbone dihedrals (`a99SBdisp_tip4pd` only).
  * `check_topology.py` – checks water, protein and ion parameters in the final topology; the build stops if a check fails.
  * `run.sh` – Runs all equilibration steps and the production runs.
  * `openmm_run.py` – OpenMM execution script.
  * `batch_submit.sh` – Helper script for submitting multiple jobs on a SLURM server.
  * `job.sh` – Example SLURM submission script for running a simulation.
  * `environment.yml` – Conda environment configuration for OpenMM.

## Build outputs

* `step3_input.parm7 / .rst7 / .pdb`, `restraints/prot_pos.txt` – simulation inputs.
* `build/topology_report.txt` – parameter check (all lines must be PASS).
* `build/build_settings.txt`, `build/build_system.leap`, `build/tleap_build.log` – exact build record.
* `build/loaded_files.txt` – every file tleap loaded, with its path.

---

## Workflow Diagram
<img src="https://github.com/symelaz/AmberSimulations/raw/HEAD/figures/workflow.png" alt="Amber Workflow" width="600">

---

## Hydrogen mass repartitioning (HMR)

* In `step5_production.inp` set `hmr = yes`, `dt = 0.004` and `integrator = langevin_middle`.
  Applied at run time by `utils/omm_hmr.py` (topology unchanged; water excluded). About 2× faster.
* Leave `hmr = no` for unfolding/kinetics runs.

| `.inp` key   | Default    | Options |
|--------------|------------|---------|
| `hmr`        | `no`       | `yes`, `no` |
| `hmr_mass`   | `3.024`    | hydrogen mass (Da) |
| `integrator` | `langevin` | `langevin`, `langevin_middle` (use with HMR / 4 fs) |

---

## Server Submission
* Pinned variables in `utils/job.sh` for server submission:

  * GPU devices: `--gres=gpu:rtx4090:2`
  * Job name: `test`
  * Partition: `gpu`
  * Conda environment: `openmm`
  * Memory: `--mem-per-cpu=8G`
  * Excluded nodes: `--exclude=gnode23`

* Use `utils/batch_submit.sh utils/job.sh 10 [dependency id]` to submit multiple dependent jobs.

---

## Notes / Best Practices

* Ensure the **PDB file is pre-optimized** (missing hydrogens may be added automatically).
* Simulations require the **same GPU and CUDA version** if using checkpoint files.
* OpenMM v7.7 was used for this but you can also try with newer versions of OpenMM.
* Systems with multiple chains are being renumbered by Amber with consecutive numbers You need to rename and renumber the different chains after md are finished to go back to your original chain/resid numbering.
* Funky ligands might fail to be parameterized by Amber automatically.
* Very large biological systems might fail to be parameterized by Amber. And for Membrane systems I also recommend to go with CHARMM force-field. DNA however is being handled better by Amber.
* Step5 production runs are controlled by the `.inp` file parameters (`dt`, `nstout`, `Target_ns`).
* Disulfide bonds are not detected automatically.
* GAFF2 ligands have not been validated with TIP4P-D water.
* Extended / disordered structures (e.g. AlphaFold models with low pLDDT): with `BOX=rect` the chain can reach its periodic image when it rotates; consider `BOX=cubic` (large) or compacting the model first.
