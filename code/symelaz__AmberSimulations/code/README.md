# Amber Simulations of Soluble Proteins

This repository provides an **automated pipeline** to prepare and run **Amber molecular dynamics simulations** of soluble proteins, with or without a ligand, using **Python’s OpenMM library**. The workflow prepares all necessary input files, solvation, ligand parameters, and simulation configurations for equilibration and production runs.

Two force-field models are available with the same workflow; only the force field, water, ions and equilibration steps change.

---

## Force-field models

| `MODEL`                    | Protein                 | Water                          | Ions (default salt)                     | Equilibration | Typical use |
|----------------------------|-------------------------|--------------------------------|-----------------------------------------|---------------|-------------|
| `ff19SB_tip3p` (default)   | ff19SB                  | TIP3P                          | Joung-Cheatham, neutralize only (Na+/Cl−) | 2 steps       | folded proteins, protein–ligand |
| `ff19SB_tip4pew`           | ff19SB                  | TIP4P-Ew                       | Joung-Cheatham, neutralize + 0.15 M NaCl  | 2 steps       | folded proteins where TIP4P-Ew is wanted for consistency with other runs. **Caution:** ff19SB was parameterised/validated against OPC, not TIP4P-Ew; this pairing is not the validated one. |
| `a99SBdisp_tip4pd`         | a99SB-disp (+ extra backbone dihedrals) | TIP4P-D (a99SB-disp variant)    | CHARMM22, neutralize + 0.15 M KCl        | 3 steps       | disordered proteins, unfolding; protocol of Sarthak et al. 2023 (JCTC 19:3721) |

The a99SB-disp/TIP4P-D files come from `pipeline_amber_tip4pd` (T. Lemmin, M. McCully; MIT licence, see `utils/forcefields/a99SBdisp_tip4pd/LICENSE_pipeline_amber_tip4pd`).

### Water models

Exact charges and geometry of every water model this repo can build, as measured from the built
topology (`utils/check_topology.py` prints this same block for every build - compare it against this
table to be sure which water you actually got).

| Model | Sites | q(H) | q(M) | d(O-M) | O sigma (A) | O epsilon (kcal/mol) | `MODEL` / `WATER` |
|---|--:|--:|--:|--:|--:|--:|---|
| TIP3P | 3 | +0.417 | - | - | 3.1507 | 0.1520 | `ff19SB_tip3p` |
| TIP4P-Ew | 4 | +0.52422 | -1.04844 | 0.1250 A | 3.1644 | 0.16275 | `ff19SB_tip4pew` |
| TIP4P-D (authors, default) | 4 | +0.59 | -1.18 | 0.1546 A | 3.1650 | 0.238764 | `a99SBdisp_tip4pd` `WATER=authors` |
| TIP4P-D (AmberTools) | 4 | +0.58 | -1.16 | 0.1546 A | 3.1650 | 0.238764 | `a99SBdisp_tip4pd` `WATER=amber` |

`M` is the massless virtual site (EPW/EP), not a real atom; O itself carries zero charge in every
4-site model. TIP4P-D's larger O epsilon (about 47% more dispersion than TIP4P-Ew/OPC-class models) is
the deliberate fix of Piana/Donchev/Robustelli/Shaw (2015, JPCB 119:5113): standard water models
under-solvate disordered/unfolded chains relative to folded ones, biasing simulated ensembles toward
being too compact; TIP4P-D's extra water-protein dispersion corrects for that. This is also why the
`a99SBdisp_tip4pd` model is recommended for disordered proteins and unfolding over `ff19SB_tip4pew`.

**`ff19SB_tip4pew` is a non-validated pairing.** ff19SB's backbone torsions were parameterised and
benchmarked against OPC, not TIP4P-Ew. Both are rigid 4-site models with similar geometry, so it runs
fine, but the secondary-structure/folding balance was not tuned for this particular water. Use it when
TIP4P-Ew is wanted for consistency with other runs in the project; prefer `ff19SB_tip3p` (its own
validated pairing) or OPC otherwise.


---

## Build options

All options are environment variables in front of `build.sh`; anything not set uses the default.

| Variable | Default | Options | Meaning |
|----------|---------|---------|---------|
| `MODEL`  | `ff19SB_tip3p` | `ff19SB_tip3p`, `ff19SB_tip4pew`, `a99SBdisp_tip4pd` | force field + water + ions (table above) |
| `WATER`  | `authors` | `authors`, `amber` | `a99SBdisp_tip4pd` only. `authors`: water of pipeline_amber_tip4pd (H +0.59 / EP −1.18). `amber`: AmberTools `leaprc.water.tip4pd-a99SBdisp` (H +0.58 / EP −1.16). Same LJ. |
| `BOX`    | `rect`  | `rect`, `cubic` | box shape |
| `PAD`    | `15.0`  | number (Å) | minimum solvent padding around the solute |
| `MINBOX` | unset   | number (Å) | target **minimum box side**; overrides `PAD` with a per-axis padding so every side is ≥ `MINBOX` (needs `BOX=rect`). For elongated solutes that tumble, what matters is the shortest box side vs. the solute's diameter, not the padding. The padding is calibrated with one probe solvate pass, not guessed. |
| `PADFLOOR` | `8.0` | number (Å) | with `MINBOX`: smallest per-axis padding allowed |
| `SS_CUTOFF` | `2.5` | number (Å) | SG-SG distance below which two CYS are bonded as a disulfide |
| `CN_CUTOFF` | `2.0` | number (Å) | C(i)-N(i+1) distance above which the backbone is treated as broken (TER, no bond) |
| `SALT`   | `0` (ff19SB_tip3p) / `0.15` (ff19SB_tip4pew, a99SBdisp_tip4pd) | number (M) | salt added on top of neutralization (0 = neutralize only) |
| `CATION` | `Na+` (ff19SB_tip3p, ff19SB_tip4pew) / `K+` (a99SBdisp_tip4pd) | `Na+`, `K+` | cation for neutralization and salt (anion is always Cl−) |

Examples:

```bash
bash utils/build.sh input.pdb                                         # ff19SB / TIP3P, as before
MODEL=ff19SB_tip4pew bash utils/build.sh input.pdb                    # ff19SB / TIP4P-Ew, 0.15 M NaCl
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

**Faster engine (optional, ~2.5x faster):** same build, then

```bash
conda activate openmm_faster
bash utils/run.sh <gpu_ids> --faster
```

See "Faster engine (`--faster`)" below; setup in `Installation_readme.md`, section 5.

---

## Faster engine (`--faster`)

`run.sh --faster` runs the same simulation with **NVIDIA cuPotential**
([github.com/NVIDIA/cuPotential](https://github.com/NVIDIA/cuPotential)) in place of three OpenMM
components: the nonbonded force (PME electrostatics + LJ), the `LangevinMiddleIntegrator` and the
isotropic `MonteCarloBarostat`. Bonded terms, CMAP, restraints, constraints, HMR, reporters and all
output files stay the same. Without `--faster` nothing changes (classic engine).

**Same model.** The same cutoff, Ewald tolerance target, LJ dispersion correction, 1-4 exceptions,
temperature, friction, pressure and barostat frequency are used, read from the System the classic
path builds. cuPotential computes long-range electrostatics with its own splitting/spreading
functions instead of standard PME.

**Measured on our systems** (A100, ff19SB/TIP4P-Ew, 4 fs HMR, NPT, `run.sh` settings), against a
converged reference (stock OpenMM, double precision, Ewald tolerance 1e-6, same LJ model):

| system | classic | `--faster` | speedup | force error classic / faster |
|---|---|---|---|---|
| nanobody, 42,855 atoms | 744 ns/day | 1,933 ns/day | 2.6x | 8.1e-4 / 3.9e-4 |
| nanobody + antigen, 148,565 atoms | 290 ns/day | 732 ns/day | 2.5x | 8.1e-4 / 4.3e-4 |

200 ps NPT from the same equilibrated state: temperature on target and box volume within 0.15 %
of the classic run. Speedups depend on the GPU (NVIDIA reports ~2-4x across GPU generations).

**What `--faster` does per step:**

| `.inp` | nonbonded | integrator | barostat |
|---|---|---|---|
| `integrator = langevin_middle` (production, 4 fs HMR) | cuPotential | cuPotential LangevinMiddle | cuPotential (isotropic) |
| `integrator = langevin` (restrained equilibration) | cuPotential | stock Langevin | stock |

Precision is `mixed` with `--faster` (cuPotential's validated configuration), `single` without.
Each `.out` file states the engine, every component that was swapped, and the precision.

**Refused, with an explanation, instead of run differently:** non-PME electrostatics, `vdw` other
than `CutoffPeriodic`, implicit solvent, anisotropic/membrane barostats, simulated annealing,
`e14scale != 1`, `-ff charmm`. The plugin also does not support LJPME or triclinic boxes.

**Checkpoints.** A `.chk` file only loads into the engine that wrote it. `run.sh` records the
engine in `run_engine.txt`; if you switch engines between production chunks it restarts from the
portable `.rst` state instead, and says so. Equilibration -> production always goes through `.rst`.

**Tests** (GPU node, `conda activate openmm_faster`):

| script | what it answers | time |
|---|---|---|
| `python utils/check_faster.py [parm7 rst7]` | Is this machine/environment ready? Driver, CUDA, GPU, accuracy vs a converged reference, 1000 NPT steps. | ~1 min |
| `bash utils/test_faster_pipeline.sh <system_dir> [gpu]` | Does the whole `run.sh --faster` pipeline work on this system? Minimisation, restrained equilibration, production chunks, checkpoint restart, engine switch, outputs. Runs in a temporary copy. | a few min |
| `python utils/benchmark_faster.py <parm7> <step5_1.rst>` | How much faster, and is it as accurate, on this system and GPU? | 5-15 min |

All three end with a clear PASS/FAIL or numbers; `check_faster.py` and `test_faster_pipeline.sh`
exit with code 0 only if everything passed.

**Before switching a project over**, run one full replica with `--faster` and compare it with a
classic replica (RMSD/RMSF, energies, density). The tests above are short (up to 200 ps).

**One GPU per run:** two GPUs (`run.sh 0,1 --faster`) work but were slower in our test
(996 vs 1,380 ns/day).

**Requirements:** OpenMM 8.4 (the classic `openmm` env has 7.7; `--faster` uses its own env
`openmm_faster`), an NVIDIA GPU of compute capability 7.5-12.0 (RTX 4090, A100, H100: yes), and a
CUDA 13 driver. `run.sh --faster` stops immediately if the plugin is missing.

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
  * `run.sh` – Runs all equilibration steps and the production runs (`--faster`: cuPotential engine).
  * `openmm_run.py` – OpenMM execution script.
  * `omm_cupotential.py` – the `--faster` engine: swaps in cuPotential, refuses unsupported settings.
  * `install_faster.sh` – builds the `openmm_faster` environment with the cuPotential plugin.
  * `check_faster.py` – checks a machine/environment/system for `--faster` (run on a GPU node).
  * `test_faster_pipeline.sh` – end-to-end test of `run.sh --faster` on a built system (temporary copy).
  * `benchmark_faster.py` – classic vs `--faster` on a built system: accuracy, ns/day, NPT stability.
  * `batch_submit.sh` – Helper script for submitting multiple jobs on a SLURM server.
  * `job.sh` – Example SLURM submission script for running a simulation.
  * `environment.yml` – Conda environment configuration for OpenMM.

## Disulfide bonds

`utils/find_disulfides.py` runs automatically inside `build.sh`, right after the protein is extracted
and before tleap: it scans for CYS pairs whose SG atoms are closer than 2.5 A (override with
`SS_CUTOFF=...`), renames the matching residues CYS -> CYX (otherwise tleap adds a free thiol
hydrogen and leaves them unbonded), and writes the `bond receptor.N.SG receptor.M.SG` command tleap
needs. The residue numbers it uses are tleap's own sequential numbering within the loaded unit (1, 2,
3, ... in file order), not the PDB's chain/resSeq - the same numbering the "Notes" section below
describes for multi-chain systems.

`utils/check_topology.py` then verifies, in the finished topology, that each detected pair really is
bonded SG-SG, both residues are CYX, and neither carries a stray HG - not merely that no error was
printed - and that **every** CYX in the topology has an SG-SG bond (a CYX without a partner is a
sulfur radical). A build with zero cysteine pairs is a normal, unflagged outcome.

**Near misses are warned about, not bonded.** A CYS pair with SG-SG between the cutoff and 6 A is
printed as `WARNING near-miss` in `build/disulfides.log`. That is almost always a native disulfide
whose geometry was damaged upstream (loop remodelling, grafting, a predicted model). Fix the input
structure or raise `SS_CUTOFF` deliberately. Do not ignore the warning.

## Chain breaks and gaps

`utils/extract_protein.tcl` writes the protein with VMD, and VMD drops every `TER` of the input.
tleap only starts a new chain at a `TER`, so without them it bonds C(i) to N(i+1) across **every**
break, including intentional gaps (unresolved loops left open) and loops that loop closure failed to
close. tleap only prints `There is a bond of X angstroms between C and N atoms` and carries on.
Minimisation and equilibration then pull the stretched bond shut and distort the structure around it.

`utils/mark_chain_breaks.py` therefore runs right after extraction and puts a `TER` back at every
chain-ID change and at every real backbone gap (C(i)-N(i+1) > 2.0 A, `CN_CUTOFF=...`). The gaps are
listed in `build/chain_breaks.log`; they stay **open** (charged termini), they are not modelled.
`check_topology.py` fails the build if any inter-residue C-N bond in the final coordinates is > 2 A.
If a gap should be closed, model the loop before building.

## Build outputs

* `step3_input.parm7 / .rst7 / .pdb`, `restraints/prot_pos.txt` – simulation inputs.
* `build/chain_breaks.log` – chain boundaries and backbone gaps that were kept open (TER).
* `build/disulfides.leap`, `build/disulfides.json` – disulfides found by `find_disulfides.py` (bond
  commands and, respectively, the same pairs with original chain/resSeq labels and SG-SG distance).
* `build/topology_report.txt` – parameter check (all lines must be PASS), including a water-model
  fingerprint block to compare against the table above.
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
* Simulations require the **same GPU and CUDA version** if using checkpoint files, and the same engine (classic or `--faster`).
* OpenMM v7.7 was used for this but you can also try with newer versions of OpenMM.
* Systems with multiple chains are being renumbered by Amber with consecutive numbers You need to rename and renumber the different chains after md are finished to go back to your original chain/resid numbering.
* Funky ligands might fail to be parameterized by Amber automatically.
* Very large biological systems might fail to be parameterized by Amber. And for Membrane systems I also recommend to go with CHARMM force-field. DNA however is being handled better by Amber.
* Step5 production runs are controlled by the `.inp` file parameters (`dt`, `nstout`, `Target_ns`).
* Disulfide bonds are detected automatically (see "Disulfide bonds" above); only genuine SG-SG contacts within the cutoff are bonded, nothing else changes. Check `build/disulfides.log` for near-miss warnings.
* Backbone gaps are kept open, not modelled (see "Chain breaks and gaps" above).
* GAFF2 ligands have not been validated with TIP4P-D water.
* Extended / disordered structures (e.g. AlphaFold models with low pLDDT): with `BOX=rect` the chain can reach its periodic image when it rotates; consider `BOX=cubic` (large) or compacting the model first.
