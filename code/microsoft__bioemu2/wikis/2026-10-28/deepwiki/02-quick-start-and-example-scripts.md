# Quick Start and Example Scripts

<details>
<summary>Relevant source files</summary>

The following files were used as context for generating this wiki page:

- [README.md](README.md)
- [examples/assets/1ten_protein_constraint.pdb](examples/assets/1ten_protein_constraint.pdb)
- [examples/boltz/protein_ligand.yaml](examples/boltz/protein_ligand.yaml)
- [examples/sample_a_glycoprotein.py](examples/sample_a_glycoprotein.py)
- [examples/sample_a_phosphorylated_complex.py](examples/sample_a_phosphorylated_complex.py)
- [examples/sample_protein_ligand_constraints.py](examples/sample_protein_ligand_constraints.py)
- [examples/sample_protein_protein_constraints.py](examples/sample_protein_protein_constraints.py)
- [examples/sample_structures.py](examples/sample_structures.py)
- [examples/simulate_a_trajectory.py](examples/simulate_a_trajectory.py)

</details>



## Purpose and Scope
This page details the command-line interface (CLI) quick start, the primary Python session manager (`EmuSession`), functional shortcuts, and a technical walkthrough of the example scripts and asset inputs provided in the `examples/` directory [README.md:98-168](). It outlines how biomolecular systems are defined, how inference and propagation are initiated, and how auxiliary assets (such as YAML configs and reference PDBs) map to runtime entities.

---

## Command Line Interface (CLI) Quick Start
The `bioemu2` command-line utility exposes subcommands for equilibrium sampling, structure prediction, and time propagation [README.md:102-127](). Executing these commands interacts directly with `src/bioemu2/cli.py` (referenced across the documentation) to parse arguments, load models, build molecular topologies, and invoke sampling loops.

### Core CLI Subcommands
* **`sample-equilibrium`**: Samples diverse conformations representing the equilibrium ensemble of a target system [README.md:102-108]().
* **`predict-structure`**: Generates high-accuracy folded static structures from sequence or PDB inputs [README.md:110-117]().
* **`propagate`**: Simulates forward time-lagged trajectory steps starting from initial structures (`x0`) using implicit transfer operator models [README.md:119-127]().

```bash
# Sample equilibrium conformations for Adenylate Kinase (PDB 1AKE)
bioemu2 sample-equilibrium \
  --pdb-id 1AKE \
  --out outputs/adk

# Predict a folded structure for Ubiquitin (PDB 1UBQ)
bioemu2 predict-structure \
  --pdb-id 1UBQ \
  --num-samples 1 \
  --out outputs/ubiquitin

# Propagate sampled conformations forward by 1 ns
bioemu2 propagate \
  --x0-cif outputs/adk/topology.cif \
  --x0-xtc outputs/adk/samples.xtc \
  --lag-time-ns 1 \
  --out outputs/adk-propagated
```

Sources: [README.md:102-129]()

---

## Python API: `EmuSession` and Functional Shortcuts
The primary object-oriented API for executing inference tasks is `EmuSession`, complemented by functional shortcuts like `sample_equilibrium` and `simulate_trajectory` [README.md:131-168](). 

### `EmuSession` Initialization and Data Flow
`EmuSession` binds a `System` object, a loaded model instance, and optional conditioning parameters (such as MSAs or spatial `Contact` constraints) into a unified context [README.md:131-156]().

```python
from pathlib import Path
from bioemu2 import EmuSession, System, load_model, save_frames

# Load the default pretrained model
model = load_model()

# 1. Structure prediction mode for Ubiquitin
ubiquitin = System.from_pdb_id("1UBQ")
prediction = EmuSession(ubiquitin, model=model, mode="structure").sample(num_samples=1)
save_frames(prediction, Path("outputs/ubiquitin"))

# 2. Equilibrium sampling mode for Adenylate Kinase
adk = System.from_pdb_id("1AKE")
emu = EmuSession(adk, model=model, mode="equilibrium")
samples = emu.sample(num_samples=8)
save_frames(samples, Path("outputs/adk"))

# 3. Time propagation mode (1 ns lag time)
endpoints = emu.propagate(samples, lag_time_ns=1.0)
save_frames(endpoints, Path("outputs/adk-propagated"))
```

Sources: [README.md:131-160](), [examples/simulate_a_trajectory.py:37-72]()

```mermaid
title "CLI and EmuSession Data Flow to Code Entities"
graph TD
    CLI["CLI Command (bioemu2 sample-equilibrium)"] --> CLIMod["src/bioemu2/cli.py"]
    UserPy["User Python Script"] --> Session["EmuSession Class"]
    CLIMod --> Session
    Session --> Model["load_model()"]
    Session --> Sys["System Object (System.from_pdb_id)"]
    Sys --> MSA["MSA Processing Module"]
    Model --> Sampler["EDMSampler Class"]
    MSA --> Sampler
    Session --> Sampler
    Sampler --> Output["save_frames() (topology.cif & samples.xtc)"]
```
Sources: [README.md:102-160]()

---

## Asset Inputs and Molecular Specification Formats
BioEmu2 accepts diverse input assets to specify molecular topologies, including PDB identifiers, custom FASTA/YAML configurations, and reference coordinate files [README.md:102-168]().

### Configuration Formats
1. **PDB Identifiers**: Parsed dynamically via `System.from_pdb_id("1AKE")` [README.md:148]().
2. **YAML System Specifications**: Complex compositions containing explicit protein chains, ligands, and sequences are defined in YAML [examples/boltz/protein_ligand.yaml:1-8]().
3. **PDB Design References**: Explicit coordinate assets (e.g., `examples/assets/1ten_protein_constraint.pdb`) supply spatial reference frames for restraint and comparison workflows [examples/sample_protein_protein_constraints.py:49-51]().

```yaml
version: 1
sequences:
  - protein:
      id: [A, B]
      sequence: RMKQLEDKVEELLSKNYHLENEVARLKKLVGER  # GCN4 leucine zipper (PDB 2ZTA)
  - ligand:
      id: L
      smiles: CCO
```
Sources: [examples/boltz/protein_ligand.yaml:1-8](), [examples/sample_protein_protein_constraints.py:49-51]()

---

## Walkthrough of Example Scripts
The `examples/` directory contains complete, executable scripts illustrating various advanced features, post-translational modifications, and steering workflows.

### Summary of Example Scripts

| Script Path | Primary Focus / Workflow | Key Classes / Functions Used |
| :--- | :--- | :--- |
| `examples/sample_structures.py` | Basic protein structure prediction and equilibrium sampling of Adenylate Kinase | `System.from_pdb_id`, `EmuSession`, `load_model`, `save_frames`, `view_trajectory` [examples/sample_structures.py:21-107]() |
| `examples/sample_a_glycoprotein.py` | Building and sampling a glycosylated Fc fragment with explicit branched glycans and disulfides | `Protein`, `Branched`, `ModifiedResidue`, `System`, `EmuSession`, `Bond` [examples/sample_a_glycoprotein.py:45-163]() |
| `examples/sample_a_phosphorylated_complex.py` | Sampling protein complexes containing post-translational modifications (pTyr / `PTR`) | `Protein`, `ModifiedResidue`, `System`, `EmuSession` [examples/sample_a_phosphorylated_complex.py:28-63]() |
| `examples/sample_protein_ligand_constraints.py` | Conditioning complex predictions using spatial protein-ligand contact constraints | `Protein`, `Ligand`, `Contact`, `ContactToken`, `EDMSampler`, `EmuSession` [examples/sample_protein_ligand_constraints.py:32-153]() |
| `examples/sample_protein_protein_constraints.py` | Biasing protein-protein binder docking interfaces via residue-residue distance bounds | `MSA`, `Contact`, `ContactToken`, `EDMSampler`, `EmuSession`, `System` [examples/sample_protein_protein_constraints.py:29-94]() |
| `examples/simulate_a_trajectory.py` | Simulating multi-step time-lagged trajectories using implicit transfer operators | `System`, `sample_equilibrium`, `simulate_trajectory`, `save_frames` [examples/simulate_a_trajectory.py:24-72]() |

Sources: [examples/sample_structures.py:21-107](), [examples/sample_a_glycoprotein.py:45-163](), [examples/sample_a_phosphorylated_complex.py:28-63](), [examples/sample_protein_ligand_constraints.py:32-153](), [examples/sample_protein_protein_constraints.py:29-94](), [examples/simulate_a_trajectory.py:24-72]()

```mermaid
title "Example Script Processing Pipeline to Code Entities"
graph TD
    YAML["examples/boltz/protein_ligand.yaml"] --> Parser["System Parsing Module"]
    PDBAsset["examples/assets/1ten_protein_constraint.pdb"] --> Parser
    Parser --> SysInst["System Class Instantiation"]
    SysInst --> Sess["EmuSession Initialization"]
    Sess --> ModLoad["load_model()"]
    ModLoad --> EDMS["EDMSampler / propagate()"]
    EDMS --> OutFiles["Output Files (topology.cif, samples.xtc)"]
```
Sources: [examples/boltz/protein_ligand.yaml:1-8](), [examples/sample_protein_protein_constraints.py:49-60](), [examples/simulate_a_trajectory.py:37-72]()

---

## Summary of Outputs and Trajectory Persistence
All execution paths across the CLI and Python API serialize output frames through `save_frames()`, which exports trajectory coordinate arrays into standard Binary Portable Trajectory (`.xtc`) files paired with Crystallographic Information File (`.cif`) topology definitions [README.md:158-160](). These saved structures can be reloaded directly using MDTraj for downstream analysis [examples/sample_structures.py:90-99]().

Sources: [README.md:158-160](), [examples/sample_structures.py:90-99]()

---