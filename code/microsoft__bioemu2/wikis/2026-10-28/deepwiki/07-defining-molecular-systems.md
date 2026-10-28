# Defining Molecular Systems

<details>
<summary>Relevant source files</summary>

The following files were used as context for generating this wiki page:

- [docs/boltz.md](docs/boltz.md)
- [docs/system.md](docs/system.md)
- [src/bioemu2/assets/examples/systems/ake.yaml](src/bioemu2/assets/examples/systems/ake.yaml)
- [src/bioemu2/assets/examples/systems/chignolin.yaml](src/bioemu2/assets/examples/systems/chignolin.yaml)
- [src/bioemu2/assets/examples/systems/meteorin_like_7tlw.yaml](src/bioemu2/assets/examples/systems/meteorin_like_7tlw.yaml)
- [src/bioemu2/assets/examples/systems/myoglobin_101m.yaml](src/bioemu2/assets/examples/systems/myoglobin_101m.yaml)
- [src/bioemu2/assets/examples/systems/shp2_3tkz.yaml](src/bioemu2/assets/examples/systems/shp2_3tkz.yaml)

</details>



## Purpose and Scope

This page details how molecular systems are defined, parsed, and converted within BioEmu2. The molecular composition of a system—comprising proteins, RNA, DNA, ligands, modified residues, branched polymers, and non-standard covalent bonds—is encapsulated by the `System` class and serialized using YAML specifications. We cover native system construction via the Python API, asset examples, Boltz-format input conversion (`convert_boltz`), and how contact conditioning and covalent bonds integrate into the data interface.

Sources: `docs/system.md:1-23]`, `docs/boltz.md:1-10`

---

## Molecular System Representations (`System` and Entities)

The `System` class represents the aggregate molecular composition passed to BioEmu2 [docs/system.md:18-23](). It manages multiple molecular entities (`Protein`, `RNA`, `DNA`, `Ligand`, `Branched`), optional non-standard covalent bonds, and canonical-to-author chain ID mappings (`asym_id_to_author_id`) [docs/system.md:18-23]().

### Natural Language to Code Entity Mapping

The following diagram maps natural language descriptions of molecular assemblies to the corresponding BioEmu2 Python classes and YAML schema keys.

```mermaid
graph TD
    subgraph "NaturalLanguageSpace"
        NLProtein[""Single/Multi-Chain Protein Sequence""]
        NLLigand[""Small Molecule / Ion (e.g. ATP, SO4)""]
        NLMod[""Phosphorylated / Modified Residue""]
        NLBranched[""Glycan / Branched Polymer""]
    end

    subgraph "CodeEntitySpace"
        ClassSystem[""System(entities, bonds, asym_id_to_author_id)""]
        ClassProtein[""Protein(sequence, asym_ids, mods)""]
        ClassLigand[""Ligand(ccd_code / smiles, asym_ids)""]
        ClassMod[""ModifiedResidue(res_id, ccd_code)""]
        ClassBranched[""Branched(mon_seq, asym_ids)""]
    end

    NLProtein --> ClassProtein
    NLLigand --> ClassLigand
    NLMod --> ClassMod
    NLBranched --> ClassBranched
    
    ClassProtein --> ClassSystem
    ClassLigand --> ClassSystem
    ClassMod --> ClassProtein
    ClassBranched --> ClassSystem
```
Sources: `docs/system.md:18-105]`, `docs/system.md:125-148]`

### Programmatic System Construction
Systems can be built directly using Python constructors or via shorthand factory methods like `System.from_description()` [docs/system.md:33-104](). For instance, defining a homodimer or a protein-ligand complex requires explicit specification of chain identifiers (`asym_ids`) [docs/system.md:35-49]():

```python
from bioemu2 import Protein, Ligand, System

# Creating a protein-ligand system
protein_with_atp = System(
    [
        Protein(sequence="VIGPPGCGKGTQCEKLAEAGVEVVAA", asym_ids=["A"]),
        Ligand(ccd_code="ATP", asym_ids=["B"]),
        Ligand(ccd_code="MG", asym_ids=["C"]),
    ]
)
```
Sources: `docs/system.md:33-104]`, `docs/system.md:188-202]`

---

## YAML System Specifications and Asset Examples

BioEmu2 includes several curated system YAML files under `src/bioemu2/assets/examples/systems/` that demonstrate various structural configurations, including disulfide bonds, multi-chain assemblies, and ligand-bound complexes.

| Asset File | Description | Key Features |
| --- | --- | --- |
| `src/bioemu2/assets/examples/systems/ake.yaml` | E. coli adenylate kinase (PDB 4AKE) | Single-chain protein [src/bioemu2/assets/examples/systems/ake.yaml:1-8]() |
| `src/bioemu2/assets/examples/systems/chignolin.yaml` | Chignolin (CLN025) | Short beta-hairpin test system [src/bioemu2/assets/examples/systems/chignolin.yaml:1-8]() |
| `src/bioemu2/assets/examples/systems/meteorin_like_7tlw.yaml` | Mouse Meteorin-like NTR domain | Explicit non-standard disulfide bonds between `CYS` residues [src/bioemu2/assets/examples/systems/meteorin_like_7tlw.yaml:1-11]() |
| `src/bioemu2/assets/examples/systems/myoglobin_101m.yaml` | Sperm whale myoglobin | Multi-ligand system (`SO4`, `HEM`, `NBN`) with contact conditioning [src/bioemu2/assets/examples/systems/myoglobin_101m.yaml:1-19]() |
| `src/bioemu2/assets/examples/systems/shp2_3tkz.yaml` | SHP-2 N-SH2 domain bound to peptides | Multi-chain protein complex with author chain mapping and contact restraints [src/bioemu2/assets/examples/systems/shp2_3tkz.yaml:1-15]() |

Sources: `src/bioemu2/assets/examples/systems/ake.yaml:1-9]`, `src/bioemu2/assets/examples/systems/chignolin.yaml:1-9]`, `src/bioemu2/assets/examples/systems/meteorin_like_7tlw.yaml:1-12]`, `src/bioemu2/assets/examples/systems/myoglobin_101m.yaml:1-20]`, `src/bioemu2/assets/examples/systems/shp2_3tkz.yaml:1-16]`

---

## Boltz Format Inputs and Conversion (`convert_boltz`)

BioEmu2 provides built-in utilities to translate Boltz-1 and Boltz-2 format input YAML files into native BioEmu2 system YAML files via the CLI command `bioemu2 convert-boltz` or the Python function `convert_boltz()` [docs/boltz.md:1-10](), [docs/boltz.md:73-83]().

### Conversion Data Flow

```mermaid
graph TD
    BoltzYAML[""BoltzInputYAML (version: 1)""] --> CLIConv[""bioemu2 convert-boltz CLI / convert_boltz()""]
    CLIConv --> FieldFilter[""Filter & Map Supported Fields (sequences, mods, ligands, bonds, contacts)""]
    FieldFilter --> WarnDiscard[""Emit Warnings for Unsupported Fields (affinity, pocket, templates, msa)""]
    FieldFilter --> BioEmuYAML[""BioEmu2 System YAML Output""]
    BioEmuYAML --> LoadSys[""load_system_and_contacts()""]
    LoadSys --> EmuSessionInit[""Initialize EmuSession / Prediction""]
```
Sources: `docs/boltz.md:1-119]`, `docs/boltz.md:131-141]`

### Conversion Behavior and Supported Fields

During conversion:
- **Sequences**: Protein, RNA, and DNA sequences are preserved, including unknown residues (`X`, `N`) [docs/boltz.md:111]().
- **Modifications**: Polymer modifications (`ccd` at 1-based `position`) are mapped to `ModifiedResidue` entries [docs/boltz.md:113]().
- **Ligands**: SMILES or CCD-defined ligands are preserved [docs/boltz.md:114-115]().
- **MSAs**: Per-chain MSA entries inside Boltz inputs are **ignored with a warning**. MSAs must be handled via `get_msa(system, source=...)` [docs/boltz.md:67-71](), [docs/boltz.md:133]().

Sources: `docs/boltz.md:106-141]`

---

## Contact Conditioning and Covalent Bonds in Systems

Systems can declare inter-entity covalent bonds and contact constraints that guide sampling or model conditioning [docs/boltz.md:116-118](), [docs/system.md:155-178]().

### Defining Non-Standard Bonds
Disulfide bonds or covalent ligand attachments are explicitly listed under the `bonds` key:
```yaml
bonds:
  - atom1: {asym_id: A, res_id: 20, res_name: CYS, name: SG}
    atom2: {asym_id: A, res_id: 113, res_name: CYS, name: SG}
```
Sources: `src/bioemu2/assets/examples/systems/meteorin_like_7tlw.yaml:7-11]`

### Contact Conditioning Constraints
Contact constraints specify maximum allowable distances between tokens (residues or ligand atoms) to condition the generative model:
```yaml
contacts:
  - token1: {asym_id: A, res_id: 96}
    token2: {asym_id: C, atom_name: FE}
    max_distance: 6.0
```
These are parsed alongside the system via `load_system_and_contacts()` and passed into sampling sessions [docs/boltz.md:82-89]().

Sources: `src/bioemu2/assets/examples/systems/myoglobin_101m.yaml:16-20]`, `docs/boltz.md:116-118]`

---