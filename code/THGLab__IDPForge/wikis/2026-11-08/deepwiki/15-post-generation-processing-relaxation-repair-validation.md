# Post-Generation Processing: Relaxation, Repair & Validation

<details>
<summary>Relevant source files</summary>

The following files were used as context for generating this wiki page:

- [data/AF-P05231_ndr.npz](data/AF-P05231_ndr.npz)
- [idpforge/misc.py](idpforge/misc.py)
- [idpforge/utils/relax.py](idpforge/utils/relax.py)
- [idpforge/utils/structure_repair.py](idpforge/utils/structure_repair.py)
- [idpforge/utils/structure_validation.py](idpforge/utils/structure_validation.py)
- [relax_raw.py](relax_raw.py)

</details>



After the diffusion model generates raw coordinates, the IDPForge pipeline performs a multi-stage post-processing sequence to ensure physical realism, stereochemical correctness, and topological integrity. This process transforms raw neural network outputs into energy-minimized, validated PDB ensembles.

### Overview of the Post-Generation Pipeline

The post-generation lifecycle is managed primarily by `output_to_pdb()` in `idpforge/misc.py` [119-135](). The pipeline consists of:
1.  **Conversion**: Transforming `atom14` internal representations to `atom37` PDB format [139-140]().
2.  **Pre-Relaxation Screening**: Checking backbone continuity and NaN coordinates before expensive minimization [192-201]().
3.  **AMBER Relaxation**: Using OpenMM to minimize the structure, typically with the folded domain residues restrained [idpforge/utils/relax.py:21-25]().
4.  **Structural Repair**: Fixing D-amino acid flips and Histidine naming artifacts introduced by the forcefield [idpforge/utils/structure_repair.py:35-178]().
5.  **Validation Suite**: A comprehensive check for clashes, chirality, bond integrity, and topological knots [idpforge/utils/structure_validation.py:206-258]().

### Pipeline Workflow Diagram

The following diagram illustrates the flow from raw tensor output to a validated PDB file, mapping the logical steps to the specific code entities responsible for them.

**Post-Generation Processing Flow**
```mermaid
graph TD
    subgraph "Raw Output Stage"
        A["Model Output Tensors"] --> B["output_to_pdb()"]
    end

    subgraph "Relaxation & Repair Stage"
        B --> C["relax_protein()"]
        C -- "AMBER Minimization" --> D["relax.AmberRelaxation"]
        D --> E{"Needs Repair?"}
        E -- "D-isomer / HIS naming" --> F["repair_chirality()"]
        E -- "D-isomer / HIS naming" --> G["fix_histidine_naming()"]
        F --> H["Re-Relaxation"]
        G --> H
    end

    subgraph "Validation Stage"
        H --> I["validate_structure_post_relax()"]
        E -- "Clean" --> I
        I --> J["check_chirality()"]
        I --> K["check_bond_integrity()"]
        I --> L["check_clashes()"]
        I --> M["AlphaKnot2 Hybrid Topology"]
    end

    subgraph "Finalization"
        M -- "PASS" --> N["*_validated.pdb"]
        M -- "FAIL" --> O["Discard Conformer"]
    end
```
**Sources:** `idpforge/misc.py` [119-300](), `idpforge/utils/relax.py` [21-92](), `idpforge/utils/structure_validation.py` [206-258]().

---

### AMBER Relaxation (relax_protein)

The relaxation stage utilizes the `AmberRelaxation` wrapper from OpenFold to perform energy minimization. A key feature of IDPForge is the use of a `viol_mask` to gate acceptance: the system calculates the fraction of residues with stereochemical violations and rejects conformers exceeding the `viol_threshold` (default 0.02) [idpforge/utils/relax.py:22-63]().

For multi-domain proteins or IDRs in context, the `exclude_residues` configuration allows the folded domains to remain rigid (via harmonic restraints) while the disordered regions are free to minimize [relax_raw.py:58-61]().

For details, see [AMBER Relaxation (relax_protein)](#4.4.1).

**Sources:** `idpforge/utils/relax.py` [21-92](), `relax_raw.py` [47-62]().

---

### Structural Repair & Re-Relaxation

Neural generation and forcefield minimization can occasionally produce stereochemical artifacts. IDPForge includes specialized repair utilities:
*   **Chirality Repair**: `repair_chirality()` detects D-amino acids by calculating the scalar triple product of the N-CA-C-CB volume. If a flip is detected, it reflects the sidechain atoms across the N-CA-C plane to restore L-isomer geometry [idpforge/utils/structure_repair.py:35-154]().
*   **Histidine Naming**: `fix_histidine_naming()` addresses a known issue where AMBER relaxation scrambles atom identities in the HIS imidazole ring. It reassigns names based on spatial geometry and hydrogen bonding patterns [idpforge/utils/structure_repair.py:156-178]().

If repairs are applied, the structure is subjected to a second round of relaxation to ensure the new coordinates are energetically stable [relax_raw.py:153-166]().

**Sources:** `idpforge/utils/structure_repair.py` [35-178](), `relax_raw.py` [123-166]().

---

### Structural Validation Suite

The final quality gate is the `validate_structure_post_relax()` function. This suite performs high-fidelity checks to ensure the ensemble is suitable for scientific analysis.

| Check | Code Entity | Description |
| :--- | :--- | :--- |
| **Chirality** | `check_chirality()` | Uses scalar triple product to ensure L-amino acid geometry [idpforge/utils/structure_validation.py:99-119](). |
| **Bonds** | `check_bond_integrity()` | Graph-based check for stretched or broken covalent bonds [idpforge/utils/structure_validation.py:130-190](). |
| **Clashes** | `check_clashes()` | Uses `cKDTree` for fast VdW overlap detection [idpforge/utils/structure_validation.py:192-204](). |
| **Topology** | `AlphaKnot2` | Hybrid Alexander/HOMFLY polynomial check to detect unwanted knots [idpforge/utils/structure_validation.py:23-27](). |

For details, see [Structural Validation Suite](#4.4.2).

**Sources:** `idpforge/utils/structure_validation.py` [1-258]().

---

### Code Mapping: Logic to Entities

This diagram maps the high-level processing requirements to the specific classes and functions in the `idpforge` codebase.

**System Entity Mapping**
```mermaid
classDiagram
    class OutputManager {
        idpforge/misc.py
        output_to_pdb()
    }
    class Minimizer {
        idpforge/utils/relax.py
        relax_protein()
        AmberRelaxation
    }
    class RepairToolkit {
        idpforge/utils/structure_repair.py
        repair_chirality()
        fix_histidine_naming()
    }
    class ValidationSuite {
        idpforge/utils/structure_validation.py
        validate_structure_post_relax()
        check_bond_integrity()
        check_clashes()
    }

    OutputManager --> Minimizer : "1. Minimizes Raw PDB"
    Minimizer --> RepairToolkit : "2. Corrects Artifacts"
    RepairToolkit --> Minimizer : "3. Re-minimizes (if repaired)"
    Minimizer --> ValidationSuite : "4. Final Quality Gate"
```
**Sources:** `idpforge/misc.py` [119-135](), `idpforge/utils/relax.py` [21-25](), `idpforge/utils/structure_repair.py` [35-178](), `idpforge/utils/structure_validation.py` [206-258]().

---