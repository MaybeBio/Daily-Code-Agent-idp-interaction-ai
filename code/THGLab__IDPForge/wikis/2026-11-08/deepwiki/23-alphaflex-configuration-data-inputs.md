# AlphaFlex Configuration & Data Inputs

<details>
<summary>Relevant source files</summary>

The following files were used as context for generating this wiki page:

- [AlphaFlex/Data_Inputs/AlphaFlex_database_Nov2025.json](AlphaFlex/Data_Inputs/AlphaFlex_database_Nov2025.json)
- [AlphaFlex/Data_Inputs/knot_screening.json](AlphaFlex/Data_Inputs/knot_screening.json)
- [AlphaFlex/README.md](AlphaFlex/README.md)
- [AlphaFlex/config.py](AlphaFlex/config.py)
- [AlphaFlex/utils/file_ops.py](AlphaFlex/utils/file_ops.py)
- [AlphaFlex/utils/smart_scoring.py](AlphaFlex/utils/smart_scoring.py)

</details>



This page documents the centralized configuration system, input databases, and reference schemas utilized by the AlphaFlex pipeline. These resources facilitate large-scale IDR ensemble generation by providing structural metadata, topological screening, and pipeline control constants.

## Pipeline Configuration (`AlphaFlex/config.py`)

The `config.py` file serves as the central source of truth for all pipeline constants, file paths, and hyperparameters across Steps 1 through 4 of the AlphaFlex workflow [AlphaFlex/config.py:1-122]().

### Global Paths and Execution
The configuration defines the environment and base directories for the pipeline.
* `PYTHON_EXEC`: The current Python interpreter used for subprocess calls [AlphaFlex/config.py:14-14]().
* `PROJECT_ROOT` and `PARENT_DIR`: Used to resolve paths relative to the `AlphaFlex` directory and the base `IDPForge` installation [AlphaFlex/config.py:15-16]().
* `INPUT_DATA_DIR`: Points to `AlphaFlex/Data_Inputs`, containing the master JSON databases [AlphaFlex/config.py:18-18]().

### Key Pipeline Hyperparameters
The pipeline behavior is controlled by several threshold and attempt constants:

| Parameter | Default Value | Purpose |
| :--- | :--- | :--- |
| `TEMPLATE_N_CONFS` | 200 | Number of template conformations in Step 2 [AlphaFlex/config.py:66-66](). |
| `SAMPLE_N_CONFS` | 10 | Target conformers per IDR in Step 3 [AlphaFlex/config.py:81-81](). |
| `SAMPLE_BATCH_SIZE` | 6 | Diffusion batch size for inference [AlphaFlex/config.py:82-82](). |
| `STITCH_N_CONFORMERS`| 10 | Final ensemble size in Step 4 [AlphaFlex/config.py:104-104](). |
| `RELAX_STIFFNESS` | 10.0 | Harmonic restraint strength for AMBER [AlphaFlex/config.py:110-110](). |
| `SAMPLE_JUNCTION_KAPPA`| 0.12 | Curvature gate for IDR-Folded junctions [AlphaFlex/config.py:89-89](). |

### Adaptive Clash Scoring
AlphaFlex utilizes an adaptive thresholding mechanism for stitching, defined by `STITCH_BASE_CLASH_THRESHOLD` (10.0) and `STITCH_CLASH_INCREMENT` (5.0) [AlphaFlex/config.py:121-122](). These are used by `get_smart_threshold` to relax constraints if valid conformers are not found within initial attempts [AlphaFlex/utils/smart_scoring.py:4-6]().

**Sources:** [AlphaFlex/config.py:1-122](), [AlphaFlex/utils/smart_scoring.py:1-6]()

---

## Master Database (`AlphaFlex_database_Nov2025.json`)

This database contains the structural metadata for the AlphaFold2 9606 Human v4 proteome, enabling the classification of IDRs based on their predicted context [AlphaFlex/README.md:9-10]().

### Schema Structure
Each entry is keyed by UniProt ID and contains:
* `idrs`: A list of `[start, end]` residue ranges (1-indexed) identifying disordered regions [AlphaFlex/Data_Inputs/AlphaFlex_database_Nov2025.json:3-8]().
* `mean_pae`: A dictionary mapping domain pairs (e.g., `D1-F1`, `F1-F2`) to their mean Predicted Aligned Error in Angstroms [AlphaFlex/Data_Inputs/AlphaFlex_database_Nov2025.json:23-27]().
* `interactions`: A list of folded domain pairs (e.g., `["F1", "F2"]`) where the mean PAE is $< 15$ Å, indicating a stable interface [AlphaFlex/Data_Inputs/AlphaFlex_database_Nov2025.json:49-54]().

**Sources:** [AlphaFlex/Data_Inputs/AlphaFlex_database_Nov2025.json:1-331](), [AlphaFlex/README.md:9-10]()

---

## Knot Screening Reference (`knot_screening.json`)

This reference file contains pre-calculated topological data for all proteins containing folded domains. It is used during structure validation to ensure that generated ensembles maintain the "native" topology predicted by AlphaFold2 [AlphaFlex/README.md:14-15]().

### Data Fields
* `label`: High-level topological classification (e.g., "None" for unknotted) [AlphaFlex/Data_Inputs/knot_screening.json:2-3]().
* `reason`: Descriptive string explaining the classification (e.g., "F1=Unknot", "NoFoldedDomain") [AlphaFlex/Data_Inputs/knot_screening.json:4-4]().
* `domains`: An array of domain-specific topological reports, including the residue `range`, the specific `knot` type, and a confidence score in the `reason` field [AlphaFlex/Data_Inputs/knot_screening.json:10-20]().

**Sources:** [AlphaFlex/Data_Inputs/knot_screening.json:1-359](), [AlphaFlex/README.md:14-15]()

---

## Data Flow: Configuration to Execution

The following diagram illustrates how the configuration constants and JSON databases are consumed by the AlphaFlex pipeline scripts.

### Configuration to Code Mapping
```mermaid
graph TD
    subgraph "Natural Language Space"
        Config["Pipeline Configuration"]
        DB["Master Database"]
        KnotRef["Knot Reference"]
    end

    subgraph "Code Entity Space"
        CFG_PY["AlphaFlex/config.py"]
        STEP1["Step_1_case_label.py"]
        STEP2["Step_2_mk_ldr_template.py"]
        STEP3["Step_3_sample_conformer.py"]
        STEP4["Step_4_ldr_stitch.py"]
        
        DB_FILE["AlphaFlex_database_Nov2025.json"]
        KNOT_FILE["knot_screening.json"]
    end

    CFG_PY -- "Provides constants to" --> STEP1
    CFG_PY -- "Provides constants to" --> STEP2
    CFG_PY -- "Provides constants to" --> STEP3
    CFG_PY -- "Provides constants to" --> STEP4

    DB_FILE -- "MASTER_DB_PATH" --> STEP1
    KNOT_FILE -- "KNOT_SCREENING_PATH" --> STEP3
    KNOT_FILE -- "KNOT_SCREENING_PATH" --> STEP4
```
**Sources:** [AlphaFlex/config.py:19-22](), [AlphaFlex/README.md:1-20]()

---

## Reference Data: Residue Counts

The file `AF2_9606_HUMAN_v4_num_residues.json` (referenced via `LENGTH_REF_PATH`) provides a lookup table for total protein lengths [AlphaFlex/config.py:20-20](). This is primarily used in Step 1B for length-based subset filtering [AlphaFlex/README.md:11-11]().

### Subset Filtering Logic
The constants `SUBSET_MIN_LENGTH` and `SUBSET_MAX_LENGTH` define the global filter applied to the UniProt IDs found in the residue count database [AlphaFlex/config.py:40-41]().

**Sources:** [AlphaFlex/config.py:34-50](), [AlphaFlex/README.md:64-80]()

---

## File Operations and Staging (`AlphaFlex/utils/file_ops.py`)

AlphaFlex employs specialized file utilities to manage the high volume of PDB files generated during sampling and stitching.

### Key Utilities
* `atomic_transfer`: Uses `shutil.copy2` and `os.replace` to ensure files are moved to final directories only when fully written, preventing partial reads by downstream steps [AlphaFlex/utils/file_ops.py:8-18]().
* `rename_and_clean_final_directory`: Sanitizes output directories by removing temporary files (starting with `.tmp_` or ending in `.tmp`/`.bak`) and renumbering valid `*_relaxed.pdb` files into a continuous sequence [AlphaFlex/utils/file_ops.py:24-87]().
* `cleanup_staging_area`: Removes intermediate `raw_att*.pdb` files from the conformer generation stage [AlphaFlex/utils/file_ops.py:92-122]().

### Staging Workflow
```mermaid
graph LR
    subgraph "Staging Area"
        Raw["raw_attN.pdb"]
        Tmp[".tmp_N_relaxed.pdb"]
    end
    
    subgraph "Final Directory"
        Final["N_relaxed.pdb"]
    end

    STEP3 -- "Generates" --> Raw
    STEP4 -- "Relaxes" --> Tmp
    ATOMIC["atomic_transfer()"] -- "Moves" --> Tmp
    Tmp -- "Renames to" --> Final
    CLEAN["cleanup_staging_area()"] -- "Deletes" --> Raw
```
**Sources:** [AlphaFlex/utils/file_ops.py:1-130]()

---