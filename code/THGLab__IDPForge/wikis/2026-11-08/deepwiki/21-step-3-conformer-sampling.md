# Step 3: Conformer Sampling

<details>
<summary>Relevant source files</summary>

The following files were used as context for generating this wiki page:

- [AlphaFlex/Step_3_sample_conformer.py](AlphaFlex/Step_3_sample_conformer.py)
- [AlphaFlex/config.py](AlphaFlex/config.py)

</details>



Step 3 of the AlphaFlex pipeline is responsible for generating large-scale, validated conformer pools for each Intrinsically Disordered Region (IDR) identified in Step 1. It orchestrates the execution of the core IDPForge sampling scripts, manages state persistence to allow for interrupted runs, and applies rigorous structural filters—including knot screening and curvature gates—before conformers are accepted into the pool.

### Workflow Orchestration

The entrypoint for this stage is `Step_3_sample_conformer.py`. It functions as a high-level manager that dispatches sampling tasks to `sample_ldr.py` (for IDRs in a folded context) or `sample_idp.py` (for fully disordered proteins).

#### Key Components
*   **State Persistence**: Uses `.step3_state.json` within each IDR's output directory to track `total_attempts`, allowing the pipeline to resume generation after crashes or HPC timeouts [AlphaFlex/Step_3_sample_conformer.py:27-57]().
*   **Parallelization**: Supports HPC batching via the `--split_index` argument, allowing the total protein list to be divided across multiple jobs [AlphaFlex/Step_3_sample_conformer.py:313-324]().
*   **Validation Gating**: Enforces quality controls defined in `AlphaFlex/config.py`, such as `SAMPLE_JUNCTION_KAPPA` (junction curvature) and `SAMPLE_FOLD_CURV_RATIO` [AlphaFlex/config.py:86-89]().

#### Data Flow: Generation to Validation
The following diagram illustrates how `Step_3_sample_conformer.py` interfaces with the core sampling logic and the validation suite.

**Diagram: Conformer Sampling Control Flow**
```mermaid
graph TD
    subgraph "AlphaFlex Orchestrator"
        A["Step_3_sample_conformer.py"] -- "calls" --> B["run_idr_workflow()"]
        B -- "manages" --> C[".step3_state.json"]
    end

    subgraph "Core Sampling Subprocess"
        B -- "executes" --> D["sample_ldr.py"]
        D -- "loads" --> E["idpforge/model.py"]
        D -- "validates" --> F["idpforge/utils/structure_validation.py"]
    end

    subgraph "Structural Gates"
        F -- "checks" --> G["Knot Screening"]
        F -- "checks" --> H["Curvature Gates"]
        F -- "checks" --> I["AMBER Relaxation"]
    end

    G --> J{"Accepted?"}
    H --> J
    I --> J
    J -- "Yes" --> K["*_validated.pdb"]
    J -- "No" --> L["Discard"]
    K -- "return count" --> B
```
Sources: [AlphaFlex/Step_3_sample_conformer.py:84-142](), [AlphaFlex/Step_3_sample_conformer.py:178-239](), [AlphaFlex/config.py:78-97]()

---

### IDR and IDP Workflows

The script differentiates between "Category 0" (fully disordered proteins) and proteins containing folded domains.

#### 1. IDR Workflow (`run_idr_workflow`)
This function handles regions that are grafted onto folded templates. It performs the following:
*   **Sidecar Management**: Copies `_truncation.json` files from Step 2 to the output directory to maintain residue numbering offsets [AlphaFlex/Step_3_sample_conformer.py:186-191]().
*   **Knot Filtering**: If the protein is listed in the `knot_screening.json` database, it filters the expected knot types based on the specific truncation window of the IDR [AlphaFlex/Step_3_sample_conformer.py:159-174]().
*   **Subprocess Execution**: Calls `sample_ldr.py` with specific flags for `--fold_curv_ratio`, `--fold_curv_window`, and `--junction_kappa` [AlphaFlex/Step_3_sample_conformer.py:94-110]().

#### 2. IDP Workflow (`run_idp_workflow`)
For proteins classified as fully disordered (Category 0), the script bypasses template-based sampling:
*   **Subprocess Execution**: Calls `sample_idp.py` directly [AlphaFlex/Step_3_sample_conformer.py:242-264]().
*   **Batching**: Uses `SAMPLE_BATCH_SIZE` to generate the requested `SAMPLE_N_CONFS` [AlphaFlex/config.py:81-82]().

Sources: [AlphaFlex/Step_3_sample_conformer.py:178-239](), [AlphaFlex/Step_3_sample_conformer.py:242-282]()

---

### Structural Validation & Knot Screening

A critical feature of Step 3 is the integration of `idpforge/utils/structure_validation.py` to ensure physical plausibility.

| Feature | Implementation | Description |
| :--- | :--- | :--- |
| **Knot Screening** | `load_knot_screening()` | Reads `knot_screening.json` to identify domains that must remain unknotted or have specific topologies [AlphaFlex/Step_3_sample_conformer.py:24](). |
| **Curvature Gates** | `--junction_kappa` | Measures the curvature ($\kappa$) at the junction between a folded domain and an IDR to prevent "kinked" transitions [AlphaFlex/Step_3_sample_conformer.py:102](). |
| **Fold Stability** | `--fold_curv_ratio` | Ensures that the folded domain template has not been significantly distorted during the diffusion or relaxation process [AlphaFlex/Step_3_sample_conformer.py:100](). |
| **Intermediate Cleanup** | `_cleanup_dir()` | Automatically deletes `*_raw.pdb` and `*_relaxed.pdb` files, retaining only the final `*_validated.pdb` to save disk space [AlphaFlex/Step_3_sample_conformer.py:59-67](). |

**Code-to-System Mapping: Validation Logic**
```mermaid
classDiagram
    class Step_3_sample_conformer {
        +run_idr_workflow()
        +generate_conformers()
        +_filter_knot_spec_for_window()
    }
    class sample_ldr {
        +main()
        +validate_structure_post_relax()
    }
    class structure_validation {
        +check_knots()
        +calc_junction_curvature()
        +check_fold_stability()
    }
    
    Step_3_sample_conformer --> sample_ldr : subprocess.run()
    sample_ldr --> structure_validation : calls
```
Sources: [AlphaFlex/Step_3_sample_conformer.py:94-124](), [AlphaFlex/Step_3_sample_conformer.py:159-174](), [AlphaFlex/config.py:86-89]()

---

### HPC Execution & Splitting

To handle large datasets (e.g., the Human Proteome), `Step_3_sample_conformer.py` includes logic to split the workload.

*   **Protein Selection**: The script reads the `idp_cases_to_run.json` generated in Step 2 [AlphaFlex/Step_3_sample_conformer.py:307-311]().
*   **Chunking**: Using the `--split_index` and `--total_splits` arguments, the script calculates a subset of proteins to process: `proteins[split_index::total_splits]` [AlphaFlex/Step_3_sample_conformer.py:321-324]().
*   **Logging**: Each protein's output is captured in a dedicated log file within the `logs/Step_3/` directory, facilitating debugging of specific failures [AlphaFlex/Step_3_sample_conformer.py:301-304]().

Sources: [AlphaFlex/Step_3_sample_conformer.py:285-330](), [AlphaFlex/config.py:56]()

---