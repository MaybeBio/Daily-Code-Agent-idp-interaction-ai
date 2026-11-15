# Entrypoints and Batch Processing

<details>
<summary>Relevant source files</summary>

The following files were used as context for generating this wiki page:

- [generate_dataset.py](generate_dataset.py)
- [scripts/build_ensemble.py](scripts/build_ensemble.py)

</details>



This page provides an overview of the two primary ways to interact with the IDP-o pipeline: the single-sequence CLI and the batch dataset generator. While the core logic of ensemble generation resides in specialized modules, these entrypoints orchestrate the flow of data, manage scratch space, and handle output serialization.

## Overview of Entrypoints

IDP-o provides two main scripts for users to generate structural ensembles:

1.  **`scripts/build_ensemble.py`**: The primary CLI for processing a single protein sequence. It orchestrates the entire four-stage pipeline (search, extraction, and assembly) for one target.
2.  **`generate_dataset.py`**: A high-level batch processor that consumes CSV or FASTA files containing multiple sequences. It functions as a wrapper that delegates work to `build_ensemble.py` via subprocesses.

### System Orchestration Diagram
The following diagram illustrates how the batch generator interfaces with the single-sequence builder and the underlying pipeline modules.

```mermaid
graph TD
    subgraph "Batch Layer (generate_dataset.py)"
        Input["Input CSV/FASTA"] --> Parser["pd.read_csv / FASTA Parser"]
        Parser --> Loop["Sequence Loop"]
        Loop --> Subproc["subprocess.check_call"]
    end

    subgraph "Execution Layer (build_ensemble.py)"
        Subproc --> BE_Main["build_ensemble.py:main"]
        BE_Main --> Search["fasta_search_in_foldcomp_database.py"]
        BE_Main --> Extract["extract_structures_from_foldcomp_database.py"]
        BE_Main --> Join["join_fragments.py"]
    end

    subgraph "Storage"
        BE_Main --> Scratch["/tmp/tmp-{name}/"]
        Join --> Output["Final Ensemble (.h5, .xtc, .pdb)"]
    end
```
**Sources:** [scripts/build_ensemble.py:25-38](), [scripts/build_ensemble.py:61-80](), [generate_dataset.py:97-118]()

---

## Single-Sequence Ensemble Builder: `build_ensemble.py`

The `build_ensemble.py` script is the central orchestrator for a single sequence. It is responsible for validating the input amino acid sequence [scripts/build_ensemble.py:47-50](), setting up the temporary directory structure for intermediate files [scripts/build_ensemble.py:52-58](), and calling the three main pipeline stages in order.

For a detailed reference on CLI arguments, output formats, and internal logic, see **[Single-Sequence Ensemble Builder: build_ensemble.py](#3.1)**.

### Key Responsibilities
*   **Pipeline Orchestration**: Sequentially executes `fasta_search`, `extract_structures`, and `join_fragments` [scripts/build_ensemble.py:61-80]().
*   **Format Support**: Manages the final assembly into various formats including HDF5 (`.h5`), XTC, DCD, and PDB [scripts/build_ensemble.py:117-119]().
*   **Intermediate Management**: Handles the creation and cleanup of `byte_starts.pkl` and fragment ensemble directories within a designated `--scratch_folder` [scripts/build_ensemble.py:52-57]().

**Sources:** [scripts/build_ensemble.py:15-81]()

---

## Batch Dataset Generator: `generate_dataset.py`

The `generate_dataset.py` script is designed for large-scale dataset production, such as the creation of the **IDRome-o** dataset. It provides robust handling for input files, duplicate removal, and error logging.

For details on input parsing, parallelization strategies, and error sidecar files, see **[Batch Dataset Generator: generate_dataset.py](#3.2)**.

### Execution Flow
The script maps sequence data from a source file to a series of CLI calls. It isolates each sequence by creating unique scratch folders under `/tmp/tmp-{name}` to prevent data collisions during concurrent execution [generate_dataset.py:101-114]().

### Data Flow Mapping
The diagram below shows how `generate_dataset.py` maps natural language/CSV concepts to code entities used by `build_ensemble.py`.

```mermaid
graph LR
    subgraph "Input Space (CSV/FASTA)"
        CSV_Row["CSV Row / FASTA Record"]
        SeqCol["'fasta' column"]
        NameCol["'seq_name' column"]
    end

    subgraph "Code Entity Space (build_ensemble.py arguments)"
        SeqCol -- maps to --> ArgSeq["--sequence"]
        NameCol -- maps to --> ArgOut["--outpath (via {name}.{format})"]
        NameCol -- maps to --> ArgScratch["--scratch_folder (via tmp-{name})"]
    end

    subgraph "Error Handling"
        CheckCall["check_call()"] -- on Exception --> Sidecar["{name}.txt (Error Log)"]
    end
```
**Sources:** [generate_dataset.py:67-81](), [generate_dataset.py:100-114](), [generate_dataset.py:120-124]()

### Key Features
*   **Input Flexibility**: Supports both FASTA headers and CSV columns defined via `--column_names` [generate_dataset.py:45-49]().
*   **Parallelization Strategy**: Includes a `--shuffle` flag to randomize the processing order, allowing multiple instances of the script to run on the same input file with a lower probability of processing the same sequence simultaneously [generate_dataset.py:50-51](), [generate_dataset.py:94-96]().
*   **Fault Tolerance**: If a sequence fails to build, the script captures the `traceback` and writes it to a `.txt` sidecar file next to the intended output path, then continues to the next sequence [generate_dataset.py:120-124]().

**Sources:** [generate_dataset.py:24-59](), [generate_dataset.py:84-88]()

---