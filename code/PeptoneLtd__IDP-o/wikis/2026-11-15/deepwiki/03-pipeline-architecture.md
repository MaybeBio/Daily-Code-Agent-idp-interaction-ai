# Pipeline Architecture

<details>
<summary>Relevant source files</summary>

The following files were used as context for generating this wiki page:

- [README.md](README.md)
- [scripts/build_ensemble.py](scripts/build_ensemble.py)

</details>



The IDP-o pipeline is a four-stage process designed to generate structural ensembles for Intrinsically Disordered Proteins (IDPs) by leveraging fragment-based assembly from the AlphaFold Database (AFDB). The pipeline transitions from raw sequence data to a fully reconstructed 3D ensemble by decomposing sequences into overlapping 6-mer fragments, searching a massive structural database, and stochastically assembling the results.

### High-Level Workflow
The orchestration of the pipeline is handled by `scripts/build_ensemble.py`, which sequentially invokes the three primary computational modules after ensuring the database is prepared [scripts/build_ensemble.py:60-80]().

#### Pipeline Overview
```mermaid
graph TD
    subgraph "0. Setup"
        A["prepare_foldcomp_fasta.py"] --> B[".fasta with Byte Offsets"]
    end

    subgraph "1. Fragment Search"
        C["fasta_search_in_foldcomp_database.py"] --> D["byte_starts.pkl"]
    end

    subgraph "2. Structure Extraction"
        E["extract_structures_from_foldcomp_database.py"] --> F["Fragment HDF5 Ensembles"]
    end

    subgraph "3. Fragment Assembly"
        G["join_fragments.py"] --> H["Final Ensemble (.h5, .pdb, etc)"]
    end

    B --> C
    D --> E
    F --> G
```
**Sources:** [scripts/build_ensemble.py:25-81](), [README.md:8-10]()

---

### 1. Database Preparation
Before running searches, the system requires a specialized FASTA index of the FoldComp database. This index includes the byte-offset for each sequence within the binary FoldComp database file, allowing for rapid random access during extraction.

*   **Primary Script:** `prepare_foldcomp_fasta.py`
*   **Key Artifacts:** `afdb_uniprot_v4.fasta` (containing `>byte_offset` headers).
*   **Role:** Downloads `afdb_uniprot_v4` and synchronizes the `.index` and `.lookup` files to map sequences to their binary locations [README.md:23-39]().

For details, see [Database Preparation: prepare_foldcomp_fasta.py](#2.1).

---

### 2. Fragment Search
The sequence of interest is cut into overlapping 6-mer fragments. A GPU-accelerated search scans the indexed FASTA file to find all occurrences of these fragments across the entire AFDB.

*   **Primary Script:** `fasta_search_in_foldcomp_database.py`
*   **Key Artifacts:** `byte_starts.pkl` [scripts/build_ensemble.py:55]().
*   **Role:** Uses CuPy-based vectorized matching to identify sequence hits and records the starting byte position of the parent protein in the FoldComp database [scripts/build_ensemble.py:61]().

For details, see [Fragment Search: fasta_search_in_foldcomp_database.py](#2.2).

#### Data Flow: Search to Extraction
```mermaid
graph LR
    subgraph "fasta_search_in_foldcomp_database.py"
        direction TB
        SEQ["Target Sequence"] -- "6-mer sliding window" --> FRAGS["Fragments"]
        DB["Offset FASTA"] -- "CuPy Stream Match" --> MATCH["Match Indices"]
    end

    MATCH -- "List[Dict]" --> PKL["byte_starts.pkl"]

    subgraph "extract_structures_from_foldcomp_database.py"
        PKL -- "Random Access" --> FCMP["FoldComp Binary"]
    end
```
**Sources:** [scripts/build_ensemble.py:60-70](), [README.md:24-32]()

---

### 3. Structure Extraction
Once the locations of matching sequences are known, the pipeline extracts the specific 3D coordinates for those fragments.

*   **Primary Script:** `extract_structures_from_foldcomp_database.py`
*   **Key Artifacts:** `fragment_ensembles/` directory containing HDF5 files for each fragment [scripts/build_ensemble.py:52-54]().
*   **Role:** Parses FoldComp binaries, de-quantizes torsion angles, and reconstructs backbone coordinates using JAX-accelerated NeRF (Natural Extension Reference Frame) [scripts/build_ensemble.py:62-70](). It can optionally filter out `cis-omega` angles [scripts/build_ensemble.py:131-135]().

For details, see [Structure Extraction: extract_structures_from_foldcomp_database.py](#2.3).

---

### 4. Fragment Assembly
The final stage stochastically joins the extracted fragments into full-length protein configurations.

*   **Primary Script:** `join_fragments.py`
*   **Key Artifacts:** Final ensemble file (e.g., `ensemble.h5`, `ensemble.pdb`) [scripts/build_ensemble.py:116-119]().
*   **Role:** Uses a hierarchical merging strategy to align overlapping residues via SVD-based affine transformations. It filters for steric clashes and can sort the resulting ensemble by RMSD for better visualization [scripts/build_ensemble.py:72-80](), [scripts/build_ensemble.py:137-141]().

For details, see [Fragment Assembly: join_fragments.py](#2.4).

---

### Filesystem and Data Flow Summary

The pipeline relies on a `scratch_folder` to pass data between stages.

| Stage | Input | Output | Code Entity |
| :--- | :--- | :--- | :--- |
| **Preparation** | AFDB FoldComp | `afdb_uniprot_v4.fasta` | `prepare_foldcomp_fasta.py` |
| **Search** | Sequence + FASTA | `byte_starts.pkl` | `fasta_search_in_foldcomp_database.py` |
| **Extraction** | `byte_starts.pkl` + FoldComp DB | `fragment_ensembles/*.h5` | `extract_structures_from_foldcomp_database.py` |
| **Assembly** | `fragment_ensembles/` | `ensemble.h5` | `join_fragments.py` |

**Sources:** [scripts/build_ensemble.py:25-80](), [README.md:43-52]()

---