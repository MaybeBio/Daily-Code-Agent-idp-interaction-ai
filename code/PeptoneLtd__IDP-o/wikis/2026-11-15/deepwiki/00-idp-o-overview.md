# IDP-o Overview

<details>
<summary>Relevant source files</summary>

The following files were used as context for generating this wiki page:

- [LICENSE](LICENSE)
- [README.md](README.md)
- [assets/idp-o.png](assets/idp-o.png)

</details>



IDP-o is a high-level computational framework designed to generate structural ensembles of **Intrinsically Disordered Proteins (IDPs)** using a fragment-based assembly approach. Unlike traditional molecular dynamics or de novo folding, IDP-o leverages the vast structural information in the AlphaFold Protein Structure Database (AFDB) by decomposing target sequences into overlapping fragments, retrieving matching structural hits, and stitching them into physically plausible ensembles.

The system is optimized for performance, utilizing GPU acceleration via **JAX** and **CuPy**, and employs the **FoldComp** format for efficient storage and retrieval of millions of protein structures.

### The IDP-o Workflow
The generation process follows a linear pipeline from raw sequence to a structural ensemble:

1.  **Sequence Fragmentation**: The input sequence is cut into 6-residue fragments with a 2-residue overlap [README.md:8-8]().
2.  **Database Search**: A GPU-accelerated search identifies matching fragments within the `afdb_uniprot_v4` database [README.md:23-39]().
3.  **Structure Extraction**: Coordinates for the identified hits are reconstructed from FoldComp binaries.
4.  **Ensemble Assembly**: Fragments are aligned and merged using a hierarchical strategy, filtered for steric clashes, and sorted by RMSD.

For a step-by-step guide on environment setup and running your first ensemble, see [Getting Started](#1.1).

### System Architecture and Data Flow

The following diagram illustrates how the natural language concepts of the pipeline map to specific code entities and scripts within the repository.

**IDP-o Pipeline Mapping**
```mermaid
graph TD
    subgraph "Input Space"
        SEQ["Protein Sequence"]
        DB["FoldComp Database (afdb_uniprot_v4)"]
    end

    subgraph "Code Entity Space (Scripts)"
        PREP["prepare_foldcomp_fasta.py"]
        SEARCH["fasta_search_in_foldcomp_database.py"]
        EXTRACT["extract_structures_from_foldcomp_database.py"]
        JOIN["join_fragments.py"]
        CLI["build_ensemble.py"]
    end

    subgraph "Output Space"
        IDX["Offset FASTA (.fasta)"]
        HITS["Byte Starts (.pkl)"]
        COORDS["Coordinates (.h5)"]
        ENS["Final Ensemble (.h5/.xtc)"]
    end

    DB --> PREP
    PREP --> IDX
    SEQ --> CLI
    CLI --> SEARCH
    IDX --> SEARCH
    SEARCH --> HITS
    HITS --> EXTRACT
    EXTRACT --> COORDS
    COORDS --> JOIN
    JOIN --> ENS

    style CLI stroke-width:4px
```
**Sources:** [README.md:36-52](), [README.md:8-10]()

### Core Components

IDP-o is organized into specialized scripts that handle distinct phases of the ensemble generation.

| Component | Code Entity | Role |
| :--- | :--- | :--- |
| **Data Preprocessor** | `prepare_foldcomp_fasta.py` | Generates the byte-offset FASTA index required for rapid database lookups [README.md:36-36](). |
| **Search Engine** | `fasta_search_in_foldcomp_database.py` | Performs vectorized 6-mer matching against the database using CuPy. |
| **Structure Reconstructor** | `extract_structures_from_foldcomp_database.py` | Uses JAX and NeRF to convert FoldComp binaries into 3D coordinates. |
| **Assembly Engine** | `join_fragments.py` | Merges fragments using SVD alignment and steric filtering. |
| **Orchestrator** | `build_ensemble.py` | The main CLI that coordinates all stages for a single sequence [README.md:44-52](). |

For a deep dive into the scientific theory, such as NeRF reconstruction and fragment overlaps, see [Key Concepts and Terminology](#1.2).

### Technical Integration

The system is designed to run within a containerized environment to manage complex dependencies like CUDA, JAX, and FoldComp.

**Component Interaction and Dependencies**
```mermaid
graph LR
    subgraph "External Dependencies"
        JC["JAX / CuPy"]
        FC["FoldComp"]
        MT["MDTraj"]
        NX["nerfax"]
    end

    subgraph "IDP-o Core"
        direction TB
        BE["build_ensemble.py"]
        GD["generate_dataset.py"]
    end

    GD -- "subprocesses" --> BE
    BE -- "Linear Algebra" --> JC
    BE -- "IO" --> FC
    BE -- "Topology" --> MT
    BE -- "Hydrogen Placement" --> NX

```
**Sources:** [README.md:55-55](), [LICENSE:1-76]()

### Licensing and Citation
IDP-o is released under the **Apache License 2.0** [README.md:72-75](). If you use this software in your research, please cite the work as described in the [README.md:58-67]().

---
**Next Steps:**
- To install and run the software, visit [Getting Started](#1.1).
- To understand the underlying fragment-based methodology, visit [Key Concepts and Terminology](#1.2).

---