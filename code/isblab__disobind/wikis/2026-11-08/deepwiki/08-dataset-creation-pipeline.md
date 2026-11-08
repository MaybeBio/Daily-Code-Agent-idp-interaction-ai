# Dataset Creation Pipeline

<details>
<summary>Relevant source files</summary>

The following files were used as context for generating this wiki page:

- [analysis/README.md](analysis/README.md)
- [dataset/1_disobind_databases.py](dataset/1_disobind_databases.py)
- [dataset/README.md](dataset/README.md)

</details>



## Purpose and Scope

This document provides a high-level overview of the Disobind dataset creation pipeline. The system transforms raw structural data from the Protein Data Bank (PDB) and various disorder databases (DisProt, IDEAL, MobiDB, DIBS, MFIB, FuzDB) into curated, non-redundant training and evaluation datasets. 

The pipeline handles the complexity of mapping PDB chains to UniProt sequences, identifying disordered regions, creating binary interaction complexes, and ensuring structural novelty in test sets.

**Sources:** [dataset/README.md:1-27](), [analysis/README.md:1-68]()

---

## Pipeline Stages

The pipeline is organized into five sequential stages, each controlled by a numbered script in the `dataset/` directory.

| Stage | Script | Description | Details |
| :--- | :--- | :--- | :--- |
| **1** | `1_disobind_database.py` | **Data Collection**: Aggregates PDB IDs from specialized disorder databases. | [Data Collection from Databases](#3.1) |
| **2** | `2_create_database_dataset_files.py` | **Complex Creation**: Downloads structures and fragments them into binary complexes. | [Creating Binary Complexes](#3.2) |
| **3** | `3_create_merged_binary_complexes.py` | **Merging**: Combines overlapping fragments into unified interaction windows. | [Creating Binary Complexes](#3.2) |
| **4** | `4_create_non_redundant_dataset.py` | **Redundancy Reduction**: Performs sequence clustering and OOD splitting. | [Non-Redundant Dataset Creation](#3.3) |
| **5** | `create_input_embeddings.py` | **Embeddings**: Generates ProtT5-XL-U50 features for the final sequences. | [Generating Embeddings](#3.4) |

**Sources:** [dataset/README.md:1-27](), [dataset/2_create_database_dataset_files.py:44-147]()

---

## System Architecture: Code to Entity Mapping

The following diagrams illustrate how high-level pipeline stages map to specific classes and modules within the codebase.

### Data Acquisition and Processing Flow
This diagram maps the data flow from external APIs to the primary processing classes.

```mermaid
graph TB
    subgraph "External Resources"
        PDB[("PDB/RCSB")]
        UniProt[("UniProt API")]
        SIFTS[("SIFTS Mapping")]
    end

    subgraph "dataset/from_APIs_with_love.py"
        API_Helper["send_request()"]
        PDB_DL["download_pdb()"]
        Uni_DL["get_uniprot_seq()"]
        SIFTS_DL["download_SIFTS_Uni_PDB_mapping()"]
    end

    subgraph "Core Processing Classes"
        Illuminati["TheIlluminati<br/>(1_disobind_database.py)"]
        Parser["parse_pdbs_for_idrs<br/>(2_create_database_dataset_files.py)"]
        DS["Dataset<br/>(3_create_merged_binary_complexes.py)"]
    end

    PDB --> PDB_DL
    UniProt --> Uni_DL
    SIFTS --> SIFTS_DL
    
    Illuminati -- "Collects IDs" --> Parser
    Parser -- "where_the_magic_happens()" --> DS
    DS -- "merge_binary_complexes()" --> Output[("Merged HDF5s")]
```
**Sources:** [dataset/1_disobind_databases.py:31-36](), [dataset/2_create_database_dataset_files.py:44-50](), [dataset/3_create_merged_binary_complexes.py:39-45](), [dataset/from_APIs_with_love.py:19-25]()

### Dataset Partitioning and Embedding Generation
This diagram bridges the transition from raw structural fragments to model-ready HDF5 files.

```mermaid
graph LR
    subgraph "Redundancy Reduction"
        NRD["NonRedundantDataset<br/>(4_create_non_redundant_dataset.py)"]
        MMSeqs["mmseqs_cluster()"]
        PDB70_Check["get_difference_PDB_set()"]
    end

    subgraph "Feature Engineering"
        Embed["Embeddings<br/>(create_input_embeddings.py)"]
        T5["ProtT5-XL-U50"]
    end

    NRD --> MMSeqs
    NRD --> PDB70_Check
    PDB70_Check -- "OOD Split" --> Embed
    Embed --> T5
    T5 -- "Final Data" --> H5[("HDF5 Train/Test Sets")]
```
**Sources:** [dataset/4_create_non_redundant_dataset.py:31-40](), [dataset/create_input_embeddings.py:1-20]()

---

## Stage Details and Child Pages

### 3.1 Data Collection from Databases
The initial stage uses the `TheIlluminati` class in `1_disobind_database.py` to parse flat files and query APIs for PDB entries containing disordered proteins. It considers databases like DIBS, MFIB, and FuzDB.
For details, see [Data Collection from Databases](#3.1).

**Sources:** [dataset/1_disobind_databases.py:31-56](), [dataset/README.md:4-7]()

### 3.2 Creating and Merging Binary Complexes
This stage involves two major steps:
1. **Fragmentation**: `parse_pdbs_for_idrs` downloads structures and identifies disordered segments, creating initial binary pairs.
2. **Merging**: `Dataset.module3()` merges these pairs into larger interaction windows (up to 100-200 residues) to capture the full context of the binding interface.
For details, see [Creating Binary Complexes](#3.2).

**Sources:** [dataset/2_create_database_dataset_files.py:1003-1201](), [dataset/3_create_merged_binary_complexes.py:736-869]()

### 3.3 Non-Redundant Dataset Creation
To ensure the model generalizes, the `NonRedundantDataset` class removes sequence redundancy using MMSeqs2. It specifically creates an **Out-Of-Distribution (OOD)** test set by excluding any PDB IDs present in the AlphaFold2 PDB70 training set.
For details, see [Non-Redundant Dataset Creation](#3.3).

**Sources:** [dataset/4_create_non_redundant_dataset.py:380-402](), [dataset/4_create_non_redundant_dataset.py:437-473]()

### 3.4 Generating Embeddings
The final stage uses the `Embeddings` class to convert protein sequences into high-dimensional vectors using a pre-trained ProtT5 model. These are stored in HDF5 format for efficient loading during training.
For details, see [Generating Embeddings](#3.4).

**Sources:** [dataset/create_input_embeddings.py:1-27]()

### 3.5 Special Dataset Preparation
Disobind includes utilities for specialized benchmarks, such as the IDPPI dataset and NMR ensemble processing (PEDS). These require specific aggregation logic for contact maps across multiple structural models.
For details, see [Special Dataset Preparation](#3.5).

**Sources:** [dataset/README.md:42-56](), [analysis/README.md:24-34]()

### 3.6 Data Processing Utilities
The pipeline relies on a robust set of utilities in `utility.py` and `from_APIs_with_love.py` for tasks like coordinate extraction (`get_coordinates`), contact map calculation (`get_contact_map`), and handling PDB obsolescence.
For details, see [Data Processing Utilities](#3.6).

**Sources:** [dataset/from_APIs_with_love.py:1-952](), [dataset/3_create_merged_binary_complexes.py:29-32]()

---