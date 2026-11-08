# Dataset Creation Classes

<details>
<summary>Relevant source files</summary>

The following files were used as context for generating this wiki page:

- [dataset/2_create_database_dataset_files.py](dataset/2_create_database_dataset_files.py)
- [dataset/create_input_embeddings.py](dataset/create_input_embeddings.py)
- [dataset/get_ped_ensembles.py](dataset/get_ped_ensembles.py)
- [dataset/prepare_entry_from_pdb.py](dataset/prepare_entry_from_pdb.py)

</details>



This page provides API reference documentation for the core classes used in the Disobind dataset creation pipeline. These classes handle downloading PDB structures, creating binary complexes, redundancy reduction, and embedding generation, as well as specialized preparation for ensemble data.

**Sources:** [dataset/2_create_database_dataset_files.py](), [dataset/3_create_merged_binary_complexes.py](), [dataset/4_create_non_redundant_dataset.py](), [dataset/create_input_embeddings.py](), [dataset/get_ped_ensembles.py](), [dataset/prepare_entry_from_pdb.py]()

---

## Overview of Dataset Creation Classes

The dataset creation system consists of several primary classes that execute sequentially to transform raw database records into training-ready tensors.

```mermaid
graph LR
    subgraph "Stage 1: Data Collection"
        A["parse_pdbs_for_idrs"]
        AP["APedsTale"]
        EP["EntryFromPDB"]
    end
    
    subgraph "Stage 2: Complex Merging"
        B["Dataset"]
    end
    
    subgraph "Stage 3: Redundancy Reduction"
        C["NonRedundantDataset"]
    end
    
    subgraph "Stage 4: Embedding Generation"
        D["Embeddings"]
    end
    
    A -->|"Binary_complexes_None/"| B
    AP -->|"ped_test_input.csv"| D
    EP -->|"misc_test_input.csv"| D
    B -->|"merged_binary_complexes/"| C
    C -->|"train/test keys"| D
    D -->|"Train/Dev/Test .npy"| E["Training Pipeline"]
```

**Sources:** [dataset/2_create_database_dataset_files.py:44-145](), [dataset/3_create_merged_binary_complexes.py:39-100](), [dataset/4_create_non_redundant_dataset.py:31-143](), [dataset/create_input_embeddings.py:18-104](), [dataset/get_ped_ensembles.py:20-52](), [dataset/prepare_entry_from_pdb.py:20-68]()

---

## Class: parse_pdbs_for_idrs

The `parse_pdbs_for_idrs` class orchestrates initial data collection, downloads structures, and splits multi-chain PDBs into binary protein complexes.

**Location:** [dataset/2_create_database_dataset_files.py:44-1420]()

### Constructor Parameters
- `cores` (int): Number of CPU cores for parallel processing [dataset/2_create_database_dataset_files.py:53]().
- `create_dataset` (bool): If True, creates Disobind training dataset; if False, creates the StrIDR database [dataset/2_create_database_dataset_files.py:57]().

### Key Methods
- `forward()`: Main execution loop that coordinates directory creation, disorder database loading, and parallelized PDB/UniProt downloads [dataset/2_create_database_dataset_files.py:147-285]().
- `download_pdb_and_sift()`: Worker function that handles superseded PDB IDs, downloads PDB/CIF files, and retrieves SIFTS mappings [dataset/2_create_database_dataset_files.py:571-777]().
- `create_binary_complexes()`: Splits non-binary complexes into all possible pairwise combinations between entities [dataset/2_create_database_dataset_files.py:1364-1420]().

**Sources:** [dataset/2_create_database_dataset_files.py:44-1420]()

---

## Class: Dataset

The `Dataset` class merges overlapping binary complexes from the same UniProt ID pair to handle structural heterogeneity and different fragments.

**Location:** [dataset/3_create_merged_binary_complexes.py:39-1100]()

### Key Methods
- `forward()`: Orchestrates the filtering of invalid PDBs, grouping of overlapping sets, and final merging of contact maps [dataset/3_create_merged_binary_complexes.py:103-232]().
- `create_contact_map()`: Aggregates contact maps across all models (e.g., NMR ensembles) using an 8Å CA-CA distance threshold [dataset/3_create_merged_binary_complexes.py:294-372]().
- `merge_binary_complexes()`: Performs logical OR operations on contact maps and unions residue positions to create a single merged representation [dataset/3_create_merged_binary_complexes.py:795-869]()

**Sources:** [dataset/3_create_merged_binary_complexes.py:39-1100]()

---

## Class: NonRedundantDataset

The `NonRedundantDataset` class creates sequence non-redundant train and out-of-distribution (OOD) test sets using MMseqs2.

**Location:** [dataset/4_create_non_redundant_dataset.py:31-1035]()

### Redundancy Logic
- **OOD Set:** Created using a strict 20% sequence identity cutoff against the AlphaFold2 PDB70 training set [dataset/4_create_non_redundant_dataset.py:547-640]().
- **Train Set:** Created from remaining entries, optionally clustered at a 40% identity cutoff [dataset/4_create_non_redundant_dataset.py:114]().

### AF2/AF3 Integration
- `create_af2_input()`: Generates FASTA files and a path list for AlphaFold2-Multimer predictions [dataset/4_create_non_redundant_dataset.py:760-804]().
- `create_af3_input()`: Generates JSON batch files compatible with the AlphaFold3 server [dataset/4_create_non_redundant_dataset.py:806-858]().

**Sources:** [dataset/4_create_non_redundant_dataset.py:31-1035]()

---

## Class: Embeddings

The `Embeddings` class generates protein language model embeddings and creates final padded `.npy` files for training.

**Location:** [dataset/create_input_embeddings.py:18-712]()

### Configuration
- `scope`: "global" (full sequence) or "local" (fragment only) [dataset/create_input_embeddings.py:34]().
- `embedding_type`: Supports "T5", "ProstT5", "ProSE", "ESM2", and "BERT" [dataset/create_input_embeddings.py:36]().
- `max_len`: Pads all sequences and contact maps to a fixed size (default 100) [dataset/create_input_embeddings.py:42]().

### Data Flow
1. `initialize()`: Calls `get_embeddings` from `utility.py` to generate vectors [dataset/create_input_embeddings.py:239-273]().
2. `split_dataset()`: Partitions keys into Train, Dev, and Test sets [dataset/create_input_embeddings.py:153]().
3. `apply_padding()`: Zero-pads embeddings to `max_len` and creates a `target_mask` [dataset/create_input_embeddings.py:371-408]().
4. `create_input()`: Concatenates padded embeddings and contact maps into a single 3D tensor [dataset/create_input_embeddings.py:412-489]().

**Sources:** [dataset/create_input_embeddings.py:18-712]()

---

## Specialized Dataset Classes

### APedsTale
Downloads entries from the Protein Ensemble Database (PED) and processes ensembles of disordered proteins.
- **Purpose:** Creates specialized test sets for conformational ensembles [dataset/get_ped_ensembles.py:21-25]().
- **Key Logic:** Aggregates contact maps across all structural models in the PED entry [dataset/get_ped_ensembles.py:140-155]().

### EntryFromPDB
A utility class for preparing Disobind input/target files for a specific user-defined list of PDB IDs.
- **Purpose:** Ad-hoc dataset creation for specific complexes [dataset/prepare_entry_from_pdb.py:20-23]().
- **Key Logic:** Maps PDB chains to UniProt IDs via SIFTS and extracts residue-level coordinates [dataset/prepare_entry_from_pdb.py:157-175]().

**Sources:** [dataset/get_ped_ensembles.py:20-210](), [dataset/prepare_entry_from_pdb.py:20-200]()

---

## Data Entity Space Mapping

**Title: Code Entities to Data Outputs**

```mermaid
graph TD
    subgraph "Code Entity Space"
        P_CLASS["parse_pdbs_for_idrs"]
        D_CLASS["Dataset"]
        N_CLASS["NonRedundantDataset"]
        E_CLASS["Embeddings"]
    end

    subgraph "System Output Space"
        B_JSON["Binary_complexes/*.json"]
        M_H5["merged_binary_complexes/*.h5"]
        K_TXT["train_binary_complexes.txt"]
        N_DATA["Train_set_global.npy"]
    end

    P_CLASS --"generates"--> B_JSON
    D_CLASS --"generates"--> M_H5
    N_CLASS --"generates"--> K_TXT
    E_CLASS --"generates"--> N_DATA
```

**Sources:** [dataset/2_create_database_dataset_files.py:94-96](), [dataset/3_create_merged_binary_complexes.py:68-72](), [dataset/4_create_non_redundant_dataset.py:97-101](), [dataset/create_input_embeddings.py:95-97]()

**Title: Class Interactions with Utility Functions**

```mermaid
graph TD
    subgraph "Dataset Creation Classes"
        Embeddings["Embeddings"]
        NonRedundantDataset["NonRedundantDataset"]
        Dataset["Dataset"]
    end

    subgraph "utility.py Functions"
        get_emb["get_embeddings"]
        mm_cluster["mmseqs_cluster"]
        get_cmap["get_contact_map"]
    end

    Embeddings --"calls"--> get_emb
    NonRedundantDataset --"calls"--> mm_cluster
    Dataset --"calls"--> get_cmap
```

**Sources:** [dataset/create_input_embeddings.py:14](), [dataset/4_create_non_redundant_dataset.py:27](), [dataset/3_create_merged_binary_complexes.py:35]()

---