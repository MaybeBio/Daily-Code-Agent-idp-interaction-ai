# Creating Binary Complexes

<details>
<summary>Relevant source files</summary>

The following files were used as context for generating this wiki page:

- [dataset/2_create_database_dataset_files.py](dataset/2_create_database_dataset_files.py)
- [dataset/3_create_merged_binary_complexes.py](dataset/3_create_merged_binary_complexes.py)
- [dataset/from_APIs_with_love.py](dataset/from_APIs_with_love.py)
- [dataset/parse_sifts.py](dataset/parse_sifts.py)

</details>



## Purpose and Scope

This page describes the process of creating binary protein complexes from PDB structures for the Disobind dataset. A binary complex consists of two protein chains: one containing intrinsically disordered regions (IDRs) and one partner protein (which may be ordered or disordered). This step follows [Data Collection from Databases](#3.1) and precedes [Non-Redundant Dataset Creation](#3.3).

The binary complex creation involves:
1. Identifying IDRs and ordered regions in protein chains from PDB structures.
2. Classifying chains based on disorder content.
3. Creating all possible binary combinations between IDR-containing chains and partner chains.
4. Validating structural coordinates and contact maps for the resulting pairs.
5. Merging overlapping fragments to maximize sequence coverage.

---

## Overview

Binary complex creation is performed by the `parse_pdbs_for_idrs` class in `dataset/2_create_database_dataset_files.py` and refined by the `Dataset` class in `dataset/3_create_merged_binary_complexes.py`.

| Script | Class | Key Stage | Purpose |
|:-------|:------|:----------|:--------|
| `2_create_database_dataset_files.py` | `parse_pdbs_for_idrs` | `parallelize_dataset_creation()` | Identify disordered and ordered chains in each PDB structure. |
| `2_create_database_dataset_files.py` | `parse_pdbs_for_idrs` | `create_binary_complexes()` | Generate all-vs-all binary combinations from identified chains. |
| `3_create_merged_binary_complexes.py` | `Dataset` | `module1()` | Validate coordinates and ensure contact existence between chains. |
| `3_create_merged_binary_complexes.py` | `Dataset` | `module3()` | Merge overlapping fragments for the same UniProt pair. |

**Sources:** [dataset/2_create_database_dataset_files.py:44-45](), [dataset/2_create_database_dataset_files.py:287-381](), [dataset/3_create_merged_binary_complexes.py:39-40](), [dataset/3_create_merged_binary_complexes.py:103-210]()

---

## IDR Identification and Fragment Classification

### Fragment Processing Workflow

The `where_the_magic_happens()` method processes each chain in a PDB structure to identify IDR-containing fragments.

**Diagram: Fragment Processing Workflow**

```mermaid
flowchart TD
    Start["PDB Chain Data<br/>(entity_id, asym_id, auth_asym_id, uniprot_ids)"] --> GetMapping["Get SIFTS Mapping<br/>get_sifts_mapping()"]
    
    GetMapping --> CheckMapping{"Mapping<br/>Exists?"}
    CheckMapping -->|No| LogNotMapped["Log: chain_not_mapped"]
    CheckMapping -->|Yes| RemoveMissing["Remove Missing Residues<br/>remove_missing_residues()"]
    
    RemoveMissing --> CheckLen{"Fragments<br/>>= min_len?"}
    CheckLen -->|No| Skip1["Skip Chain"]
    CheckLen -->|Yes| FragmentChain{"Fragments<br/>> max_len?"}
    
    FragmentChain -->|Yes| SplitFrags["Fragment Chain<br/>fragment_chain()"]
    FragmentChain -->|No| GetDisorder["Get Disorder Regions<br/>find_disorder_regions()"]
    SplitFrags --> GetDisorder
    
    GetDisorder --> CheckDisorder{"Disorder<br/>Regions Found?"}
    CheckDisorder -->|No| LogNoDisorder["Log: no_disordered_regions"]
    CheckDisorder -->|Yes| CalcOverlap["Calculate Overlap<br/>get_overlap()"]
    
    CalcOverlap --> CheckInSeq{"disorder_in_seq<br/>!= []?"}
    CheckInSeq -->|No| LogNoDisorder
    CheckInSeq -->|Yes| CheckInStruct{"disorder_in_struct<br/>!= []?"}
    
    CheckInStruct -->|No| Skip2["Skip Chain"]
    CheckInStruct -->|Yes| ClassifyFrag["Classify Fragment<br/>IDR >= min_disorder_percent<br/>→ prot1<br/>else → prot2"]
    
    ClassifyFrag --> StoreData["Store in idr_pdbs_selected<br/>{PDB: {UniID: {prot1/prot2: data}}}"]
    
    LogNotMapped --> End([End])
    Skip1 --> End
    LogNoDisorder --> End
    Skip2 --> End
    StoreData --> End
```

**Sources:** [dataset/2_create_database_dataset_files.py:1003-1249]()

### Key Processing Steps

#### 1. SIFTS Mapping and Missing Residues
The pipeline uses `get_sifts_mapping` from `from_APIs_with_love.py` to align PDB residues to UniProt accessions. The `remove_missing_residues()` method removes residues that lack structural coordinates. Missing residues split a chain into separate fragments.

**Sources:** [dataset/2_create_database_dataset_files.py:470-511](), [dataset/from_APIs_with_love.py:287-320]()

#### 2. Fragment Length Control
Fragments are filtered and split based on `min_len` (default 20) and `max_len` (default 200). `fragment_chain()` recursively calls `split_frag()` to ensure no fragment exceeds the length limit.

**Sources:** [dataset/2_create_database_dataset_files.py:514-554]()

#### 3. Disorder Classification
`find_disorder_regions()` retrieves annotations from DisProt, IDEAL, and MobiDB. A fragment is classified as **prot1** (the disordered query) if its `disorder_percentage` (disordered residues in structure / total fragment length) is $\ge$ `min_disorder_percent` (default 0.2). Otherwise, it is **prot2** (the partner).

**Sources:** [dataset/2_create_database_dataset_files.py:1102-1169](), [dataset/utility.py:911-967]()

---

## Binary Complex Generation and Validation

### All-vs-All Combination
`create_binary_complexes()` generates pairs for every `prot1` and `prot2` fragment found in a PDB. If `uni_id1 == uni_id2`, the system checks for positional overlap using `check_for_overlap()` to prevent self-interaction artifacts.

**Sources:** [dataset/2_create_database_dataset_files.py:1254-1450](), [dataset/utility.py:821-837]()

### Structural Validation (Module 1)
In `3_create_merged_binary_complexes.py`, the `module1()` method performs strict structural checks:
1. **Coordinate Retrieval**: Uses `get_coordinates()` to extract $C\alpha$ positions.
2. **Contact Map Calculation**: Uses `get_contact_map()` with a default `contact_threshold` of 8Å.
3. **Contact Validation**: Complexes with zero contacts between the two chains are discarded.

**Sources:** [dataset/3_create_merged_binary_complexes.py:133-153](), [dataset/3_create_merged_binary_complexes.py:280-450](), [dataset/utility.py:100-147]()

---

## Merging Overlapping Complexes

To reduce redundancy and increase the context for the model, overlapping fragments of the same UniProt pair are merged.

**Diagram: Merging and Validation Pipeline**

```mermaid
flowchart TD
    Start["Valid Binary Complexes<br/>(from Module 1)"] --> Sort["Sort by Residue Positions<br/>sort_by_residue_positions()"]
    Sort --> Module2["Module 2: Grouping<br/>Group overlapping sets by UniProt Pair"]
    Module2 --> Module3["Module 3: Merging<br/>module3()"]
    
    subgraph MergingLogic["Merging Logic (Dataset.merge_complexes)"]
        FindOverlap["Check for Overlap<br/>check_for_overlap()"]
        MergePos["Merge Positions<br/>merge_residue_positions()"]
        CheckMax["Check Max Length<br/>merged_seq_exceeds_maxlen()"]
    end
    
    Module3 --> MergingLogic
    MergingLogic --> FinalValid["Final Merged Binary Complexes"]
```

**Sources:** [dataset/3_create_merged_binary_complexes.py:155-210](), [dataset/3_create_merged_binary_complexes.py:644-800](), [dataset/utility.py:840-863]()

### Merging Rules
1. **Overlap Requirement**: Two fragments are merged if their UniProt position ranges overlap.
2. **Length Constraint**: The resulting merged sequence must not exceed `max_len` (100 for this stage).
3. **Coordinate Handling**: If multiple PDB structures cover the same region, the system retains the coordinates from the first encountered structure in the sorted list.

**Sources:** [dataset/3_create_merged_binary_complexes.py:80-82](), [dataset/3_create_merged_binary_complexes.py:720-780]()

---

## Data Structures and Storage

### Directory Structure
```
database/
├── Disobind_dataset_200_20_0.2/
│   ├── IDR_PDBs/                # Classified chains (.json)
│   └── Binary_complexes_None/    # Initial all-vs-all pairs (.json)
└── v_{version}/
    ├── Valid_Binary_Complexes/   # Structurally validated pairs
    └── merged_binary_complexes/  # Final merged pairs (.json)
```

**Sources:** [dataset/2_create_database_dataset_files.py:83-105](), [dataset/3_create_merged_binary_complexes.py:46-72]()

### Binary Complex Schema
The final complexes are stored as JSON dictionaries indexed by a key created via `create_key()` (format: `uni_id1--uni_id2`).

| Field | Description |
|:---|:---|
| `PDB ID` | Source PDB accession |
| `Asym ID1 / 2` | Internal PDB chain identifiers |
| `Uniprot positions1 / 2` | Residue indices in UniProt sequence |
| `PDB positions1 / 2` | Residue indices in PDB file |

**Sources:** [dataset/2_create_database_dataset_files.py:557-568](), [dataset/3_create_merged_binary_complexes.py:74-78]()

---