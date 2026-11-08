# Data Collection from Databases

<details>
<summary>Relevant source files</summary>

The following files were used as context for generating this wiki page:

- [database/input_files/DIBS_MFIB.csv](database/input_files/DIBS_MFIB.csv)
- [database/input_files/DisProt.csv](database/input_files/DisProt.csv)
- [database/input_files/IDEAL.csv](database/input_files/IDEAL.csv)
- [database/input_files/Merged_PDB_IDs.txt](database/input_files/Merged_PDB_IDs.txt)
- [database/input_files/MobiDB.csv](database/input_files/MobiDB.csv)
- [database/raw/DIBS_complete_17Apr24.txt](database/raw/DIBS_complete_17Apr24.txt)
- [dataset/1_disobind_databases.py](dataset/1_disobind_databases.py)

</details>



## Purpose and Scope

This page describes the initial data collection phase of the Disobind dataset creation pipeline. It covers how to collect PDB IDs from multiple disorder protein databases, download PDB structures, obtain PDB-to-UniProt mappings via SIFTS, and download UniProt sequences. This is the first step in creating training data for Disobind models.

The primary entry point for this stage is the script `1_disobind_database.py`, which utilizes the `TheIlluminati` class to parse raw database files and the `parse_pdbs_for_idrs` class to orchestrate the retrieval of structural data.

**Sources:** [dataset/1_disobind_databases.py:1-36](), [dataset/2_create_database_dataset_files.py:1-150]()

---

## Overview of Disorder Databases

Disobind collects PDB IDs and disorder annotations from several specialized databases. These databases provide the ground truth for intrinsically disordered regions (IDRs) and their binding interfaces.

| Database | Full Name | Code Reference / Implementation |
|----------|-----------|---------------------------------|
| **DisProt** | Database of Protein Disorder | Parsed from `DisProt.csv` [dataset/1_disobind_databases.py:471-475]() |
| **IDEAL** | Intrinsically Disordered proteins with Extensive Annotations and Literature | Parsed from `IDEAL.csv` [dataset/1_disobind_databases.py:477-481]() |
| **MobiDB** | Database of Protein Disorder and Mobility Annotations | Parsed from `MobiDB.csv` [dataset/1_disobind_databases.py:483-487]() |
| **DIBS** | Disordered Binding Sites | `dibs_mfib()` method [dataset/1_disobind_databases.py:175-214]() |
| **MFIB** | Mutual Folding Induced by Binding | `dibs_mfib()` method [dataset/1_disobind_databases.py:175-214]() |
| **FuzDB** | Fuzzy Complexes Database | `fuzdb()` method [dataset/1_disobind_databases.py:273-315]() |
| **PDBtot** | PDB Total | `pdb_tot_cdr()` method [dataset/1_disobind_databases.py:318-364]() |
| **PDBcdr** | PDB Circular Dichroism | `pdb_tot_cdr()` method [dataset/1_disobind_databases.py:318-364]() |

**Sources:** [dataset/1_disobind_databases.py:38-56](), [dataset/1_disobind_databases.py:175-364]()

---

## Data Collection Workflow

The following diagram bridges the high-level data flow with the specific code entities responsible for each step.

### System Data Flow and Code Entities

```mermaid
flowchart TD
    Start["Start: Run 1_disobind_database.py"] --> Illuminati["Class: TheIlluminati<br/>Parse raw DB txt/tsv/xlsx"]
    
    subgraph "Database Parsing"
        Illuminati --> DIBS_MFIB["dibs_mfib()<br/>Extracts PDB IDs & UniProt boundaries"]
        Illuminati --> FuzDB["fuzdb()<br/>Parses browse_fuzdb.tsv"]
        Illuminati --> PDB_CDR["pdb_tot_cdr()<br/>Parses FP_pdbcdr_modified.xlsx"]
    end

    DIBS_MFIB --> Merge["merge_all_db_pdbs()<br/>Creates Merged_PDB_IDs.txt"]
    FuzDB --> Merge
    PDB_CDR --> Merge

    Merge --> Parser["Class: parse_pdbs_for_idrs<br/>forward()"]
    
    subgraph "Structure Retrieval"
        Parser --> PDB_REST["from_pdb_rest_api_with_love()<br/>Fetch entity/chain metadata"]
        Parser --> SIFTS["download_SIFTS_Uni_PDB_mapping()<br/>Get PDB-UniProt alignment"]
        Parser --> Download["download_pdb()<br/>Fetch .pdb or .cif structure"]
    end
    
    Download --> UniProt["parallelize_uniprot_seq_download()<br/>Fetch FASTA via get_uniprot_seq()"]
    UniProt --> Final["Filtered_PDB_entries.h5<br/>Uniprot_seq.json"]
```

**Workflow Summary:**

1.  **Parsing**: `TheIlluminati` class reads raw database files (e.g., `DIBS_complete_17Apr24.txt`) and converts them into structured dictionaries using `convert_txt_to_dict` [dataset/1_disobind_databases.py:58-172]().
2.  **Merging**: PDB IDs from all sources are aggregated and written to `Merged_PDB_IDs.txt` [dataset/1_disobind_databases.py:429-456]().
3.  **Validation**: `parse_pdbs_for_idrs` checks for obsolete PDB IDs using `get_superseding_pdb_id` [dataset/from_APIs_with_love.py:726-749]().
4.  **Mapping**: SIFTS mappings are downloaded to align PDB residues to UniProt positions [dataset/from_APIs_with_love.py:586-612]().
5.  **Sequencing**: UniProt sequences are fetched and filtered by length (default `uni_max_len = 10000`) [dataset/2_create_database_dataset_files.py:779-818]().

**Sources:** [dataset/1_disobind_databases.py:429-467](), [dataset/2_create_database_dataset_files.py:147-285](), [dataset/from_APIs_with_love.py:586-612]()

---

## Key Classes and Implementation

### `TheIlluminati` Class
This class is responsible for the initial extraction of PDB IDs and UniProt metadata from external disorder databases.

*   **`convert_txt_to_dict(txt)`**: A manual parser for the DIBS/MFIB text format. It identifies entries via `[Entry]` tags and extracts PDB IDs, chain IDs, and UniProt boundaries [dataset/1_disobind_databases.py:58-172]().
*   **`merge_all_db_pdbs()`**: Aggregates IDs from all loaded DataFrames, ensuring a unique set of PDB identifiers for the download pipeline [dataset/1_disobind_databases.py:429-456]().

**Sources:** [dataset/1_disobind_databases.py:31-56](), [dataset/1_disobind_databases.py:58-172]()

### `parse_pdbs_for_idrs` Class
The main orchestrator for downloading and validating structural data.

*   **`parallelize_PDB_download()`**: Uses a `multiprocessing.Pool` to execute `download_pdb_and_sift` across multiple cores [dataset/2_create_database_dataset_files.py:571-610]().
*   **`download_pdb_and_sift(pdb_id)`**: The atomic unit of work for structure collection. It performs the following sequence:
    1.  Check for superseding IDs [dataset/2_create_database_dataset_files.py:642-646]().
    2.  Fetch metadata via `from_pdb_rest_api_with_love` [dataset/2_create_database_dataset_files.py:657-658]().
    3.  Verify the entry is not a monomer (if `create_dataset=True`) [dataset/2_create_database_dataset_files.py:699-703]().
    4.  Download SIFTS mapping [dataset/2_create_database_dataset_files.py:734-738]().
    5.  Download PDB/CIF file [dataset/2_create_database_dataset_files.py:749-753]().

**Sources:** [dataset/2_create_database_dataset_files.py:613-776]()

---

## Detailed Data Flow: PDB to UniProt Mapping

One of the most critical steps is ensuring that structural coordinates can be mapped back to the canonical UniProt sequence.

### Mapping Entity Relationships

```mermaid
erDiagram
    "PDB ID" ||--o{ "Entity ID" : contains
    "Entity ID" ||--o{ "Asym ID (Chain)" : instance_of
    "Asym ID (Chain)" }|--|| "UniProt ID" : maps_via_SIFTS
    "UniProt ID" ||--o{ "Disorder Region" : has_annotation
    
    "parse_pdbs_for_idrs" {
        string pdb_id
        list uniprot_ids
        dict sifts_mapping
    }
```

**Implementation Details:**
- **SIFTS Parsing**: The `parse_sifts.py` script (called via shell) converts XML mappings into a TSV format containing `PDB_pos` and `Uniprot_pos` [dataset/from_APIs_with_love.py:527-582]().
- **Chimeric Chains**: The pipeline explicitly checks if a single PDB chain maps to multiple UniProt IDs using `is_chimera()`. These are typically excluded to avoid ambiguity in interaction labeling [dataset/from_APIs_with_love.py:843-853]().
- **Non-Protein Entities**: Entities like DNA, RNA, or small molecules are filtered out by checking the `rcsb_entity_polymer_info.type` field in the PDB REST API response [dataset/from_APIs_with_love.py:826-834]().

**Sources:** [dataset/from_APIs_with_love.py:527-582](), [dataset/from_APIs_with_love.py:826-853]()

---

## Output Files and Formats

The collection process produces several intermediate and final files in the `database/input_files/` and project-specific dataset directories.

| File Name | Format | Content Description |
|-----------|--------|---------------------|
| `Merged_PDB_IDs.txt` | Text (CSV) | A flat list of all PDB IDs collected from disorder databases [dataset/1_disobind_databases.py:456](). |
| `Selected_PDBs_info.h5` | HDF5 (DataFrame) | Metadata including PDB ID, Entity ID, Asym ID, and associated UniProt IDs [dataset/2_create_database_dataset_files.py:221](). |
| `Uniprot_seq.json` | JSON | Dictionary mapping UniProt Accession to `[sequence, protein_name]` [dataset/2_create_database_dataset_files.py:248](). |
| `Filtered_PDB_entries.h5` | HDF5 (DataFrame) | The final set of PDB-UniProt mappings after excluding obsolete IDs and long sequences [dataset/2_create_database_dataset_files.py:273](). |

**Sources:** [dataset/1_disobind_databases.py:456](), [dataset/2_create_database_dataset_files.py:221-273]()

---

## Parallelization and API Handling

To handle the thousands of PDB entries, the pipeline employs robust parallelization and error handling.

*   **CPU Parallelization**: The number of cores is set via the `-c` argument in `1_disobind_database.py` [dataset/1_disobind_databases.py:524-528]().
*   **API Resilience**: The `from_APIs_with_love.py` module wraps all `requests.get` calls in retry logic with exponential backoff.
    *   `max_trials` (default 10-25) [dataset/from_APIs_with_love.py:10-11]()
    *   `wait_time` (default 5-20 seconds) [dataset/from_APIs_with_love.py:10-11]()
*   **Large File Handling**: The system checks for the existence of files before downloading, allowing for interrupted runs to be resumed without redundant API calls [dataset/2_create_database_dataset_files.py:734-753]().

**Sources:** [dataset/1_disobind_databases.py:524-528](), [dataset/from_APIs_with_love.py:10-20](), [dataset/2_create_database_dataset_files.py:734-753]()

---