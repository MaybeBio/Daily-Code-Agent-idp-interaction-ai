# pinder.data

<details>
<summary>Relevant source files</summary>

The following files were used as context for generating this wiki page:

- [pyproject.toml](pyproject.toml)
- [src/pinder-data/pinder/data/pipeline/tasks.py](src/pinder-data/pinder/data/pipeline/tasks.py)
- [tests/data/test_get_data.py](tests/data/test_get_data.py)
- [tests/data/test_graphql.py](tests/data/test_graphql.py)
- [tests/data/test_tasks.py](tests/data/test_tasks.py)

</details>



The `pinder.data` package is responsible for generating the PINDER dataset from raw RCSB PDB data. It orchestrates the download of mmCIF files, processes biological assemblies, identifies interacting chains, annotates structures with domain and functional information, and produces the curated index and metadata files that power the core PINDER system.

For information about using the generated dataset, see [pinder.core](#9.1). For evaluation tools, see [pinder.eval](#9.4).

## Purpose and Scope

The `pinder.data` pipeline transforms raw PDB structures into a structured dataset suitable for protein-protein docking research. Key responsibilities include:

- **Data acquisition** from RCSB PDB via rsync
- **Biological assembly generation** from mmCIF coordinate files
- **Interacting chain identification** based on spatial proximity and buried surface area
- **Metadata extraction** including resolution, method, oligomeric state
- **Annotation enrichment** with Pfam, ECOD, SCOP, and CATH domain classifications
- **Task orchestration** for parallel processing of thousands of structures
- **Quality control** and dataset curation

This package is typically used by PINDER maintainers to regenerate the dataset rather than by end users who consume the pre-generated data.

Sources: [pyproject.toml:38-47](), [tests/data/test_get_data.py:1-368]()

## Architecture Overview

The pipeline follows a multi-stage architecture where each stage processes inputs from the previous stage and produces outputs consumed by the next.

```mermaid
graph TB
    subgraph "Stage 1: Download"
        RSYNC["rcsb_rsync.download_two_char_codes()"]
        MMCIF["mmCIF Files<br/>*_xyz-enrich.cif.gz"]
    end
    
    subgraph "Stage 2: Process"
        PROCESS["get_data.process_mmcif()"]
        BIOASM["generate_bio_assembly()"]
        ENTITIES["get_entities()"]
        CHAINS["get_interacting_chains()"]
        METADATA["get_metadata()"]
    end
    
    subgraph "Stage 3: Annotate"
        GRAPHQL["graphql.fetch_entry_annotations()"]
        PFAM["Pfam domains"]
        FEATURES["Structural features"]
        ANNOTS["ECOD/SCOP/CATH"]
    end
    
    subgraph "Stage 4: Aggregate"
        COLLECT["graphql.collect_rcsb_annotations()"]
        INDEX["index.parquet"]
        META["metadata.parquet"]
    end
    
    subgraph "Task Orchestration"
        TASKS["tasks.get_stage_tasks()"]
        CACHE["cache.get_uningested_mmcif()"]
        SCATTER["scatter.chunk_collection()"]
        RUN["tasks.run_task()"]
    end
    
    RSYNC --> MMCIF
    MMCIF --> PROCESS
    PROCESS --> BIOASM
    BIOASM --> ENTITIES
    BIOASM --> CHAINS
    PROCESS --> METADATA
    
    ENTITIES --> GRAPHQL
    GRAPHQL --> PFAM
    GRAPHQL --> FEATURES
    GRAPHQL --> ANNOTS
    
    PFAM --> COLLECT
    FEATURES --> COLLECT
    ANNOTS --> COLLECT
    METADATA --> COLLECT
    
    COLLECT --> INDEX
    COLLECT --> META
    
    TASKS --> CACHE
    CACHE --> SCATTER
    SCATTER --> RUN
    RUN --> PROCESS
```

Sources: [tests/data/test_get_data.py:12-368](), [tests/data/test_graphql.py:1-100]()

## Data Acquisition from RCSB

The pipeline downloads mmCIF files from RCSB PDB using rsync, organizing them by two-character codes derived from PDB IDs.

### Two-Character Code Organization

RCSB organizes structures using the middle two characters of the PDB ID. For example, PDB entry `7cm8` is stored under directory `cm/pdb_00007cm8/`.

```mermaid
graph LR
    RCSB["RCSB PDB<br/>rsync server"]
    DIRS["rsync directories<br/>get_rsync_directories()"]
    CODES["two_char_codes<br/>['1a', 'cm', 'fa', ...]"]
    DOWNLOAD["download_two_char_codes()"]
    LOCAL["local filesystem<br/>nextgen_rcsb/cm/pdb_00007cm8/"]
    
    RCSB --> DIRS
    DIRS --> CODES
    CODES --> DOWNLOAD
    DOWNLOAD --> LOCAL
```

### Key Functions

| Function | Purpose | Location |
|----------|---------|----------|
| `get_rsync_directories()` | Lists available two-character codes on RCSB server | [rcsb_rsync module]() |
| `get_rsync_two_char_pdb_entries()` | Lists PDB entries for a specific two-character code | [rcsb_rsync module]() |
| `download_two_char_codes()` | Downloads all mmCIF files for specified codes | [rcsb_rsync module]() |
| `get_two_char_codes_not_downloaded()` | Identifies missing two-character codes | [rcsb_rsync module]() |

Sources: [tests/data/test_get_data.py:12-38](), [src/pinder-data/pinder/data/pipeline/tasks.py:36-41]()

## mmCIF Processing Pipeline

The core processing logic extracts biological assemblies, identifies entities, and detects interacting chains from mmCIF coordinate files.

### Biological Assembly Generation

The `generate_bio_assembly()` function reads the `_pdbx_struct_assembly*` categories from mmCIF files and applies transformation matrices to create the biological assembly.

```mermaid
graph TB
    MMCIF["mmCIF file<br/>_xyz-enrich.cif.gz"]
    READ["read_mmcif_file()"]
    PDBX["pdbx_file object"]
    BIOASM["generate_bio_assembly()"]
    ASSEMBLY["bio_asm: AtomArray"]
    ENTITY_MAP["entity_map: DataFrame"]
    
    MMCIF --> READ
    READ --> PDBX
    PDBX --> BIOASM
    BIOASM --> ASSEMBLY
    BIOASM --> ENTITY_MAP
    
    ENTITY_MAP --> |"chain_id -> asym_id"| MAPPING["chain/asym_id mapping"]
    ENTITY_MAP --> |"asym_id -> entity_id"| ENTITIES["entity assignments"]
```

The `entity_map` DataFrame maps assembly chains to asymmetric unit IDs and entity IDs, which is crucial for tracking chains across transformations.

Sources: [tests/data/test_get_data.py:242-362]()

### Entity Extraction

The `get_entities()` function extracts information about molecular entities from the mmCIF file, including chain assignments, entity types, and sequence information.

```mermaid
graph LR
    PDBX["pdbx_file"]
    ENTITIES["get_entities()"]
    DF["DataFrame"]
    
    PDBX --> ENTITIES
    ENTITIES --> DF
    
    DF --> COLS["Columns:<br/>entity_id<br/>asym_id<br/>pdb_strand_id<br/>chain<br/>entry_id"]
```

Sources: [tests/data/test_get_data.py:107-149]()

### Interacting Chain Detection

The `get_interacting_chains()` function identifies pairs of chains that interact based on spatial proximity and calculates buried surface area (SASA).

```mermaid
graph TB
    BIOASM["bio_asm: AtomArray"]
    ENTITIES["entities: DataFrame"]
    PARAMS["Parameters:<br/>contact_threshold=5.0<br/>backbone_only=False"]
    
    INTERACTING["get_interacting_chains()"]
    
    CONTACTS["Contact Detection<br/>atom pairs within threshold"]
    SASA["SASA Calculation<br/>buried_sasa"]
    RESIDUES["Residue Counting<br/>n_residue_pairs"]
    
    OUTPUT["DataFrame:<br/>chain_1, chain_2<br/>asym_id_1, asym_id_2<br/>pdb_strand_id_1/2<br/>n_residue_pairs<br/>n_residues<br/>buried_sasa<br/>chain_1/2_residues"]
    
    BIOASM --> INTERACTING
    ENTITIES --> INTERACTING
    PARAMS --> INTERACTING
    
    INTERACTING --> CONTACTS
    INTERACTING --> SASA
    INTERACTING --> RESIDUES
    
    CONTACTS --> OUTPUT
    SASA --> OUTPUT
    RESIDUES --> OUTPUT
```

The function supports two modes:
- **`backbone_only=False`**: Considers all atoms for contact detection
- **`backbone_only=True`**: Only considers backbone atoms (CA, C, N, O)

Sources: [tests/data/test_get_data.py:40-163]()

### Complete Processing Workflow

The `process_mmcif()` function orchestrates the complete processing of a single mmCIF file.

| Step | Function | Output |
|------|----------|--------|
| 1. Read file | `read_mmcif_file()` | `pdbx_file` object |
| 2. Generate assembly | `generate_bio_assembly()` | `bio_asm` AtomArray, `entity_map` |
| 3. Extract entities | `get_entities()` | Entity DataFrame |
| 4. Map sequences | `sequence_mapping()` | Residue-level mapping |
| 5. Find interactions | `get_interacting_chains()` | Interacting chain pairs |
| 6. Extract metadata | `get_metadata()` | Resolution, method, dates |
| 7. Write outputs | File I/O | `{pdb_id}-metadata.tsv`, `{pdb_id}-assembly.cif` |

Sources: [tests/data/test_get_data.py:217-240](), [tests/data/test_get_data.py:364-368]()

## Annotation Pipeline

The annotation system enriches structures with domain classifications and functional information from RCSB using GraphQL queries.

### GraphQL Data Fetching

The `fetch_entry_annotations()` function queries the RCSB GraphQL API to retrieve comprehensive annotations for a PDB entry.

```mermaid
graph TB
    PDB_ID["pdb_id: '6Q0R'"]
    GRAPHQL["graphql.fetch_entry_annotations()"]
    API["RCSB GraphQL API<br/>data.rcsb.org/graphql"]
    JSON["annotation JSON file"]
    
    PARSE["parse_annotation_data()"]
    
    PFAM["pfam_df<br/>Domain annotations<br/>19 columns"]
    FEATURES["feature_df<br/>Structural features<br/>14 columns<br/>ASA, HELIX_P, SHEET, BINDING_SITE"]
    ANNOTS["annotation_df<br/>ECOD/SCOP/CATH<br/>12 columns"]
    EC["ec_df<br/>EC numbers"]
    
    PDB_ID --> GRAPHQL
    GRAPHQL --> API
    API --> JSON
    JSON --> PARSE
    
    PARSE --> PFAM
    PARSE --> FEATURES
    PARSE --> ANNOTS
    PARSE --> EC
```

### Annotation Types

| Annotation Type | Description | Column Count | Example Fields |
|-----------------|-------------|--------------|----------------|
| Pfam | Protein family domains | 19 | `rcsb_id`, `rcsb_pfam_identifier`, `pfam_accession` |
| Features | Structural features | 14 | `type`, `asym_id`, `beg_seq_id`, `end_seq_id` |
| Annotations | Domain classifications | 12 | `annotation_id`, `name`, `type`, `provenance_source` |
| EC | Enzyme classifications | Variable | `ec_number`, `enzyme_name` |

**Feature Types** include:
- `ASA` - Accessible surface area
- `HELIX_P`, `SHEET` - Secondary structure
- `BINDING_SITE` - Ligand binding sites
- `SCOP`, `SCOP2B_SUPERFAMILY`, `CATH`, `ECOD` - Domain classifications
- `UNOBSERVED_RESIDUE_XYZ`, `UNOBSERVED_ATOM_XYZ` - Missing coordinates
- `CIS-PEPTIDE` - Unusual peptide bonds

Sources: [tests/data/test_graphql.py:13-64]()

### Batch Annotation Processing

The `populate_rcsb_annotations()` function processes multiple PDB entries in parallel.

```mermaid
graph LR
    PDB_IDS["pdb_ids: list<br/>['6Q0R', '2A79', ...]"]
    POPULATE["populate_rcsb_annotations()"]
    WORKERS["parallel=True<br/>max_workers=2"]
    
    DIRS["Output directories:<br/>rcsb_annotations/pfam/<br/>rcsb_annotations/features/<br/>rcsb_annotations/annotations/"]
    
    FILES["Individual CSV files:<br/>{pdb_id}.csv.gz"]
    
    COLLECT["collect_rcsb_annotations()"]
    
    AGGREGATED["Aggregated files:<br/>pfam.csv.gz<br/>annotations_ecod.csv.gz<br/>annotations_scop.csv.gz<br/>annotations_cath.csv.gz<br/>features_*.csv.gz"]
    
    PDB_IDS --> POPULATE
    WORKERS --> POPULATE
    POPULATE --> DIRS
    DIRS --> FILES
    FILES --> COLLECT
    COLLECT --> AGGREGATED
```

Sources: [tests/data/test_graphql.py:66-100]()

## Task Orchestration System

The task orchestration system manages the parallel processing of thousands of structures through a batching and caching mechanism.

### Task Flow

```mermaid
graph TB
    INPUT["Input Type<br/>'cif', 'dimers', 'pdb_ids'"]
    
    GLOB["get_stage_inputs()<br/>DATA_GLOBS[input_type]"]
    
    GLOBS["Available glob functions:<br/>cif_glob()<br/>dimer_glob()<br/>entry_glob()<br/>pdb_id_glob()<br/>foldseek_pdb_glob()<br/>graph_type_glob()"]
    
    INPUTS["inputs: list<br/>Paths or strings"]
    
    CACHE["cache_func()<br/>Filter already processed"]
    
    REMAINING["remaining: list<br/>Unprocessed inputs"]
    
    SCATTER["scatter_method()<br/>chunk_collection()"]
    
    BATCHES["batches: list[list]<br/>batch_size groups"]
    
    RUN["run_task()<br/>task_func(**task_input)"]
    
    INPUT --> GLOB
    GLOB --> GLOBS
    GLOBS --> INPUTS
    INPUTS --> CACHE
    CACHE --> REMAINING
    REMAINING --> SCATTER
    SCATTER --> BATCHES
    BATCHES --> RUN
```

### Stage Task Generation

The `get_stage_tasks()` function orchestrates the creation of task batches for a pipeline stage.

| Parameter | Type | Purpose |
|-----------|------|---------|
| `data_dir` | `Path` | Root directory for data files |
| `input_type` | `str` | Type of input: 'cif', 'dimers', 'pdb_ids', etc. |
| `batch_size` | `int` | Number of inputs per batch |
| `cache_func` | `Callable` | Function to filter processed inputs |
| `cache_kwargs` | `dict` | Arguments for cache function |
| `scatter_method` | `Callable` | Function to split into batches (default: `chunk_collection`) |
| `scatter_kwargs` | `dict` | Arguments for scatter function |

Sources: [src/pinder-data/pinder/data/pipeline/tasks.py:117-138]()

### Task Execution

The `run_task()` function validates and executes a task function with proper error handling.

```mermaid
graph TB
    FUNC["task_func: Callable"]
    INPUT["task_input: dict"]
    ITER["iterable_kwarg: str | None"]
    
    INSPECT["Inspect function signature"]
    VALIDATE["Validate required arguments"]
    
    CHECK["Check iterable_kwarg length"]
    SKIP["Skip if empty"]
    
    EXECUTE["Execute task_func(**task_input)"]
    
    FUNC --> INSPECT
    INPUT --> INSPECT
    INSPECT --> VALIDATE
    VALIDATE --> CHECK
    ITER --> CHECK
    CHECK --> SKIP
    CHECK --> EXECUTE
```

The function performs signature validation to ensure all required parameters are provided before execution.

Sources: [src/pinder-data/pinder/data/pipeline/tasks.py:141-182](), [tests/data/test_tasks.py:95-104]()

### Available Input Types

The `DATA_GLOBS` dictionary maps input type names to glob functions:

```python
DATA_GLOBS = {
    "cif": cif_glob,                    # mmCIF files: */*/*enrich.cif.gz
    "dimers": dimer_glob,               # Dimer PDB files: */*/*__*--*__*.pdb
    "entries": entry_glob,              # Entry directories
    "two_char_codes": two_char_code_rsync,  # RCSB two-char codes to download
    "pdb_ids": pdb_id_glob,            # PDB IDs from entry directories
    "foldseek": foldseek_pdb_glob,     # Foldseek input PDBs: *.pdb
    "graph_types": graph_type_glob,    # Alignment graph types
    "putative_apo_pairings": putative_apo_pairings,  # Apo/holo pairings
}
```

Sources: [src/pinder-data/pinder/data/pipeline/tasks.py:87-96]()

### Caching Mechanism

Cache functions filter out already-processed inputs to avoid redundant work. Example cache functions:

| Cache Function | Purpose | Module |
|----------------|---------|--------|
| `get_uningested_mmcif()` | Find mmCIF files not yet processed | `cache` module |
| `get_pisa_unannotated()` | Find structures missing PISA annotations | `cache` module |
| `get_rcsb_unannotated()` | Find structures missing RCSB annotations | `cache` module |

Sources: [tests/data/test_tasks.py:22-93]()

## Command-Line Interface

The `pinder.data` package provides several CLI commands for pipeline execution and data management.

### Main Pipeline Command

The `pinder_data` command is the entry point for running the complete data generation pipeline.

```bash
pinder_data --help
```

Entry point: [pyproject.toml:39]()

### Related Commands

| Command | Purpose | Entry Point |
|---------|---------|-------------|
| `pinder_qc` | Quality control and validation | [pyproject.toml:46]() |
| `pinder_plot` | Visualization and plotting | [pyproject.toml:47]() |

## Module Reference

The `pinder.data` package contains several sub-modules organized by functionality:

| Module | Primary Purpose | Key Components |
|--------|-----------------|----------------|
| `rcsb_rsync` | Download data from RCSB PDB | `download_two_char_codes()`, `get_rsync_directories()` |
| `get_data` | Process mmCIF files | `process_mmcif()`, `generate_bio_assembly()`, `get_interacting_chains()` |
| `annotation.graphql` | Fetch RCSB annotations | `fetch_entry_annotations()`, `populate_rcsb_annotations()` |
| `pipeline.tasks` | Task orchestration | `get_stage_tasks()`, `run_task()` |
| `pipeline.cache` | Caching processed results | Cache validation functions |
| `pipeline.scatter` | Batch distribution | `chunk_collection()` |
| `foldseek_utils` | FoldSeek clustering utilities | `create_fasta_from_systems()` |
| `get_apo` | Apo structure identification | `get_putative_pairings()` |
| `get_dimers` | Dimer processing | Entry directory utilities |
| `qc` | Quality control | Validation and statistics |
| `plot` | Visualization | Plotting utilities |

### Dependencies

The `pinder.data` package requires additional dependencies beyond the core package:

```python
# From pyproject.toml [project.optional-dependencies.data]
data = [
    "fire",                      # CLI framework
    "gemmi",                     # mmCIF parsing
    "networkx",                  # Graph algorithms
    "scikit-learn <= 1.2.2",    # Clustering (PRODIGY-cryst)
    "python-graphql-client",     # GraphQL queries
    "mpire",                     # Parallel processing
    "matplotlib",                # Plotting
    "seaborn",                   # Statistical visualization
    "colormath",                 # Color manipulation
]
```

Sources: [pyproject.toml:96-106]()

## Integration with Core Dataset

The outputs of `pinder.data` are consumed by `pinder.core` to provide the user-facing dataset interface.

```mermaid
graph TB
    subgraph "pinder.data Output"
        PDBS["pdbs/<br/>Processed PDB files"]
        INDEX["index.parquet<br/>2.4M+ entries"]
        META["metadata.parquet<br/>Detailed annotations"]
        MAPPINGS["mappings/<br/>UniProt mappings"]
    end
    
    subgraph "pinder.core Consumption"
        LOADER["PinderLoader"]
        SYSTEM["PinderSystem"]
        DATASET["PinderDataset"]
    end
    
    INDEX --> LOADER
    META --> LOADER
    LOADER --> SYSTEM
    PDBS --> SYSTEM
    MAPPINGS --> SYSTEM
    SYSTEM --> DATASET
```

For details on how to use the generated dataset, see [pinder.core](#9.1) and [Data Loading System](#4).

Sources: [pyproject.toml:138-147]()

---