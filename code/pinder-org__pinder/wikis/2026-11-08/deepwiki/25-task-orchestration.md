# Task Orchestration

<details>
<summary>Relevant source files</summary>

The following files were used as context for generating this wiki page:

- [src/pinder-data/pinder/data/pipeline/tasks.py](src/pinder-data/pinder/data/pipeline/tasks.py)
- [tests/data/test_graphql.py](tests/data/test_graphql.py)
- [tests/data/test_tasks.py](tests/data/test_tasks.py)

</details>



## Purpose and Scope

This document describes the task orchestration system in the PINDER data generation pipeline. The orchestration layer handles parallel processing of large-scale data transformations by batching inputs, caching processed results, and distributing work across multiple workers. This system sits between the data acquisition layer (see [Data Acquisition from RCSB](#5.1)) and the specific processing tasks like mmCIF parsing (see [mmCIF Processing](#5.3)) and annotation fetching (see [RCSB Annotations](#5.4)).

The orchestration system is implemented in [`src/pinder-data/pinder/data/pipeline/tasks.py`](src/pinder-data/pinder/data/pipeline/tasks.py:1-182)() and provides a unified interface for executing different pipeline stages with automatic progress tracking, failure handling, and incremental processing.

## Task Orchestration Flow

The task orchestration system follows a three-stage pattern: (1) collect inputs from the filesystem or external sources, (2) filter out already-processed items using cache functions, and (3) batch remaining work into parallel tasks.

```mermaid
graph TB
    subgraph "Stage 1: Input Collection"
        DATA_DIR["data_dir (Path)"]
        INPUT_TYPE["input_type (str)<br/>e.g., 'cif', 'pdb_ids'"]
        DATA_GLOBS["DATA_GLOBS dict<br/>Maps type → glob function"]
        GLOB_FUNC["glob function<br/>e.g., cif_glob(), pdb_id_glob()"]
        RAW_INPUTS["raw inputs<br/>list[Path] or list[str]"]
    end
    
    subgraph "Stage 2: Cache Filtering"
        CACHE_FUNC["cache_func<br/>e.g., get_uningested_mmcif()"]
        CACHE_KWARGS["cache_kwargs<br/>Additional parameters"]
        FILTERED_INPUTS["filtered inputs<br/>Items not yet processed"]
    end
    
    subgraph "Stage 3: Task Batching"
        SCATTER_METHOD["scatter_method<br/>Default: chunk_collection()"]
        BATCH_SIZE["batch_size (int)<br/>Items per batch"]
        TASK_BATCHES["task batches<br/>list[list[...]]"]
    end
    
    subgraph "Execution"
        RUN_TASK["run_task()<br/>Execute individual batch"]
        TASK_FUNC["task_func<br/>Actual processing function"]
        VALIDATION["Signature validation"]
    end
    
    DATA_DIR --> GLOB_FUNC
    INPUT_TYPE --> DATA_GLOBS
    DATA_GLOBS --> GLOB_FUNC
    GLOB_FUNC --> RAW_INPUTS
    
    RAW_INPUTS --> CACHE_FUNC
    CACHE_KWARGS --> CACHE_FUNC
    CACHE_FUNC --> FILTERED_INPUTS
    
    FILTERED_INPUTS --> SCATTER_METHOD
    BATCH_SIZE --> SCATTER_METHOD
    SCATTER_METHOD --> TASK_BATCHES
    
    TASK_BATCHES --> RUN_TASK
    RUN_TASK --> VALIDATION
    VALIDATION --> TASK_FUNC
    
    style DATA_GLOBS fill:#f9f9f9
    style CACHE_FUNC fill:#f9f9f9
    style SCATTER_METHOD fill:#f9f9f9
    style RUN_TASK fill:#f9f9f9
```

**Sources:** [src/pinder-data/pinder/data/pipeline/tasks.py:99-139]()

## Input Collection System

The orchestration system uses a registry of input collection functions to gather work items for different pipeline stages. Each input type has a corresponding "glob" function that scans the filesystem or queries external sources.

### Input Types and Glob Functions

| Input Type | Glob Function | Return Type | Description |
|------------|---------------|-------------|-------------|
| `cif` | `cif_glob()` | `list[Path]` | mmCIF files matching `*/*/*enrich.cif.gz` |
| `dimers` | `dimer_glob()` | `list[Path]` | Dimer PDB files matching `*/*/*__*--*__*.pdb` |
| `entries` | `entry_glob()` | `list[Path]` | PDB entry directories |
| `two_char_codes` | `two_char_code_rsync()` | `list[str]` | Two-character codes from RCSB rsync |
| `pdb_ids` | `pdb_id_glob()` | `list[str]` | PDB IDs extracted from entry directories |
| `foldseek` | `foldseek_pdb_glob()` | `list[Path]` | PDB files for FoldSeek clustering |
| `graph_types` | `graph_type_glob()` | `list[tuple[str, bool]]` | Graph type specifications |
| `putative_apo_pairings` | `putative_apo_pairings()` | `pd.DataFrame` | Apo-holo pairing candidates |

**Sources:** [src/pinder-data/pinder/data/pipeline/tasks.py:23-96]()

### DATA_GLOBS Registry

The `DATA_GLOBS` dictionary maps input type strings to their corresponding glob functions, enabling dynamic dispatch based on stage requirements:

```python
DATA_GLOBS: dict[str, Callable[[Path], INPUT_TYPES]] = {
    "cif": cif_glob,
    "dimers": dimer_glob,
    "entries": entry_glob,
    "two_char_codes": two_char_code_rsync,
    "pdb_ids": pdb_id_glob,
    "foldseek": foldseek_pdb_glob,
    "graph_types": graph_type_glob,
    "putative_apo_pairings": putative_apo_pairings,
}
```

The `get_stage_inputs()` function uses this registry to collect inputs for a given pipeline stage:

**Sources:** [src/pinder-data/pinder/data/pipeline/tasks.py:87-105]()

## Caching Strategy

The orchestration system implements incremental processing by filtering out items that have already been processed. This is accomplished through cache functions that determine which inputs still require processing.

```mermaid
graph LR
    subgraph "Cache Functions"
        GET_UNINGESTED["get_uningested_mmcif()<br/>Returns unprocessed mmCIF files"]
        GET_PISA["get_pisa_unannotated()<br/>Returns entries without PISA data"]
        GET_RCSB["get_rcsb_unannotated()<br/>Returns entries without RCSB annotations"]
    end
    
    subgraph "Cache Delta Calculation"
        ALL_INPUTS["all_inputs<br/>Complete input set"]
        CACHE_FUNC["cache_func(batches, **kwargs)<br/>Filters out processed items"]
        REMAINING["remaining<br/>Items needing processing"]
    end
    
    subgraph "Cache Artifacts"
        INGESTED_DIR["data_dir/ingested/<br/>Processed mmCIF markers"]
        PISA_DIR["data_dir/pisa/<br/>PISA annotation files"]
        RCSB_DIR["pinder_dir/rcsb_annotations/<br/>RCSB annotation files"]
    end
    
    ALL_INPUTS --> CACHE_FUNC
    
    GET_UNINGESTED -.->|checks| INGESTED_DIR
    GET_PISA -.->|checks| PISA_DIR
    GET_RCSB -.->|checks| RCSB_DIR
    
    CACHE_FUNC --> REMAINING
    
    style CACHE_FUNC fill:#f9f9f9
```

The `get_cache_delta()` function applies a cache function to filter inputs:

**Sources:** [src/pinder-data/pinder/data/pipeline/tasks.py:108-114](), [tests/data/test_tasks.py:22-92]()

### Example: Apo Pairing Cache

The `putative_apo_pairings()` function demonstrates a more complex caching strategy that checks for existing metrics to avoid recomputing pairing evaluations:

**Sources:** [src/pinder-data/pinder/data/pipeline/tasks.py:68-84]()

## Task Batching and Scatter Methods

After filtering cached items, the orchestration system batches remaining work into parallel tasks. The default scatter method is `chunk_collection()`, which divides inputs into fixed-size batches.

```mermaid
graph TB
    subgraph "Scatter Input"
        FILTERED["filtered_inputs<br/>1000 items"]
        BATCH_SIZE["batch_size: 50"]
    end
    
    subgraph "chunk_collection()"
        CHUNK["Divide into chunks<br/>of size batch_size"]
    end
    
    subgraph "Scatter Output"
        BATCH1["Batch 0<br/>Items 0-49"]
        BATCH2["Batch 1<br/>Items 50-99"]
        BATCH3["Batch 2<br/>Items 100-149"]
        BATCHN["...<br/>Batch 19<br/>Items 950-999"]
    end
    
    FILTERED --> CHUNK
    BATCH_SIZE --> CHUNK
    
    CHUNK --> BATCH1
    CHUNK --> BATCH2
    CHUNK --> BATCH3
    CHUNK --> BATCHN
    
    style CHUNK fill:#f9f9f9
```

The `get_stage_tasks()` function orchestrates the complete flow from input collection to task batching:

**Sources:** [src/pinder-data/pinder/data/pipeline/tasks.py:117-138]()

## Task Execution and Validation

The `run_task()` function executes individual task batches with built-in validation and error handling. It inspects the task function's signature to ensure all required arguments are provided before execution.

### Task Execution Flow

```mermaid
graph TB
    subgraph "Task Input"
        TASK_FUNC["task_func<br/>Callable to execute"]
        TASK_INPUT["task_input<br/>dict of parameters"]
        ITERABLE_KWARG["iterable_kwarg<br/>Name of batch parameter"]
    end
    
    subgraph "Validation"
        INSPECT["inspect.signature(task_func)<br/>Get required parameters"]
        CHECK_MISSING["Check for missing<br/>required arguments"]
        CHECK_EMPTY["Check if iterable_kwarg<br/>is empty"]
    end
    
    subgraph "Execution"
        LOG_START["Log task start<br/>with input details"]
        EXECUTE["task_func(**task_input)<br/>Execute with unpacked kwargs"]
    end
    
    subgraph "Error Handling"
        LOG_ERROR["Log missing arguments"]
        LOG_SKIP["Log skipped empty batch"]
        RETURN["Early return"]
    end
    
    TASK_FUNC --> INSPECT
    TASK_INPUT --> CHECK_MISSING
    ITERABLE_KWARG --> CHECK_EMPTY
    
    INSPECT --> CHECK_MISSING
    CHECK_MISSING -->|"Missing args"| LOG_ERROR
    LOG_ERROR --> RETURN
    
    CHECK_MISSING -->|"All present"| CHECK_EMPTY
    CHECK_EMPTY -->|"Empty batch"| LOG_SKIP
    LOG_SKIP --> RETURN
    
    CHECK_EMPTY -->|"Has items"| LOG_START
    LOG_START --> EXECUTE
    
    style INSPECT fill:#f9f9f9
    style EXECUTE fill:#f9f9f9
```

**Sources:** [src/pinder-data/pinder/data/pipeline/tasks.py:141-181]()

### Signature Validation Logic

The task execution system validates function signatures before execution to catch configuration errors early:

1. **Extract required parameters** using `inspect.signature(task_func)`
2. **Check for missing arguments** by comparing signature parameters against `task_input` keys
3. **Validate iterable argument** by checking if the designated batch parameter is non-empty
4. **Log and skip** if validation fails, preventing wasted computation

**Sources:** [src/pinder-data/pinder/data/pipeline/tasks.py:154-176]()

## Pipeline Stage Integration

The orchestration system integrates with multiple pipeline stages, each with specific input types, cache functions, and processing tasks.

### Example Pipeline Stages

```mermaid
graph TB
    subgraph "Stage: ingest_rcsb_files"
        INGEST_INPUT["input_type: 'cif'"]
        INGEST_CACHE["cache_func:<br/>get_uningested_mmcif()"]
        INGEST_SCATTER["scatter_method:<br/>chunk_collection()"]
        INGEST_BATCH["batch_size: 2"]
        INGEST_TASKS["Task batches:<br/>~4 batches for 8 files"]
    end
    
    subgraph "Stage: get_pisa_annotations"
        PISA_INPUT["input_type: 'cif'"]
        PISA_CACHE["cache_func:<br/>get_pisa_unannotated()"]
        PISA_SCATTER["scatter_method:<br/>chunk_collection()"]
        PISA_BATCH["batch_size: 2"]
        PISA_TASKS["Task batches:<br/>~6 batches for 12 files"]
    end
    
    subgraph "Stage: get_rcsb_annotations"
        RCSB_INPUT["input_type: 'pdb_ids'"]
        RCSB_CACHE["cache_func:<br/>get_rcsb_unannotated()"]
        RCSB_SCATTER["scatter_method:<br/>chunk_collection()"]
        RCSB_BATCH["batch_size: 1"]
        RCSB_TASKS["Task batches:<br/>~13 batches for 13 IDs"]
    end
    
    INGEST_INPUT --> INGEST_CACHE
    INGEST_CACHE --> INGEST_SCATTER
    INGEST_SCATTER --> INGEST_BATCH
    INGEST_BATCH --> INGEST_TASKS
    
    PISA_INPUT --> PISA_CACHE
    PISA_CACHE --> PISA_SCATTER
    PISA_SCATTER --> PISA_BATCH
    PISA_BATCH --> PISA_TASKS
    
    RCSB_INPUT --> RCSB_CACHE
    RCSB_CACHE --> RCSB_SCATTER
    RCSB_SCATTER --> RCSB_BATCH
    RCSB_BATCH --> RCSB_TASKS
    
    style INGEST_CACHE fill:#f9f9f9
    style PISA_CACHE fill:#f9f9f9
    style RCSB_CACHE fill:#f9f9f9
```

**Sources:** [tests/data/test_tasks.py:22-92]()

### Stage Configuration Parameters

Each pipeline stage is configured with the following parameters passed to `get_stage_tasks()`:

| Parameter | Type | Description |
|-----------|------|-------------|
| `data_dir` | `Path` | Root directory for pipeline data |
| `input_type` | `str` | Key into `DATA_GLOBS` registry |
| `batch_size` | `int` | Number of items per task batch |
| `cache_func` | `Callable` or `None` | Function to filter processed items |
| `cache_kwargs` | `dict` | Additional arguments for cache function |
| `scatter_method` | `Callable` | Batching function (default: `chunk_collection`) |
| `scatter_kwargs` | `dict` | Additional arguments for scatter method |

**Sources:** [src/pinder-data/pinder/data/pipeline/tasks.py:117-125]()

## Complete Orchestration Example

The following diagram shows how the orchestration system processes a complete pipeline stage from input collection to task execution:

```mermaid
sequenceDiagram
    participant Pipeline
    participant get_stage_tasks
    participant get_stage_inputs
    participant DATA_GLOBS
    participant cache_func
    participant scatter_method
    participant run_task
    participant task_func
    
    Pipeline->>get_stage_tasks: data_dir, input_type="cif",<br/>batch_size=2, cache_func=...
    
    get_stage_tasks->>get_stage_inputs: data_dir, input_type="cif"
    get_stage_inputs->>DATA_GLOBS: Lookup "cif"
    DATA_GLOBS-->>get_stage_inputs: cif_glob function
    get_stage_inputs->>get_stage_inputs: Execute cif_glob(data_dir)
    get_stage_inputs-->>get_stage_tasks: [file1.cif, file2.cif, ..., file8.cif]
    
    get_stage_tasks->>cache_func: Filter processed files
    cache_func-->>get_stage_tasks: [file1.cif, file2.cif, ..., file8.cif]<br/>(all need processing)
    
    get_stage_tasks->>scatter_method: batches=[...], batch_size=2
    scatter_method-->>get_stage_tasks: [[file1, file2], [file3, file4],<br/>[file5, file6], [file7, file8]]
    
    get_stage_tasks-->>Pipeline: 4 task batches
    
    loop For each batch
        Pipeline->>run_task: task_func, task_input={cifs: [file1, file2]}
        run_task->>run_task: Validate signature
        run_task->>run_task: Check for empty batch
        run_task->>task_func: Execute with cifs=[file1, file2]
        task_func-->>run_task: Completed
        run_task-->>Pipeline: Success
    end
```

**Sources:** [src/pinder-data/pinder/data/pipeline/tasks.py:99-181](), [tests/data/test_tasks.py:95-104]()

## Task Dependencies and Imports

The orchestration system imports processing functions from other modules in the data pipeline:

| Module | Imported Functions | Purpose |
|--------|-------------------|---------|
| `pinder.data.rcsb_rsync` | `get_rsync_directories()`, `get_two_char_codes_not_downloaded()`, `download_two_char_codes()` | RCSB data synchronization |
| `pinder.data.get_dimers` | `get_pdb_entry_dirs()` | Entry directory discovery |
| `pinder.data.get_apo` | `get_putative_pairings()` | Apo-holo pairing candidates |
| `pinder.data.pipeline.scatter` | `chunk_collection()` | Default batching strategy |

**Sources:** [src/pinder-data/pinder/data/pipeline/tasks.py:1-12](), [tests/data/test_tasks.py:95-104]()

## Logging and Progress Tracking

The orchestration system uses structured logging to track pipeline progress. Key log messages include:

- **Input summary**: `"{input_type}: {len(batches)} inputs -> {len(tasks)} task batches"`
- **Task start**: `"Running {task_name}: {arg_msg}"`
- **Empty batch skip**: `"Skipping {task_name}: {iterable_kwarg} inputs are empty..."`
- **Missing arguments**: `"Missing required arguments for {task_name}: {missing_args}"`

**Sources:** [src/pinder-data/pinder/data/pipeline/tasks.py:137-179]()

---