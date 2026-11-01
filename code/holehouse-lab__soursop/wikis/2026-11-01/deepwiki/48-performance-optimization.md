# Performance Optimization

<details>
<summary>Relevant source files</summary>

The following files were used as context for generating this wiki page:

- [.github/workflows/soursop-ci.yml](.github/workflows/soursop-ci.yml)
- [soursop/sssampling.py](soursop/sssampling.py)
- [soursop/sstrajectory.py](soursop/sstrajectory.py)
- [soursop/ssutils.py](soursop/ssutils.py)
- [soursop/tests/test_sssampling.py](soursop/tests/test_sssampling.py)
- [soursop/tests/test_ssutils.py](soursop/tests/test_ssutils.py)

</details>



## Purpose and Scope

This page documents performance optimization strategies and features available in SOURSOP for improving analysis speed, reducing memory footprint, and controlling computational resources. Topics covered include thread management for NumPy operations, parallel trajectory loading, lazy loading mechanisms, and caching strategies.

For information about general trajectory loading, see [Loading Trajectories](#3.1). For information about the plugin system architecture, see [Plugin Architecture](#8.1).

---

## Thread Control for NumPy Operations

SOURSOP provides explicit control over the number of threads used by NumPy's underlying BLAS libraries (MKL or OpenBLAS). This is critical when running multiple SOURSOP analyses concurrently or when working in environments with limited CPU resources.

### Setting Thread Count

The `ssutils.set_numpy_threads()` function allows you to limit NumPy's thread usage:

```python
from soursop import ssutils

# Limit NumPy to 4 threads
num_threads_set, library = ssutils.set_numpy_threads(4)
print(f"Set {num_threads_set} threads using {library}")
```

### Architecture and Implementation

The thread control system automatically detects the BLAS backend (MKL on Windows, MKL or OpenBLAS on Unix-like systems) and applies the appropriate configuration:

```mermaid
flowchart TD
    User["User calls set_numpy_threads(n)"]
    
    subgraph Detection["Library Detection"]
        Platform["platform.system()"]
        Windows{"Windows?"}
        LibSearch["_identify_library_paths()"]
        Locate["_locate_libraries()"]
    end
    
    subgraph Config["Thread Configuration"]
        MKL["_set_mkl_numpy_threads()"]
        OpenBLAS["_set_openblas_numpy_threads()"]
        CDLL["ctypes.CDLL / LoadLibrary"]
    end
    
    subgraph Paths["Path Resolution"]
        Env["Check CONDA_PREFIX / VIRTUAL_ENV"]
        Walk["os.walk() site-packages"]
        Filter["Filter by library name"]
    end
    
    User --> Platform
    Platform --> Windows
    Windows -->|Yes| MKL
    Windows -->|No| LibSearch
    LibSearch --> Locate
    Locate --> Env
    Env --> Walk
    Walk --> Filter
    Filter --> MKL
    Filter --> OpenBLAS
    
    MKL --> CDLL
    OpenBLAS --> CDLL
    CDLL --> Return["Return (threads_set, library)"]
```

**Sources:** [soursop/ssutils.py:137-151](), [soursop/ssutils.py:29-49](), [soursop/ssutils.py:89-109]()

### Platform-Specific Behavior

| Platform | Library Search | File Extension | Notes |
|----------|---------------|----------------|-------|
| Windows | Direct `mkl` import | `.dll` | Uses conda's bundled MKL |
| Linux | Environment scanning | `.so*` | Searches CONDA_PREFIX or VIRTUAL_ENV |
| macOS | Environment scanning | `.dylib*` | Searches for Intel MKL or OpenBLAS |

### Why Thread Control Matters

NumPy operations (especially linear algebra routines like matrix multiplications, SVD, and eigenvalue computations) will by default attempt to use all available CPU cores. This can cause issues when:

1. **Running multiple SOURSOP analyses in parallel** - Each process may spawn many threads, leading to oversubscription
2. **Working on shared compute resources** - Respecting resource limits is important
3. **Memory bandwidth limitations** - More threads don't always mean faster execution

**Sources:** [soursop/ssutils.py:131-136](), [soursop/ssutils.py:52-86]()

---

## Parallel Trajectory Loading

SOURSOP provides parallel trajectory loading capabilities through the `parallel_load_trjs()` function, which significantly reduces I/O time when loading multiple trajectory files.

### Architecture

```mermaid
flowchart LR
    subgraph Input["Input"]
        TrajList["traj_list: List[str]"]
        Topology["topology: str"]
        NCPUs["n_procs: int"]
    end
    
    subgraph ParallelLoad["parallel_load_trjs()"]
        Pool["multiprocessing.Pool(n_procs)"]
        Partial["functools.partial(SSTrajectory)"]
        Map["pool.map()"]
    end
    
    subgraph Workers["Worker Processes"]
        W1["Worker 1<br/>SSTrajectory(traj1)"]
        W2["Worker 2<br/>SSTrajectory(traj2)"]
        W3["Worker 3<br/>SSTrajectory(traj3)"]
        WN["Worker N<br/>SSTrajectory(trajN)"]
    end
    
    subgraph Output["Output"]
        SSTList["List[SSTrajectory]"]
    end
    
    TrajList --> ParallelLoad
    Topology --> ParallelLoad
    NCPUs --> Pool
    Pool --> Map
    Partial --> Map
    
    Map --> W1
    Map --> W2
    Map --> W3
    Map --> WN
    
    W1 --> SSTList
    W2 --> SSTList
    W3 --> SSTList
    WN --> SSTList
```

**Sources:** [soursop/sstrajectory.py:31-33](), [soursop/sssampling.py:297-310]()

### Usage in SamplingQuality

The `SamplingQuality` class automatically uses parallel loading when multiple trajectories are provided:

```mermaid
flowchart TD
    Init["SamplingQuality.__init__()"]
    CheckLen{"len(traj_list)"}
    CheckForce{"force_sequential?"}
    
    subgraph Sequential["Sequential Loading"]
        SeqLoop["for traj in traj_list"]
        SeqLoad["SSTrajectory([traj])"]
    end
    
    subgraph Parallel["Parallel Loading"]
        ParLoad["parallel_load_trjs()"]
        Pool2["Pool(n_cpus)"]
        MultiLoad["Concurrent SSTrajectory creation"]
    end
    
    Init --> CheckLen
    CheckLen -->|== 1| SeqLoad
    CheckLen -->|> 1| CheckForce
    CheckForce -->|True| Sequential
    CheckForce -->|False| Parallel
    
    SeqLoop --> SeqLoad
    ParLoad --> Pool2
    Pool2 --> MultiLoad
```

**Key Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `n_cpus` | int | `os.cpu_count()` | Number of parallel processes |
| `force_sequential` | bool | False | Override parallel loading |
| `truncate` | bool | False | Match trajectory lengths |

**Sources:** [soursop/sssampling.py:255-310](), [soursop/sssampling.py:107-178]()

### When to Use Sequential vs Parallel Loading

**Use Parallel Loading (Default) When:**
- Loading multiple trajectory files (>1)
- Files are on fast storage (SSD, local disk)
- Sufficient RAM is available (each process loads a full trajectory)

**Use Sequential Loading (`force_sequential=True`) When:**
- Memory is limited (loads one trajectory at a time)
- Files are on slow network storage (parallel I/O may hurt performance)
- Debugging trajectory loading issues
- Only 1-2 trajectories are being loaded

**Sources:** [soursop/sssampling.py:280-310]()

---

## Lazy Loading Mechanisms

SOURSOP implements lazy loading for certain expensive operations that may not always be needed.

### Single Protein Trajectory Lazy Loading

The `SSTrajectory` class uses a decorator-based lazy loading pattern for the combined protein trajectory:

```mermaid
flowchart TD
    subgraph SSTrajectory["SSTrajectory Instance"]
        Init["__init__()"]
        SingleProtTraj["__single_protein_traj = None"]
    end
    
    subgraph LazyDecorator["@lazy_loading_single_protein_trajectory"]
        Check{"__single_protein_traj<br/>is None?"}
        Load["__get_all_proteins()"]
        Cache["Set __single_protein_traj"]
    end
    
    subgraph Methods["Decorated Methods"]
        GetRg["get_overall_radius_of_gyration()"]
        GetAsph["get_overall_asphericity()"]
        GetRh["get_overall_hydrodynamic_radius()"]
    end
    
    Init --> SingleProtTraj
    GetRg --> Check
    GetAsph --> Check
    GetRh --> Check
    
    Check -->|Yes| Load
    Load --> Cache
    Cache --> Execute["Execute method"]
    Check -->|No| Execute
```

**Implementation Details:**

The `lazy_loading_single_protein_trajectory` decorator checks if the full protein trajectory has been loaded. If not, it extracts all protein atoms and creates a single `SSProtein` object:

- **Trigger:** First call to `get_overall_*()` methods
- **Cached:** Once loaded, reused for all subsequent calls
- **Memory Impact:** Single merged trajectory object instead of per-chain objects

**Sources:** [soursop/sstrajectory.py:34-58](), [soursop/sstrajectory.py:233-235](), [soursop/sstrajectory.py:665-689]()

---

## Caching Mechanisms

SOURSOP implements internal caching to avoid redundant computations.

### SamplingQuality Precomputed Cache

The `SamplingQuality` class maintains a `__precomputed` dictionary for expensive calculations:

```mermaid
flowchart LR
    subgraph Cache["__precomputed Dict"]
        direction TB
        TrjHel["'trj_helicity'"]
        RefHel["'ref_helicity'"]
        Custom["Custom keys"]
    end
    
    subgraph Method["compute_frac_helicity()"]
        CheckCache{"recompute==False<br/>AND in cache?"}
        ReturnCache["Return cached values"]
        Compute["Compute helicity"]
        StoreCache["Store in __precomputed"]
    end
    
    Call["Method Call"] --> CheckCache
    CheckCache -->|Yes| ReturnCache
    CheckCache -->|No| Compute
    Compute --> StoreCache
    StoreCache --> Cache
```

**Cached Computations:**

| Key | Computation | Cost | Recompute Flag |
|-----|-------------|------|----------------|
| `trj_helicity` | DSSP secondary structure | High | `recompute` parameter |
| `ref_helicity` | Reference DSSP | High | `recompute` parameter |

**Usage Pattern:**

```python
# First call - computes and caches
helicity_1, ref_helicity_1 = sq.compute_frac_helicity(proteinID=0)

# Second call - returns cached values
helicity_2, ref_helicity_2 = sq.compute_frac_helicity(proteinID=0)

# Force recomputation
helicity_3, ref_helicity_3 = sq.compute_frac_helicity(proteinID=0, recompute=True)
```

**Sources:** [soursop/sssampling.py:180](), [soursop/sssampling.py:434-479]()

---

## Memory Management Strategies

### Trajectory Truncation

When analyzing ongoing simulations or comparing trajectories of different lengths, the `truncate` parameter enables memory-efficient processing:

```mermaid
flowchart TD
    Truncate{"truncate=True?"}
    GetLengths["Get all trajectory lengths"]
    FindMin["min_length = np.min(lengths)"]
    
    subgraph TruncateLoop["For each trajectory"]
        Slice["traj[0:min_length]"]
        NewSST["SSTrajectory(TRJ=sliced)"]
    end
    
    Truncate -->|Yes| GetLengths
    GetLengths --> FindMin
    FindMin --> TruncateLoop
    Slice --> NewSST
    
    Truncate -->|No| UseOriginal["Use full trajectories"]
```

**When to Truncate:**

- **Ongoing simulations:** Different replicas may have different lengths
- **Memory constraints:** Analyzing a consistent subset reduces memory
- **Fair comparisons:** Ensures all analyses use same amount of data

**Memory Savings Example:**

```
Without truncation:
  Traj1: 10,000 frames × 500 atoms × 8 bytes = 40 MB
  Traj2: 50,000 frames × 500 atoms × 8 bytes = 200 MB
  Total: 240 MB

With truncation to 10,000 frames:
  Traj1: 10,000 frames × 500 atoms × 8 bytes = 40 MB
  Traj2: 10,000 frames × 500 atoms × 8 bytes = 40 MB
  Total: 80 MB (67% reduction)
```

**Sources:** [soursop/sssampling.py:312-378]()

### Stride Parameter

All trajectory loading functions accept a `stride` parameter (passed via `**kwargs`) to load every nth frame:

```python
# Load every 10th frame
sq = SamplingQuality(
    traj_list,
    reference_list,
    top_file='topology.pdb',
    stride=10  # Passed to SSTrajectory
)
```

**Memory Impact:** Loading with `stride=10` reduces memory by ~90%.

**Sources:** [soursop/sssampling.py:156-159]()

---

## Precomputed Reference Data

For sampling quality analysis without explicit reference trajectories, SOURSOP uses precomputed dihedral distributions from excluded volume (EV) simulations.

### PrecomputedDihedralInterface

```mermaid
flowchart TD
    subgraph Decision["Reference Selection"]
        RefList{"reference_list<br/>provided?"}
    end
    
    subgraph ExplicitRef["Explicit Reference"]
        LoadRef["Load reference trajectories"]
        CompDih["Compute dihedrals"]
        RefAngles["ref_psi_angles<br/>ref_phi_angles"]
    end
    
    subgraph PrecomputedRef["Precomputed Reference"]
        GetSeq["Extract sequence"]
        CleanSeq["Remove caps <, >"]
        Interface["PrecomputedDihedralInterface"]
        LoadEV["Load from EV_RESIDUE_MAPPER"]
        EVData["PHI_EV_ANGLES_DICT<br/>PSI_EV_ANGLES_DICT"]
    end
    
    RefList -->|Yes| ExplicitRef
    RefList -->|No| PrecomputedRef
    
    LoadRef --> CompDih
    CompDih --> RefAngles
    
    GetSeq --> CleanSeq
    CleanSeq --> Interface
    Interface --> LoadEV
    LoadEV --> EVData
    EVData --> RefAngles
```

**Benefits:**

- **No reference trajectory needed:** Saves disk I/O and memory
- **Consistent baseline:** All comparisons use same EV model
- **Fast initialization:** Lookup is faster than computing dihedrals

**Data Source:** Precomputed from extensive EV simulations stored in `ssdata.py`

**Sources:** [soursop/sssampling.py:211-233](), [soursop/ssdata.py:26-29]()

---

## Performance Best Practices

### Optimal Configuration Matrix

| Scenario | `n_cpus` | `force_sequential` | `stride` | `truncate` | `set_numpy_threads` |
|----------|----------|-------------------|----------|-----------|---------------------|
| Single trajectory, full analysis | 1 | N/A | 1 | False | 4-8 |
| Multiple trajectories, HPC | `cpu_count()` | False | 1 | False | 1-2 |
| Memory-limited system | 2-4 | True | 5-10 | True | 2 |
| Quick exploratory analysis | 4-8 | False | 10 | True | 2-4 |
| Production analysis, many replicates | `cpu_count()` | False | 1 | False | 1 |

### Workflow Optimization Strategy

```mermaid
flowchart TD
    Start["Analysis Task"]
    
    subgraph Planning["Planning Phase"]
        CountFiles{"Number of<br/>trajectory files?"}
        CheckMem{"Memory<br/>sufficient?"}
        NeedSpeed{"Need full<br/>resolution?"}
    end
    
    subgraph Config["Configuration"]
        ParallelLoad["n_cpus=cpu_count()<br/>force_sequential=False"]
        SeqLoad["force_sequential=True"]
        UseStride["stride=5-10"]
        UseTrunc["truncate=True"]
        LimitThreads["set_numpy_threads(2)"]
    end
    
    Start --> CountFiles
    CountFiles -->|> 5| ParallelLoad
    CountFiles -->|1-5| CheckMem
    CheckMem -->|Low| SeqLoad
    CheckMem -->|OK| ParallelLoad
    
    ParallelLoad --> LimitThreads
    SeqLoad --> NeedSpeed
    NeedSpeed -->|No| UseStride
    NeedSpeed -->|Yes| UseTrunc
```

**Sources:** [soursop/sssampling.py:107-184](), [soursop/ssutils.py:137-151]()

### Profiling and Monitoring

**Monitoring Thread Usage:**

```python
from threadpoolctl import threadpool_info

# Check current BLAS configuration
info = threadpool_info()
for library in info:
    print(f"{library['user_api']}: {library['num_threads']} threads")
```

**Monitoring Memory:**

```python
import psutil
import os

process = psutil.Process(os.getpid())
mem_mb = process.memory_info().rss / 1024**2
print(f"Current memory usage: {mem_mb:.1f} MB")
```

---

## Common Performance Issues

### Issue 1: Slow Parallel Loading

**Symptoms:** `parallel_load_trjs()` is slower than sequential loading

**Causes:**
- Network file system (NFS, SMB) with poor parallel I/O
- Too many processes competing for I/O bandwidth
- Small trajectory files where process overhead dominates

**Solutions:**
1. Set `force_sequential=True`
2. Copy trajectories to local storage first
3. Reduce `n_cpus` to 2-4

**Sources:** [soursop/sssampling.py:280-310]()

### Issue 2: High Memory Usage

**Symptoms:** System runs out of memory during analysis

**Causes:**
- Loading many large trajectories simultaneously
- NumPy operations creating large intermediate arrays

**Solutions:**
1. Use `stride` to reduce loaded frames
2. Enable `truncate=True` to match trajectory lengths
3. Set `force_sequential=True` to load one at a time
4. Process trajectories in batches

**Sources:** [soursop/sssampling.py:312-378]()

### Issue 3: NumPy Using Too Many Threads

**Symptoms:** System becomes unresponsive, CPU usage at 100% across all cores

**Causes:**
- NumPy defaulting to all available cores
- Multiple SOURSOP processes each spawning many threads

**Solutions:**
```python
from soursop import ssutils

# Limit NumPy to 2 threads per process
ssutils.set_numpy_threads(2)
```

**Sources:** [soursop/ssutils.py:137-151]()

---

## Testing Performance Optimizations

The test suite includes verification of thread control functionality:

```mermaid
flowchart LR
    subgraph TestSuite["test_ssutils.py"]
        TestThreads["test_set_numpy_threads()"]
        SetTwo["Set threads to 2"]
        Verify["Assert set_threads == 2"]
        CheckLib["Assert library != 'unknown'"]
    end
    
    TestThreads --> SetTwo
    SetTwo --> Verify
    Verify --> CheckLib
```

**Test Execution:**
```bash
pytest soursop/tests/test_ssutils.py::test_set_numpy_threads -v
```

**Sources:** [soursop/tests/test_ssutils.py:15-19](), [.github/workflows/soursop-ci.yml:46-50]()

---

## Summary

SOURSOP provides comprehensive performance optimization capabilities:

1. **Thread Control:** Explicit management of NumPy/BLAS threads via `ssutils.set_numpy_threads()`
2. **Parallel I/O:** Automatic parallel trajectory loading with `parallel_load_trjs()` 
3. **Lazy Loading:** Deferred initialization of expensive trajectory operations
4. **Caching:** Internal result caching in `SamplingQuality.__precomputed`
5. **Memory Management:** Trajectory truncation and stride options
6. **Precomputed Data:** EV reference data eliminates need for explicit reference trajectories

Optimal performance requires balancing parallelism, memory usage, and computational resources based on the specific analysis workflow and system constraints.