# System Architecture

<details>
<summary>Relevant source files</summary>

The following files were used as context for generating this wiki page:

- [README.md](README.md)
- [docs/index.rst](docs/index.rst)
- [docs/modules/sssampling.rst](docs/modules/sssampling.rst)
- [soursop/ssprotein.py](soursop/ssprotein.py)
- [soursop/sssampling.py](soursop/sssampling.py)
- [soursop/sstrajectory.py](soursop/sstrajectory.py)
- [soursop/tests/test_sssampling.py](soursop/tests/test_sssampling.py)

</details>



## Purpose and Scope

This document describes the overall system architecture of SOURSOP, including the design philosophy, core components, and how they interact. It covers the three primary classes (`SSTrajectory`, `SSProtein`, `SamplingQuality`), the data flow pipeline, supporting utilities, and the layered module structure.

For detailed API documentation of individual classes, see:
- [SSTrajectory: Loading and Multi-Chain Analysis](#3)
- [SSProtein: Single Protein Analysis](#4)
- [SamplingQuality: PENGUIN Pipeline](#5)

For information about extending SOURSOP, see [Plugin Extension System](#8).

## Core Architecture Overview

SOURSOP is built as a three-layer system with a clear separation of concerns:

```mermaid
graph TB
    subgraph "User Interface Layer"
        SST["SSTrajectory<br/>Multi-chain trajectory loading<br/>System-level analysis"]
        SSP["SSProtein<br/>Single-chain analysis<br/>50+ methods"]
        SSQ["SamplingQuality<br/>PENGUIN conformational<br/>quality assessment"]
    end
    
    subgraph "Data Backend Layer"
        MDTraj["mdtraj.Trajectory<br/>Core trajectory representation"]
        Topology["mdtraj.Topology<br/>System structure"]
    end
    
    subgraph "Analysis Layer"
        SpecMod["Specialized Modules<br/>ssnmr, sspre, ssmutualinformation"]
        UtilMod["Utility Modules<br/>ssdata, sstools, ssutils"]
    end
    
    subgraph "Infrastructure Layer"
        IO["ssio<br/>Message handling"]
        Exc["ssexceptions<br/>Error types"]
        Configs["configs<br/>System settings"]
    end
    
    SST --> MDTraj
    SST --> Topology
    SST -->|"extracts chains"| SSP
    SSP --> MDTraj
    SSQ --> SST
    SSQ --> SSP
    
    SSP --> SpecMod
    SSQ --> SpecMod
    
    SST --> UtilMod
    SSP --> UtilMod
    SSQ --> UtilMod
    SpecMod --> UtilMod
    
    UtilMod --> IO
    UtilMod --> Exc
    UtilMod --> Configs
```

**Sources:** [soursop/sstrajectory.py:62-661](), [soursop/ssprotein.py:41-2506](), [soursop/sssampling.py:106-1200]()

## Core Classes

### SSTrajectory: Multi-Chain Trajectory Loading

`SSTrajectory` is the entry point for loading simulation data. It wraps `mdtraj.Trajectory` and provides protein-specific functionality.

| Property | Type | Description |
|----------|------|-------------|
| `traj` | `mdtraj.Trajectory` | Underlying trajectory object |
| `proteinTrajectoryList` | `list[SSProtein]` | List of extracted protein chains |
| `n_frames` | `int` | Number of trajectory frames |
| `n_proteins` | `int` | Number of protein chains detected |
| `unitcell` | `np.ndarray` | Unit cell dimensions (Angstroms) |

**Key Initialization Parameters:**
- `trajectory_filename`: Path to trajectory file (`.xtc`, `.dcd`)
- `pdb_filename`: Topology file (`.pdb`, `.gro`)
- `protein_grouping`: Manual residue grouping (optional)
- `explicit_residue_checking`: Validate each residue individually

**Automatic Protein Extraction:**
The class automatically identifies protein chains by scanning the topology for valid residue names defined in `ALL_VALID_RESIDUE_NAMES` from `ssdata`. Chain extraction occurs in `__get_proteins()` [soursop/sstrajectory.py:439-559](), which cycles through each chain and creates `SSProtein` objects via atom slicing.

**Sources:** [soursop/sstrajectory.py:62-243](), [soursop/ssdata.py:1-50]()

### SSProtein: Single Protein Analysis

`SSProtein` represents a single protein chain and provides over 50 analysis methods. Each method operates on the chain's trajectory frames.

| Property | Type | Description |
|----------|------|-------------|
| `traj` | `mdtraj.Trajectory` | Single-chain trajectory |
| `topology` | `mdtraj.Topology` | Chain topology |
| `n_frames` | `int` | Number of frames |
| `n_residues` | `int` | Number of residues (including caps) |
| `resid_with_CA` | `list[int]` | Residues with CA atoms |
| `ncap` | `bool` | N-terminal cap present |
| `ccap` | `bool` | C-terminal cap present |

**Caching and Performance:**
The class implements aggressive memoization through private dictionaries (`__CA_residue_atom`, `__residue_atom_table`, `__SASA_saved`, `__all_angles`) to avoid expensive topology lookups. The `__residue_atom_lookup()` method [soursop/ssprotein.py:684-749]() converts O(n) topology selections to O(1) dictionary lookups after first access.

**Coarse-Grained Optimization:**
For one-bead-per-residue systems, `__check_cg_onebead()` [soursop/ssprotein.py:531-547]() detects this and triggers fast initialization via `__initialize_cg_atoms()` [soursop/ssprotein.py:665-679](), providing ~30x speedup.

**Sources:** [soursop/ssprotein.py:41-218](), [soursop/ssprotein.py:531-679]()

### SamplingQuality: PENGUIN Pipeline

`SamplingQuality` implements the PENGUIN (Pipeline for Evaluating coNformational heteroGeneity in Unstructured proteINs) methodology for assessing conformational sampling quality.

| Attribute | Type | Description |
|-----------|------|-------------|
| `trajs` | `list[SSTrajectory]` | Simulated trajectories |
| `ref_trajs` | `list[SSTrajectory]` | Reference trajectories (or None) |
| `psi_angles` | `np.ndarray` | Psi dihedral angles |
| `phi_angles` | `np.ndarray` | Phi dihedral angles |
| `method` | `str` | `"1D angle distributions"` or `"2D angle distributions"` |
| `bwidth` | `float` | Histogram bin width (radians) |

**Parallel Loading:**
Trajectories are loaded in parallel using `parallel_load_trjs()` [soursop/sstrajectory.py:1382-1413]() unless `force_sequential=True`. The `n_cpus` parameter controls parallelization.

**Reference Models:**
When `reference_list=None`, precomputed excluded volume (EV) dihedral distributions are loaded via `PrecomputedDihedralInterface` [soursop/sssampling.py:1203-1339](). These EV angle distributions are stored in `PHI_EV_ANGLES_DICT` and `PSI_EV_ANGLES_DICT` from `ssdata`.

**Sources:** [soursop/sssampling.py:106-234](), [soursop/ssdata.py:300-450]()

## Data Flow Pipeline

```mermaid
flowchart TD
    Files["Trajectory Files<br/>xtc/dcd + pdb/gro"]
    
    Load["SSTrajectory.__init__<br/>__readTrajectory<br/>mdtraj.load"]
    
    Parse["__get_proteins<br/>topology.chains iteration<br/>valid residue checking"]
    
    Slice["trajectory.atom_slice<br/>per-chain extraction"]
    
    Create["SSProtein.__init__<br/>__get_resid_with_CA<br/>cache initialization"]
    
    Store["proteinTrajectoryList<br/>list[SSProtein]"]
    
    Analyze1["SSProtein methods<br/>get_radius_of_gyration<br/>get_distance_map<br/>get_secondary_structure_DSSP"]
    
    Analyze2["SSTrajectory methods<br/>get_interchain_distance_map<br/>get_overall_radius_of_gyration"]
    
    Analyze3["SamplingQuality<br/>compute_dihedral_hellingers<br/>quality_plot"]
    
    Files --> Load
    Load --> Parse
    Parse --> Slice
    Slice --> Create
    Create --> Store
    Store --> Analyze1
    Store --> Analyze2
    Store --> Analyze3
```

**Loading Phase** [soursop/sstrajectory.py:283-364](): `mdtraj.load()` reads trajectory and topology files. Optional `pdblead=True` prepends the PDB structure as the first frame.

**Parsing Phase** [soursop/sstrajectory.py:439-559](): For each chain in `topology.chains`, the first residue name is checked against `ALL_VALID_RESIDUE_NAMES`. If `explicit_residue_checking=True`, every residue is validated individually [soursop/sstrajectory.py:498-511]().

**Extraction Phase**: `trajectory.atom_slice()` creates per-chain sub-trajectories. Each is wrapped in an `SSProtein` object which initializes CA atom lookup via `__get_resid_with_CA()` [soursop/ssprotein.py:613-663]().

**Analysis Phase**: Users access methods on `SSProtein` objects for single-chain analysis or `SSTrajectory` methods for multi-chain analysis.

**Sources:** [soursop/sstrajectory.py:67-243](), [soursop/sstrajectory.py:439-559](), [soursop/ssprotein.py:59-171]()

## Supporting Modules

### Data and Constants (ssdata)

Provides amino acid mappings, residue validation lists, and precomputed reference data.

| Symbol | Type | Purpose |
|--------|------|---------|
| `ALL_VALID_RESIDUE_NAMES` | `list[str]` | Residues recognized as protein |
| `THREE_TO_ONE` | `dict` | 3-letter to 1-letter AA codes |
| `ONE_TO_THREE` | `dict` | 1-letter to 3-letter AA codes |
| `PHI_EV_ANGLES_DICT` | `dict` | Precomputed phi angles for EV model |
| `PSI_EV_ANGLES_DICT` | `dict` | Precomputed psi angles for EV model |
| `EV_RESIDUE_MAPPER` | `dict` | Maps residues to EV parameters |

**Sources:** [soursop/ssdata.py:1-450]()

### Numerical Utilities (sstools)

Helper functions for numerical operations and file discovery.

| Function | Purpose |
|----------|---------|
| `find_trajectory_files()` | Discover trajectory/topology pairs in directories |
| `apply_PBC_correction()` | Periodic boundary condition correction |
| `angle_mean()` | Circular mean for angular data |
| `angle_std()` | Circular standard deviation |

**Sources:** [soursop/sstools.py:1-300]()

### Thread Control (ssutils)

Functions for parallel execution and keyword validation.

| Function | Purpose |
|----------|---------|
| `validate_keyword_option()` | Validate keyword arguments against allowed values |
| `check_thread_count()` | Ensure thread count is reasonable |
| `get_cores()` | Determine available CPU cores |

Used extensively in `SamplingQuality.__validate_arguments()` [soursop/sssampling.py:235-253]() to ensure valid parameters.

**Sources:** [soursop/ssutils.py:1-150]()

### Error Handling

| Module | Components | Purpose |
|--------|-----------|---------|
| `ssexceptions` | `SSException` class | Custom exception type |
| `ssio` | `warning_message()`, `exception_message()`, `debug_message()` | Formatted user messages |

**Sources:** [soursop/ssexceptions.py:1-50](), [soursop/ssio.py:1-100]()

## Module Dependency Architecture

```mermaid
graph TD
    subgraph External["External Dependencies"]
        mdtraj["mdtraj 1.9.5-1.9.7"]
        numpy["numpy >= 1.20.0"]
        scipy["scipy >= 1.5.0"]
        pandas["pandas >= 0.23.0"]
    end
    
    subgraph Core["Core Analysis Classes"]
        SSTrajectory["sstrajectory.py<br/>SSTrajectory<br/>parallel_load_trjs"]
        SSProtein["ssprotein.py<br/>SSProtein<br/>2000+ lines"]
        SamplingQuality["sssampling.py<br/>SamplingQuality<br/>hellinger_distance<br/>rel_entropy"]
    end
    
    subgraph Specialized["Specialized Modules"]
        SSNMR["ssnmr.py<br/>get_offset_predicted_shifts"]
        SSPRE["sspre.py<br/>SSPRE class"]
        SSMI["ssmutualinformation.py<br/>get_mutual_information"]
        SSPoly["sspolymer.py<br/>get_overlap_concentration"]
    end
    
    subgraph Utilities["Utilities"]
        ssdata["ssdata.py<br/>ALL_VALID_RESIDUE_NAMES<br/>PHI_EV_ANGLES_DICT<br/>PSI_EV_ANGLES_DICT"]
        sstools["sstools.py<br/>find_trajectory_files<br/>apply_PBC_correction"]
        ssutils["ssutils.py<br/>validate_keyword_option<br/>check_thread_count"]
        ssio["ssio.py<br/>warning_message<br/>exception_message"]
        ssexceptions["ssexceptions.py<br/>SSException"]
        internal_data["_internal_data.py<br/>BBSEG2<br/>MAX_SASA"]
    end
    
    mdtraj --> SSTrajectory
    numpy --> SSTrajectory
    numpy --> SSProtein
    numpy --> SamplingQuality
    scipy --> SSProtein
    pandas --> SSProtein
    
    SSTrajectory --> SSProtein
    SSTrajectory --> SamplingQuality
    SSProtein --> SamplingQuality
    
    SSProtein --> SSNMR
    SSProtein --> SSPRE
    SSProtein --> SSMI
    SSProtein --> SSPoly
    
    ssdata --> SSProtein
    ssdata --> SamplingQuality
    sstools --> SSProtein
    sstools --> SamplingQuality
    ssutils --> SSProtein
    ssutils --> SamplingQuality
    ssio --> SSProtein
    ssio --> SamplingQuality
    ssexceptions --> SSProtein
    ssexceptions --> SamplingQuality
    ssexceptions --> SSNMR
    ssexceptions --> SSPRE
    internal_data --> SSProtein
```

**Dependency Notes:**

1. **External Layer**: Core scientific Python libraries (mdtraj, numpy, scipy, pandas) provide trajectory I/O and numerical computation.

2. **Core Layer**: The three main classes have minimal inter-dependencies. `SSTrajectory` creates `SSProtein` objects. `SamplingQuality` consumes both but doesn't inherit from either.

3. **Utilities Layer**: All core and specialized modules depend on utilities for validation (`ssutils`), data lookups (`ssdata`), and error handling (`ssexceptions`, `ssio`).

4. **Lazy Loading**: `SSTrajectory` uses the `@lazy_loading_single_protein_trajectory` decorator [soursop/sstrajectory.py:34-58]() to defer full protein trajectory construction until needed by methods like `get_overall_radius_of_gyration()` [soursop/sstrajectory.py:665-689]().

**Sources:** [soursop/sstrajectory.py:1-30](), [soursop/ssprotein.py:1-30](), [soursop/sssampling.py:1-34]()

## Initialization and Object Lifecycle

### SSTrajectory Initialization Sequence

```mermaid
sequenceDiagram
    participant User
    participant SSTrajectory
    participant mdtraj
    participant SSProtein
    participant ssdata
    
    User->>SSTrajectory: __init__(traj_file, pdb_file)
    SSTrajectory->>ssdata: Load ALL_VALID_RESIDUE_NAMES
    SSTrajectory->>mdtraj: md.load(traj_file, top=pdb_file)
    mdtraj-->>SSTrajectory: mdtraj.Trajectory object
    SSTrajectory->>SSTrajectory: __get_proteins(traj, debug)
    loop For each chain
        SSTrajectory->>ssdata: Check if residue in ALL_VALID_RESIDUE_NAMES
        ssdata-->>SSTrajectory: Valid/Invalid
        SSTrajectory->>mdtraj: trajectory.atom_slice(chain_atoms)
        mdtraj-->>SSTrajectory: Chain sub-trajectory
        SSTrajectory->>SSProtein: __init__(sub_trajectory)
        SSProtein->>SSProtein: __get_resid_with_CA()
        SSProtein->>SSProtein: Initialize cache dictionaries
        SSProtein-->>SSTrajectory: SSProtein instance
    end
    SSTrajectory-->>User: SSTrajectory with proteinTrajectoryList
```

**Key Steps:**

1. **File Loading** [soursop/sstrajectory.py:283-364](): `mdtraj.load()` creates the base trajectory object.

2. **Residue Validation** [soursop/sstrajectory.py:481-493](): Each chain's first residue is checked against `ALL_VALID_RESIDUE_NAMES` (or all residues if `explicit_residue_checking=True`).

3. **Atom Slicing** [soursop/sstrajectory.py:534-544](): Valid chains are extracted via `trajectory.atom_slice()`.

4. **SSProtein Creation** [soursop/ssprotein.py:59-171](): Each chain becomes an `SSProtein`. The expensive `__get_resid_with_CA()` operation identifies CA atoms.

5. **Storage**: `SSProtein` objects are stored in `proteinTrajectoryList` [soursop/sstrajectory.py:553]().

**Sources:** [soursop/sstrajectory.py:67-243](), [soursop/sstrajectory.py:439-559](), [soursop/ssprotein.py:59-171]()

### SamplingQuality Initialization Sequence

```mermaid
sequenceDiagram
    participant User
    participant SamplingQuality
    participant parallel_load_trjs
    participant SSTrajectory
    participant PrecomputedDihedralInterface
    participant ssdata
    
    User->>SamplingQuality: __init__(traj_list, reference_list, top_file)
    SamplingQuality->>SamplingQuality: __validate_arguments()
    alt reference_list provided
        SamplingQuality->>parallel_load_trjs: Load traj_list
        parallel_load_trjs->>SSTrajectory: Create SSTrajectory objects in parallel
        SSTrajectory-->>parallel_load_trjs: List of SSTrajectory
        parallel_load_trjs-->>SamplingQuality: self.trajs
        SamplingQuality->>parallel_load_trjs: Load reference_list
        parallel_load_trjs-->>SamplingQuality: self.ref_trajs
        SamplingQuality->>SamplingQuality: __compute_dihedrals(precomputed=False)
    else No reference_list
        SamplingQuality->>parallel_load_trjs: Load traj_list
        parallel_load_trjs-->>SamplingQuality: self.trajs
        SamplingQuality->>SamplingQuality: __compute_dihedrals(precomputed=True)
        SamplingQuality->>SamplingQuality: Extract sequence from trajs[0]
        SamplingQuality->>PrecomputedDihedralInterface: __init__(sequence, bins)
        PrecomputedDihedralInterface->>ssdata: Load PHI_EV_ANGLES_DICT, PSI_EV_ANGLES_DICT
        ssdata-->>PrecomputedDihedralInterface: EV angle data
        PrecomputedDihedralInterface->>PrecomputedDihedralInterface: sample_angles()
        PrecomputedDihedralInterface-->>SamplingQuality: ref_phi_angles, ref_psi_angles
    end
    SamplingQuality-->>User: Initialized SamplingQuality
```

**Sources:** [soursop/sssampling.py:106-234](), [soursop/sssampling.py:1203-1339](), [soursop/sstrajectory.py:1382-1413]()

## Extension Mechanisms

### Plugin System

SOURSOP supports custom analysis plugins in `soursop/plugins/` directory. Plugins can access core `SSProtein` and `SSTrajectory` objects.

**Example Structure:**
```
soursop/
└── plugins/
    ├── __init__.py
    └── sparrow_plugin.py
```

Plugins can define standalone functions or classes that operate on `SSProtein` objects. See [Plugin Extension System](#8) for details.

**Sources:** [soursop/plugins/sparrow_plugin.py:1-100]()

### Custom SSProtein Methods

Users can extend `SSProtein` by accessing the underlying `traj` and `topology` attributes, which are standard `mdtraj` objects. Custom analyses can leverage the caching infrastructure by using `__residue_atom_lookup()` patterns.

**Sources:** [soursop/ssprotein.py:684-749]()

## Performance Considerations

### Caching Strategy

`SSProtein` caches expensive operations in private dictionaries:

- `__CA_residue_atom`: CA atom indices per residue
- `__residue_atom_table`: General atom lookups per residue
- `__residue_COM`: Center-of-mass calculations
- `__SASA_saved`: Solvent-accessible surface area
- `__all_angles`: Dihedral angle calculations

Cache can be manually cleared with `reset_cache()` [soursop/ssprotein.py:176-218]().

**Sources:** [soursop/ssprotein.py:136-147]()

### Coarse-Grained Optimization

For one-bead-per-residue systems, `__check_cg_onebead()` [soursop/ssprotein.py:531-547]() detects when `n_atoms == n_residues`. This triggers `__initialize_cg_atoms()` [soursop/ssprotein.py:665-679]() which pre-populates `__residue_atom_table` for all residues in a single pass, achieving ~30x initialization speedup.

**Sources:** [soursop/ssprotein.py:149-153](), [soursop/ssprotein.py:531-679]()

### Parallel Loading

`parallel_load_trjs()` [soursop/sstrajectory.py:1382-1413]() uses Python's `multiprocessing.Pool` to load multiple trajectories concurrently. The `n_procs` parameter controls parallelization (defaults to `os.cpu_count()`).

**Sources:** [soursop/sstrajectory.py:1382-1413]()

---