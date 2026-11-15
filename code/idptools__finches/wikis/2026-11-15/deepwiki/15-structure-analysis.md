# Structure Analysis

<details>
<summary>Relevant source files</summary>

The following files were used as context for generating this wiki page:

- [Makefile](Makefile)
- [docs/requirements.txt](docs/requirements.txt)
- [finches/PDB_structure_tools.py](finches/PDB_structure_tools.py)
- [finches/data/reference_sequence_info.py](finches/data/reference_sequence_info.py)
- [finches/multicomponent_evolution_at_sequence.py](finches/multicomponent_evolution_at_sequence.py)
- [finches/utils/folded_domain_utils.py](finches/utils/folded_domain_utils.py)

</details>



This page provides API documentation for FINCHES structure analysis components, which enable the analysis of interactions between intrinsically disordered regions (IDRs) and folded protein domains using PDB structure files. The structure analysis system processes PDB files to identify surface-accessible residues, calculate solvent accessible surface areas (SASA), and compute epsilon values for IDR-folded domain interactions.

For general usage guidance on IDR-folded domain analysis, see [IDR-Folded Domain Analysis](#3.4). For core epsilon calculation concepts, see [Matrix Calculations](#4.3).

## System Architecture

The structure analysis system consists of two main components: the `FoldedDomain` class for comprehensive structural analysis and a collection of utility functions for PDB processing.

### Structure Analysis Component Overview

```mermaid
graph TB
    subgraph "Main API Classes"
        FD["FoldedDomain<br/>finches/utils/folded_domain_utils.py"]
    end
    
    subgraph "Utility Functions"
        PST["PDB_structure_tools<br/>finches/PDB_structure_tools.py"]
        SASA["SASA Calculation<br/>_accessibility_parse()"]
        EXTRACT["Domain Extraction<br/>pdb_to_SDFDresidues_and_xyzs()"]
    end
    
    subgraph "External Dependencies" 
        MDTRAJ["MDTraj<br/>PDB parsing & SASA"]
        SOURSOP["Soursop<br/>Alternative PDB parser"]
        METAPREDICT["MetaPredict<br/>Domain boundary prediction"]
    end
    
    subgraph "Data Processing"
        SURFACE["Surface Detection<br/>probe_radius, surface_thresh"]
        GRAPH["Neighborhood Graph<br/>NetworkX graph construction"]
        COORDS["Coordinate Extraction<br/>Center of mass calculation"]
    end
    
    subgraph "Analysis Outputs"
        EPS["Surface Epsilon Values<br/>calculate_surface_epsilon()"]
        IWD["Inverse Weighted Distance<br/>calculate_IWD()"]
        VIS["Visualization Files<br/>write_*_vis_file()"]
    end
    
    FD --> SURFACE
    FD --> GRAPH  
    FD --> COORDS
    FD --> EPS
    FD --> IWD
    FD --> VIS
    
    PST --> SASA
    PST --> EXTRACT
    
    FD --> MDTRAJ
    PST --> SOURSOP
    PST --> METAPREDICT
    
    SURFACE --> EPS
    GRAPH --> IWD
    COORDS --> GRAPH
```

Sources: [finches/utils/folded_domain_utils.py:1-1000](), [finches/PDB_structure_tools.py:1-630]()

## FoldedDomain Class

The `FoldedDomain` class is the primary interface for folded domain structural analysis, providing comprehensive functionality for surface residue identification, epsilon calculations, and structural characterization.

### Class Initialization and Core Properties

```mermaid
graph LR
    subgraph "FoldedDomain.__init__()"
        INIT["FoldedDomain.__init__()<br/>pdbfilename, start, end<br/>probe_radius, surface_thresh"]
    end
    
    subgraph "Core Properties"
        SEQ["self.sequence<br/>Full AA sequence"]
        TRAJ["self.traj<br/>MDTraj trajectory"]
        SASA_PROP["self.sasa<br/>Per-residue SASA values"]
        SURF_VEC["self.surface_vector<br/>Binary surface indicators"]
        SURF_IDX["self.surface_indices<br/>Surface residue indices"]
        SURF_POS["self.surface_positions<br/>3D coordinates dict"]
    end
    
    subgraph "Lazy Properties"
        GRAPH["self.surface_graph<br/>@property NetworkX graph"]
        NEIGHBORS["self.surface_neighbours<br/>@property Distance-based neighbors"]
        PATHS["self.all_shortest_paths<br/>@property Graph shortest paths"]
    end
    
    INIT --> SEQ
    INIT --> TRAJ
    INIT --> SASA_PROP
    INIT --> SURF_VEC
    INIT --> SURF_IDX
    INIT --> SURF_POS
    
    GRAPH --> NEIGHBORS
    NEIGHBORS --> PATHS
```

The `FoldedDomain` constructor accepts several key parameters:

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `pdbfilename` | str | Required | Path to PDB file |
| `start` | int | None | Start residue index (0-based) |
| `end` | int | None | End residue index (0-based) |
| `probe_radius` | float | 1.4 | Probe radius for SASA calculation (Angstroms) |
| `surface_thresh` | float | 0.10 | Threshold for surface residue classification |
| `sasa_mode` | str | 'v1' | SASA comparison mode ('v1' or 'v2') |
| `residue_overide_mapping` | dict | {} | Non-standard residue name mapping |
| `ignore_warnings` | bool | False | Suppress start/end warnings |
| `SASA_ONLY` | bool | False | Calculate only SASA, skip other analysis |

Sources: [finches/utils/folded_domain_utils.py:66-196]()

### Surface Residue Detection

The surface detection algorithm uses solvent accessible surface area (SASA) calculations to identify residues exposed to solvent. The system supports two SASA comparison modes:

```mermaid
flowchart TD
    subgraph "SASA Calculation"
        PDB["PDB Structure"] --> MDTRAJ["MDTraj.shrake_rupley()"]
        MDTRAJ --> SASA_VALUES["Per-residue SASA values<br/>(Angstrom²)"]
    end
    
    subgraph "Surface Classification"
        SASA_VALUES --> MODE_CHECK{"sasa_mode"}
        MODE_CHECK -->|"v1"| V1["Compare to sidechain<br/>MAX_SASA_DATA[aa][0]"]
        MODE_CHECK -->|"v2"| V2["Compare to sidechain + backbone<br/>MAX_SASA_DATA[aa][0] + [1]"]
        
        V1 --> THRESHOLD1["SASA > surface_thresh * max_sidechain"]
        V2 --> THRESHOLD2["SASA > surface_thresh * (sidechain + backbone)"]
        
        THRESHOLD1 --> SURFACE["Surface residue classification"]
        THRESHOLD2 --> SURFACE
    end
    
    subgraph "Output Properties"
        SURFACE --> SURF_VECTOR["self.surface_vector<br/>[0,1,1,0,1,...]"]
        SURFACE --> SURF_INDICES["self.surface_indices<br/>[2,3,5,...]"]
        SURFACE --> SURF_FRACTION["self.surface_fraction<br/>0.35"]
    end
```

The `MAX_SASA_DATA` dictionary contains reference SASA values for each amino acid type, calculated from GXG tripeptide simulations with excluded volume interactions.

Sources: [finches/utils/folded_domain_utils.py:14-36](), [finches/utils/folded_domain_utils.py:197-237]()

### Neighborhood Analysis and Graph Construction

The `FoldedDomain` class builds spatial relationship graphs of surface residues for distance-based analysis:

```mermaid
graph TB
    subgraph "get_nearest_neighbour_res()"
        POSITIONS["Surface residue positions<br/>Center of mass coordinates"]
        DISTANCE_MATRIX["Euclidean distance matrix<br/>scipy.spatial.distance.cdist()"]
        THRESHOLD["Distance threshold filter<br/>default: 9.0 Angstroms"]
    end
    
    subgraph "Graph Construction"
        NEIGHBORS["Neighbor identification<br/>residues within threshold"]
        NETWORKX["NetworkX graph creation<br/>Weighted by distances"]
        DIJKSTRA["Shortest path calculation<br/>nx.all_pairs_dijkstra_path()"]
    end
    
    subgraph "Output Properties"
        SURF_GRAPH["self.surface_graph<br/>NetworkX.Graph object"]
        SURF_NEIGHBORS["self.surface_neighbours<br/>Dict[int, List[Tuple]]"]
        SURFACE_DIST["self.surface_distance_surface<br/>Over-surface distances"]
        STRAIGHT_DIST["self.surface_distance_straight_line<br/>Euclidean distances"]
    end
    
    POSITIONS --> DISTANCE_MATRIX
    DISTANCE_MATRIX --> THRESHOLD
    THRESHOLD --> NEIGHBORS
    NEIGHBORS --> NETWORKX
    NETWORKX --> DIJKSTRA
    
    NETWORKX --> SURF_GRAPH
    NEIGHBORS --> SURF_NEIGHBORS
    DIJKSTRA --> SURFACE_DIST
    DISTANCE_MATRIX --> STRAIGHT_DIST
```

Sources: [finches/utils/folded_domain_utils.py:324-440]()

### Epsilon Calculation Methods

The `FoldedDomain` class provides several methods for calculating epsilon values between surface residues and input sequences:

| Method | Return Type | Description |
|--------|-------------|-------------|
| `calculate_surface_epsilon()` | Dict[int, List] | Per-residue epsilon values with local context |
| `calculate_mean_surface_epsilon()` | float | Mean epsilon across all surface residues |
| `calculate_attractive_surface_epsilon()` | List[float] | Epsilon values below threshold (attractive) |
| `calculate_repulsive_surface_epsilon()` | List[float] | Epsilon values above threshold (repulsive) |

The surface epsilon calculation algorithm considers local chemical environment:

```mermaid
flowchart TD
    subgraph "Surface Epsilon Algorithm"
        SURFACE_RES["For each surface residue"]
        CENTER_TYPE["Identify center residue type<br/>hydrophobe, positive, negative"]
        NEIGHBORS["Get neighboring residues<br/>within distance threshold"]
        REORDER["Reorder neighbor sequence<br/>by chemical similarity"]
        EPSILON["Calculate epsilon with<br/>IMCObject.calculate_epsilon_value()"]
        NORMALIZE["Normalize by neighbor count"]
    end
    
    subgraph "Chemical Grouping"
        HYDROPHOBE["Hydrophobes: I,V,L,A,M<br/>Group with center residue"]
        POSITIVE["Positive: K,R<br/>Group with center residue"]  
        NEGATIVE["Negative: D,E<br/>Group with center residue"]
        OTHER["Other residues<br/>No special grouping"]
    end
    
    SURFACE_RES --> CENTER_TYPE
    CENTER_TYPE --> NEIGHBORS
    NEIGHBORS --> REORDER
    REORDER --> EPSILON
    EPSILON --> NORMALIZE
    
    CENTER_TYPE --> HYDROPHOBE
    CENTER_TYPE --> POSITIVE
    CENTER_TYPE --> NEGATIVE
    CENTER_TYPE --> OTHER
```

Sources: [finches/utils/folded_domain_utils.py:486-696]()

### Structural Analysis Methods

The class provides specialized structural analysis methods:

#### Inverse Weighted Distance (IWD) Calculation

The `calculate_IWD()` method computes clustering measures for specific residue types:

```python
def calculate_IWD(self, residues, positions=None, 
                  distance_mode='surface', calculate_null=False, 
                  number_null_iterations=100)
```

**Parameters:**
- `residues`: List of amino acid types to analyze
- `positions`: Specific residue positions (alternative to residue types)
- `distance_mode`: 'surface' or 'straight_line' distance calculation
- `calculate_null`: Whether to compute null distribution
- `number_null_iterations`: Number of random samplings for null

**Returns:** `[iwd_value, residue_count, null_distribution]`

Sources: [finches/utils/folded_domain_utils.py:703-821]()

#### Patch-Based IDR Interaction Analysis

The `calculate_idr_surface_patch_interactions()` method performs sliding window analysis:

```python
def calculate_idr_surface_patch_interactions(self, interacting_sequence, 
                                           IMCObject, idr_tile_size=31, 
                                           patch_radius=12)
```

This method:
1. Defines surface patches around each surface residue
2. Slides a window along the IDR sequence
3. Calculates epsilon values between patches and IDR windows
4. Returns interaction matrices for visualization

Sources: [finches/utils/folded_domain_utils.py:891-985]()

### Visualization Output Methods

The class provides methods to generate files for molecular visualization software:

| Method | Output Format | Description |
|--------|---------------|-------------|
| `write_SASA_vis_file()` | Text | Binary surface/buried assignment |
| `write_epsilon_vis_file()` | Text | Per-residue epsilon values |

The output format follows the pattern: `<residue_index> A <value>`

Sources: [finches/utils/folded_domain_utils.py:827-889]()

## PDB Structure Tools

The `PDB_structure_tools` module provides utility functions for PDB file processing and domain analysis using the soursop library.

### Core PDB Processing Functions

```mermaid
graph TB
    subgraph "PDB Processing Pipeline"
        PDB_FILE["PDB File Input"]
        SOURSOP["SSTrajectory(pdb)<br/>Soursop parsing"]
        SEQUENCE["get_amino_acid_sequence()<br/>Extract sequence"]
        SASA_CALC["get_all_SASA()<br/>Calculate accessibility"]
        COORDS["get_multiple_CA_index()<br/>Extract coordinates"]
    end
    
    subgraph "Analysis Functions"
        EXTRACT_FUNC["pdb_to_SDFDresidues_and_xyzs()<br/>Main extraction function"]
        ACCESS_FUNC["_get_accessible_indicies()<br/>Threshold-based filtering"]
        FILTER_FUNC["_filter_indicies()<br/>Domain boundary filtering"]
        TRACK_FUNC["pdb_to_xyz_cordinate_track()<br/>Full coordinate extraction"]
    end
    
    subgraph "Domain Analysis"
        META["extract_flanking_domain_combinations()<br/>MetaPredict integration"]
        DISORDER["metapredict.predict_disorder_domains()"]
        BOUNDARIES["IDR/FD boundary detection"]
    end
    
    PDB_FILE --> SOURSOP
    SOURSOP --> SEQUENCE
    SOURSOP --> SASA_CALC
    SOURSOP --> COORDS
    
    SASA_CALC --> ACCESS_FUNC
    ACCESS_FUNC --> FILTER_FUNC
    FILTER_FUNC --> EXTRACT_FUNC
    COORDS --> EXTRACT_FUNC
    
    SEQUENCE --> META
    META --> DISORDER
    DISORDER --> BOUNDARIES
```

### Key Function Specifications

#### `pdb_to_SDFDresidues_and_xyzs()`

Primary function for extracting surface-accessible folded domain residues:

```python
def pdb_to_SDFDresidues_and_xyzs(pdb, FD_start, FD_end, 
                                 issolate_domain=False)
```

**Parameters:**
- `pdb`: Path to PDB file
- `FD_start`: Domain start residue index  
- `FD_end`: Domain end residue index
- `issolate_domain`: Whether to isolate domain for SASA calculation

**Returns:**
- `SAFD_seq`: Surface-accessible folded domain sequence
- `SAFD_idxs`: Indices of surface residues
- `SAFD_cords`: 3D coordinates of surface residues

Sources: [finches/PDB_structure_tools.py:249-340]()

#### `extract_flanking_domain_combinations()`

Automated domain boundary detection using MetaPredict:

```python
def extract_flanking_domain_combinations(pdb, return_domain_lists=False)
```

This function:
1. Parses PDB sequence with soursop
2. Predicts disorder with MetaPredict
3. Identifies adjacent IDR-folded domain pairs
4. Returns boundary combinations for analysis

Sources: [finches/PDB_structure_tools.py:572-628]()

### Utility Functions

| Function | Purpose | Key Parameters |
|----------|---------|----------------|
| `_accessibility_parse()` | Extract SASA from soursop | `PO`: Protein object |
| `_get_accessible_indicies()` | Filter by accessibility threshold | `threshold=10` |
| `calculate_distance()` | 3D Euclidean distance | `coord1`, `coord2` |
| `map_SAFD_vector_to_full_folded_domain()` | Map partial to full vectors | `partial_vector`, `SAFD_idxs` |

Sources: [finches/PDB_structure_tools.py:22-568]()

## Integration with Computation Engine

The structure analysis components integrate with FINCHES computation engine through the `IMCObject` parameter:

```mermaid
sequenceDiagram
    participant User
    participant FoldedDomain
    participant IMCObject
    participant EpsilonCalc
    
    User->>FoldedDomain: calculate_surface_epsilon(sequence, IMCObject)
    FoldedDomain->>FoldedDomain: Identify surface residues
    FoldedDomain->>FoldedDomain: Build neighbor sequences
    loop For each surface residue
        FoldedDomain->>IMCObject: calculate_epsilon_value(neighbor_seq, input_seq)
        IMCObject->>EpsilonCalc: Core epsilon calculation
        EpsilonCalc-->>IMCObject: Epsilon value
        IMCObject-->>FoldedDomain: Epsilon value
    end
    FoldedDomain-->>User: Surface epsilon dictionary
```

The `IMCObject` is typically obtained from frontend classes:

```python
# From Mpipi_frontend or CALVADOS_frontend
imc_object = frontend.IMC_object
surface_eps = folded_domain.calculate_surface_epsilon(sequence, imc_object)
```

Sources: [finches/utils/folded_domain_utils.py:524-528](), [finches/utils/folded_domain_utils.py:636-638]()

## Common Usage Patterns

### Basic Folded Domain Analysis

```python
# Initialize folded domain
fd = FoldedDomain('protein.pdb', probe_radius=1.4, surface_thresh=0.10)

# Analyze surface properties  
surface_fraction = fd.surface_fraction
surface_residues = fd.surface_indices

# Calculate epsilon with IDR sequence
surface_eps = fd.calculate_surface_epsilon('DEDEDEDE', imc_object)
mean_eps = fd.calculate_mean_surface_epsilon('DEDEDEDE', imc_object)
```

### Domain Boundary Detection

```python
# Automatic domain detection
combinations = extract_flanking_domain_combinations('protein.pdb')
for idr_bounds, fd_bounds in combinations:
    idr_start, idr_end = idr_bounds
    fd_start, fd_end = fd_bounds
    # Process each IDR-FD pair
```

### Patch-Based Analysis

```python
# Sliding window analysis
interaction_dict, matrix, mean_vector = fd.calculate_idr_surface_patch_interactions(
    'KEKEKEKE' * 10, imc_object, idr_tile_size=31, patch_radius=12
)
```

Sources: [finches/utils/folded_domain_utils.py:66-196](), [finches/PDB_structure_tools.py:572-628]()

---