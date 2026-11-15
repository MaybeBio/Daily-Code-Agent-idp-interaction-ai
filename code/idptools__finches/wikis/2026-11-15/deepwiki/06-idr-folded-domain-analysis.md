# IDR-Folded Domain Analysis

<details>
<summary>Relevant source files</summary>

The following files were used as context for generating this wiki page:

- [Makefile](Makefile)
- [docs/idr_fd.rst](docs/idr_fd.rst)
- [docs/requirements.txt](docs/requirements.txt)
- [finches/PDB_structure_tools.py](finches/PDB_structure_tools.py)
- [finches/data/reference_sequence_info.py](finches/data/reference_sequence_info.py)
- [finches/multicomponent_evolution_at_sequence.py](finches/multicomponent_evolution_at_sequence.py)
- [finches/utils/folded_domain_utils.py](finches/utils/folded_domain_utils.py)

</details>



This page covers FINCHES capabilities for analyzing interactions between intrinsically disordered regions (IDRs) and folded protein domains using PDB structure files. This analysis enables quantification of binding specificity between IDRs and specific surface patches on folded domains.

For general epsilon calculations between sequences, see [Epsilon Calculations](#3.2). For generating 2D interaction maps between IDR sequences, see [Interaction Maps](#3.3).

## Overview

IDR-folded domain analysis in FINCHES involves:

1. **Structure Processing**: Loading PDB files and identifying surface-exposed residues using SASA calculations
2. **Surface Characterization**: Building graphs of surface residues and defining local patches
3. **Epsilon Calculations**: Computing interaction strengths between IDR sequences and surface patches
4. **Interaction Mapping**: Generating 2D matrices showing IDR position vs. folded domain position interactions

The core workflow uses the `FoldedDomain` class to process PDB structures and integrate with FINCHES epsilon calculation engines.

**IDR-Folded Domain Analysis Workflow**

```mermaid
flowchart TD
    PDB["PDB Structure File"] --> FD["FoldedDomain Class"]
    FD --> SASA["SASA Calculation<br/>identify_surface_residues()"]
    SASA --> GRAPH["Surface Graph Construction<br/>get_nearest_neighbour_res()"]
    GRAPH --> PATCHES["Surface Patch Definition<br/>patch communities"]
    
    IDR_SEQ["IDR Sequence"] --> WINDOW["Sliding Window<br/>n-residue fragments"]
    
    PATCHES --> EPSILON["Epsilon Calculations<br/>calculate_surface_epsilon()"]
    WINDOW --> EPSILON
    
    EPSILON --> INTERMAP["Interaction Matrix<br/>calculate_idr_surface_patch_interactions()"]
    EPSILON --> MEAN["Mean Surface Epsilon<br/>calculate_mean_surface_epsilon()"]
    EPSILON --> ATTRACTIVE["Attractive/Repulsive<br/>Regions Analysis"]
```

Sources: [finches/utils/folded_domain_utils.py:66-984](), [docs/idr_fd.rst:1-87]()

## FoldedDomain Class Architecture

The `FoldedDomain` class serves as the primary interface for structure-based analysis. It integrates MDTraj for structure processing, NetworkX for graph analysis, and SciPy for distance calculations.

**FoldedDomain Class Structure**

```mermaid
classDiagram
    class FoldedDomain {
        +str sequence
        +mdtraj.Trajectory traj
        +numpy.array sasa
        +list surface_indices
        +numpy.array surface_vector
        +dict surface_positions
        +networkx.Graph surface_graph
        +dict surface_neighbours
        +dict surface_distance_surface
        +dict surface_distance_straight_line
        
        +__init__(pdbfilename, start, end, probe_radius, surface_thresh, sasa_mode)
        +get_nearest_neighbour_res(distance_thresh)
        +calculate_surface_epsilon(input_sequence, IMCObject)
        +calculate_mean_surface_epsilon(input_sequence, IMCObject)
        +calculate_idr_surface_patch_interactions(interacting_sequence, IMCObject, idr_tile_size, patch_radius)
        +calculate_IWD(residues, positions, distance_mode)
        +write_SASA_vis_file(filename)
        +write_epsilon_vis_file(surface_epsilon, filename)
    }
    
    class MDTraj {
        +load_pdb()
        +shrake_rupley()
        +compute_center_of_mass()
    }
    
    class NetworkX {
        +Graph()
        +all_pairs_dijkstra_path()
        +all_pairs_dijkstra_path_length()
    }
    
    class IMCObject {
        +calculate_epsilon_value()
    }
    
    FoldedDomain --> MDTraj : "structure processing"
    FoldedDomain --> NetworkX : "graph analysis"
    FoldedDomain --> IMCObject : "epsilon calculations"
```

Sources: [finches/utils/folded_domain_utils.py:66-248](), [finches/utils/folded_domain_utils.py:256-306]()

## Initialization and Surface Identification

### Basic Initialization

The `FoldedDomain` class initializes with a PDB file and automatically performs SASA calculations to identify surface residues:

```python
from finches.utils import folded_domain_utils

# Initialize with PDB file
fd = folded_domain_utils.FoldedDomain('protein.pdb')

# Access computed properties
print(f"Sequence length: {len(fd.sequence)}")
print(f"Surface fraction: {fd.surface_fraction}")
print(f"Surface residue indices: {fd.surface_indices}")
```

### SASA-Based Surface Detection

Surface residues are identified using the Shrake-Rupley algorithm with configurable parameters:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `probe_radius` | 1.4 Å | Probe radius for SASA calculation |
| `surface_thresh` | 0.10 | Threshold fraction of max SASA for surface classification |
| `sasa_mode` | 'v1' | SASA comparison mode ('v1' or 'v2') |

The `sasa_mode` parameter controls surface residue identification:
- **v1**: Compares residue SASA to sidechain-only max SASA
- **v2**: Compares residue SASA to combined sidechain + backbone max SASA (more stringent)

Sources: [finches/utils/folded_domain_utils.py:70-196](), [finches/utils/folded_domain_utils.py:17-36]()

## Surface Graph Construction and Neighbor Analysis

### Neighbor Detection

The `get_nearest_neighbour_res()` method builds a spatial graph of surface residues within a distance threshold:

```python
# Build surface graph with 9.0 Å distance threshold
fd.get_nearest_neighbour_res(distance_thresh=9.0)

# Access neighbor information
for residue_idx in fd.surface_indices:
    neighbors = fd.surface_neighbours[residue_idx]
    print(f"Residue {residue_idx}: {len(neighbors)} neighbors")
```

### Distance Calculations

The class maintains two distance metrics between surface residues:

- **`surface_distance_straight_line`**: Euclidean distances between residue positions
- **`surface_distance_surface`**: Graph-based distances following surface connectivity

**Surface Graph and Distance Analysis**

```mermaid
graph TB
    subgraph "Surface Residue Detection"
        PDB_LOAD["md.load_pdb()"] --> SASA_CALC["md.shrake_rupley()"]
        SASA_CALC --> SURFACE_ID["surface_vector identification"]
        SURFACE_ID --> COORDS["surface_positions extraction"]
    end
    
    subgraph "Graph Construction"
        COORDS --> DIST_MATRIX["scipy.spatial.distance.cdist()"]
        DIST_MATRIX --> NEIGHBORS["neighbor identification<br/>distance_thresh=9.0Å"]
        NEIGHBORS --> NX_GRAPH["nx.Graph() creation"]
        NX_GRAPH --> SHORTEST_PATHS["nx.all_pairs_dijkstra_path()"]
    end
    
    subgraph "Distance Metrics"
        SHORTEST_PATHS --> SURFACE_DIST["surface_distance_surface<br/>graph-based distances"]
        DIST_MATRIX --> EUCLIDEAN_DIST["surface_distance_straight_line<br/>euclidean distances"]
    end
```

Sources: [finches/utils/folded_domain_utils.py:324-440](), [finches/utils/folded_domain_utils.py:360-431]()

## Epsilon Calculations

### Surface Epsilon Analysis

The `calculate_surface_epsilon()` method computes interaction strengths between an input IDR sequence and each surface residue's local environment:

```python
# Calculate surface epsilon values
surface_eps = fd.calculate_surface_epsilon(idr_sequence, IMC_object)

# surface_eps is a dictionary mapping residue_index -> [center_residue, neighbor_sequence, epsilon_value]
for idx, (center_res, neighbor_seq, eps_val) in surface_eps.items():
    print(f"Residue {idx} ({center_res}): ε = {eps_val:.3f}")
```

### Local Environment Reordering

The method implements chemical environment-aware reordering of neighbor residues:

| Center Residue Type | Reordering Strategy |
|-------------------|-------------------|
| Hydrophobic (I,V,L,A,M) | Place hydrophobic neighbors first |
| Negative (D,E) | Place negative neighbors first |
| Positive (K,R) | Place positive neighbors first |
| Other | No reordering |

This ensures that epsilon calculations reflect the local chemical context around each surface residue.

### Aggregate Metrics

```python
# Mean surface epsilon
mean_eps = fd.calculate_mean_surface_epsilon(idr_sequence, IMC_object)

# Attractive surface residues (ε < 0)
attractive_eps = fd.calculate_attractive_surface_epsilon(idr_sequence, IMC_object)

# Repulsive surface residues (ε > 0)  
repulsive_eps = fd.calculate_repulsive_surface_epsilon(idr_sequence, IMC_object)
```

Sources: [finches/utils/folded_domain_utils.py:486-584](), [finches/utils/folded_domain_utils.py:543-575]()

## IDR-Surface Interaction Mapping

### Patch-Based Interaction Analysis

The `calculate_idr_surface_patch_interactions()` method generates comprehensive interaction matrices between IDR sequences and folded domain surface patches:

```python
# Generate IDR-surface interaction matrix
interaction_data = fd.calculate_idr_surface_patch_interactions(
    interacting_sequence=idr_sequence,
    IMCObject=IMC_object,
    idr_tile_size=31,      # Sliding window size for IDR
    patch_radius=12        # Surface patch radius in Å
)

# Returns: [interaction_dict, epsilon_matrix, mean_vector]
interaction_dict, eps_matrix, mean_vector = interaction_data
```

### Interaction Matrix Structure

The output epsilon matrix has dimensions `(folded_domain_length, idr_length - window_size + 1)`:
- **Rows**: Folded domain residue positions (surface residues have values, buried residues are NaN)
- **Columns**: IDR sliding window positions
- **Values**: Epsilon interactions between IDR fragments and surface patches

**IDR-Surface Interaction Mapping Pipeline**

```mermaid
flowchart LR
    subgraph "IDR Processing"
        IDR_SEQ["IDR Sequence<br/>length N"] --> SLIDING["Sliding Window<br/>size = idr_tile_size"]
        SLIDING --> IDR_FRAGS["N - window_size + 1<br/>IDR fragments"]
    end
    
    subgraph "Surface Processing"
        SURFACE_RES["Surface Residues"] --> PATCH_DEF["Patch Definition<br/>patch_radius = 12Å"]
        PATCH_DEF --> PATCHES["Surface Patches<br/>neighbor communities"]
    end
    
    subgraph "Epsilon Calculation"
        IDR_FRAGS --> EPS_CALC["calculate_epsilon_value()"]
        PATCHES --> EPS_CALC
        EPS_CALC --> MATRIX["Interaction Matrix<br/>(FD_length × IDR_fragments)"]
    end
    
    subgraph "Output Products"
        MATRIX --> INTERMAP["2D Interaction Map"]
        MATRIX --> MEAN_VEC["Mean Interaction Vector"]
        MATRIX --> PATCH_DICT["Per-patch Interaction Dict"]
    end
```

Sources: [finches/utils/folded_domain_utils.py:891-984](), [docs/idr_fd.rst:8-48]()

## Advanced Analysis Methods

### Inverse Weighted Distance (IWD)

The `calculate_IWD()` method quantifies spatial clustering of specific residue types on the protein surface:

```python
# Calculate IWD for positively charged residues
iwd_result = fd.calculate_IWD(
    residues=['K', 'R'],           # Target residue types
    distance_mode='surface',       # 'surface' or 'straight_line'
    calculate_null=True,           # Generate null distribution
    number_null_iterations=100
)

iwd_value, num_residues, null_distribution = iwd_result
```

### Visualization Output

The class provides methods for generating visualization files compatible with molecular graphics software:

```python
# Binary surface/buried visualization
fd.write_SASA_vis_file('surface_residues.txt')

# Epsilon value coloring
surface_eps = fd.calculate_surface_epsilon(idr_sequence, IMC_object)
fd.write_epsilon_vis_file(surface_eps, 'epsilon_values.txt')
```

### Integration with Computational Engines

The analysis integrates with FINCHES frontend objects through the `IMCObject`:

| Frontend Class | IMC Object Access | Usage |
|---------------|------------------|-------|
| `Mpipi_frontend` | `mf.IMC_object` | Mpipi forcefield calculations |
| `CALVADOS_frontend` | `cf.IMC_object` | CALVADOS forcefield calculations |

Sources: [finches/utils/folded_domain_utils.py:703-821](), [finches/utils/folded_domain_utils.py:827-889]()

## Data Structures and Properties

### Key Object Properties

| Property | Type | Description |
|----------|------|-------------|
| `sequence` | str | Full amino acid sequence |
| `sasa` | numpy.array | Per-residue SASA values |
| `surface_vector` | list | Binary surface classification (0/1) |
| `surface_indices` | list | Indices of surface residues |
| `surface_positions` | dict | 3D coordinates of surface residues |
| `surface_graph` | networkx.Graph | Graph connecting surface residues |
| `surface_neighbours` | dict | Neighbor lists with distances |

### Lazy Property Evaluation

Several computationally expensive properties use lazy evaluation through the `@property` decorator:

- `all_shortest_paths`: Computed on first access via NetworkX
- `surface_graph`: Built during first neighbor analysis
- `surface_distance_surface`: Generated with graph construction
- `surface_distance_straight_line`: Calculated with distance matrix

This design pattern ensures efficient memory usage and computation only when needed.

Sources: [finches/utils/folded_domain_utils.py:86-102](), [finches/utils/folded_domain_utils.py:256-306]()

---