# Angles and Dynamics

<details>
<summary>Relevant source files</summary>

The following files were used as context for generating this wiki page:

- [soursop/data/test_data/gs6_distance_map_mean.npy](soursop/data/test_data/gs6_distance_map_mean.npy)
- [soursop/data/test_data/gs6_distance_map_std.npy](soursop/data/test_data/gs6_distance_map_std.npy)
- [soursop/ssprotein.py](soursop/ssprotein.py)
- [soursop/tests/test_ssproteins.py](soursop/tests/test_ssproteins.py)

</details>



This page documents the angle and dynamics analysis methods in the `SSProtein` class. These methods compute various angular measurements from protein conformations including backbone and sidechain dihedral angles, angular correlations between residue pairs, and sidechain orientation analysis.

For information about secondary structure analysis using angles, see [Secondary Structure Analysis](#4.4). For dihedral angle distributions used in sampling quality assessment, see [Dihedral Analysis](#5.2).

---

## Overview of Angle Analysis Methods

The `SSProtein` class provides four main categories of angle-related analysis:

| Method | Purpose | Returns |
|--------|---------|---------|
| `get_angles()` | Extract backbone (φ, ψ, ω) and sidechain (χ1-χ5) dihedral angles | Atom names and angle values |
| `get_angle_decay()` | Measure angular correlation between CA-CA vectors at different sequence separations | Mean angles by separation distance |
| `get_D_vector()` | Compute local heterogeneity in structural dynamics | Per-residue D-vector values |
| `get_sidechain_alignment_angle()` | Calculate angle between sidechain vectors of two residues | Alignment angles across frames |

**Diagram: Angle Analysis System Architecture**

```mermaid
graph TD
    subgraph Input["SSProtein Trajectory"]
        Traj["self.traj<br/>mdtraj.Trajectory"]
        Topo["self.topology<br/>Topology object"]
    end
    
    subgraph Cache["Angle Caching System"]
        AngleCache["self.__all_angles<br/>Dictionary cache"]
    end
    
    subgraph Methods["Angle Analysis Methods"]
        GetAngles["get_angles(angle_type)<br/>Dihedral extraction"]
        AngleDecay["get_angle_decay()<br/>CA-CA angular correlation"]
        DVector["get_D_vector()<br/>Local dynamics measure"]
        SCAlign["get_sidechain_alignment_angle()<br/>Sidechain orientation"]
    end
    
    subgraph MDTraj["MDTraj Backend"]
        ComputePhi["md.compute_phi()"]
        ComputePsi["md.compute_psi()"]
        ComputeOmega["md.compute_omega()"]
        ComputeChi["md.compute_chi1/2/3/4/5()"]
    end
    
    subgraph Output["Analysis Output"]
        AtomList["Atom name lists<br/>[4 atoms per angle]"]
        AngleValues["Angle arrays<br/>[n_frames x n_angles]"]
        DecayMatrix["Decay matrix<br/>[separation x angle]"]
        PairDict["All-pairs dictionary<br/>{i-j: angles}"]
        DVals["D-vector array<br/>[n_residues]"]
    end
    
    Traj --> GetAngles
    Topo --> GetAngles
    Traj --> AngleDecay
    Traj --> DVector
    Traj --> SCAlign
    
    GetAngles --> ComputePhi
    GetAngles --> ComputePsi
    GetAngles --> ComputeOmega
    GetAngles --> ComputeChi
    
    ComputePhi --> AngleCache
    ComputePsi --> AngleCache
    ComputeOmega --> AngleCache
    ComputeChi --> AngleCache
    
    AngleCache --> AtomList
    AngleCache --> AngleValues
    AngleDecay --> DecayMatrix
    AngleDecay --> PairDict
    DVector --> DVals
    SCAlign --> AngleValues
```

**Sources:** [soursop/ssprotein.py:1-200](), [soursop/ssprotein.py:146](), [soursop/tests/test_ssproteins.py:1019-1045]()

---

## Dihedral Angles

Dihedral angles describe the geometry of the protein backbone and sidechains. The `SSProtein` class provides access to both backbone dihedrals (φ, ψ, ω) and up to five sidechain dihedrals (χ1-χ5) per residue.

### Angle Types

**Diagram: Dihedral Angle Types**

```mermaid
graph LR
    subgraph Backbone["Backbone Dihedrals"]
        Phi["φ (phi)<br/>C-N-CA-C"]
        Psi["ψ (psi)<br/>N-CA-C-N"]
        Omega["ω (omega)<br/>CA-C-N-CA"]
    end
    
    subgraph Sidechain["Sidechain Dihedrals"]
        Chi1["χ1 (chi1)<br/>N-CA-CB-XG"]
        Chi2["χ2 (chi2)<br/>CA-CB-XG-XD"]
        Chi3["χ3 (chi3)<br/>CB-XG-XD-XE"]
        Chi4["χ4 (chi4)<br/>XG-XD-XE-XZ"]
        Chi5["χ5 (chi5)<br/>XD-XE-XZ-XH"]
    end
    
    Backbone -.->|"20 standard residues"| All["All residues except<br/>ACE/NME caps"]
    Sidechain -.->|"Residue-dependent"| Chi1
    Chi1 -.->|"If CB exists"| Chi2
    Chi2 -.->|"Long sidechains"| Chi3
    Chi3 -.->|"Arg, Lys, Met"| Chi4
    Chi4 -.->|"Arg only"| Chi5
```

**Sources:** [soursop/tests/test_ssproteins.py:1019-1045]()

### The `get_angles()` Method

The `get_angles()` method extracts dihedral angles from the trajectory using MDTraj's angle computation functions. Results are cached in `self.__all_angles` to avoid redundant computation.

**Method Signature:**
```python
get_angles(angle_type, stride=1, verbose=True)
```

**Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `angle_type` | `str` | One of: `'phi'`, `'psi'`, `'omega'`, `'chi1'`, `'chi2'`, `'chi3'`, `'chi4'`, `'chi5'` |
| `stride` | `int` | Frame spacing for analysis (default: 1) |
| `verbose` | `bool` | Print status messages (default: True) |

**Returns:**

A tuple `(atom_names, angle_values)` where:
- `atom_names`: List of lists, each containing 4 atom identifier strings defining the dihedral
- `angle_values`: NumPy array of shape `[n_frames, n_angles]` with angles in degrees

**Example Output Structure:**

For a peptide with 6 CA-containing residues, `get_angles('phi')` returns:

```python
atom_names = [
    ['ACE1-C', 'GLY2-N', 'GLY2-CA', 'GLY2-C'],
    ['GLY2-C', 'SER3-N', 'SER3-CA', 'SER3-C'],
    # ... one list per phi angle
]

angle_values = np.array([
    [-157.64, -126.71, ...],  # Frame 0
    [-155.32, -128.45, ...],  # Frame 1
    # ... one row per frame
])
```

**Diagram: Angle Caching Mechanism**

```mermaid
graph TD
    Request["get_angles(angle_type)"]
    CheckCache{"angle_type in<br/>self.__all_angles?"}
    Cached["Return cached result"]
    Compute["Compute angles"]
    
    subgraph Computation["Angle Computation"]
        SelectMD["Select MDTraj function<br/>md.compute_phi/psi/omega/chi"]
        CallMD["Call MDTraj function<br/>with self.traj"]
        Process["Process results<br/>Convert to degrees"]
        Store["Store in self.__all_angles[angle_type]"]
    end
    
    Request --> CheckCache
    CheckCache -->|"Yes"| Cached
    CheckCache -->|"No"| Compute
    Compute --> SelectMD
    SelectMD --> CallMD
    CallMD --> Process
    Process --> Store
    Store --> Cached
```

**Sources:** [soursop/ssprotein.py:146](), [soursop/tests/test_ssproteins.py:1019-1045]()

### Angle Availability by Residue Type

Not all angles are defined for all residues:

| Residue Property | Available Angles |
|-----------------|------------------|
| ACE/NME caps | None (no CA atom) |
| N-terminal residue | ψ, ω only (no φ) |
| C-terminal residue | φ, ω only (no ψ) |
| Glycine | φ, ψ, ω only (no CB, therefore no χ) |
| Standard residues with CB | φ, ψ, ω, χ1 |
| Long sidechains (Ile, Leu, etc.) | Up to χ2 or χ3 |
| Arg, Lys, Met | Up to χ4 |
| Arginine only | χ5 |

**Sources:** [soursop/tests/test_ssproteins.py:1037-1044]()

---

## Angle Decay Analysis

Angle decay measures the angular correlation between CA-CA vectors as a function of sequence separation. This provides insight into chain flexibility and persistence length.

### Conceptual Overview

The angle decay analysis computes the angle between CA-CA vectors for pairs of residues separated by varying numbers of residues along the sequence. High angles (close to 180°) indicate extended conformations, while lower angles indicate chain bending.

**Diagram: Angle Decay Calculation**

```mermaid
graph TD
    subgraph Chain["Protein Chain"]
        R1["Residue i"]
        R2["Residue i+1"]
        R3["Residue i+2"]
        R4["Residue i+3"]
        Dots["..."]
        RN["Residue i+n"]
    end
    
    subgraph Vectors["CA-CA Vectors"]
        V1["Vector i→i+1"]
        V2["Vector i+1→i+2"]
        V3["Vector i+2→i+3"]
        VN["Vector i+n-1→i+n"]
    end
    
    subgraph Angles["Angular Measurements"]
        A1["Angle(V1, V2)<br/>Separation = 1"]
        A2["Angle(V1, V3)<br/>Separation = 2"]
        AN["Angle(V1, VN)<br/>Separation = n"]
    end
    
    R1 -.-> V1
    R2 -.-> V1
    R2 -.-> V2
    R3 -.-> V2
    R3 -.-> V3
    
    V1 --> A1
    V2 --> A1
    V1 --> A2
    V3 --> A2
    V1 --> AN
    VN --> AN
    
    subgraph Output["Output Matrix"]
        Matrix["return_matrix[separation]<br/>= mean(all angles)"]
    end
    
    A1 --> Matrix
    A2 --> Matrix
    AN --> Matrix
```

**Sources:** [soursop/tests/test_ssproteins.py:672-732]()

### The `get_angle_decay()` Method

**Method Signature:**
```python
get_angle_decay(R1=None, R2=None, stride=1, return_all_pairs=False, verbose=True)
```

**Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `R1` | `int` or `None` | First residue in region (default: first CA residue) |
| `R2` | `int` or `None` | Last residue in region (default: last CA residue) |
| `stride` | `int` | Frame spacing for analysis (default: 1) |
| `return_all_pairs` | `bool` | If True, return dictionary of all individual pairs (default: False) |
| `verbose` | `bool` | Print status messages (default: True) |

**Returns:**

When `return_all_pairs=False`:
- `return_matrix`: NumPy array where `return_matrix[i]` = `[separation, mean_angle]` for separation distance `i`

When `return_all_pairs=True`:
- `(return_matrix, all_pairs)` tuple where `all_pairs` is a dictionary with keys like `"i-j"` mapping to angle arrays

**Return Matrix Structure:**

```python
return_matrix = np.array([
    [1, 145.2],  # Separation 1: mean angle 145.2°
    [2, 132.4],  # Separation 2: mean angle 132.4°
    [3, 125.1],  # Separation 3: mean angle 125.1°
    # ... decreasing angles indicate chain bending
])
```

**All Pairs Dictionary:**

```python
all_pairs = {
    "1-2": np.array([148.3, 142.1, ...]),  # Angles for residues 1→2 across frames
    "1-3": np.array([135.7, 129.2, ...]),  # Angles for residues 1→3 across frames
    "2-3": np.array([147.1, 143.8, ...]),  # Angles for residues 2→3 across frames
    # ... one entry per unique pair
}
```

The number of entries in `all_pairs` equals `∑(n-1, n-2, ..., 1)` where `n` is the number of CA-containing residues.

**Diagram: Angle Decay Data Flow**

```mermaid
graph TD
    Input["get_angle_decay(R1, R2, stride,<br/>return_all_pairs)"]
    
    subgraph Selection["Region Selection"]
        GetRegion["Determine R1, R2<br/>Default: all CA residues"]
        GetCA["Get CA indices<br/>via get_multiple_CA_index()"]
    end
    
    subgraph Computation["Vector and Angle Computation"]
        Pairs["Generate all unique pairs<br/>(i, j) where j > i"]
        Vectors["Compute CA-CA vectors<br/>for each pair"]
        Angles["Compute angles between<br/>consecutive vectors"]
        Group["Group by separation<br/>distance |j-i|"]
    end
    
    subgraph Output["Output Generation"]
        Average["Compute mean angle<br/>per separation"]
        Matrix["Build return_matrix<br/>[separation, mean_angle]"]
        Dict["Build all_pairs dict<br/>{i-j: angle_array}"]
    end
    
    Decision{"return_all_pairs?"}
    ReturnBoth["Return (matrix, dict)"]
    ReturnMatrix["Return matrix only"]
    
    Input --> GetRegion
    GetRegion --> GetCA
    GetCA --> Pairs
    Pairs --> Vectors
    Vectors --> Angles
    Angles --> Group
    Group --> Average
    Average --> Matrix
    Group --> Dict
    Matrix --> Decision
    Decision -->|"True"| ReturnBoth
    Decision -->|"False"| ReturnMatrix
    Dict --> ReturnBoth
```

**Sources:** [soursop/tests/test_ssproteins.py:672-732]()

---

## D-Vector Analysis

The D-vector provides a per-residue measure of local structural heterogeneity. It quantifies how much a residue's local environment varies across the ensemble.

### The `get_D_vector()` Method

The D-vector is computed by analyzing the distribution of local structural features for each residue.

**Method Signature:**
```python
get_D_vector(fragment_size=10, stride=1, backbone=True, weights=False, verbose=True)
```

**Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `fragment_size` | `int` | Size of local fragment around each residue (default: 10) |
| `stride` | `int` | Frame spacing for analysis (default: 1) |
| `backbone` | `bool` | Use backbone atoms only (default: True) |
| `weights` | `array-like` or `False` | Frame weights for reweighting (default: False) |
| `verbose` | `bool` | Print status messages (default: True) |

**Returns:**

- `d_vector`: NumPy array of length `n_residues` containing D-vector values for each residue

Higher D-vector values indicate greater local heterogeneity (more structural variability in the local environment).

**Diagram: D-Vector Computation Process**

```mermaid
graph TD
    Input["get_D_vector(fragment_size, stride)"]
    
    subgraph Iteration["Per-Residue Iteration"]
        SelectRes["For each residue i"]
        DefineFragment["Define fragment<br/>i ± fragment_size/2"]
        ExtractAtoms["Extract atoms<br/>in fragment region"]
    end
    
    subgraph Analysis["Local Structure Analysis"]
        ComputeDist["Compute pairwise distances<br/>within fragment"]
        BuildMatrix["Build distance matrix<br/>per frame"]
        Flatten["Flatten to 1D<br/>distance vector"]
    end
    
    subgraph Heterogeneity["Heterogeneity Measure"]
        CrossFrame["Compare distance vectors<br/>across frames"]
        CalcVar["Calculate variance or<br/>distribution width"]
        DValue["Assign D-value<br/>for residue i"]
    end
    
    Collect["Collect D-values<br/>for all residues"]
    Output["Return d_vector array"]
    
    Input --> SelectRes
    SelectRes --> DefineFragment
    DefineFragment --> ExtractAtoms
    ExtractAtoms --> ComputeDist
    ComputeDist --> BuildMatrix
    BuildMatrix --> Flatten
    Flatten --> CrossFrame
    CrossFrame --> CalcVar
    CalcVar --> DValue
    DValue --> Collect
    Collect --> Output
```

**Sources:** [soursop/tests/test_ssproteins.py:79]()

---

## Sidechain Alignment Angles

Sidechain alignment angles measure the relative orientation between sidechain vectors of two residues, providing information about packing and interactions.

### The `get_sidechain_alignment_angle()` Method

This method computes the angle between sidechain vectors defined by CA→sidechain_atom for two residues.

**Method Signature:**
```python
get_sidechain_alignment_angle(R1, R2, sidechain_atom_1='default', sidechain_atom_2='default', stride=1, verbose=True)
```

**Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `R1` | `int` | First residue index |
| `R2` | `int` | Second residue index |
| `sidechain_atom_1` | `str` | Atom name for R1 sidechain vector (default: residue-specific) |
| `sidechain_atom_2` | `str` | Atom name for R2 sidechain vector (default: residue-specific) |
| `stride` | `int` | Frame spacing for analysis (default: 1) |
| `verbose` | `bool` | Print status messages (default: True) |

**Returns:**

- NumPy array of alignment angles (degrees) for each frame

**Default Sidechain Atoms:**

The method uses residue-specific default atoms defined in `DEFAULT_SIDECHAIN_VECTOR_ATOMS`:

```python
# Examples of default sidechain vector atoms
'ALA': 'CB'
'CYS': 'SG'
'ASP': 'CG'
'GLU': 'CD'
'PHE': 'CZ'
'TRP': 'CH2'
# ... etc.
```

**Limitations:**

- Does not work for residues without sidechains (GLY)
- Does not work for caps (ACE, NME)
- Both residues must have the specified sidechain atoms

**Diagram: Sidechain Alignment Angle Geometry**

```mermaid
graph LR
    subgraph Residue1["Residue R1"]
        CA1["CA atom"]
        SC1["Sidechain atom<br/>(e.g., CB, CG)"]
    end
    
    subgraph Residue2["Residue R2"]
        CA2["CA atom"]
        SC2["Sidechain atom<br/>(e.g., CB, CG)"]
    end
    
    subgraph Vectors["Sidechain Vectors"]
        V1["Vector 1<br/>CA1 → SC1"]
        V2["Vector 2<br/>CA2 → SC2"]
    end
    
    subgraph Angle["Alignment Angle"]
        Theta["θ = angle(V1, V2)<br/>Range: 0° to 180°"]
    end
    
    CA1 -.-> V1
    SC1 -.-> V1
    CA2 -.-> V2
    SC2 -.-> V2
    V1 --> Theta
    V2 --> Theta
    
    Interpretation["θ ≈ 0°: Parallel sidechains<br/>θ ≈ 90°: Perpendicular<br/>θ ≈ 180°: Anti-parallel"]
```

**Sources:** [soursop/ssprotein.py:23](), [soursop/tests/test_ssproteins.py:782-882]()

---

## Common Patterns and Best Practices

### Caching Behavior

All angle methods use caching to improve performance:

```python
# First call computes and caches
phi_data = protein.get_angles('phi')  # Computes

# Subsequent calls return cached results
phi_data_again = protein.get_angles('phi')  # Instant
```

To clear the cache:
```python
protein.reset_cache()  # Clears self.__all_angles and other cached data
```

### Working with Stride

Using `stride` reduces computational cost and memory usage:

```python
# Analyze every frame (slow for large trajectories)
angles = protein.get_angles('phi', stride=1)

# Analyze every 10th frame (10x faster)
angles = protein.get_angles('phi', stride=10)
```

### Handling Missing Angles

Some residues may not have all angle types:

```python
phi_atoms, phi_values = protein.get_angles('phi')

# Check which residues have phi angles
print(f"Number of phi angles: {len(phi_atoms)}")
print(f"Total CA residues: {len(protein.resid_with_CA)}")

# N-terminal residue typically lacks phi
# C-terminal residue typically lacks psi
```

### Angle Decay Interpretation

Typical angle decay patterns:

| Mean Angle at Separation | Interpretation |
|--------------------------|----------------|
| > 150° | Highly extended, stiff chain |
| 120-150° | Moderately extended |
| 90-120° | Partially collapsed |
| < 90° | Highly collapsed, compact |

**Sources:** [soursop/ssprotein.py:146](), [soursop/ssprotein.py:176-218](), [soursop/tests/test_ssproteins.py:672-732]()

---