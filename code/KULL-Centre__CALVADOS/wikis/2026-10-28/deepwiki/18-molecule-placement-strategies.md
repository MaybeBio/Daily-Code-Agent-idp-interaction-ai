# Molecule Placement Strategies

<details>
<summary>Relevant source files</summary>

The following files were used as context for generating this wiki page:

- [calvados/build.py](calvados/build.py)
- [calvados/components.py](calvados/components.py)
- [calvados/sequence.py](calvados/sequence.py)
- [calvados/sim.py](calvados/sim.py)
- [examples/single_dsRNA/input/domains.yaml](examples/single_dsRNA/input/domains.yaml)
- [examples/single_dsRNA/input/dspolyR12.pdb](examples/single_dsRNA/input/dspolyR12.pdb)
- [examples/single_dsRNA/input/fastalib.fasta](examples/single_dsRNA/input/fastalib.fasta)
- [examples/single_dsRNA/input/residues_C2RNA.csv](examples/single_dsRNA/input/residues_C2RNA.csv)
- [examples/single_dsRNA/prepare.py](examples/single_dsRNA/prepare.py)

</details>



## Purpose and Scope

This page documents the algorithms and strategies used by CALVADOS to position molecules in the simulation box during system initialization. Molecule placement occurs during the `build_system` phase and determines the initial spatial configuration of all components before simulation begins. 

For information about system building workflow, see [Sim Class & System Building](#4.1). For details on simulation execution after placement, see [Equilibration & Production](#4.3). For component-specific properties and structure building, see [Component Class Hierarchy](#3.1).

---

## Placement Workflow Overview

Molecule placement in CALVADOS follows a two-stage process: first, each component builds its internal structure (molecular conformation), then the entire molecule is positioned in the simulation box according to the chosen topology.

```mermaid
graph TB
    subgraph "Stage 1: Molecular Conformation"
        COMP["Component<br/>calc_x_setup()"]
        SPIRAL["build_spiral()<br/>Archimedes spiral"]
        COMPACT["build_compact()<br/>Cubic lattice"]
        LINEAR["build_linear()<br/>Extended chain"]
        PDB["geometry_from_pdb()<br/>From structure"]
        
        COMP --> SPIRAL
        COMP --> COMPACT
        COMP --> LINEAR
        COMP --> PDB
        
        SPIRAL --> XINIT["xinit<br/>Internal coordinates"]
        COMPACT --> XINIT
        LINEAR --> XINIT
        PDB --> XINIT
    end
    
    subgraph "Stage 2: Box Placement"
        TOPOL["topol parameter"]
        SLAB["Slab placement<br/>build_xyzgrid()"]
        GRID["Grid placement<br/>build_xyzgrid()"]
        CENTER["Center placement<br/>box * 0.5"]
        SHIFT["Shift ref bead<br/>center - xinit[ref_bead]"]
        RANDOM["Random placement<br/>random_placement()"]
        BILAYER["Bilayer placement<br/>build_xybilayer()"]
        
        TOPOL --> SLAB
        TOPOL --> GRID
        TOPOL --> CENTER
        TOPOL --> SHIFT
        TOPOL --> RANDOM
        TOPOL --> BILAYER
    end
    
    XINIT --> SLAB
    XINIT --> GRID
    XINIT --> CENTER
    XINIT --> SHIFT
    XINIT --> RANDOM
    XINIT --> BILAYER
    
    SLAB --> POS["pos array<br/>System positions"]
    GRID --> POS
    CENTER --> POS
    SHIFT --> POS
    RANDOM --> POS
    BILAYER --> POS
    
    style XINIT fill:#e1f5ff
    style POS fill:#ffe1e1
    style TOPOL fill:#fff4e1
```

**Sources:** [calvados/sim.py:292-336](), [calvados/components.py:63-70](), [calvados/build.py:84-294]()

---

## Topology Modes

The `topol` configuration parameter controls the overall placement strategy. Each mode is designed for specific simulation scenarios.

| Topology | Description | Use Case | Grid Function | Constraints |
|----------|-------------|----------|---------------|-------------|
| `slab` | Molecules in central slab region | Phase separation studies | `build_xyzgrid` | Crowders in outer regions |
| `grid` | Regular 3D grid placement | Multi-component systems | `build_xyzgrid` | All molecules on grid |
| `center` | Single molecule at box center | Single IDR simulations | None | Only 1 molecule |
| `shift_ref_bead` | Center with reference bead aligned | Folded protein analysis | None | Only 1 molecule |
| `random` | Random non-overlapping placement | Default fallback | None | Clash detection required |

**Sources:** [calvados/sim.py:121-122](), [calvados/sim.py:292-314]()

### Slab Topology

Designed for liquid-liquid phase separation simulations, the slab topology concentrates proteins/RNA in a central slab region while placing crowders in outer regions.

```mermaid
graph TB
    subgraph "Slab Geometry Setup"
        PARAMS["Parameters:<br/>box, slab_width, slab_outer"]
        XYZGRID1["build_xyzgrid(nproteins+nrnas,<br/>[Lx, Ly, slab_width])"]
        SHIFT1["Shift to box center:<br/>z += box[2]/2 - slab_width/2"]
        
        XYZGRID2["build_xyzgrid(ncrowders/2,<br/>[Lx, Ly, box[2]/2 - slab_outer])"]
        APPEND1["Append to xyzgrid"]
        SHIFT2["Duplicate + shift:<br/>z += box[2]/2 + slab_outer"]
        APPEND2["Append to xyzgrid"]
        
        PARAMS --> XYZGRID1
        XYZGRID1 --> SHIFT1
        
        PARAMS --> XYZGRID2
        XYZGRID2 --> APPEND1
        APPEND1 --> SHIFT2
        SHIFT2 --> APPEND2
        
        SHIFT1 --> FINAL["Final xyzgrid:<br/>proteins/RNA in center,<br/>crowders in outer regions"]
        APPEND2 --> FINAL
    end
    
    FINAL --> PLACE["place_molecule()<br/>assigns grid positions"]
    
    style FINAL fill:#e1ffe1
```

**Sources:** [calvados/sim.py:172-178]()

### Grid Topology

Standard 3D lattice placement for general multi-component systems.

```mermaid
graph LR
    NMOL["nmolecules"] --> BUILD["build_xyzgrid(nmolecules, box)"]
    BUILD --> GRID["xyzgrid array"]
    GRID --> COUNTER["grid_counter<br/>(increments per molecule)"]
    COUNTER --> ASSIGN["x0 = xyzgrid[grid_counter]<br/>xs = x0 + xinit"]
    
    style GRID fill:#e1f5ff
```

The `build_xyzgrid` function distributes molecules as uniformly as possible across the box volume, accounting for different box aspect ratios.

**Sources:** [calvados/sim.py:179-180](), [calvados/sim.py:302-305]()

---

## Initial Molecular Conformations

Before placement in the box, each component builds its internal coordinate array (`xinit`). The conformation builder is selected automatically based on molecule type, or can be controlled via the `comp_setup` parameter.

```mermaid
graph TB
    subgraph "Conformation Selection Logic"
        MOLTYPE["molecule_type"]
        
        MOLTYPE -->|protein, cyclic, seastar, ptm_protein| COMPACT_S["comp_setup = 'compact'"]
        MOLTYPE -->|rna| SPIRAL_S["comp_setup = 'spiral'"]
        MOLTYPE -->|lipid, cooke_lipid, crowder| LINEAR_S["comp_setup = 'linear'"]
        
        COMPACT_S --> CALCX["calc_x_setup(comp_setup)"]
        SPIRAL_S --> CALCX
        LINEAR_S --> CALCX
    end
    
    subgraph "calc_x_setup dispatcher"
        CALCX --> COND1{comp_setup?}
        COND1 -->|'spiral'| SPIRAL_B["build_spiral(bondlengths, arc, n_per_res)"]
        COND1 -->|'compact'| COMPACT_B["build_compact(nbeads, d)"]
        COND1 -->|else| LINEAR_B["build_linear(bondlengths, n_per_res, ys)"]
        
        SPIRAL_B --> XINIT["xinit array<br/>(nbeads × 3)"]
        COMPACT_B --> XINIT
        LINEAR_B --> XINIT
    end
    
    style XINIT fill:#ffe1e1
```

**Sources:** [calvados/sim.py:50-87](), [calvados/components.py:63-70]()

### Spiral Conformation

Creates an Archimedes spiral, commonly used for unfolded proteins and RNA.

```mermaid
graph TB
    SPIRAL["build_spiral(bondlengths, arc=0.38, separation=0.7, n_per_res=1)"]
    
    INIT["Initialize:<br/>r = arc<br/>b = separation / (2π)<br/>φ = r / b"]
    
    LOOP["For each residue i:"]
    CONVERT["(x, y) = (r·cos(φ), r·sin(φ))"]
    BEADS["For j in range(n_per_res):<br/>coords[i·n_per_res + j] = [x, y, j·z]"]
    UPDATE["φ += arc / r<br/>r = b·φ"]
    
    SPIRAL --> INIT
    INIT --> LOOP
    LOOP --> CONVERT
    CONVERT --> BEADS
    BEADS --> UPDATE
    UPDATE --> LOOP
    
    style SPIRAL fill:#e1f5ff
```

Parameters:
- `arc`: Distance between consecutive points along the spiral (default 0.38 nm)
- `separation`: Distance between consecutive turnings (default 0.7 nm)
- `n_per_res`: Number of beads per residue (1 for proteins, 2 for RNA)

**Sources:** [calvados/build.py:108-124]()

### Compact Conformation

Generates a cubic lattice arrangement, suitable for globular proteins or crowders.

The algorithm fills a 3D grid by alternating directions to create a space-filling path:

```mermaid
graph TB
    INPUT["build_compact(nbeads, d=0.38)"]
    
    CALC["N = ⌈∛nbeads⌉ - 1<br/>(grid size N+1 per dimension)"]
    
    INIT["Initialize position:<br/>i, j, k = 0, 0, 0<br/>directions: di, dj, dk = 1, 1, 1<br/>counters: cti, ctj, ctk = 0, 0, 0"]
    
    LOOP["For bead 0 to nbeads-1:"]
    APPEND["xs.append([i, j, k])"]
    
    CHECK1{ctk == N?}
    CHECK2{ctj == N?}
    
    MOVEI["i += di<br/>cti += 1<br/>ctj = 0<br/>dj *= -1"]
    MOVEJ["j += dj<br/>ctj += 1"]
    MOVEK["k += dk<br/>ctk += 1"]
    
    RESET1["ctk = 0<br/>dk *= -1"]
    RESET2["cti = 0<br/>cty = 0<br/>x = 0<br/>y = 0<br/>ctz += 1<br/>z += dz<br/>zplane *= -1"]
    
    SCALE["xs = (xs - 0.5·N) · d<br/>(center at origin)"]
    
    INPUT --> CALC
    CALC --> INIT
    INIT --> LOOP
    LOOP --> APPEND
    APPEND --> CHECK1
    
    CHECK1 -->|Yes| CHECK2
    CHECK1 -->|No| MOVEK
    
    CHECK2 -->|Yes| MOVEI
    CHECK2 -->|No| MOVEJ
    
    MOVEI --> RESET1
    MOVEJ --> RESET1
    MOVEK --> LOOP
    
    MOVEI --> RESET2
    RESET1 --> LOOP
    RESET2 --> LOOP
    
    LOOP -->|Done| SCALE
    
    style SCALE fill:#ffe1e1
```

**Sources:** [calvados/build.py:126-152]()

### Linear Conformation

Builds an extended chain growing along the z-axis, used for lipids and membrane systems.

**Sources:** [calvados/build.py:84-100]()

---

## Grid Generation Algorithms

### 3D Grid (`build_xyzgrid`)

The 3D grid algorithm distributes N molecules uniformly across a box with dimensions `[Lx, Ly, Lz]`, accounting for non-cubic boxes.

```mermaid
graph TB
    INPUT["build_xyzgrid(N, box)"]
    
    RATIO["r = box / sum(box)<br/>(relative proportions)"]
    CALC["a = ∛(N / product(r))<br/>nxyz = floor(a · r)"]
    
    ADJ1["While product(nxyz) < N:<br/>increase dimension with max deviation"]
    ADJ2["While product(nxyz) > N:<br/>decrease dimension with max count"]
    
    SPACING["dx = box[0] / nxyz[0]<br/>dy = box[1] / nxyz[1]<br/>dz = box[2] / nxyz[2]"]
    
    FILL["Fill grid with alternating offsets:<br/>zplane alternation (dx/2, dy/2)<br/>xyplane alternation (dz/2)"]
    
    INPUT --> RATIO
    RATIO --> CALC
    CALC --> ADJ1
    ADJ1 --> ADJ2
    ADJ2 --> SPACING
    SPACING --> FILL
    
    style FILL fill:#e1ffe1
```

Key features:
- Maintains aspect ratio of the box
- Uses alternating offsets between layers to improve uniformity
- Iteratively adjusts grid dimensions to match exact molecule count

**Sources:** [calvados/build.py:222-294]()

### 2D Grid (`build_xygrid`)

For bilayer and slab systems, the 2D grid distributes molecules across the xy-plane at a fixed z-coordinate.

**Sources:** [calvados/build.py:195-220]()

---

## Random Placement with Clash Detection

When no specific topology is set, molecules are placed randomly while avoiding overlaps.

```mermaid
graph TB
    START["random_placement(box, xs_others, xinit, ntries=10000)"]
    
    TRY["ntry = 0"]
    LOOP["ntry += 1"]
    CHECK{ntry > ntries?}
    
    DRAW["x0 = draw_starting_vec(box)<br/>(random position)"]
    ADD["xs = x0 + xinit"]
    
    WALLS["check_walls(xs, box)"]
    CLASH["check_clash(xs, xs_others, box, cutoff=0.7)"]
    
    WALLCHECK{Outside box?}
    CLASHCHECK{Overlap detected?}
    
    SUCCESS["Return xs"]
    FAIL["Raise ValueError:<br/>'Tried 10000x to add molecule'"]
    
    START --> TRY
    TRY --> LOOP
    LOOP --> CHECK
    CHECK -->|Yes| FAIL
    CHECK -->|No| DRAW
    DRAW --> ADD
    ADD --> WALLS
    WALLS --> WALLCHECK
    WALLCHECK -->|Yes| LOOP
    WALLCHECK -->|No| CLASH
    CLASH --> CLASHCHECK
    CLASHCHECK -->|Yes| LOOP
    CLASHCHECK -->|No| SUCCESS
    
    style SUCCESS fill:#e1ffe1
    style FAIL fill:#ffe1e1
```

### Clash Detection Algorithm

The `check_clash` function uses MDAnalysis distance calculations to detect particle overlaps:

```python
# Simplified logic from calvados/build.py:52-63
def check_clash(x, pos, box, cutoff=0.7):
    boxfull = np.append(box, [90,90,90])  # box + angles
    xothers = np.array(pos)
    if len(xothers) == 0:
        return False  # no other particles
    d = distances.distance_array(x, xothers, boxfull)
    if np.amin(d) < cutoff:
        return True  # clash detected
    else:
        return False  # no clash
```

**Sources:** [calvados/build.py:154-170](), [calvados/build.py:42-63]()

---

## Bilayer Systems

Lipid bilayers require specialized placement to create membrane structures.

```mermaid
graph TB
    subgraph "Bilayer Setup"
        GRID["build_xygrid(nlipids * 1.05, box)<br/>(2D grid with 5% excess)"]
        PROTEIN["If proteins/RNA present:<br/>build_xyzgrid for proteins<br/>in regions above/below membrane"]
    end
    
    subgraph "Bilayer Insertion"
        SELECT["Select grid point: bilayergrid[0]"]
        ORIENT["Try upward orientation first"]
        ALIGN["Align lipid to z = box[2]/2 ± 1.5"]
        CHECK["check_walls() and check_clash(cutoff=0.5)"]
        
        SUCCESS{Inserted?}
        FLIP["Try downward orientation"]
        SHUFFLE["Randomly shuffle bilayergrid"]
        
        SELECT --> ORIENT
        ORIENT --> ALIGN
        ALIGN --> CHECK
        CHECK --> SUCCESS
        SUCCESS -->|No| FLIP
        FLIP --> ALIGN
        SUCCESS -->|Still No| SHUFFLE
        SHUFFLE --> SELECT
        SUCCESS -->|Yes| NEXT["Next lipid"]
    end
    
    GRID --> SELECT
    
    style NEXT fill:#e1ffe1
```

The bilayer is constructed at the center of the box (z = box[2]/2) with lipids oriented perpendicular to the membrane plane. The 5% excess grid points provide flexibility when placement attempts fail due to steric clashes.

**Sources:** [calvados/sim.py:181-190](), [calvados/sim.py:320-336](), [calvados/build.py:172-193]()

---

## Placement Method Routing in Sim Class

The `Sim.build_system` method routes each component to the appropriate placement function based on molecule type:

```mermaid
graph TB
    BUILD["For each component in system:"]
    
    CHECK{molecule_type?}
    
    PMOL["place_molecule(comp)"]
    PBIL["place_bilayer(comp)"]
    
    BUILD --> CHECK
    CHECK -->|"protein, crowder,<br/>cyclic, seastar,<br/>ptm_protein, rna"| PMOL
    CHECK -->|"lipid,<br/>cooke_lipid"| PBIL
    
    PMOL --> TOPOL{topol?}
    TOPOL -->|slab| SLAB["xs = xyzgrid[grid_counter] + xinit<br/>grid_counter++"]
    TOPOL -->|grid| GRID["xs = xyzgrid[grid_counter] + xinit<br/>grid_counter++"]
    TOPOL -->|center| CENTER["xs = box * 0.5 + xinit"]
    TOPOL -->|shift_ref_bead| SHIFT["xs = box * 0.5 + xinit<br/>xs -= xinit[ref_bead]"]
    TOPOL -->|else| RANDOM["xs = random_placement(box, pos, xinit)"]
    
    SLAB --> POS["Append xs to pos array"]
    GRID --> POS
    CENTER --> POS
    SHIFT --> POS
    RANDOM --> POS
    PBIL --> POS
    
    POS --> COUNTER["nparticles += comp.nbeads"]
    
    style POS fill:#ffe1e1
```

**Sources:** [calvados/sim.py:192-207]()

---

## Configuration Example

Example showing how placement strategies are configured:

```yaml
# config.yaml
topol: 'slab'          # Topology mode
box: [25, 25, 100]     # Box dimensions [nm]
slab_width: 20.0       # Slab region width [nm]
slab_outer: 5.0        # Distance from slab to crowder regions [nm]
```

```yaml
# components.yaml
system:
  protein_A:
    molecule_type: 'protein'
    nmol: 50
    # Will be placed in central slab region
  
  crowder_PEG:
    molecule_type: 'crowder'
    nmol: 200
    # Will be placed in outer regions for slab topology
```

For a single molecule in the center:

```yaml
# config.yaml
topol: 'center'
box: [20, 20, 20]

# components.yaml
system:
  my_IDR:
    molecule_type: 'protein'
    nmol: 1
```

**Sources:** [calvados/sim.py:44](), [calvados/cfg.py]()

---

## Summary Table: Placement Strategy Selection

| Component Type | Default comp_setup | Recommended topol | Grid Used | Special Handling |
|----------------|-------------------|-------------------|-----------|------------------|
| Protein (unfolded) | `compact` | `slab`, `grid`, `random` | xyzgrid | None |
| Protein (restraints) | from PDB | `center`, `shift_ref_bead` | None | Uses PDB coordinates |
| RNA | `spiral` | `slab`, `grid` | xyzgrid | Two-bead model (n_per_res=2) |
| Lipid | `linear` | N/A | xygrid (bilayer) | Oriented perpendicular to membrane |
| Crowder | `compact` | `slab` (outer), `grid` | xyzgrid | Excluded volume particles |
| Cyclic | `compact` | `grid`, `random` | xyzgrid | Closed ring topology |

**Sources:** [calvados/sim.py:56-84](), [calvados/components.py:168-188]()

---