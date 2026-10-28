# Spatial Calculations and Boundaries

<details>
<summary>Relevant source files</summary>

The following files were used as context for generating this wiki page:

- [docs/src/documentation.md](docs/src/documentation.md)
- [src/coupling.jl](src/coupling.jl)
- [src/spatial.jl](src/spatial.jl)

</details>



This page covers the spatial infrastructure of Molly.jl, including periodic boundary conditions, coordinate transformations, and spatial utility functions. This system provides the fundamental geometric framework that underlies all molecular dynamics simulations, handling how atoms are positioned in space and how they interact across periodic boundaries.

For information about force calculations that use these spatial utilities, see [Force and Energy Calculations](#4.3). For details about the core System data structure that contains boundaries, see [System and Atoms](#3.1).

## Boundary Types

Molly.jl supports three types of boundary conditions for different simulation geometries. Each boundary type implements specific algorithms for coordinate wrapping and distance calculations.

### Boundary Type Overview

```mermaid
graph TB
    subgraph "Boundary Types"
        CB["CubicBoundary<br/>3D orthogonal box"]
        RB["RectangularBoundary<br/>2D orthogonal box"]  
        TB["TriclinicBoundary<br/>3D non-orthogonal box"]
    end
    
    subgraph "Core Operations"
        VO["vector()<br/>Distance calculation"]
        WC["wrap_coords()<br/>Coordinate wrapping"]
        RC["random_coord()<br/>Random placement"]
        VOL["volume()<br/>Box volume/area"]
    end
    
    subgraph "Advanced Operations"
        SB["scale_boundary()<br/>Resize boundaries"]
        BC["box_center()<br/>Geometric center"]
        BBL["bounding_box_lines()<br/>Visualization"]
    end
    
    CB --> VO
    RB --> VO
    TB --> VO
    
    CB --> WC
    RB --> WC
    TB --> WC
    
    CB --> RC
    RB --> RC
    TB --> RC
    
    CB --> VOL
    RB --> VOL
    TB --> VOL
    
    VO --> SB
    VO --> BC
    VO --> BBL
```

Sources: [src/spatial.jl:38-209](), [src/spatial.jl:210-254]()

### CubicBoundary

The `CubicBoundary` represents a 3D orthogonal simulation box with independent side lengths. It supports both periodic and infinite boundary conditions by setting dimensions to `Inf`.

**Key Properties:**
- Three independent side lengths stored as `SVector{3, T}`
- Supports infinite boundaries in any dimension
- Optimized distance calculations using `vector_1D()` for each axis
- Compatible with Chemfiles.jl for trajectory output

**Usage Examples:**
```julia
# Cubic periodic box
boundary = CubicBoundary(2.0u"nm")
# Rectangular box  
boundary = CubicBoundary(2.0u"nm", 3.0u"nm", 4.0u"nm")
# Infinite in z-direction
boundary = CubicBoundary(2.0u"nm", 2.0u"nm", Inf*u"nm")
```

Sources: [src/spatial.jl:38-70](), [src/spatial.jl:412-418]()

### RectangularBoundary

The `RectangularBoundary` provides 2D simulation capabilities with two orthogonal side lengths.

**Key Properties:**
- Two side lengths for x and y dimensions
- Supports 2D molecular dynamics simulations
- Infinite boundaries supported via `Inf` values
- Vector operations padded to 3D internally for compatibility

Sources: [src/spatial.jl:81-105](), [src/spatial.jl:420-425]()

### TriclinicBoundary

The `TriclinicBoundary` supports non-orthogonal 3D simulation boxes defined by three basis vectors, enabling simulations of crystal lattices and complex geometries.

**Key Properties:**
- Three basis vectors defining the parallelpiped
- Approximate or exact minimum image convention
- Precomputed trigonometric values for performance
- Not compatible with infinite boundaries

```mermaid
graph LR
    subgraph "TriclinicBoundary Construction"
        BV["Basis Vectors<br/>v1, v2, v3"]
        AL["Lengths + Angles<br/>α, β, γ"]
    end
    
    subgraph "Internal Storage"
        BVS["basis_vectors<br/>SVector{3, SVector{3, D}}"]
        ANG["α, β, γ<br/>Angles in radians"]
        RS["reciprocal_size<br/>Inverse scaling"]
        TRIG["Precomputed<br/>Trigonometric values"]
    end
    
    subgraph "Operations"
        VEC["vector()<br/>Min image distance"]
        WRAP["wrap_coords()<br/>Complex wrapping"]
        VOL["volume()<br/>Determinant-based"]
    end
    
    BV --> BVS
    AL --> BVS
    BVS --> ANG
    BVS --> RS
    BVS --> TRIG
    
    BVS --> VEC
    ANG --> VEC
    RS --> VEC
    TRIG --> VEC
    
    VEC --> WRAP
    VEC --> VOL
```

Sources: [src/spatial.jl:127-209](), [src/spatial.jl:427-450]()

## Coordinate Operations

The spatial system provides fundamental operations for working with coordinates within different boundary conditions.

### Distance Calculations

The `vector()` function calculates the displacement between two coordinates using the minimum image convention, which finds the shortest distance accounting for periodic boundaries.

```mermaid
graph TB
    subgraph "Distance Calculation Flow"
        C1["coord_1"]
        C2["coord_2"]
        B["boundary"]
        V1D["vector_1D()<br/>1D displacement"]
        VEC["vector()<br/>Full displacement"]
        VFAST["Fast path<br/>CubicBoundary"]
        VSLOW["Complex path<br/>TriclinicBoundary"]
    end
    
    C1 --> VEC
    C2 --> VEC
    B --> VEC
    
    VEC --> VFAST
    VEC --> VSLOW
    
    VFAST --> V1D
    VSLOW --> V1D
    
    subgraph "Implementation Details"
        APPROX["Approximate images<br/>(fast)"]
        EXACT["Exact images<br/>(slow but correct)"]
    end
    
    VSLOW --> APPROX
    VSLOW --> EXACT
```

Sources: [src/spatial.jl:382-461]()

### Coordinate Wrapping

The `wrap_coords()` function ensures coordinates remain within the simulation box by applying periodic boundary conditions.

**Algorithm Variations:**
- **Cubic/Rectangular:** Simple modulo operation per dimension
- **Triclinic:** Complex sequential wrapping accounting for basis vector relationships

Sources: [src/spatial.jl:481-500]()

### Random Coordinate Generation

The `random_coord()` function generates uniformly distributed random coordinates within the boundary volume.

Sources: [src/spatial.jl:365-379]()

## Physical Property Calculations

The spatial system implements several important physical calculations that depend on the geometric configuration.

### Pressure and Virial Calculations

```mermaid
graph TB
    subgraph "Pressure Calculation"
        SYS["System"]
        NEIGHBORS["Neighbors"]
        VIR["virial()<br/>Force-distance product"]
        NKT["NkT<br/>Kinetic contribution"]
        PRES["pressure()<br/>Final pressure"]
    end
    
    subgraph "Virial Components"
        PAIR["Pairwise interactions"]
        GEN["General interactions"]
        DOT["dot(force, distance)<br/>per atom pair"]
    end
    
    SYS --> VIR
    NEIGHBORS --> VIR
    VIR --> PAIR
    VIR --> GEN
    PAIR --> DOT
    GEN --> DOT
    
    VIR --> PRES
    NKT --> PRES
    
    subgraph "Formula"
        FORM["P = (NkT - 2Ξ/D) / V<br/>Ξ = virial, D = dimensions"]
    end
    
    PRES --> FORM
```

Sources: [src/spatial.jl:726-855]()

### Density and Volume

The system calculates density using the total mass and boundary volume:

```julia
density(sys) = total_mass(sys) / volume(sys.boundary)
```

Volume calculations are boundary-type specific:
- **Cubic/Rectangular:** Product of side lengths
- **Triclinic:** Determinant of basis vectors

Sources: [src/spatial.jl:268-292]()

### Molecular Properties

The spatial system provides calculations for molecular-scale properties:

- `molecule_centers()` - Center of mass for each molecule using circular mean for periodic boundaries
- `dipole_moment()` - Electric dipole moment from charges and coordinates
- `scale_coords!()` - Coordinated scaling of coordinates and boundary

Sources: [src/spatial.jl:857-960]()

## Velocity Generation

The spatial system includes Maxwell-Boltzmann velocity generation for thermal equilibration.

### Maxwell-Boltzmann Distribution

```mermaid
graph LR
    subgraph "Velocity Generation"
        MASS["Atom mass"]
        TEMP["Temperature"]
        K["Boltzmann constant"]
        SIGMA["σ = sqrt(kT/m)"]
        NORM["Normal(0, σ)"]
        VEL["Final velocity"]
    end
    
    MASS --> SIGMA
    TEMP --> SIGMA
    K --> SIGMA
    SIGMA --> NORM
    NORM --> VEL
    
    subgraph "System-level Generation"
        SYS["System"]
        DIMS["Dimensions (2D/3D)"]
        RVELS["random_velocities()"]
        GPU["GPU compatibility"]
    end
    
    VEL --> RVELS
    SYS --> DIMS
    DIMS --> RVELS
    RVELS --> GPU
```

The `maxwell_boltzmann()` function generates individual velocity components, while `random_velocities()` generates complete velocity vectors for all atoms in a system.

Sources: [src/spatial.jl:502-656]()

## Integration with Simulation System

The spatial calculations integrate deeply with the broader Molly.jl simulation infrastructure.

### System Integration

```mermaid
graph TB
    subgraph "Core System Components"
        SYS["System"]
        COORDS["coordinates"]
        BOUNDARY["boundary"]
        ATOMS["atoms"]
    end
    
    subgraph "Spatial Operations"
        FORCES["Force calculations"]
        NEIGHBORS["Neighbor finding"]
        COUPLING["Coupling methods"]
        CONSTRAINTS["Constraints"]
    end
    
    subgraph "Physical Properties"
        ENERGY["Energy calculations"]
        PRESSURE["Pressure calculation"]
        TEMPERATURE["Temperature measurement"]
        ANALYSIS["Analysis functions"]
    end
    
    SYS --> COORDS
    SYS --> BOUNDARY
    SYS --> ATOMS
    
    COORDS --> FORCES
    BOUNDARY --> FORCES
    COORDS --> NEIGHBORS
    BOUNDARY --> NEIGHBORS
    
    FORCES --> COUPLING
    NEIGHBORS --> COUPLING
    BOUNDARY --> COUPLING
    
    FORCES --> CONSTRAINTS
    COORDS --> CONSTRAINTS
    BOUNDARY --> CONSTRAINTS
    
    COORDS --> ENERGY
    BOUNDARY --> ENERGY
    COORDS --> PRESSURE
    BOUNDARY --> PRESSURE
    COORDS --> TEMPERATURE
    
    COORDS --> ANALYSIS
    BOUNDARY --> ANALYSIS
```

Sources: [src/spatial.jl:1-28](), [docs/src/documentation.md:41-50]()

### GPU Compatibility

The spatial calculation functions are designed to work transparently with both CPU and GPU array types. The boundary types are bits types, enabling efficient GPU kernels, while coordinate operations use broadcasting for GPU compatibility.

Sources: [src/spatial.jl:633-639](), [docs/src/documentation.md:142-194]()

---