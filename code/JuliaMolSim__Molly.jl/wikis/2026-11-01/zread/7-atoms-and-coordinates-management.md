

This page explores how Molly.jl represents and manipulates atomic data and coordinate systems—the foundational structures for molecular dynamics simulations. Understanding these core data structures is essential for effective simulation design and analysis, whether you're building simple gas phase systems or complex biomolecular assemblies.

## Atom Representation Architecture

The atomic representation in Molly.jl separates computational efficiency from descriptive information through a dual-structure design. This architecture enables high-performance GPU computation while maintaining rich metadata for analysis and visualization.

The `Atom` struct serves as the primary computational entity, storing essential physical properties required for force calculations and integration. It is designed as a bits type to enable GPU acceleration, containing only numeric fields needed during simulation steps.

Sources: [src/types.jl](src/types.jl#L229-L257), [src/types.jl](src/types.jl#L250-L257)

```julia
@kwdef struct Atom{T, M, C, S, E}
    index::Int = 1
    atom_type::T = 1
    mass::M = 1.0u"g/mol"
    charge::C = 0.0
    σ::S = 0.0u"nm"
    ϵ::E = 0.0u"kJ * mol^-1"
end
```

The `AtomData` struct complements this by storing descriptive metadata that doesn't participate in force calculations. This separation allows the `Atom` type to remain GPU-compatible while preserving detailed molecular information for analysis.

Sources: [src/types.jl](src/types.jl#L341-L358)

```julia
@kwdef struct AtomData
    atom_type::String = "?"
    atom_name::String = "?"
    res_number::Int = 1
    res_name::String = "???"
    chain_id::String = "A"
    element::String = "?"
    hetero_atom::Bool = false
end
```

> [!TIP]
> The dual-structure design is critical for GPU performance: always store `Atom` instances in GPU arrays, while keeping `AtomData` on CPU for analysis. This pattern ensures force calculations remain on-device without transferring unnecessary metadata.

## Atom Field Reference

| Field | Type | Purpose | Units | Default |
|-------|------|---------|-------|---------|
| `index` | `Int` | Atom position in the system vector | N/A | 1 |
| `atom_type` | `T` | Type identifier for force field lookup | N/A | 1 |
| `mass` | `M` | Atomic mass for kinetic energy and integration | g/mol | 1.0u"g/mol" |
| `charge` | `C` | Partial charge for electrostatics | e (dimensionless) | 0.0 |
| `σ` | `S` | Lennard-Jones size parameter | nm | 0.0u"nm" |
| `ϵ` | `E` | Lennard-Jones well depth | kJ/mol | 0.0u"kJ * mol^-1" |

Sources: [src/types.jl](src/types.jl#L239-L248)

## AtomData Field Reference

| Field | Type | Description |
|-------|------|-------------|
| `atom_type` | `String` | Force field atom type identifier |
| `atom_name` | `String` | PDB-format atom name |
| `res_number` | `Int` | Residue sequence number |
| `res_name` | `String` | Three-letter residue name |
| `chain_id` | `String` | Chain identifier in macromolecules |
| `element` | `String` | Chemical element symbol |
| `hetero_atom` | `Bool` | Non-standard atom flag |

Sources: [src/types.jl](src/types.jl#L342-L357)

## The System Container

The `System` struct serves as the central container for all simulation data, organizing atomic coordinates, velocities, boundary conditions, and interactions into a single coherent entity. It implements the `AtomsBase.AbstractSystem` interface for interoperability with Julia's broader ecosystem.

Sources: [src/types.jl](src/types.jl#L469-L538)

```mermaid
flowchart TD
    A[System] --> B[atoms<br/>Vector{Atom}]
    A --> C[coords<br/>Vector{SVector}]
    A --> D[boundary<br/>AbstractBoundary]
    A --> E[velocities<br/>Vector{SVector}]
    A --> F[atoms_data<br/>Vector{AtomData}]
    A --> G[topology<br/>MolecularTopology]
    A --> H[interactions<br/>pairwise/specific/general]
    A --> I[constraints<br/>Tuple]
    A --> J[loggers<br/>Tuple]
    A --> K[units<br/>force/energy]
    
    B --> B1[GPU-compatible<br/>bits type]
    C --> C1[SVector{3,<:Quantity}]
    D --> D1[Cubic/Rectangular/<br/>TriclinicBoundary]
    E --> E1[Optional: zero by default]
    F --> F1[CPU metadata only]
```

The System type is heavily parameterized to enable compile-time optimization, encoding dimensionality, array types, and unit specifications directly in the type signature. This design allows the compiler to generate specialized code for specific simulation configurations.

Sources: [src/types.jl](src/types.jl#L517-L538)

## Core System Properties

| Property | Type | Requirement | Description |
|----------|------|-------------|-------------|
| `atoms` | `A` | Required | Vector of Atom structures |
| `coords` | `C` | Required | Atomic positions (typically `SVector{3}`) |
| `boundary` | `B` | Required | Simulation boundary conditions |
| `velocities` | `V` | Optional | Atomic velocities (defaults to zero) |
| `atoms_data` | `AD` | Optional | Descriptive atom metadata |
| `topology` | `TO` | Optional | Molecular connectivity information |

Sources: [src/types.jl](src/types.jl#L483-L493)

## Coordinate Systems and Boundary Conditions

Molly.jl provides three boundary condition types, each optimized for different simulation scenarios. The boundary determines how coordinates are interpreted, how periodic images are applied, and how distances are calculated.

### Cubic Boundary

The `CubicBoundary` represents an orthogonal 3D box with equal side lengths. It's the most common choice for simple fluids and isotropic systems. Setting side lengths to `Inf` creates an unbounded dimension for non-periodic simulations.

Sources: [src/spatial.jl](src/spatial.jl#L31-L69)

```julia
# Create a 5 nm cubic box
boundary = CubicBoundary(5.0u"nm")

# Create a rectangular box from individual dimensions
boundary = CubicBoundary(4.0u"nm", 5.0u"nm", 6.0u"nm")

# Create a box with infinite z-dimension (slab geometry)
boundary = CubicBoundary(5.0u"nm", 5.0u"nm", Inf)
```

### Rectangular Boundary

The `RectangularBoundary` provides a 2D simulation box, ideal for membrane systems or surface studies where one dimension is irrelevant to the physics of interest.

Sources: [src/spatial.jl](src/spatial.jl#L84-L127)

```julia
# Create a 10x8 nm 2D box
boundary = RectangularBoundary(10.0u"nm", 8.0u"nm")
```

### Triclinic Boundary

The `TriclinicBoundary` represents the most general periodic cell with three arbitrary basis vectors. This is essential for crystalline systems and pressure-controlled simulations where box shape can change.

Sources: [src/spatial.jl](src/spatial.jl#L128-L148)

```julia
# Define by basis vectors
v1 = SVector(3.0u"nm", 0.0u"nm", 0.0u"nm")
v2 = SVector(1.5u"nm", 3.0u"nm", 0.0u"nm")
v3 = SVector(0.5u"nm", 1.0u"nm", 4.0u"nm")
boundary = TriclinicBoundary(SVector(v1, v2, v3))

# Define by lengths and angles (in radians)
L = SVector(3.0u"nm", 3.0u"nm", 4.0u"nm")
angles = SVector(π/2, π/2, π/3)  # 90°, 90°, 60°
boundary = TriclinicBoundary(L, angles)
```

## Boundary Type Comparison

| Feature | CubicBoundary | RectangularBoundary | TriclinicBoundary |
|---------|---------------|---------------------|-------------------|
| Dimensions | 3D | 2D | 3D |
| Flexibility | Orthogonal only | Orthogonal only | Arbitrary angles |
| Memory | Minimal | Minimal | Higher |
| Speed | Fastest | Fast | Moderate |
| Use Cases | Fluids, gases | Membranes, surfaces | Crystals, NPT |

Sources: [src/spatial.jl](src/spatial.jl#L31-L148), [src/spatial.jl](src/spatial.jl#L282-L304)

## Coordinate Vector Operations

Molly.jl provides comprehensive coordinate manipulation functions that respect periodic boundary conditions. These operations are fundamental for distance calculations, neighbor finding, and integration schemes.

### Vector Calculation

The `vector` function computes the minimum image vector between two coordinates, automatically applying periodic boundary conditions to find the shortest path across periodic images.

Sources: [src/spatial.jl](src/spatial.jl#L12-L13), [src/spatial.jl](src/spatial.jl#L350-L460)

```julia
# Compute minimum image vector
vec_ij = vector(coords[i], coords[j], boundary)

# For specific dimensions with 1D boundaries
vec_ij_1d = vector_1D(coords[i], coords[j], boundary, dim)
```

### Coordinate Wrapping

The `wrap_coords` function ensures all coordinates remain within the primary simulation box, which is essential for analysis and visualization.

Sources: [src/spatial.jl](src/spatial.jl#L14-L15)

```julia
# Wrap all coordinates into primary box
coords = wrap_coords(coords, boundary)

# Wrap single coordinate
coord_i = wrap_coord_1D(coord_i, boundary, dim)
```

> [!TIP]
> When using TriclinicBoundary with `approx_images=true` (default), minimum image calculations are approximate but fast. For exact calculations with highly skewed boxes, set `approx_images=false` at the cost of performance.

## Accessor Functions

Molly.jl provides standardized accessor functions for retrieving atom properties, following the AtomsBase interface for ecosystem compatibility.

### Mass and Charge Access

The `mass` and `charge` functions provide consistent access to these critical properties, regardless of the underlying atom type implementation.

Sources: [src/types.jl](src/types.jl#L283-L298)

```julia
# Get atom mass
atom_mass = mass(atoms[i])

# Get atom charge  
atom_charge = charge(atoms[i])

# These also work with custom atom types
mass(custom_atom)  # Must be implemented for custom types
```

### System-Level Accessors

Convenient functions retrieve aggregated atom properties directly from the System struct.

Sources: [src/types.jl](src/types.jl#L20-L21)

```julia
# Get all masses as a vector
masses(system)  # Returns Vector of masses

# Get all charges as a vector  
charges(system)  # Returns Vector of charges
```

## Molecular Topology

The `MolecularTopology` struct encodes connectivity information, enabling efficient determination of molecular boundaries and bonded interactions. This is essential for molecular dynamics with constraints and for computing molecular properties.

Sources: [src/types.jl](src/types.jl#L360-L397)

```julia
# Create topology from bond lists
bond_is = [1, 2, 3]  # Atom indices
bond_js = [2, 3, 4]  # Bonded partners
topology = MolecularTopology(bond_is, bond_js, n_atoms)

# Access topology
atom_molecule_inds = topology.atom_molecule_inds  # Which molecule each atom belongs to
molecule_atom_counts = topology.molecule_atom_counts  # Size of each molecule
bonded_atoms = topology.bonded_atoms  # List of bonded pairs
```

The topology automatically computes connected components from bond lists, providing both atom-to-molecule mappings and the complete set of bonded pairs for efficient constraint application.

Sources: [src/types.jl](src/types.jl#L375-L397)

## Creating a Complete System

Constructing a properly configured system requires coordinating atoms, coordinates, boundaries, and auxiliary data. The following example demonstrates the complete workflow.

```julia
using Molly
using Unitful

# Define boundary
boundary = CubicBoundary(10.0u"nm")

# Create atoms
atoms = [
    Atom(index=1, atom_type=1, mass=12.011u"g/mol", charge=0.0),
    Atom(index=2, atom_type=2, mass=1.008u"g/mol", charge=0.0),
    Atom(index=3, atom_type=2, mass=1.008u"g/mol", charge=0.0),
]

# Define coordinates (centered in box)
n_atoms = length(atoms)
coords = [
    SVector(5.0, 5.0, 5.0)u"nm",
    SVector(5.1, 5.0, 5.0)u"nm", 
    SVector(4.9, 5.0, 5.0)u"nm",
]

# Define topology (bonds)
bond_is = [1, 1]
bond_js = [2, 3]
topology = MolecularTopology(bond_is, bond_js, n_atoms)

# Create atom data for analysis
atoms_data = [
    AtomData(atom_type="CT", atom_name="C", res_name="MET"),
    AtomData(atom_type="H", atom_name="H1", res_name="MET"),
    AtomData(atom_type="H", atom_name="H2", res_name="MET"),
]

# Assemble system
system = System(
    atoms=atoms,
    coords=coords,
    boundary=boundary,
    atoms_data=atoms_data,
    topology=topology,
    pairwise_inters=(LennardJonesInteraction(cutoff=1.0u"nm"),),
)
```

Sources: [src/types.jl](src/types.jl#L540-L699)

## Initialization and Randomization

Molly.jl provides utilities for generating random coordinates and velocities following physical distributions, essential for system initialization.

### Random Coordinates

The `random_coord` function generates positions uniformly within the boundary, useful for creating disordered initial configurations.

Sources: [src/spatial.jl](src/spatial.jl#L11)

```julia
# Generate random coordinates for all atoms
coords = [random_coord(boundary) for _ in 1:n_atoms]
```

### Maxwell-Boltzmann Velocities

The `random_velocities` function generates velocities following the Maxwell-Boltzmann distribution at a specified temperature, initializing the system in thermal equilibrium.

Sources: [src/spatial.jl](src/spatial.jl#L16-L19)

```julia
# Initialize velocities at 300 K
velocities = random_velocities(system, 300.0u"K")

# Modify velocities in-place
random_velocities!(system, 300.0u"K")

# Remove center-of-mass motion to prevent drift
remove_CM_motion!(system)
```

Sources: [src/spatial.jl](src/spatial.jl#L22)

## Coordinate Scaling and Manipulation

Various coordinate manipulation functions support system preparation and analysis, including scaling for density adjustments and center-of-mass calculations.

### Coordinate Scaling

The `scale_coords!` function uniformly scales coordinates around the box center, useful for adjusting system density during equilibration.

Sources: [src/spatial.jl](src/spatial.jl#L25)

```julia
# Scale coordinates by factor 1.1 (increase volume)
scale_coords!(system, 1.1)

# Compress system (decrease volume)
scale_coords!(system, 0.95)
```

### Box Center

The `box_center` function computes the geometric center of the boundary, useful for positioning atoms or calculating distances from the box center.

Sources: [src/spatial.jl](src/spatial.jl#L9)

```julia
center = box_center(boundary)  # Returns SVector{3}
```

## Integration with Unitful

Molly.jl deeply integrates with Unitful.jl for dimensional analysis and unit conversion. Physical quantities can be expressed with explicit units or as dimensionless numbers for performance-critical code.

Sources: [src/types.jl](src/types.jl#L509-L514)

```julia
# Using units (default)
system = System(
    atoms=atoms,
    coords=coords,
    boundary=boundary,
    force_units=u"kJ * mol^-1 * nm^-1",
    energy_units=u"kJ * mol^-1",
)

# Using dimensionless numbers for performance
system = System(
    atoms=atoms,
    coords=ustrip.(coords),  # Remove units
    boundary=ustrip(boundary),
    force_units=NoUnits,
    energy_units=NoUnits,
)
```

> [!TIP]
> When switching between unitful and dimensionless modes, ensure all quantities (coordinates, velocities, forces) are consistently stripped of units or all maintain units. Mixing these approaches can lead to subtle numerical errors.

## Performance Considerations

The atom and coordinate management system is designed for high-performance molecular dynamics, with several key optimizations:

### GPU Compatibility

Atom types must be bits types to work on GPU arrays. The standard `Atom` struct satisfies this requirement, while `AtomData` remains on CPU to avoid unnecessary data transfer.

Sources: [src/types.jl](src/types.jl#L596-L599)

```julia
# GPU-compatible atom array
atoms_gpu = CuArray(atoms)  # Works because Atom is bits type

# AtomsData stays on CPU
atoms_data_cpu = atoms_data  # No GPU transfer needed
```

### Coordinate Vectorization

Coordinates are stored as `SVector` objects for stack-allocated, SIMD-friendly vector operations. This enables efficient force calculations without dynamic memory allocation.

Sources: [src/types.jl](src/types.jl#L485-L486)

```julia
# Coordinates use static vectors for performance
coords = Vector{SVector{3, typeof(1.0u"nm")}}(undef, n_atoms)

# Access is fast and allocation-free
coord_i = coords[i]
x_coord = coord_i[1]
```

## Type System for Array Operations

The System type encodes array backend information in its type parameters, enabling compile-time dispatch for CPU or GPU operations.

Sources: [src/types.jl](src/types.jl#L17-L18)

```julia
# Detect array backend
array_type(system)  # Returns Vector, CuArray, etc.
is_on_gpu(system)   # Returns true if on GPU

# Detect floating point precision
float_type(system)  # Returns Float64, Float32, etc.
```

This type-level information allows kernels to be compiled specifically for the target hardware, eliminating runtime branching in performance-critical loops.

## Error Handling and Validation

The System constructor performs comprehensive validation to ensure consistency between atom counts, coordinates, velocities, and auxiliary data.

Sources: [src/types.jl](src/types.jl#L586-L599)

```julia
# These checks are performed automatically:
# - atoms and coords must have same length
# - atoms and velocities must have same length  
# - atoms_data must match atoms length if non-empty
# - GPU atoms must be bits types
```

Explicit error messages guide users to resolve configuration issues before simulation begins, preventing subtle runtime errors.

## Next Steps

Understanding atoms and coordinates is foundational for molecular dynamics in Molly.jl. Continue exploring the system architecture with these related topics:

- [Boundary Conditions and Periodic Systems](8-boundary-conditions-and-periodic-systems) - Deep dive into periodic boundary implementations
- [System Object and Data Structures](6-system-object-and-data-structures) - Complete System architecture overview
- [Unitful Integration and Physical Units](9-unitful-integration-and-physical-units) - Advanced unit handling techniques
- [Pairwise Interactions: Lennard-Jones and Coulomb](10-pairwise-interactions-lennard-jones-and-coulomb) - How atomic properties drive force calculations