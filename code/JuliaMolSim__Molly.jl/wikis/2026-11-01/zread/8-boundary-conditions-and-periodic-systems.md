

Boundary conditions define the spatial environment in which molecular dynamics simulations occur, determining how atoms interact when they approach system boundaries. Molly.jl provides flexible support for periodic boundary conditions across multiple dimensionality and box geometries, essential for simulating bulk materials, liquids, and biomolecules while minimizing surface effects.

## Architectural Overview

Molly.jl's boundary condition system centers on three core boundary types derived from `AbstractBoundary`, each implementing the AtomsBase.jl interface for ecosystem compatibility. The architecture separates boundary representation from coordinate operations, enabling clean composition with simulation algorithms and neighbor finding strategies.

```mermaid
classDiagram
    AbstractBoundary{D,T,C} <|-- CubicBoundary
    AbstractBoundary{D,T,C} <|-- RectangularBoundary
    AbstractBoundary{D,T,C} <|-- TriclinicBoundary
    
    class AbstractBoundary{D,T,C} {
        <<abstract>>
        D: Dimension
        T: Float type
        C: Coordinate type
        +volume()
        +vector(c1, c2)
        +wrap_coords(coords)
    }
    
    class CubicBoundary {
        +side_lengths: SVector{3, C}
        +vector_1D(c1, c2, side_length)
    }
    
    class RectangularBoundary {
        +side_lengths: SVector{2, C}
        +vector_1D(c1, c2, side_length)
    }
    
    class TriclinicBoundary {
        +basis_vectors: SVector{3, SVector{3, C}}
        +α, β, γ: T
        +reciprocal_size: SVector{3}
        +approx_images: Bool
    }
    
    System --> AbstractBoundary
    NeighborFinder --> AbstractBoundary
```

## Boundary Types and Construction

### Cubic Boundary

The `CubicBoundary` defines a three-dimensional orthogonal box with independent side lengths, supporting both finite and infinite dimensions [src/spatial.jl#L30-L60](src/spatial.jl#L30-L60). This is the most commonly used boundary type for bulk simulations and is compatible with all neighbor finding algorithms.

```julia
# Equal side lengths
boundary = CubicBoundary(2.0u"nm")  # 2nm × 2nm × 2nm cube

# Different side lengths
boundary = CubicBoundary(2.0u"nm", 3.0u"nm", 4.0u"nm")

# Infinite in z-direction (slab geometry)
boundary = CubicBoundary(2.0u"nm", 2.0u"nm", Inf)

# From box matrix
boundary = CubicBoundary(SMatrix{3,3}([2.0 0 0; 0 3.0 0; 0 0 4.0]))
```

> [!TIP]
> Infinite dimensions are implemented using Julia's `Inf` type and prevent periodic wrapping in that direction. This is useful for interface simulations where one dimension should remain non-periodic.

### Rectangular Boundary

The `RectangularBoundary` provides a two-dimensional orthogonal box for planar simulations [src/spatial.jl#L62-L95](src/spatial.jl#L62-L95). This boundary type automatically extends to three dimensions when padded for vector operations, with the third dimension treated as infinite.

```julia
# Square boundary
boundary = RectangularBoundary(3.0u"nm")  # 3nm × 3nm

# Rectangular boundary
boundary = RectangularBoundary(3.0u"nm", 4.0u"nm")

# Infinite dimension
boundary = RectangularBoundary(3.0u"nm", Inf)
```

### Triclinic Boundary

The `TriclinicBoundary` implements a general parallelepiped defined by three basis vectors [src/spatial.jl#L97-L200](src/spatial.jl#L97-L200). This boundary type is essential for crystallographic systems, anisotropic materials, and simulations requiring non-orthogonal periodic boundaries. The first basis vector must align with the x-axis, and the second must lie in the xy-plane.

```julia
# From basis vectors
v1 = SVector(2.2, 0.0, 0.0)u"nm"
v2 = SVector(1.0, 1.7, 0.0)u"nm"
v3 = SVector(1.4, 0.5, 1.0)u"nm"
boundary = TriclinicBoundary(v1, v2, v3)

# From lengths and angles (α, β, γ in radians)
lengths = SVector(2.2, 2.0, 1.8)u"nm"
angles = deg2rad.(SVector(50.0, 40.0, 60.0))
boundary = TriclinicBoundary(lengths, angles)

# With exact image finding (slower but more accurate)
boundary = TriclinicBoundary(lengths, angles; approx_images=false)
```

The `approx_images` parameter controls the minimum image convention strategy. When `true` (default), an efficient approximation is used that is exact for distances shorter than half the shortest box dimension. When `false`, all 26 periodic images are explicitly checked for the closest one [src/spatial.jl#L487-L500](src/spatial.jl#L487-L500).

| Property | CubicBoundary | RectangularBoundary | TriclinicBoundary |
|----------|---------------|---------------------|-------------------|
| **Dimensions** | 3D | 2D | 3D |
| **Basis vectors** | Orthogonal | Orthogonal | General parallelepiped |
| **Infinite dimensions** | Supported | Supported | Not supported |
| **Tree neighbor finder** | Compatible | Compatible | Incompatible |
| **Construction** | Side lengths | Side lengths | Basis vectors or (lengths, angles) |

## Coordinate Operations with Periodic Boundaries

### Minimum Image Convention

The minimum image convention ensures that distances between atoms are computed to the closest periodic image, critical for force calculations in periodic systems. Molly.jl implements this through the `vector()` and `vector_1D()` functions [src/spatial.jl#L397-L485](src/spatial.jl#L397-L485).

For 1D coordinates (used internally for each axis in orthogonal boundaries):

```julia
# Displacement with periodic wrapping
v12 = vector_1D(4.0u"nm", 9.0u"nm", 10.0u"nm")  # Returns -2.0u"nm"
# The path through the boundary (9→10→0→4) is shorter than direct path
```

For multidimensional coordinates:

```julia
# 3D displacement with minimum image convention
coord1 = SVector(4.0, 1.0, 6.0)u"nm"
coord2 = SVector(6.0, 9.0, 4.0)u"nm"
boundary = CubicBoundary(SVector(10.0, 10.0, 10.0)u"nm")
dr = vector(coord1, coord2, boundary)
# Returns SVector(2.0, -2.0, -2.0)u"nm"
```

The triclinic boundary implementation uses a sequential wrapping strategy that accounts for the non-orthogonal geometry, applying transformations in the z, y, then x directions to ensure correct minimum image identification [src/spatial.jl#L465-L485](src/spatial.jl#L465-L485).

> [!TIP]
> The minimum image convention assumes the largest non-bonded interaction distance is less than half the box size. If this assumption is violated, you may encounter artifacts. In such cases, increase the box size or use a larger cutoff distance.

### Coordinate Wrapping

Coordinate wrapping ensures atoms remain within the simulation box, applied automatically during integration steps [src/spatial.jl#L540-L560](src/spatial.jl#L540-L560). This is essential for maintaining correct periodic boundary conditions over long simulations.

```julia
# Wrap individual 1D coordinate
wrapped = wrap_coord_1D(12.5u"nm", 10.0u"nm")  # Returns 2.5u"nm"

# Wrap 2D coordinates
coord = SVector(12.0, -3.0)u"nm"
boundary = RectangularBoundary(10.0u"nm", 5.0u"nm")
wrapped = wrap_coords(coord, boundary)  # Returns SVector(2.0, 2.0)u"nm"

# Wrap 3D coordinates
coord = SVector(11.0, -2.0, 6.0)u"nm"
boundary = CubicBoundary(SVector(10.0, 5.0, 4.0)u"nm")
wrapped = wrap_coords(coord, boundary)  # Returns SVector(1.0, 3.0, 2.0)u"nm"
```

Simulators automatically apply coordinate wrapping at each integration step, as seen in the VelocityVerlet implementation where coordinates are wrapped after position updates [src/simulators.jl#L93-L94](src/simulators.jl#L93-L94):

```julia
sys.coords .+= sys.velocities .* sim.dt .+ accels_t .* dt_sq_div2
sys.coords .= wrap_coords.(sys.coords, (sys.boundary,))
```

### Coordinate Unwrapping

For trajectory analysis and visualization, you may need to reconstruct continuous paths that cross periodic boundaries. The `unwrap_global()` function uses molecular topology and neighbor information to produce unwrapped coordinates [src/spatial.jl#L600-L650](src/spatial.jl#L600-L650):

```julia
# Unwrap coordinates using bond topology
unwrapped = unwrap_global(system.coords, system.boundary, system.topology)

# Include neighbor list for non-bonded continuity
neighbors = find_neighbors(system)
unwrapped = unwrap_global(system.coords, system.boundary, 
                          system.topology; neighbors=neighbors)
```

## Boundary Properties and Analysis

### Geometric Properties

Molly.jl provides functions for computing key geometric properties of boundary boxes [src/spatial.jl#L319-L405](src/spatial.jl#L319-L405):

```julia
boundary = CubicBoundary(2.0u"nm", 3.0u"nm", 4.0u"nm")

# Volume (or area for 2D boundaries)
vol = volume(boundary)  # Returns 24.0u"nm^3"

# Box center
center = box_center(boundary)  # Returns SVector(1.0, 1.5, 2.0)u"nm"

# Density (requires system with mass information)
sys = System(atoms=atoms, coords=coords, boundary=boundary)
rho = density(sys)  # Returns density in kg/m³

# Scale boundary (isotropic)
scaled = scale_boundary(boundary, 1.1)  # 10% larger in all dimensions

# Scale boundary (anisotropic)
scaled = scale_boundary(boundary, SVector(1.1, 1.2, 1.3))
```

For infinite boundaries, `volume()` returns `Inf` and `density()` returns `0`, reflecting that these quantities are undefined in non-periodic directions [src/spatial.jl#L330-L355](src/spatial.jl#L330-L355).

### Box Matrix and Cell Vectors

The `boxmatrix()` function provides a 3×3 matrix representation of the boundary vectors, useful for interfacing with external tools and analysis [src/spatial.jl#L71-L79](src/spatial.jl#L71-L79):

```julia
# Cubic boundary
boundary = CubicBoundary(2.0u"nm", 3.0u"nm", 4.0u"nm")
box_matrix = boxmatrix(boundary)
# Returns SMatrix{3,3}([2.0 0 0; 0 3.0 0; 0 0 4.0]u"nm")

# Triclinic boundary
boundary = TriclinicBoundary(lengths, angles)
box_matrix = boxmatrix(boundary)
# Returns matrix where columns are basis vectors

# AtomsBase interface
cell_vectors = AtomsBase.cell_vectors(boundary)
# Returns tuple of basis vectors
```

The invert_box_vectors function computes the reciprocal box vectors, essential for converting between Cartesian and fractional coordinates [src/spatial.jl#L293-L315](src/spatial.jl#L293-L315):

```julia
reciprocal = invert_box_vectors(boundary)
# Used internally for coordinate transformations
```

## Integration with Simulation Components

### System Integration

The `System` struct stores the boundary as a required component [src/types.jl#L444-L500](src/types.jl#L444-L500):

```julia
system = System(
    atoms=atoms,
    coords=coords,
    boundary=CubicBoundary(2.0u"nm"),
    velocities=velocities,
    # ... other parameters
)

# Access boundary
boundary = system.boundary

# Check for infinite boundaries
has_infinite = has_infinite_boundary(system)
n_infinite = n_infinite_dims(system)
```

### Neighbor Finding

All neighbor finding algorithms in Molly.jl use the boundary for distance calculations [src/neighbors.jl#L60-L110](src/neighbors.jl#L60-L110):

```julia
# Distance-based neighbor finder
finder = DistanceNeighborFinder(
    eligible=eligible_matrix,
    dist_cutoff=1.0u"nm",
    special=special_matrix,
    n_steps=10
)

# The find_neighbors function uses vector() with boundary
neighbors = find_neighbors(system, finder)
# Internally computes: vector(coord_i, coord_j, system.boundary)
```

The `TreeNeighborFinder` explicitly requires finite boundaries and does not support `TriclinicBoundary` [src/neighbors.jl#L180-L185](src/neighbors.jl#L180-L185):

```julia
# Tree neighbor finder (CPU, optimized for large systems)
finder = TreeNeighborFinder(
    eligible=eligible_matrix,
    dist_cutoff=1.0u"nm"
)
# Only works with CubicBoundary or RectangularBoundary (finite)
```

### Force Calculations

Pairwise interactions use the boundary through the `vector()` function for computing interatomic distances and directions. For example, the Lennard-Jones interaction [src/interactions/lennard_jones.jl#L60-L80](src/interactions/lennard_jones.jl#L60-L80):

```julia
# In force calculation kernels
dr = vector(coords[i], coords[j], boundary)
r = norm(dr)
f = pairwise_force(interaction, r, params)
```

This ensures that all force calculations respect periodic boundary conditions automatically, regardless of the interaction type.

## Practical Examples

### Setting Up a Periodic Liquid Simulation

```julia
using Molly

# Define cubic boundary for liquid simulation
boundary = CubicBoundary(3.0u"nm")  # 3nm cube

# Create water molecules (simplified)
n_water = 100
coords = [random_coord(boundary) for _ in 1:n_water]
velocities = random_velocities(coords, 300u"K")
atoms = [Atom(atom_type=1, mass=18.0u"g/mol", 
              charge=0.0u"e_au", σ=0.3u"nm", ϵ=0.65u"kJ * mol^-1") 
         for _ in 1:n_water]

# Define interactions with cutoff
lj = LennardJones(use_neighbors=true, 
                  cutoff=DistanceCutoff(1.0u"nm"))

# Create system
system = System(
    atoms=atoms,
    coords=coords,
    boundary=boundary,
    velocities=velocities,
    pairwise_inters=(lj,),
    neighbor_finder=DistanceNeighborFinder(
        eligible=trues(n_water, n_water),
        dist_cutoff=1.2u"nm",
        n_steps=10
    )
)

# Run simulation with periodic boundary conditions
simulator = VelocityVerlet(dt=0.002u"ps", coupling=AndersenThermostat(300u"K"))
simulate!(system, simulator, 1000)
```

### Slab Geometry with Mixed Boundary Conditions

```julia
# Slab with periodic in x,y but infinite in z
boundary = CubicBoundary(5.0u"nm", 5.0u"nm", Inf)

# Verify infinite dimension
println("Has infinite boundary: ", has_infinite_boundary(system))
println("Number of infinite dimensions: ", n_infinite_dims(system))

# This is useful for interface simulations
# Atoms can move freely in z-direction without wrapping
```

### Triclinic Crystal Simulation

```julia
# Triclinic box for crystal simulation
a, b, c = 2.5u"nm", 3.0u"nm", 3.5u"nm"
α, β, γ = 70°, 80°, 90°
boundary = TriclinicBoundary(SVector(a, b, c), SVector(α, β, γ))

# Verify geometry
println("Volume: ", volume(boundary))
println("Center: ", box_center(boundary))
println("Basis vectors: ", AtomsBase.cell_vectors(boundary))

# Create crystal coordinates (simplified)
coords = [SVector(x, y, z) .* boundary.basis_vectors 
          for x in 0:0.2:1, y in 0:0.2:1, z in 0:0.2:1]

system = System(
    atoms=atoms,
    coords=coords,
    boundary=boundary,
    # ... interactions
)
```

## Best Practices and Considerations

### Box Size Selection

When selecting boundary dimensions, ensure the box is large enough to avoid artifacts from periodic images:

- The cutoff distance should be ≤ half the smallest box dimension
- For liquid simulations, use at least 3× the correlation length of the property of interest
- For crystal simulations, use integer multiples of the unit cell dimensions

```julia
# Good practice: verify cutoff constraint
cutoff = 1.0u"nm"
min_box_side = minimum(boundary.side_lengths)
@assert cutoff < min_box_side / 2 "Cutoff too large for periodic box"
```

### Temperature and Pressure Control

Boundary conditions interact with temperature and pressure coupling algorithms:

```julia
# Berendsen barostat scales the boundary
barostat = BerendsenBarostat(1u"bar", 1.0u"ps", 
                             thermostat=AndersenThermostat(300u"K"))

# The barostat automatically scales coordinates and boundary
simulator = VelocityVerlet(dt=0.002u"ps", coupling=barostat)
```

### GPU Considerations

When using GPU acceleration, ensure boundary types are compatible with your chosen neighbor finder [src/neighbors.jl#L25-L55](src/neighbors.jl#L25-L55):

```julia
using CUDA

# Transfer system to GPU
atoms_gpu = cu(atoms)
coords_gpu = cu(coords)

# GPUNeighborFinder is recommended for NVIDIA GPUs
finder = GPUNeighborFinder(
    eligible=cu(trues(n_atoms, n_atoms)),
    dist_cutoff=1.2u"nm",
    special=cu(zeros(Bool, n_atoms, n_atoms))
)

# Works with any boundary type on GPU
system_gpu = System(
    atoms=atoms_gpu,
    coords=coords_gpu,
    boundary=boundary,  # Boundary stays on CPU
    neighbor_finder=finder
)
```

### Visualization

When visualizing periodic systems, consider whether to show periodic images:

```julia
using MollyGLMakie

# Visualize with periodic images
visualize(system, boundary_lines=true)

# Visualize single box
visualize(system, show_box=true, boundary_lines=false)
```

The `bounding_box_lines()` function provides line coordinates for drawing box outlines in 3D plots [src/spatial.jl#L368-L395](src/spatial.jl#L368-L395).

## Next Steps

Understanding boundary conditions is essential for building correct molecular simulations. For comprehensive simulation setup, continue to:

- **[Unitful Integration and Physical Units](9-unitful-integration-and-physical-units)** - Learn about physical units and unit handling in Molly.jl
- **[Pairwise Interactions: Lennard-Jones and Coulomb](10-pairwise-interactions-lennard-jones-and-coulomb)** - Understand how forces are computed across periodic boundaries
- **[Integration Methods: Verlet, Velocity Verlet, and Langevin](15-integration-methods-verlet-velocity-verlet-and-langevin)** - See how boundary conditions are applied during time integration
- **[Neighbor Finding Algorithms and Cell Lists](20-neighbor-finding-algorithms-and-cell-lists)** - Explore efficient neighbor search in periodic systems