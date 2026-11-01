

This guide explains how to implement custom interactions in Molly.jl, from simple pairwise potentials to complex multi-body specific interactions. Molly's interaction system follows a dual-path architecture: **pairwise interactions** for generalized atom pair calculations and **specific interactions** for pre-defined atom groups like bonds and angles.

## Interaction Type Architecture

Molly's interaction hierarchy uses Julia's type system to enforce structure and enable efficient dispatch. The foundation consists of abstract base types that all interactions extend:

```mermaid
classDiagram
    class AbstractInteraction {
        <<abstract>>
    }
    
    class NBodyInteraction~N~ {
        <<abstract>>
        N: Int
    }
    
    class SpecificInteractionList~N~ {
        <<abstract>>
        is: Vector{Int}
        inters: Vector
        types: Vector{String}
    }
    
    class PairwiseInteraction
    class InteractionList2Atoms
    class InteractionList3Atoms
    class InteractionList4Atoms
    
    AbstractInteraction <|-- NBodyInteraction
    AbstractInteraction <|-- SpecificInteractionList
    NBodyInteraction <|-- PairwiseInteraction
    SpecificInteractionList <|-- InteractionList2Atoms
    SpecificInteractionList <|-- InteractionList3Atoms
    SpecificInteractionList <|-- InteractionList4Atoms
    
    note for PairwiseInteraction "Alias for NBodyInteraction{2}"
```

The `AbstractInteraction` type serves as the root, with two primary branches: `NBodyInteraction{N}` for interactions involving N atoms (where N=2 is the currently supported value for pairwise), and `SpecificInteractionList{N}` for managing lists of interactions between specific atom sets [src/types.jl#L20-L42](src/types.jl#L20-L42).

> [!TIP]
> Custom pairwise interactions must subtype `PairwiseInteraction` (an alias for `NBodyInteraction{2}`). This enables automatic integration with Molly's neighbor finding, cutoff strategies, and GPU acceleration systems.

## Pairwise Interaction Implementation

Pairwise interactions calculate forces and energies between all or most atom pairs in the system. Examples include Lennard-Jones van der Waals forces and Coulomb electrostatics.

### Defining the Interaction Structure

A pairwise interaction struct should include parameters that don't depend on distance, along with optional cutoff and neighbor-finding configuration:

```julia
using Molly: PairwiseInteraction, NoCutoff

@kwdef struct MyPairwiseInteraction{C, T} <: PairwiseInteraction
    cutoff::C = NoCutoff()           # Cutoff strategy
    use_neighbors::Bool = false      # Enable neighbor list optimization
    strength::T                     # Custom interaction parameters
end

# Indicate whether this interaction uses neighbor lists
Molly.use_neighbors(inter::MyPairwiseInteraction) = inter.use_neighbors
```

The key insight is separating **distance-dependent computations** (handled in `pairwise_force`/`pairwise_pe`) from **configuration parameters** (stored in the struct) [src/interactions/lennard_jones.jl#L20-L35](src/interactions/lennard_jones.jl#L20-L35).

### Implementing Core Functions

For pairwise interactions, you must implement four functions:

1. **`force()`**: Main force calculation with all parameters
2. **`pairwise_force()`**: Distance-based force magnitude (for cutoff support)
3. **`potential_energy()`**: Main energy calculation
4. **`pairwise_pe()`**: Distance-based potential energy (for cutoff support)

Here's a complete example implementing a simple inverse-square force:

```julia
using Molly: force, pairwise_force, potential_energy, pairwise_pe, force_cutoff, pe_cutoff
using Unitful: u"kJ", u"mol^-1", u"nm", u"nm^-1"
using StaticArrays: SVector

@kwdef struct InverseSquare{C, T} <: PairwiseInteraction
    cutoff::C = NoCutoff()
    use_neighbors::Bool = false
    coefficient::T = 1.0u"kJ * mol^-1 * nm^2"
end

Molly.use_neighbors(inter::InverseSquare) = inter.use_neighbors

# Main force calculation - called during simulation
function Molly.force(inter::InverseSquare, dr, atom_i, atom_j, 
                    force_units=u"kJ * mol^-1 * nm^-1", args...)
    r = norm(dr)
    params = (inter.coefficient,)
    f = force_cutoff(inter.cutoff, inter, r, params)  # Apply cutoff strategy
    return (f / r) * dr  # Convert magnitude to vector
end

# Distance-based force magnitude - used by cutoff system
function Molly.pairwise_force(::InverseSquare, r, (coeff,))
    return coeff / r^2
end

# Main energy calculation
function Molly.potential_energy(inter::InverseSquare, dr, atom_i, atom_j,
                                energy_units=u"kJ * mol^-1", args...)
    r = norm(dr)
    params = (inter.coefficient,)
    return pe_cutoff(inter.cutoff, inter, r, params)  # Apply cutoff strategy
end

# Distance-based potential - used by cutoff system
function Molly.pairwise_pe(::InverseSquare, r, (coeff,))
    return -coeff / r
end
```

The `force_cutoff` and `pe_cutoff` functions from `cutoffs.jl` automatically handle different cutoff strategies (distance, shifted potential, shifted force, spline) by evaluating your `pairwise_force`/`pairwise_pe` functions at appropriate distances [src/cutoffs.jl#L16-L48](src/cutoffs.jl#L16-L48).

> [!TIP]
> The `pairwise_force` and `pairwise_pe` functions must NOT contain terms that depend on distance beyond what's passed in the `r` parameter. This is critical for correct cutoff switching behavior.

### Supporting Additional Features

For production-quality interactions, you should also implement:

```julia
# Zero value for gradient accumulation
function Base.zero(is::InverseSquare{C, T}) where {C, T}
    return InverseSquare(is.cutoff, is.use_neighbors, zero(T))
end

# Addition for replica exchange
function Base.:+(is1::InverseSquare, is2::InverseSquare)
    return InverseSquare(is1.cutoff, is1.use_neighbors, 
                        is1.coefficient + is2.coefficient)
end
```

## Specific Interaction Implementation

Specific interactions operate on pre-defined sets of atoms, such as bonds (2 atoms), angles (3 atoms), or torsions (4 atoms). These use `InteractionListNAtoms` structures to organize the atom indices and interaction parameters [src/types.jl#L44-L85](src/types.jl#L44-L85).

### Two-Atom Specific Interactions (Bonds)

Here's an example implementing a Morse bond potential:

```julia
using Molly: SpecificForce2Atoms, InteractionList2Atoms

@kwdef struct MorseBond{K, R, A}
    k::K   # Dissociation energy
    r0::R  # Equilibrium distance
    a::A   # Width parameter
end

# Force calculation returns SpecificForce2Atoms with forces on both atoms
function Molly.force(mb::MorseBond, coord_i, coord_j, boundary, args...)
    ab = vector(coord_i, coord_j, boundary)
    r = norm(ab)
    
    # Morse potential: V(r) = k * (1 - exp(-a*(r-r0)))^2
    # Force magnitude: F(r) = 2*k*a*(1 - exp(-a*(r-r0)))*exp(-a*(r-r0))
    exp_term = exp(-mb.a * (r - mb.r0))
    f_mag = 2 * mb.k * mb.a * (1 - exp_term) * exp_term
    
    f = f_mag * normalize(ab)
    return SpecificForce2Atoms(f, -f)  # Newton's third law
end

# Energy calculation
function Molly.potential_energy(mb::MorseBond, coord_i, coord_j, boundary, args...)
    dr = vector(coord_i, coord_j, boundary)
    r = norm(dr)
    exp_term = exp(-mb.a * (r - mb.r0))
    return mb.k * (1 - exp_term)^2
end

# Integration with system via InteractionList2Atoms
# Usage in System construction:
# specific_inter_lists = (
#     bonds = InteractionList2Atoms(
#         [1, 2, 3],           # Atom indices i
#         [2, 3, 4],           # Atom indices j
#         [MorseBond(k=..., r0=..., a=...), ...],  # Interaction parameters
#         ["C-C", "C-C", "C-H"]  # Type labels (optional)
#     ),
# )
```

The force calculation returns a `SpecificForce2Atoms` struct containing force vectors for both atoms [src/force.jl#L75-L85](src/force.jl#L75-L85). The first force applies to atom i, the second to atom j, with the sign convention ensuring Newton's third law (equal and opposite forces).

### Three-Atom Specific Interactions (Angles)

Angle potentials involve three atoms, with the middle atom (index j) as the vertex:

```julia
using Molly: SpecificForce3Atoms

@kwdef struct CosineAngle{K, N, P}
    k::K   # Force constant
    n::N   # Periodicity
    p::P   # Phase shift (radians)
end

function Molly.force(ca::CosineAngle, coords_i, coords_j, coords_k, boundary, args...)
    # Compute angle using boundary-aware vectors
    ba = vector_pad3D(coords_j, coords_i, boundary)
    bc = vector_pad3D(coords_j, coords_k, boundary)
    
    # Bond angle: θ = acos((ba·bc) / (|ba||bc|))
    θ = acosbound(dot(ba, bc) / (norm(ba) * norm(bc)))
    
    # Force magnitude derivative
    # V(θ) = k*(1 + cos(n*θ - p))
    # ∂V/∂θ = -k*n*sin(n*θ - p)
    angle_deriv = -ca.k * ca.n * sin(ca.n * θ - ca.p)
    
    # Project force onto perpendicular directions
    cross_ba_bc = ba × bc
    if iszero_value(cross_ba_bc)
        zf = zero(angle_deriv ./ ba)
        return SpecificForce3Atoms(zf, zf, zf)
    end
    
    pa = normalize(trim3D(ba × cross_ba_bc, boundary))
    pc = normalize(trim3D(-bc × cross_ba_bc, boundary))
    
    fa = (angle_deriv / norm(ba)) * pa
    fc = (angle_deriv / norm(bc)) * pc
    fb = -fa - fc  # Force balance
    
    return SpecificForce3Atoms(fa, fb, fc)
end

function Molly.potential_energy(ca::CosineAngle, coords_i, coords_j, coords_k, 
                                boundary, args...)
    θ = bond_angle(coords_i, coords_j, coords_k, boundary)
    return ca.k * (1 + cos(ca.n * θ - ca.p))
end
```

Note the use of `vector_pad3D` and `trim3D` for handling both 2D and 3D systems uniformly. The force projection uses cross products to compute perpendicular directions to the bonds [src/interactions/harmonic_angle.jl#L45-L60](src/interactions/harmonic_angle.jl#L45-L60).

### Four-Atom Specific Interactions (Torsions)

Torsion potentials involve four atoms defining a dihedral angle:

```julia
using Molly: SpecificForce4Atoms

@kwdef struct RyckertBellemansTorsion{C, K}
    c0::C
    c1::C
    c2::C
    c3::C
    c4::C
    c5::C
    k::K
end

function Molly.force(rbt::RyckertBellemansTorsion, coords_i, coords_j, 
                    coords_k, coords_l, boundary, args...)
    # Compute dihedral angle φ
    # Force distribution follows chain rule through angle derivatives
    # (Full implementation would compute φ and project onto 4 perpendicular directions)
    
    # Pseudo-code for structure:
    # φ = dihedral(coords_i, coords_j, coords_k, coords_l, boundary)
    # dV_dφ = polynomial_derivative(rbt, φ)
    # fi, fj, fk, fl = project_torsion_forces(coords, boundary, dV_dφ)
    
    return SpecificForce4Atoms(fi, fj, fk, fl)
end
```

The force distribution for torsions is mathematically complex, requiring careful handling of the chain rule through angle derivatives. See `PeriodicTorsion` in the source for a complete example [src/interactions/periodic_torsion.jl#L52-L70](src/interactions/periodic_torsion.jl#L52-L70).

## Cutoff Strategy Integration

Molly provides four built-in cutoff strategies that automatically apply to your interactions when you use `force_cutoff` and `pe_cutoff`:

| Cutoff Type | Potential Treatment | Force Treatment | Continuity |
|-------------|-------------------|-----------------|------------|
| `NoCutoff` | Full potential | Full force | Not applicable |
| `DistanceCutoff` | Zero beyond cutoff | Zero beyond cutoff | Discontinuous |
| `ShiftedPotentialCutoff` | V - V(rc) | Full force | Continuous potential |
| `ShiftedForceCutoff` | V - (r-rc)V'(rc) - V(rc) | F - F(rc) | Continuous potential & force |

The cutoff is automatically applied in your force/energy functions through the `force_cutoff`/`pe_cutoff` dispatch [src/cutoffs.jl#L50-L145](src/cutoffs.jl#L50-L145). You only need to implement the physics; the cutoff system handles the rest.

## GPU Acceleration Support

For GPU acceleration, define a `force_gpu` function that mirrors your CPU `force` implementation:

```julia
# GPU-specific force calculation (fallbacks to CPU version by default)
Molly.force_gpu(inter::MyPairwiseInteraction, dr, ai, aj, fu, sp, ci, cj, bnd, vi, vj, sn) = 
    Molly.force(inter, dr, ai, aj, fu, sp, ci, cj, bnd, vi, vj, sn)
```

The `force_gpu` functions are defined in `force.jl` with default fallbacks to the CPU implementations [src/force.jl#L56-L61](src/force.jl#L56-L61). For most interactions, the default fallback is sufficient. Custom GPU implementations are only needed for specialized optimizations.

## Integration with System

Once your interaction is defined, integrate it into a `System`:

```julia
using Molly

atoms = [Atom(charge=1.0, σ=0.3u"nm", ϵ=0.2u"kJ * mol^-1") for _ in 1:100]
coords = [SVector(0.1*i, 0.1*j, 0.1*k)u"nm" for i in 1:10, j in 1:10, k in 1:1 |> vec]
boundary = CubicBoundary(2.0u"nm")

system = System(
    atoms=atoms,
    coords=coords,
    boundary=boundary,
    pairwise_inters=(lj=LennardJones(cutoff=DistanceCutoff(1.0u"nm")),
                    custom=MyPairwiseInteraction(strength=10.0u"kJ * mol^-1 * nm^2")),
    specific_inter_lists=(bonds=InteractionList2Atoms(
        Int32[1, 2], Int32[2, 3], 
        [MorseBond(k=100.0u"kJ/mol", r0=0.15u"nm", a=10.0u"nm^-1"),
         MorseBond(k=120.0u"kJ/mol", r0=0.14u"nm", a=12.0u"nm^-1")],
    ),)
)

# The custom interactions are now automatically used in forces and energy calculations
forces(system)
potential_energy(system)
```

## Implementation Checklist

Use this checklist when implementing custom interactions:

**Pairwise Interactions:**
- [ ] Subtype `PairwiseInteraction`
- [ ] Define struct with cutoff and parameters
- [ ] Implement `use_neighbors` method
- [ ] Implement `force()` calling `force_cutoff`
- [ ] Implement `pairwise_force()` with distance-based calculation
- [ ] Implement `potential_energy()` calling `pe_cutoff`
- [ ] Implement `pairwise_pe()` with distance-based calculation
- [ ] Implement `Base.zero()` for replica exchange
- [ ] Implement `Base.:+()` for gradient accumulation

**Specific Interactions (N atoms):**
- [ ] Define parameter struct (not a subtype of interaction base)
- [ ] Implement `force()` returning `SpecificForceNAtoms`
- [ ] Implement `potential_energy()`
- [ ] Implement `Base.zero()` for parameter struct
- [ ] Implement `Base.:+()` for parameter struct
- [ ] Create `InteractionListNAtoms` with atom indices and parameters
- [ ] Use `vector`, `vector_pad3D`, `bond_angle` helper functions

**Advanced Features:**
- [ ] Implement `inject_interaction` for force field file reading
- [ ] Implement `extract_parameters!` for parameter export
- [ ] Add GPU optimization if needed via `force_gpu`
- [ ] Add special weight handling for 1-4 interactions
- [ ] Write tests verifying energy conservation and force correctness

## Testing Custom Interactions

Molly's test suite provides patterns for validating interaction implementations:

```julia
using Molly: SVector, u"nm", u"kJ", u"mol^-1"
using Test

# Setup test coordinates
c1 = SVector(1.0, 1.0, 1.0)u"nm"
c2 = SVector(1.3, 1.0, 1.0)u"nm"
a1 = Atom(charge=1.0, σ=0.3u"nm", ϵ=0.2u"kJ * mol^-1")
boundary = CubicBoundary(2.0u"nm")

# Test force calculation
inter = MyPairwiseInteraction(strength=10.0u"kJ * mol^-1 * nm^2")
dr = vector(c1, c2, boundary)
f = force(inter, dr, a1, a1)

@test isapprox(f, expected_force, atol=1e-9u"kJ * mol^-1 * nm^-1")

# Test energy calculation
pe = potential_energy(inter, dr, a1, a1)
@test isapprox(pe, expected_energy, atol=1e-9u"kJ * mol^-1")
```

The test patterns in `test/interactions.jl` demonstrate comprehensive validation including cutoff switching, special interaction weights, and energy conservation [test/interactions.jl#L1-L100](test/interactions.jl#L1-L100).

## Next Steps

After implementing your custom interaction, you may want to explore:

- **[Pairwise Interactions: Lennard-Jones and Coulomb](10-pairwise-interactions-lennard-jones-and-coulomb)** - For deeper understanding of standard pairwise implementations
- **[Specific Interactions: Bonds, Angles, and Torsions](11-specific-interactions-bonds-angles-and-torsions)** - For patterns in multi-body force calculations
- **[GPU Acceleration with CUDA](21-gpu-acceleration-with-cuda)** - For optimizing interactions for large-scale simulations
- **[Cutoff Strategies and Force Calculation](23-cutoff-strategies-and-force-calculation)** - For implementing custom cutoff schemes

The interaction system is designed to be modular and extensible. By following the patterns established in the built-in interactions and properly implementing the required interface functions, your custom interactions will integrate seamlessly with Molly's simulation, analysis, and optimization infrastructure.