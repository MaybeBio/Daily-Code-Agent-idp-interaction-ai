

Specific interactions form the foundation of molecular mechanics force fields by describing bonded interactions between atoms. Unlike pairwise interactions that act between all atom pairs, specific interactions operate on predefined subsets of atoms: bonds (2 atoms), angles (3 atoms), and torsions (4 atoms). These interactions maintain molecular geometry, reproduce vibrational spectra, and encode conformational preferences essential for accurate biomolecular simulations. Molly.jl implements a comprehensive suite of specific interaction types through a unified interface that supports both CPU and GPU execution, enabling efficient simulation of systems ranging from simple molecules to large proteins.

## Architecture Overview

The specific interaction architecture in Molly.jl follows a hierarchical type system organized by the number of participating atoms. Each interaction type implements a common interface with force and potential energy methods, enabling uniform treatment across different interaction categories. The system employs specialized container types—`InteractionList2Atoms`, `InteractionList3Atoms`, and `InteractionList4Atoms`—that store atom indices and interaction parameters, facilitating batch processing and efficient memory access patterns. This design allows the simulation engine to handle thousands of interactions efficiently while maintaining type safety and enabling automatic differentiation.

```mermaid
classDiagram
    AbstractInteraction <|-- SpecificInteractionList
    SpecificInteractionList <|-- InteractionList2Atoms
    SpecificInteractionList <|-- InteractionList3Atoms
    SpecificInteractionList <|-- InteractionList4Atoms
    
    InteractionList2Atoms "1" <-- "*" HarmonicBond
    InteractionList2Atoms "1" <-- "*" FENEBond
    InteractionList2Atoms "1" <-- "*" MorseBond
    
    InteractionList3Atoms "1" <-- "*" HarmonicAngle
    InteractionList3Atoms "1" <-- "*" CosineAngle
    InteractionList3Atoms "1" <-- "*" UreyBradley
    
    InteractionList4Atoms "1" <-- "*" PeriodicTorsion
    InteractionList4Atoms "1" <-- "*" RBTorsion
    
    class AbstractInteraction {
        <<abstract>>
    }
    
    class SpecificInteractionList {
        <<abstract>>
        +N : Int
        +is : Vector{Int32}
        +inters : Vector
        +types : Vector{String}
    }
    
    class InteractionList2Atoms {
        +is : Vector{Int32}
        +js : Vector{Int32}
        +inters : Vector
    }
    
    class InteractionList3Atoms {
        +is : Vector{Int32}
        +js : Vector{Int32}
        +ks : Vector{Int32}
        +inters : Vector
    }
    
    class InteractionList4Atoms {
        +is : Vector{Int32}
        +js : Vector{Int32}
        +ks : Vector{Int32}
        +ls : Vector{Int32}
        +inters : Vector
    }
    
    class HarmonicBond {
        +k : ForceConstant
        +r0 : Distance
    }
    
    class HarmonicAngle {
        +k : ForceConstant
        +θ0 : Angle
    }
    
    class PeriodicTorsion {
        +periodicities : NTuple{N, Int}
        +phases : NTuple{N, Float}
        +ks : NTuple{N, Energy}
        +proper : Bool
    }
```

The interaction type system distinguishes between proper torsions (dihedral angles around single bonds) and improper torsions (out-of-plane bending constraints), with `PeriodicTorsion` supporting both through the `proper` boolean field. All interaction types support parameter injection and extraction for force field customization, with the `inject_interaction` and `extract_parameters!` functions enabling dynamic parameter modification during simulation setup or replica exchange calculations. The zero and addition operators are overloaded for all interaction types, supporting automatic differentiation and parameter optimization workflows.

Sources: [src/types.jl](src/types.jl#L1-L1573), [src/interactions/harmonic_bond.jl](src/interactions/harmonic_bond.jl#L1-L56)

## Bond Interactions

Bond interactions constrain distances between covalently bonded atoms, maintaining molecular structure and governing bond vibrations. Molly.jl provides three bond types: HarmonicBond for standard covalent bonds, FENEBond for polymer chains with limited extensibility, and MorseBond for bonds with realistic dissociation behavior. Each bond type implements both force and potential energy functions that operate on atom coordinates while respecting periodic boundary conditions through the vector function for minimum image calculations.

### Harmonic Bond

The HarmonicBond potential represents the most common bond interaction, modeling bond stretching as a simple harmonic oscillator around an equilibrium distance. The potential energy follows the classical formulation V(r) = (1/2)k(r - r₀)², where k is the force constant and r₀ is the equilibrium bond length. The force calculation computes the derivative analytically, yielding a force proportional to the displacement from equilibrium: f = -k(r - r₀) along the bond direction. This simple quadratic approximation works well for small bond displacements typical in equilibrium molecular dynamics simulations.

The HarmonicBond struct stores force constant k and equilibrium distance r₀ as typed parameters, enabling compile-time optimization and unit-aware calculations. The force implementation uses vector normalization to ensure proper directionality while maintaining numerical stability near equilibrium distances. Special handling for periodic boundaries ensures that bonded atoms across periodic images are treated correctly, with the vector function computing the minimum image displacement between atom positions.

Sources: [src/interactions/harmonic_bond.jl](src/interactions/harmonic_bond.jl#L1-L56), [src/spatial.jl](src/spatial.jl#L500-L600)

### FENE Bond

The FENEBond (Finitely Extensible Nonlinear Elastic) potential addresses limitations of harmonic bonds for polymer simulations by preventing unphysical bond extension. The potential combines a logarithmic entropic spring term with a WCA (Weeks-Chandler-Andersen) repulsive core, creating bonds that are soft at equilibrium but become infinitely stiff as they approach their maximum extension. This model is particularly suited for coarse-grained polymer simulations like bead-spring models, where it captures the finite extensibility of real polymer chains while preventing bond crossing.

The mathematical formulation V(r) = -(1/2)kr₀²ln(1 - (r/r₀)²) + V_WCA(r) features two components: a logarithmic term that dominates at small extensions and a purely repulsive WCA term that prevents bond dissociation at short distances. The WCA contribution switches off at r ≥ 2^(1/6)σ, creating a purely repulsive interaction that maintains excluded volume. The force calculation combines derivatives of both terms, with special conditional logic to handle the WCA cutoff. The bond becomes infinitely stiff as r approaches r₀, providing natural bond length constraints without numerical instabilities.

Sources: [src/interactions/fene_bond.jl](src/interactions/farmonic_bond.jl#L1-L65)

### Morse Bond

The MorseBond potential provides a more realistic representation of bond dissociation by incorporating anharmonicity and finite bond energy. The potential V(r) = D(1 - e^(-a(r - r₀)))² captures three key physical features: the equilibrium distance r₀, the well depth D (bond dissociation energy), and the width parameter a controlling the potential steepness. Unlike harmonic bonds that have infinite energy at all displacements, Morse bonds asymptotically approach zero energy at large separations, enabling realistic bond breaking and formation in reactive molecular dynamics simulations.

The exponential form of the Morse potential enables efficient force calculation through analytic derivatives, with the force expression f = 2Da(1 - e^(-a(r - r₀)))e^(-a(r - r₀)) computed in a single exponential evaluation. This makes Morse bonds computationally efficient despite their increased realism compared to harmonic bonds. The implementation supports parameter injection and extraction for force field compatibility, with typed parameters D (energy), a (inverse length), and r₀ (length) enabling automatic differentiation for parameter optimization workflows.

Sources: [src/interactions/morse_bond.jl](src/interactions/morse_bond.jl#L1-L39)

**Bond Interaction Comparison:**

| Interaction | Potential Form | Key Parameters | Best Use Cases | Dissociation |
|-------------|---------------|----------------|----------------|--------------|
| HarmonicBond | (1/2)k(r - r₀)² | k, r₀ | Standard biomolecules | No |
| FENEBond | -(1/2)kr₀²ln(1 - (r/r₀)²) + V_WCA | k, r₀, σ, ε | Coarse-grained polymers | No |
| MorseBond | D(1 - e^(-a(r - r₀)))² | D, a, r₀ | Reactive simulations | Yes |

## Angle Interactions

Angle interactions govern the geometry of three-atom systems, controlling bond angles that determine molecular shape and flexibility. Molly.jl implements three angle types: HarmonicAngle (standard harmonic potential), CosineAngle (trigonometric potential), and UreyBradley (combined angle plus 1-3 bond interaction). Angle forces act perpendicular to bond vectors, restoring angles toward equilibrium values while conserving the total force on the system. The implementation uses 3D vector operations with cross products for directionality, with special handling for 2D systems through vector padding and trimming functions.

### Harmonic Angle

The HarmonicAngle potential extends the harmonic oscillator concept to bond angles, with potential energy V(θ) = (1/2)k(θ - θ₀)² where θ₀ is the equilibrium angle in radians and k is the force constant. The force calculation involves computing the angle derivative with respect to atomic coordinates, yielding forces perpendicular to the bond vectors that restore the angle toward equilibrium. Implementation uses the acosbound function for numerically stable angle calculation, handling degenerate cases where atoms are collinear through cross product zero checks.

The force derivation involves projecting the angle derivative onto directions perpendicular to the bonds. The code computes perpendicular vectors pa and pc using cross products, then scales these by the angle derivative divided by bond lengths. This geometric approach ensures that forces are correctly distributed among all three atoms while maintaining the angle constraint. The middle atom receives the negative sum of the outer atom forces, conserving total momentum. Periodic boundary conditions are handled through the vector_pad3D function, which extends 2D coordinates to 3D for cross product calculations while preserving minimum image conventions.

Sources: [src/interactions/harmonic_angle.jl](src/interactions/harmonic_angle.jl#L1-L69), [src/spatial.jl](src/spatial.jl#L500-L600)

### Cosine Angle

The CosineAngle potential provides an alternative harmonic formulation using trigonometric functions: V(θ) = k(1 + cos(θ - θ₀)). This formulation is mathematically equivalent to the harmonic angle for small displacements but differs in its asymptotic behavior, making it suitable for specific force field parameterizations. The force calculation uses sin(θ - θ₀) instead of (θ - θ₀), which can provide better numerical behavior for certain angle ranges. The implementation follows the same geometric approach as HarmonicAngle but substitutes the sine function for the angle deviation term.

Cosine angles are commonly used in certain force fields (including some versions of AMBER) where the trigonometric form provides better parameter transferability across different molecular environments. The force computation requires the same cross-product-based geometry, with perpendicular forces calculated using normalized vectors in the plane of the three atoms. This potential is particularly useful for maintaining specific angle preferences in ring systems or constrained geometries where the trigonometric form captures the desired energy landscape more accurately.

Sources: [src/interactions/cosine_angle.jl](src/interactions/cosine_angle.jl#L1-L43)

### Urey-Bradley Interaction

The UreyBradley potential combines a harmonic angle term with a harmonic 1-3 bond term, providing simultaneous control over both angle geometry and distance between outer atoms. The potential V(θ, r) = (1/2)k_a(θ - θ₀)² + (1/2)k_b(r - r₀)² addresses correlations between angle bending and bond stretching that are important in certain molecular systems, particularly in force fields like CHARMM. This interaction is particularly useful for maintaining proper ring conformations and describing coupled vibrations in rigid molecular fragments.

The implementation computes angle forces identically to HarmonicAngle, then adds 1-3 bond forces between atoms i and k. The bond force follows the harmonic bond formulation, creating a direct interaction between outer atoms that constrains their distance independently of the central atom position. This additional term can help reproduce experimental vibrational frequencies more accurately than pure angle potentials, especially in systems where angle-bond coupling is significant. The force summation ensures that all three atoms receive appropriate contributions from both the angle and bond terms while conserving total force.

Sources: [src/interactions/urey_bradley.jl](src/interactions/urey_bradley.jl#L1-L62)

> [!TIP]
> Urey-Bradley interactions are particularly valuable for accurate reproduction of vibrational spectra in force fields like CHARMM. When implementing custom force fields, consider whether angle-bond coupling is significant in your molecular system—if it is, UreyBradley provides a straightforward way to capture this effect without implementing complex cross terms.

## Torsion Interactions

Torsion (dihedral angle) interactions control rotation around bonds, determining conformational preferences and barriers to internal rotation. Molly.jl provides PeriodicTorsion for standard dihedral potentials and RBTorsion for Ryckaert-Bellemans torsions commonly used in liquid crystal and lipid simulations. Torsion forces involve four atoms and depend on the signed dihedral angle φ, requiring careful handling of periodic discontinuities in angle calculations. The implementation uses atan2 with cross and dot products for robust angle computation, with GPU-specific implementations addressing automatic differentiation requirements on different architectures.

### Periodic Torsion

The PeriodicTorsion potential represents the most widely used dihedral formulation in biomolecular force fields, combining multiple Fourier terms to capture complex conformational preferences: V(φ) = Σ(k_n(1 + cos(nφ - φ_s,n))) for n = 1 to N. The parameters include periodicities n (multiplicity), phases φ_s,n, and force constants k_n for each term, with the proper boolean flag distinguishing proper (single bond) from improper (out-of-plane) torsions. Multi-term torsions enable representation of energy landscapes with multiple minima, essential for accurate protein sidechain rotamer libraries and backbone φ/ψ preferences.

The force calculation computes the derivative dV/dφ = -Σ(k_n·n·sin(nφ - φ_s,n)) and projects this onto atomic coordinates through geometric relationships involving cross products of bond vectors. The implementation uses the periodic_torsion_vectors helper function to compute intermediate quantities (bond vectors, cross products, and the dihedral angle θ), then applies periodic_torsion_force to compute individual term contributions. GPU execution uses an explicit loop summation rather than zip/sum to ensure Enzyme compatibility, demonstrating platform-specific optimization for automatic differentiation frameworks. The force distribution follows standard torsion force formulas, with the central bond atoms (j and k) receiving opposing forces balanced by contributions to the outer atoms (i and l).

Sources: [src/interactions/periodic_torsion.jl](src/interactions/periodic_torsion.jl#L1-L157), [src/spatial.jl](src/spatial.jl#L950-L1100)

### Ryckaert-Bellemans Torsion

The RBTorsion potential provides an alternative formulation based on Legendre polynomials: V(φ) = (f₁(1 + cos φ) + f₂(1 - cos 2φ) + f₃(1 + cos 3φ) + f₄)/2. This formulation is commonly used in simulations of liquid crystals and lipid bilayers where the trigonometric form provides better parameter transferability across different chain environments. The four force constants f₁-f₄ allow representation of complex torsional profiles with up to three barrier maxima and minima. Unlike PeriodicTorsion which uses integer multiples of the dihedral angle, RBTorsion uses fixed cos(nφ) terms with n = 1, 2, 3.

The force calculation computes dV/dφ = (f₁sin φ - 2f₂sin 2φ + 3f₃sin 3φ)/2 using analytic derivatives of the trigonometric terms. The geometric projection onto atomic coordinates follows the same pattern as PeriodicTorsion, using cross products to compute perpendicular directions and distributing forces according to the internal coordinate derivative chain rule. This implementation is limited to 3D systems due to the requirement for torsion angle calculations involving cross products. The simpler parameter set compared to multi-term PeriodicTorsion can make force field parameterization more straightforward for certain applications, particularly in hydrocarbon chain modeling where the fixed cosine series provides adequate conformational control.

Sources: [src/interactions/rb_torsion.jl](src/interactions/rb_torsion.jl#L1-L44)

**Torsion Interaction Comparison:**

| Interaction | Potential Form | Parameters | Flexibility | Typical Applications |
|-------------|---------------|------------|------------|---------------------|
| PeriodicTorsion | Σkₙ(1 + cos(nφ - φ_s,n)) | n terms × (periodicity, phase, k) | High (multi-term) | Proteins, nucleic acids, general molecules |
| RBTorsion | (f₁(1+cosφ) + f₂(1-cos2φ) + f₃(1+cos3φ) + f₄)/2 | f₁, f₂, f₃, f₄ | Medium | Lipids, liquid crystals, hydrocarbons |

## Force Calculation and Spatial Operations

The force calculation for specific interactions relies on careful geometric operations that respect periodic boundary conditions and numerical stability requirements. All interactions use the vector function for computing displacements between atoms, which implements the minimum image convention for periodic systems. For angle and torsion calculations, intermediate vectors and cross products are computed using 3D operations, with 2D systems handled through vector_pad3D and trim3D functions that extend coordinates temporarily for cross product calculations then project back to 2D.

Angle calculations require computing the bond angle θ = arccos((b₁·b₂)/(|b₁||b₂|)) where b₁ and b₂ are bond vectors from the central atom. The acosbound function provides numerical stability near the angle limits (0° and 180°) where the dot product approaches ±1. Torsion angles use a more complex calculation involving cross products: θ = atan2(n₁·(b₂/|b₂|), n₁·n₂) where n₁ = b₁ × b₂ and n₂ = b₂ × b₃ are normals to the planes formed by adjacent bond triples. The atan2 formulation ensures correct quadrant determination and handles the periodic nature of dihedral angles (ranging from -π to π).

Force projections for angles and torsions involve computing directions perpendicular to bond vectors, typically using cross products with intermediate normal vectors. For angle forces, the perpendicular directions lie in the plane of the three atoms, while torsion force directions involve out-of-plane components. The geometric derivations ensure that forces sum to zero across all participating atoms, conserving momentum. Special handling for degenerate geometries (collinear atoms for angles, planar geometries for torsions) prevents numerical singularities by detecting zero cross products and returning zero forces.

Sources: [src/spatial.jl](src/spatial.jl#L500-L600), [src/interactions/harmonic_angle.jl](src/interactions/harmonic_angle.jl#L32-L52)

## Parameter Management and Force Field Integration

Molly.jl provides systematic parameter management through inject_interaction and extract_parameters! functions that enable dynamic parameter modification and force field integration. These functions use string-based key prefixes (e.g., "inter_HB_bondtype_k" for HarmonicBond force constants) to map between interaction types and parameter dictionaries, supporting force field parsing from external formats and parameter optimization workflows. The extract_parameters! function populates parameter dictionaries from force field objects based on atom type combinations, enabling automatic parameter assignment during system setup.

The type-based parameter system uses atom type tuples as keys (e.g., ("C", "H", "N") for angle types) to retrieve appropriate interaction parameters from force field libraries. This allows efficient parameter lookup during simulation while maintaining flexibility for different force field conventions. The parameter injection system supports per-interaction parameter overrides, enabling replica exchange simulations with scaled parameters or customized interaction strengths for specific molecular subsets. All interaction types support Base.zero and Base.+ operations, enabling gradient computation through automatic differentiation frameworks like Enzyme.

Sources: [src/interactions/harmonic_bond.jl](src/interactions/harmonic_bond.jl#L30-L56), [src/interactions/periodic_torsion.jl](src/interactions/periodic_torsion.jl#L35-L78)

> [!TIP]
> When implementing custom interactions, ensure that your force calculation correctly handles periodic boundary conditions by always using the vector function instead of direct coordinate subtraction. This is critical for maintaining physical correctness in simulations with non-infinite boundaries and prevents artifacts from atoms interacting across periodic images.

## Performance Considerations

Specific interaction performance is dominated by memory access patterns and vectorized arithmetic operations. The interaction list structure stores atom indices contiguously, enabling cache-efficient access during force calculations. For GPU execution, interaction lists are transferred to device memory with from_device and to_device functions, with specialized force implementations (e.g., force_gpu for PeriodicTorsion) addressing platform-specific automatic differentiation requirements. The use of inline (@inline) annotations on force and energy functions enables aggressive compiler optimization and kernel fusion.

Performance optimization strategies include: (1) batching force calculations across multiple interactions to maximize arithmetic intensity, (2) using typed parameters to enable compile-time constant propagation, (3) minimizing temporary allocations in force calculations by reusing vector variables, and (4) selecting appropriate interaction types for the physical system (e.g., harmonic vs. Morse bonds) based on required accuracy vs. computational cost. GPU acceleration benefits particularly from the embarrassingly parallel nature of specific interaction calculations, where each interaction can be processed independently once atom coordinates are loaded.

Sources: [src/interactions/periodic_torsion.jl](src/interactions/periodic_torsion.jl#L100-L125), [src/types.jl](src/types.jl#L100-L200)

## Next Steps

Understanding specific interactions provides the foundation for building complete molecular systems. Continue your exploration with [Custom Interaction Implementation](12-custom-interaction-implementation) to learn how to extend Molly.jl with your own interaction types. For broader context on force fields, review [Pairwise Interactions: Lennard-Jones and Coulomb](10-pairwise-interactions-lennard-jones-and-coulomb) to understand how specific and non-bonded interactions combine in complete simulations. When ready to build full systems, consult [System Object and Data Structures](6-system-object-and-data-structures) for guidance on organizing interaction lists within simulation systems.