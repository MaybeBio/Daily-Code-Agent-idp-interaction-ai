

Energy minimization is a critical preparatory step in molecular dynamics simulations that removes high-energy configurations and atomic clashes from initial structures before productive dynamics can begin. The steepest descent algorithm implemented in Molly provides a robust first-order optimization method suitable for relieving strained molecular geometries, resolving bad contacts in complex systems like proteins, and preparing systems for subsequent MD integration schemes.

## Algorithm Architecture

The `SteepestDescentMinimizer` implements an adaptive gradient descent algorithm that navigates the potential energy landscape by moving atoms in the direction of the negative force gradient. The algorithm balances step size adaptation with energy-based acceptance criteria to efficiently locate local minima while avoiding numerical instabilities.

```mermaid
flowchart TD
    A[Initialize System] --> B[Wrap Coordinates]
    B --> C[Find Neighbors]
    C --> D[Calculate Initial Energy]
    D --> E[Initialize Step Size h₀]
    
    E --> F{Step < max_steps?}
    F -->|No| G[Return Minimized System]
    
    F -->|Yes| H[Calculate Forces]
    H --> I[Normalize Forces by Max]
    I --> J[Propose New Coordinates<br/>coords' = coords + h × F/|F|max]
    J --> K[Apply Constraints if Present]
    K --> L[Wrap Coordinates]
    L --> M[Calculate Trial Energy E']
    
    M --> N{E' < E?}
    N -->|Yes| O[Accept Step<br/>h ← 6h/5<br/>E ← E']
    N -->|No| P[Reject Step<br/>Restore coords<br/>h ← h/5]
    
    O --> Q{max force < tol?}
    P --> Q
    
    Q -->|Yes| G
    Q -->|No| F
```

### Core Implementation

The minimizer structure encapsulates four configurable parameters that control optimization behavior:

Sources: [src/simulators.jl](src/simulators.jl#L20-L35)

| Parameter | Type | Default | Purpose |
|-----------|------|---------|---------|
| `step_size` | Distance units | 0.01 nm | Initial maximum displacement per iteration |
| `max_steps` | Integer | 1,000 | Upper bound on optimization iterations |
| `tol` | Force units | 1,000 kJ mol⁻¹ nm⁻¹ | Convergence threshold based on maximum force |
| `log_stream` | IO stream | `devnull` | Output destination for progress logging |

The algorithm operates through a loop that iteratively refines atomic coordinates by:

Sources: [src/simulators.jl](src/simulators.jl#L63-L119)

1. **Force Computation**: Calculates the complete force vector on all atoms using the `forces!` function, which aggregates contributions from pairwise interactions (Lennard-Jones, Coulomb), specific interactions (bonds, angles, torsions), and general interactions
2. **Coordinate Proposals**: Proposes new coordinates by displacing atoms along the normalized force direction scaled by the current step size: `coords' = coords + h × F/|F|max`
3. **Energy Evaluation**: Computes the trial potential energy at the proposed configuration
4. **Adaptive Acceptance**: Accepts moves that reduce energy (increasing step size by 20%) and rejects moves that increase energy (decreasing step size by 80% and restoring previous coordinates)
5. **Convergence Check**: Terminates when the maximum force magnitude falls below the tolerance threshold

The adaptive step size mechanism is crucial for efficiency—it aggressively explores promising directions while being conservative when encountering energy barriers. This prevents oscillation while maintaining reasonable progress toward minima.

Sources: [src/simulators.jl](src/simulators.jl#L95-L106)

> [!TIP]
> The adaptive step size strategy (multiply by 6/5 on success, divide by 5 on failure) provides robustness across diverse molecular systems. Larger steps accelerate descent in smooth energy basins, while aggressive reduction prevents overshooting sharp potential wells common in bonded interactions.

## System Integration

The steepest descent minimizer integrates deeply with Molly's system architecture, leveraging existing force calculation infrastructure and neighbor finding algorithms.

### Force and Energy Calculations

The minimizer uses the same force calculation pipeline as MD simulations, ensuring consistent behavior across simulation types:

Sources: [src/force.jl](src/force.jl#L200-L260)

```julia
forces!(fs, sys, neighbors, buffers, Val(needs_vir), step_n; n_threads=n_threads)
```

This function computes forces from:
- **Pairwise interactions**: Lennard-Jones, Coulomb electrostatics, Buckingham potentials
- **Specific interactions**: Harmonic bonds, angles, periodic torsions, improper torsions
- **General interactions**: Implicit solvent models, position restraints

The energy calculation uses the `potential_energy` function which aggregates the same interaction types, providing the scalar energy landscape value used for step acceptance decisions.

Sources: [src/energy.jl](src/energy.jl#L140-L180)

### Boundary Conditions and Constraints

The minimizer respects periodic boundary conditions through coordinate wrapping:

Sources: [src/simulators.jl](src/simulators.jl#L71-L73)

```julia
sys.coords .+= hn .* F ./ max_force
sys.coords .= wrap_coords.(sys.coords, (sys.boundary,))
```

This ensures atoms crossing periodic boundaries remain properly mapped into the simulation box, preventing artificial energy discontinuities.

Constraint handling occurs after coordinate proposals but before energy evaluation:

Sources: [src/simulators.jl](src/simulators.jl#L72-L74)

```julia
using_constraints && apply_position_constraints!(sys, coords_copy; n_threads=n_threads)
```

Constraints enforce fixed distances (bonds) and fixed angles, maintaining molecular geometry during minimization. The constraint system uses specialized algorithms to satisfy geometric constraints while minimizing overall system energy.

Sources: [src/constraints/constraints.jl](src/constraints/constraints.jl#L1-L60)

> [!TIP]
> Constraints applied during minimization can lead to convergence difficulties in systems with tight geometric constraints. If minimization stalls near the tolerance threshold but atoms are still moving significantly, consider temporarily relaxing constraints or using a smaller initial step size.

### Neighbor Finding

The minimizer utilizes Molly's flexible neighbor finding infrastructure to efficiently identify atom pairs for non-bonded interactions:

Sources: [src/neighbors.jl](src/neighbors.jl#L1-L60)

Supported neighbor finders include:
- `NoNeighborFinder`: Direct N² calculation for small systems
- `DistanceNeighborFinder`: Simple distance-based cutoff lists
- `TreeNeighborFinder`: Octree-based spatial partitioning for medium systems
- `CellListMapNeighborFinder`: Cell-based lists for large systems
- `GPUNeighborFinder`: GPU-optimized algorithm using the Eastman-Pande method

The minimizer updates neighbor lists each iteration to account for atomic displacements, ensuring force calculations remain accurate as the system evolves.

Sources: [src/simulators.jl](src/simulators.jl#L76-L78)

## Usage Patterns and Best Practices

### Basic Minimization

The simplest invocation uses default parameters, suitable for most organic molecules and small proteins:

Sources: [test/minimization.jl](test/minimization.jl#L5-L15)

```julia
using Molly

sys = System(
    atoms=[Atom(σ=0.34u"nm", ϵ=0.65u"kJ * mol^-1") for i in 1:100],
    coords=random_coords(100),
    boundary=CubicBoundary(10.0u"nm"),
    pairwise_inters=(LennardJones(),),
)

minimizer = SteepestDescentMinimizer()
simulate!(sys, minimizer)
```

This configuration converges typical biomolecular systems within a few hundred iterations, automatically stopping when forces fall below the 1,000 kJ mol⁻¹ nm⁻¹ tolerance.

### Parameter Tuning for Challenging Systems

Systems with severe initial clashes or strained geometries may require parameter adjustment:

| Scenario | Recommended Settings | Rationale |
|----------|---------------------|-----------|
| Severe atomic overlaps | `step_size=0.001u"nm"`, `tol=5000.0` | Smaller initial steps prevent divergence; higher tolerance allows more iterations |
| Tightly constrained proteins | `step_size=0.005u"nm"`, `max_steps=5000` | Conservative steps respect constraints; extended iterations for convergence |
| GPU-accelerated systems | `tol=2000.0`, `log_stream=stdout` | Faster force evaluation justifies higher tolerance; logging monitors progress |

### Integration with Simulation Pipelines

Energy minimization typically precedes equilibration and production MD runs:

Sources: [docs/src/examples.md](docs/src/examples.md#L17-L28)

```julia
# Minimize to relieve bad contacts
minimizer = SteepestDescentMinimizer()
simulate!(sys, minimizer)

# Add thermal energy
random_velocities!(sys, 300.0u"K")

# Equilibrate with gentle thermostat
simulator = Langevin(
    dt=0.001u"ps",
    temperature=300.0u"K",
    friction=1.0u"ps^-1",
)
simulate!(sys, simulator, 10_000)

# Production run
simulate!(sys, simulator, 100_000)
```

This pattern ensures the system starts dynamics from a physically reasonable configuration, preventing instabilities from high initial forces.

### Monitoring Convergence

Progress logging provides insight into minimization behavior:

```julia
minimizer = SteepestDescentMinimizer(log_stream=stdout)
simulate!(sys, minimizer)
# Output:
# Step 0 - potential energy 1250.3 kJ mol^-1 - max force N/A - N/A
# Step 1 - potential energy 980.1 kJ mol^-1 - max force 150000.0 kJ mol^-1 nm^-1 - accepted
# Step 2 - potential energy 745.8 kJ mol^-1 - max force 120000.0 kJ mol^-1 nm^-1 - accepted
# ...
# Step 45 - potential energy -1240.5 kJ mol^-1 - max force 850.2 kJ mol^-1 nm^-1 - accepted
```

The log tracks energy reduction, maximum force magnitude, and step acceptance, enabling diagnosis of convergence issues.

## Advanced Considerations

### GPU Acceleration

The minimizer supports GPU computation through CUDA arrays:

Sources: [test/minimization.jl](test/minimization.jl#L33-L58)

```julia
using CUDA

sys = System(
    atoms=cu(atoms),
    coords=cu(coords),
    boundary=CubicBoundary(10.0u"nm"),
    pairwise_inters=(LennardJones(),),
    neighbor_finder=GPUNeighborFinder(; eligible, dist_cutoff=1.0u"nm"),
)

minimizer = SteepestDescentMinimizer()
simulate!(sys, minimizer)  # Executes on GPU
```

GPU acceleration provides significant speedup for large systems (>10,000 atoms) where neighbor finding and force calculations dominate runtime.

### Unitless Operation

For maximum performance or compatibility with external codes, the minimizer operates without units:

Sources: [test/minimization.jl](test/minimization.jl#L19-L34)

```julia
sys = System(
    atoms=[Atom(σ=0.34, ϵ=0.65, mass=1.0) for i in 1:100],
    coords=random_coords(100),
    boundary=CubicBoundary(10.0),
    pairwise_inters=(LennardJones(),),
    force_units=NoUnits,
    energy_units=NoUnits,
)

minimizer = SteepestDescentMinimizer(step_size=0.01, tol=1.0)
simulate!(sys, minimizer)
```

Unitless operation eliminates unit conversion overhead, beneficial for large-scale minimization campaigns.

### Interaction with Enzyme Autodifferentiation

The minimizer is compatible with automatic differentiation through the Enzyme extension, enabling gradient-based force field optimization:

Sources: [src/simulators.jl](src/simulators.jl#L50-L51)

```julia
using Enzyme

function loss(params)
    sys = build_system_with_params(params)
    minimizer = SteepestDescentMinimizer()
    simulate!(sys, minimizer)
    return potential_energy(sys)
end

grad = autodiff(loss, Active, params)
```

The `@inline` annotation on the simulate! method ensures Enzyme compatibility, allowing differentiable minimization workflows.

## Limitations and Alternatives

The steepest descent algorithm has inherent limitations that users should understand:

1. **Convergence Rate**: First-order methods converge slowly near minima due to diminishing gradients. Conjugate gradient or L-BFGS methods typically converge faster in smooth energy basins.

2. **Local Minima**: The algorithm converges to the nearest local minimum, which may not represent the global minimum for complex multi-well potentials. Basin hopping or simulated annealing can escape local minima.

3. **Constraint Conflicts**: Tight geometric constraints can create incompatible force directions, causing the algorithm to stall. Constraint relaxation or alternative minimization methods may be necessary.

4. **Saddle Point Detection**: Steepest descent naturally descends from saddle points, but cannot distinguish minima from saddle points during convergence. Eigenvalue analysis of the Hessian provides definitive classification.

For systems where steepest descent proves inadequate, consider:
- **Conjugate Gradient**: Faster convergence for smooth landscapes
- **L-BFGS**: Quasi-Newton method with memory efficiency
- **Simulated Annealing**: Global optimization for challenging energy landscapes
- **Replica Exchange**: Parallel exploration of multiple minima

## Summary

The `SteepestDescentMinimizer` provides a robust, integrated approach to energy minimization in Molly. Its adaptive step size strategy, comprehensive interaction support, and seamless integration with Molly's system architecture make it suitable for preparing diverse molecular systems for dynamics simulations. While first-order convergence limits performance near minima, the algorithm's reliability and simplicity ensure predictable behavior across a wide range of molecular systems.

For next steps in simulation workflows, consult [Integration Methods: Verlet, Velocity Verlet, and Langevin](15-integration-methods-verlet-velocity-verlet-and-langevin) to continue with dynamics after minimization, or explore [Thermostats: Andersen, Berendsen, and Velocity Rescaling](16-thermostats-andersen-berendsen-and-velocity-rescaling) for temperature control strategies.