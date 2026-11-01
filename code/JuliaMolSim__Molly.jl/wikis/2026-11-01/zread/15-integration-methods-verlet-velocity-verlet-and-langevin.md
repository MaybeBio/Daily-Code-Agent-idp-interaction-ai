

Integration methods are the core algorithms that advance molecular systems through time by solving Newton's equations of motion. Molly.jl provides several deterministic and stochastic integrators suited for different ensemble types and physical scenarios, each with distinct numerical properties and performance characteristics. This guide explains when and how to use Velocity Verlet, Verlet (leapfrog), and Langevin integrators for accurate molecular dynamics simulations.

## Integration Methods Overview

Molecular dynamics simulations integrate the classical equations of motion to evolve atomic positions and velocities over discrete time steps. The choice of integrator determines energy conservation properties, ensemble sampling accuracy, and numerical stability. Molly.jl's simulator architecture supports both symplectic (energy-conserving) integrators for microcanonical (NVE) simulations and stochastic integrators for canonical (NVT) and other ensembles.

The three primary integration methods covered here form a hierarchy of sophistication:

- **Velocity Verlet**: The workhorse symplectic integrator with excellent energy conservation and time-reversibility
- **Verlet (Leapfrog)**: An alternative symplectic formulation with velocities offset from positions
- **Langevin**: A stochastic integrator that samples the canonical ensemble through friction and thermal noise

## Velocity Verlet Integrator

Velocity Verlet is the most widely used symplectic integrator in molecular dynamics due to its superior numerical stability and straightforward implementation. It provides third-order accuracy in positions and second-order accuracy in velocities while preserving the Hamiltonian structure of the system.

### Algorithm and Implementation

The Velocity Verlet algorithm advances the system through two half-steps and one full-step per iteration. At each time step, positions are updated using current velocities and accelerations, then forces are recomputed at the new positions, and finally velocities are updated using the average of old and new accelerations. The implementation in Molly.jl stores forces at the current and next time steps to minimize recomputation [simulators.jl](src/simulators.jl#L152-L156).

```julia
simulator = VelocityVerlet(dt=0.002u"ps", coupling=NoCoupling(), remove_CM_motion=1)
simulate!(system, simulator, n_steps=10_000)
```

The core integration loop executes four key operations per step [simulators.jl](src/simulators.jl#L166-L198):

1. **Position update**: `sys.coords .+= sys.velocities .* sim.dt .+ accels_t .* dt_sq_div2`
2. **Force computation**: Calculate forces at new positions via `forces!()`
3. **Velocity update**: `sys.velocities .+= (accels_t .+ accels_t_dt) .* dt_div2`
4. **Coupling and constraints**: Apply thermostats, barostats, and positional constraints

### Parameter Configuration

Velocity Verlet accepts three parameters [simulators.jl](src/simulators.jl#L132-L140):

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `dt` | Time unit | Required | Integration time step (typically 0.5-2 fs for atomistic simulations) |
| `coupling` | Coupler | `NoCoupling()` | Thermostat or barostat for ensemble control |
| `remove_CM_motion` | Integer | `1` | Frequency of center-of-mass motion removal |

### Performance Characteristics

Velocity Verlet exhibits excellent energy conservation in microcanonical simulations, with energy drift typically below 0.5 kJ/mol over 10,000 steps for well-behaved systems [energy_conservation.jl](test/energy_conservation.jl#L71). The symplectic nature ensures phase space volume preservation and time-reversibility, making it ideal for long NVE simulations.

> [!TIP]
> Velocity Verlet requires only one force evaluation per time step, making it computationally efficient. The pre-computation of `dt_sq_div2` and `dt_div2` avoids redundant operations in the tight integration loop.

### Constraint Handling

Velocity Verlet fully supports holonomic constraints such as bond length constraints (SHAKE/RATTLE). When constraints are present, the integrator:
- Stores coordinate copies before position updates [simulators.jl](src/simulators.jl#L160)
- Applies position constraints after coordinate updates [simulators.jl](src/simulators.jl#L173-L174)
- Applies velocity constraints after velocity updates [simulators.jl](src/simulators.jl#L181)

## Verlet (Leapfrog) Integrator

The leapfrog Verlet integrator positions positions and velocities at alternating time points, with velocities offset by half a time step from positions. This formulation is mathematically equivalent to Velocity Verlet but has a different computational structure that can be advantageous in certain contexts.

### Algorithm and Implementation

The leapfrog scheme advances velocities by a full time step, then positions by a full time step using the new velocities. The velocity update occurs before the position update, and velocities represent the average motion between coordinate steps [simulators.jl](src/simulators.jl#L230-L282).

```julia
simulator = Verlet(dt=0.002u"ps", coupling=NoCoupling(), remove_CM_motion=1)
simulate!(system, simulator, n_steps=10_000)
```

The leapfrog integration loop proceeds as follows [simulators.jl](src/simulators.jl#L249-L280):

1. **Force and acceleration computation**: Calculate forces and derive accelerations
2. **Velocity half-step**: `sys.velocities .+= accels_t .* sim.dt`
3. **Position update**: `sys.coords .+= sys.velocities .* sim.dt`
4. **Constraint correction**: Adjust velocities based on constrained position changes
5. **Coupling and wrapping**: Apply thermostats and periodic boundary conditions

When constraints are active, velocities are recomputed from the constrained coordinate differences to ensure consistency with the constraint manifold [simulators.jl](src/simulators.jl#L263-L265).

### Parameter Configuration

Verlet uses the same parameter interface as Velocity Verlet [simulators.jl](src/simulators.jl#L220-L228):

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `dt` | Time unit | Required | Integration time step |
| `coupling` | Coupler | `NoCoupling()` | Coupling method for ensemble control |
| `remove_CM_motion` | Integer | `1` | Center-of-mass motion removal frequency |

### When to Use Leapfrog vs Velocity Verlet

The choice between leapfrog and Velocity Verlet involves practical considerations:

- **Leapfrog advantages**: Natural fit for force calculations that require velocities at the midpoint, simpler implementation when force evaluation is expensive
- **Velocity Verlet advantages**: Synchronized positions and velocities facilitate analysis and logging, direct velocity access at coordinate time points

Both methods are symplectic and exhibit identical energy conservation properties in the limit of small time steps.

> [!TIP]
> Leapfrog Verlet stores only one force array, reducing memory overhead compared to Velocity Verlet which maintains forces at both current and next time steps. This can be significant for large systems.

## Langevin Integrator

Langevin dynamics provides a route to the canonical (NVT) ensemble by coupling particles to a heat bath through friction and random forces. The Langevin integrator in Molly.jl implements the "Langevin Middle" scheme from OpenMM, which offers better numerical accuracy than naive implementations [simulators.jl](src/simulators.jl#L355-L360).

### Physical Model and Algorithm

Langevin dynamics augments Newton's equations with dissipative and stochastic terms:

```
m·dv/dt = F(r) - γ·v + √(2γkBT)·R(t)
```

where γ is the friction coefficient, kB is Boltzmann's constant, T is temperature, and R(t) is Gaussian white noise. The Molly implementation pre-computes scaling factors for efficient integration [simulators.jl](src/simulators.jl#L381-L386):

```julia
vel_scale = exp(-dt * friction)
noise_scale = sqrt(1 - vel_scale^2)
```

These factors apply velocity damping and thermal noise in a single operation, avoiding explicit friction force calculation.

### Integration Scheme

The Langevin Middle scheme applies velocity updates in three phases [simulators.jl](src/simulators.jl#L410-L430):

1. **Deterministic half-step**: `sys.velocities .+= accels_t .* sim.dt`
2. **Stochastic full-step**: Apply friction and noise `sys.velocities .= sys.velocities .* sim.vel_scale .+ noise .* sim.noise_scale`
3. **Position full-step**: `sys.coords .+= sys.velocities .* sim.dt`

This ordering preserves the symplectic structure while correctly sampling the canonical ensemble.

### Parameter Configuration

Langevin requires specification of both thermal and friction parameters [simulators.jl](src/simulators.jl#L371-L386):

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `dt` | Time unit | Required | Integration time step |
| `temperature` | Temperature | Required | Target equilibrium temperature |
| `friction` | Inverse time | Required | Friction coefficient (γ) coupling to heat bath |
| `coupling` | Coupler | `NoCoupling()` | Additional coupling methods |
| `remove_CM_motion` | Integer | `1` | Center-of-mass motion removal frequency |

```julia
simulator = Langevin(dt=0.002u"ps", temperature=300.0u"K", friction=1.0u"ps^-1")
simulate!(system, simulator, n_steps=100_000)
```

### Choosing the Friction Coefficient

The friction coefficient controls coupling strength to the thermal bath and impacts dynamics:

- **Low friction (γ < 1 ps⁻¹)**: Weak coupling, near-deterministic dynamics, slower temperature equilibration
- **Medium friction (γ = 1-10 ps⁻¹)**: Balanced dynamics and temperature control, typical range for biomolecular simulations
- **High friction (γ > 10 ps⁻¹)**: Strong coupling, overdamped regime, rapid temperature equilibration but suppressed dynamics

For protein simulations, friction coefficients of 1-5 ps⁻¹ are common [examples.md](docs/src/examples.md#L38-L43).

### Constraint Handling

Unlike deterministic integrators, Langevin dynamics with constraints requires careful velocity updates to preserve the canonical ensemble distribution. Molly.jl applies velocity constraints after the deterministic half-step to ensure constraints are satisfied before the stochastic update [simulators.jl](src/simulators.jl#L416).

## Advanced Integration Schemes

### LangevinSplitting for Custom Schemes

For researchers requiring fine-grained control over stochastic integration, Molly.jl provides the LangevinSplitting simulator that allows custom operator splitting schemes [simulators.jl](src/simulators.jl#L449-L560). This framework implements general splitting schemes composed of three operators:

- **A step**: Position update `r → r + v·dt`
- **B step**: Velocity update from potential forces `v → v + (F/m)·dt`
- **O step**: Velocity update from thermal bath `v → α·v + σ·noise`

Different splitting strings define the integration scheme:

| Splitting | Name | Description |
|-----------|------|-------------|
| `"BAOA"` | Langevin Middle | Default Langevin scheme |
| `"BAB"` | Velocity Verlet | Reduces to deterministic Velocity Verlet |
| `"ABOBA"` | BAOAB variant | Alternative ordering |

The splitting framework computes effective friction and noise scales per O-step to ensure correct ensemble sampling [simulators.jl](src/simulators.jl#L499-L501):

```julia
α_eff = exp.(-sim.friction * sim.dt .* M_inv / count('O', sim.splitting))
σ_eff = sqrt.((1 * unit(eltype(α_eff))) .- (α_eff .^ 2))
```

> [!TIP]
> Custom splitting schemes can improve accuracy or phase space volume preservation for specific systems. However, LangevinSplitting does not support constraints due to the complexity of enforcing constraints across arbitrary operator orderings.

### Overdamped Langevin for Brownian Dynamics

For systems where inertia is negligible compared to viscous forces (e.g., colloids, large biomolecules in viscous solvent), the OverdampedLangevin simulator implements Brownian dynamics using the Euler-Maruyama method [simulators.jl](src/simulators.jl#L585-L649).

```julia
simulator = OverdampedLangevin(dt=0.01u"ps", temperature=300.0u"K", friction=100.0u"ps^-1")
simulate!(system, simulator, n_steps=1_000_000)
```

The overdamped algorithm updates positions directly from forces and noise [simulators.jl](src/simulators.jl#L634-L635):

```julia
sys.coords .+= (accels_t ./ sim.friction) .* sim.dt .+ noise_prefac .* noise
```

This bypasses velocity integration entirely, which can be more efficient for highly damped systems.

## Integration Method Comparison

Understanding the strengths and limitations of each integrator is essential for selecting the appropriate method for your simulation needs.

| Property | Velocity Verlet | Verlet (Leapfrog) | Langevin | Overdamped Langevin |
|----------|----------------|-------------------|----------|---------------------|
| **Ensemble** | NVE (+ thermostats) | NVE (+ thermostats) | NVT | NVT |
| **Order** | O(dt²) | O(dt²) | O(dt²) | O(dt¹) |
| **Energy conservation** | Excellent | Excellent | N/A | N/A |
| **Time-reversible** | Yes | Yes | No | No |
| **Constraints supported** | Yes | Yes | Yes | No |
| **Force evaluations/step** | 1 | 1 | 1 | 1 |
| **Memory overhead** | Medium (2 force arrays) | Low (1 force array) | Low (1 force array) | Low (1 force array) |
| **Typical use case** | Standard MD | Alternative to VV | Constant T simulations | Brownian dynamics |

### Energy Conservation Performance

Velocity Verlet demonstrates superior energy conservation in microcanonical tests. Energy conservation benchmarks show maximum energy differences below 0.0005 kJ/mol over 10,000 steps for Lennard-Jones fluids [energy_conservation.jl](test/energy_conservation.jl#L67-L71):

```julia
max_ΔE = maximum(abs.(Es .- E0))
@test max_ΔE < 5e-4u"kJ * mol^-1"
```

This performance holds across multiple cutoff schemes (distance, shifted potential, shifted force, cubic spline) and hardware platforms (CPU multithreaded, CUDA, AMDGPU) [energy_conservation.jl](test/energy_conservation.jl#L21-L26).

## Practical Usage Patterns

### Running Energy Conservation Simulations

For microcanonical simulations testing energy conservation:

```julia
simulator = VelocityVerlet(dt=0.002u"ps", remove_CM_motion=false)
sys.loggers = (energy=TotalEnergyLogger(100),)
simulate!(sys, simulator, n_steps=10_000)
```

Setting `remove_CM_motion=false` prevents artificial energy drift from center-of-mass motion corrections [energy_conservation.jl](test/energy_conservation.jl#L17).

### Temperature-Controlled Simulations

Combine Velocity Verlet with a thermostat for NVT simulations:

```julia
simulator = VelocityVerlet(
    dt=0.002u"ps",
    coupling=VelocityRescaleThermostat(300.0u"K", 1.0u"ps"),
    remove_CM_motion=10
)
simulate!(sys, simulator, n_steps=100_000)
```

The VelocityRescaleThermostat provides canonical ensemble sampling with correct dynamics properties [coupling.jl](src/coupling.jl#L75-L95).

### Simulated Annealing with Langevin

For temperature-dependent processes like simulated annealing:

```julia
temps = [300.0, 200.0, 100.0, 0.0]u"K"
for temp in temps
    simulator = Langevin(
        dt=0.001u"ps",
        temperature=temp,
        friction=1.0u"ps^-1",
    )
    simulate!(sys, simulator, 5_000; run_loggers=:skipzero)
end
```

This pattern gradually reduces temperature to anneal into low-energy configurations [examples.md](docs/src/examples.md#L34-L44).

### Large-Scale Production Simulations

For large-scale simulations, leverage multithreading and GPU acceleration:

```julia
using CUDA

simulator = VelocityVerlet(dt=0.002u"ps")
n_threads = Threads.nthreads()
simulate!(sys, simulator, n_steps=1_000_000; n_threads=n_threads)
```

The integration loop supports parallel force computation via the `n_threads` parameter [simulators.jl](src/simulators.jl#L68-L70).

## Performance Considerations

### Time Step Selection

The maximum stable time step depends on the highest frequency motions in the system:

- **Atomistic MD with bonds constrained**: 2.0-2.5 fs
- **All-atom MD without constraints**: 0.5-1.0 fs
- **Coarse-grained systems**: 10-50 fs
- **Implicit solvent**: 2.0-5.0 fs

Larger time steps require more frequent neighbor list updates and may reduce accuracy for fast degrees of freedom.

### Coupling to Thermostats and Barostats

Integrators interact with coupling methods through the `apply_coupling!` interface. Thermostats modify velocities to control temperature, while barostats may adjust coordinates and box vectors to control pressure [coupling.jl](src/coupling.jl#L14-L36).

The integration loop checks whether coupling methods invalidate forces and recomputes forces if necessary [simulators.jl](src/simulators.jl#L186-L194):

```julia
recompute_forces = apply_coupling!(sys, buffers, sim.coupling, sim, neighbors, step_n;
                                   n_threads=n_threads, rng=rng)
if recompute_forces
    forces!(forces_t_dt, sys, neighbors, buffers, Val(needs_vir), step_n; n_threads=n_threads)
    forces_t .= forces_t_dt
    accels_t .= forces_t ./ masses(sys)
end
```

### Memory and Force Calculation Overhead

Velocity Verlet maintains two force arrays (`forces_t` and `forces_t_dt`) for efficient acceleration computation [simulators.jl](src/simulators.jl#L152). For memory-constrained large systems, the leapfrog Verlet or Langevin integrators may be preferable as they require only one force array.

## Integration with Analysis and Logging

Molly.jl's logger system works seamlessly with all integration methods. For energy conservation analysis:

```julia
sys.loggers = (
    coords=CoordinatesLogger(100),
    energy=TotalEnergyLogger(100),
    temp=TemperatureLogger(100),
)
simulate!(sys, simulator, n_steps=10_000)
```

Loggers are called at specified intervals, and the `:skipzero` option avoids logging the initial state when running sequential simulations with different parameters [simulators.jl](src/simulators.jl#L60-L62).

## Common Issues and Troubleshooting

### Energy Drift in Microcanonical Simulations

If you observe significant energy drift with Velocity Verlet:

1. Check for periodic boundary condition wrapping errors
2. Verify neighbor list cutoffs are appropriate for the system
3. Ensure constraints are correctly applied
4. Reduce the time step if high-frequency modes exist

### Temperature Not Reaching Target in Langevin

For Langevin simulations that don't equilibrate to the target temperature:

1. Increase the friction coefficient for stronger coupling
2. Verify temperature logger is sampling correctly after equilibration
3. Check for artifacts from constraints that may affect dynamics

### Instability with Large Time Steps

If simulations become unstable:

1. Reduce the time step by 50%
2. Verify constraint algorithms are compatible with the integrator
3. Check for stiff interactions (e.g., harmonic bonds with high force constants)

## Next Steps

Understanding integration methods is foundational for molecular dynamics simulations. Continue your learning journey with:

- **[Thermostats: Andersen, Berendsen, and Velocity Rescaling](16-thermostats-andersen-berendsen-and-velocity-rescaling)**: Detailed exploration of temperature coupling methods
- **[Barostats and Pressure Control](17-barostats-and-pressure-control)**: Pressure control techniques for isothermal-isobaric simulations
- **[Energy Minimization with Steepest Descent](18-energy-minimization-with-steepest-descent)**: System preparation before dynamics

For comprehensive integration with Molly's simulation framework, explore the complete [API documentation](api.md) and additional examples in the [examples gallery](examples.md).