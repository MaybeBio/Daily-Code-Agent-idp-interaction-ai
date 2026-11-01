

This documentation explores Molly.jl's capability for gradient-based force field optimization using automatic differentiation. This powerful approach enables systematic refinement of force field parameters by computing exact gradients through molecular dynamics simulations, supporting optimization of atom properties, interaction parameters, and even neural network potentials through end-to-end differentiable simulations.

## Architectural Foundation

The gradient-based optimization framework in Molly.jl builds on Enzyme's source-to-source automatic differentiation, enabling reverse and forward mode AD across complex molecular simulations. The architecture supports differentiable computation through force calculations, energy evaluations, and time integration, even in the presence of sophisticated interactions like PME electrostatic summation and implicit solvation models.

```mermaid
flowchart TB
    A[Force Field Parameters] --> B[System Construction]
    B --> C[Differentiable Simulation]
    C --> D[Loss Function Evaluation]
    D --> E[Gradient Computation via Enzyme]
    E --> F[Parameter Update]
    F --> B
    
    subgraph Parameter Types
        A1[Atom Properties<br/>σ, ϵ, charge]
        A2[Bond Parameters<br/>k, r0]
        A3[Angle Parameters<br/>k, θ0]
        A4[Torsion Parameters<br/>k, phases, periodicities]
        A5[Implicit Solvent Parameters<br/>GB coefficients]
        A6[Neural Network Weights]
    end
    
    subgrad Loss Functions
        D1[Distance-based objectives]
        D2[RMSD/dRMSD]
        D3[RDF matching]
        D4[Radius of gyration]
        D5[Energy/temperature targets]
    end
    
    A --> A1 & A2 & A3 & A4 & A5 & A6
    D --> D1 & D2 & D3 & D4 & D5
```

The core integration with Enzyme is handled through `MollyEnzymeExt.jl`, which defines augmented primal and reverse rules for key operations including FFT operations used in PME calculations. This extension marks certain functions as inactive for differentiation purposes, such as `check_units`, `random_coord`, and `neighbor_finder` initialization, ensuring efficient gradient computation while maintaining numerical stability [MollyEnzymeExt.jl](ext/MollyEnzymeExt.jl#L1-L59).

> [!TIP]
> Differentiable simulation requires running without units (`units=false` in System construction) and demands significant memory for reverse-mode AD, especially on GPUs. Consider using gradient checkpointing with Checkpointing.jl for larger simulations to manage memory requirements.

## Parameter Injection System

Molly.jl provides a systematic approach for managing force field parameters through the `inject_gradients` and `extract_parameters!` functions. This system enables separate parameter storage while maintaining differentiability through the entire simulation stack. The parameter dictionary uses structured keys that follow the pattern `"inter_<InteractionType>_<AtomType1>/<AtomType2>_<ParameterName>"` for pairwise interactions and `"inter_<InteractionType>_<ParameterName>"` for general interactions.

For biomolecular force fields, the parameter dictionary includes atom-level properties and interaction-specific parameters. Each interaction type defines its `inject_interaction` and `extract_parameters!` methods to handle parameter injection and extraction. For example, `HarmonicBond` defines injection for bond force constant `k` and equilibrium distance `r0` [harmonic_bond.jl](src/interactions/harmonic_bond.jl#L18-L35):

```julia
function inject_interaction(inter::HarmonicBond, inter_type, params_dic)
    key_prefix = "inter_HB_$(inter_type)_"
    return HarmonicBond(
        dict_get(params_dic, key_prefix * "k" , inter.k ),
        dict_get(params_dic, key_prefix * "r0", inter.r0),
    )
end
```

The `LennardJones` interaction handles the `weight_14` parameter for special 1-4 interactions, demonstrating the flexibility of the injection system [lennard_jones.jl](src/interactions/lennard_jones.jl#L55-L63). This architecture enables simultaneous optimization of hundreds of parameters across different interaction types, as demonstrated in the GB99dms force field development [differentiable.md](docs/src/differentiable.md#L423-L428).

## Pairwise Interaction Gradient Optimization

The simplest entry point to force field optimization involves tuning pairwise interaction parameters. A canonical example demonstrates optimizing the Lennard-Jones σ parameter to achieve a target minimum separation distance in a fluid simulation. The loss function computes the mean minimum separation after simulation and measures deviation from the desired distance [differentiable.md](docs/src/differentiable.md#L22-L95).

```julia
function loss(σ, coords, velocities, boundary, pairwise_inters,
              neighbor_finder, simulator, n_steps, n_atoms, atom_mass, dist_true)
    atoms = [Atom(i, 1, atom_mass, 0.0, σ, 0.2) for i in 1:n_atoms]

    sys = System(
        atoms=atoms,
        coords=coords,
        boundary=boundary,
        velocities=velocities,
        pairwise_inters=pairwise_inters,
        neighbor_finder=neighbor_finder,
        force_units=NoUnits,
        energy_units=NoUnits,
    )

    simulate!(sys, simulator, n_steps)
    mms_end = mean_min_separation(sys.coords, boundary)
    loss_val = abs(mms_end - dist_true)

    return loss_val
end
```

Gradient computation uses Enzyme's reverse mode AD with the `set_runtime_activity` function to mark active inputs for differentiation. The training loop iteratively updates the σ parameter using the computed gradients, converging toward the theoretical value that produces the desired interatomic separation [differentiable.md](docs/src/differentiable.md#L96-L156).

> [!TIP]
> When optimizing pairwise parameters, always validate gradients against finite differencing using central_fdm or similar methods, as shown in the test suite. The stochastic nature of molecular dynamics can cause gradient variability—consider averaging gradients across multiple independent simulations or initial conditions.

The test suite validates gradient correctness by comparing Enzyme-computed gradients against central finite difference approximations across various configurations including CPU/GPU execution, forward/reverse modes, and different numerical precisions [gradients.jl](test/gradients.jl#L1-L200). This extensive testing framework ensures gradient accuracy across the full range of supported interaction types and computational backends.

## Specific Interaction Gradient Optimization

Specific interactions such as bonds, angles, and torsions provide rich parameter optimization opportunities. A representative example demonstrates optimizing the equilibrium angle θ in a `HarmonicAngle` interaction to achieve a target end-to-end distance in triatomic molecules [differentiable.md](docs/src/differentiable.md#L168-L291).

The loss function constructs the `InteractionList3Angles` with the trainable angle parameter, simulates the system dynamics, and computes the loss as the absolute difference between the final end-to-end distance and the target distance. This pattern extends naturally to multiple parameters—for instance, simultaneously optimizing bond lengths, angle parameters, and torsion barriers to match reference structural data.

Differentiable simulation supports all specific interaction types implemented in Molly.jl:

| Interaction Type | Optimizable Parameters | Gradient Support | Use Cases |
|------------------|----------------------|------------------|-----------|
| HarmonicBond | k (force constant), r0 (equilibrium distance) | ✓ | Bond length refinement, stiffness tuning |
| HarmonicAngle | k, θ0 (equilibrium angle) | ✓ | Angular geometry optimization |
| PeriodicTorsion | k, phases, periodicities | ✓ | Rotational barrier fitting |
| FENEBond | k, r0, nonlinear parameters | ✓ | Polymer force field development |
| UreyBradley | k, r0 | ✓ | Nonbonded 1-3 interactions |

The test suite demonstrates optimization of bond and angle parameters in a multi-component system with pairwise interactions, bonds, angles, torsions, and even implicit solvent models [gradients.jl](test/gradients.jl#L201-L350). This comprehensive validation confirms that gradients propagate correctly through the entire force calculation pipeline including neighbor finding, cutoff handling, and thermostat coupling.

## Neural Network Potentials

Molly.jl integrates seamlessly with Flux neural networks, enabling machine-learned force fields and potential energy surfaces. Neural networks can be embedded directly into interaction functions, providing automatic gradients for network weights through the simulation. A simple example demonstrates training a network to control interatomic distances in a three-atom system [differentiable.md](docs/src/differentiable.md#L292-L422).

```julia
struct NNBonds{T}
    model::T
end

function AtomsCalculators.forces!(fs, sys, inter::NNBonds; kwargs...)
    vec_ij = vector(sys.coords[1], sys.coords[3], sys.boundary)
    dist = norm(vec_ij)
    f = inter.model([dist])[1] * normalize(vec_ij)
    fs .+= [f, zero(f), -f]
    return fs
end
```

The training loop uses Flux's optimization framework (ADAM in this example) to update network weights based on gradients computed through the full molecular dynamics simulation. The loss function targets a specific triangular geometry, demonstrating how learned potentials can encode complex geometric constraints that would be difficult to parameterize with traditional force fields.

This capability extends to more sophisticated architectures including graph neural networks that operate on molecular graphs, attention-based models for long-range interactions, and equivariant networks that respect physical symmetries. The key requirement is that custom interactions must define appropriate `zero` and `+` methods to support gradient accumulation.

## Biomolecular Force Field Optimization

Molly.jl has been used to train the GB99dms force field for implicit solvent molecular dynamics of proteins, involving differentiable simulations of one million steps with loss functions based on residue-residue distance matching to explicit solvent simulations [differentiable.md](docs/src/differentiable.md#L423-L428). This real-world application demonstrates the scalability of the gradient-based optimization framework to complex biological systems.

The protein test case shows how the parameter injection system handles the full complexity of a biomolecular force field [gradients.jl](test/gradients.jl#L351-L517):

```julia
params_dic = Dict(
    "atom_C8_σ" => 0.33996695084235345,
    "atom_C8_ϵ" => 0.4577296,
    "atom_N_ϵ" => 0.71128,
    "inter_PT_C/N/CT/C_k_1" => -0.142256,
    "inter_GB_screen_O" => 1.061039,
    # ... hundreds of parameters
)
```

The test suite validates gradients for multiple parameter types including atom properties (`atom_*_σ`, `atom_*_ϵ`), periodic torsion parameters (`inter_PT_*/k_*`), and implicit solvent parameters (`inter_GB_screen_*`) in a realistic protein system. Gradients are verified against finite difference approximations across CPU, GPU, and parallel execution modes.

## Molecular Loss Functions

Choosing an appropriate loss function is critical for successful force field optimization. The loss function defines the physical target and influences optimization convergence. Molly.jl supports various loss function patterns for different applications [differentiable.md](docs/src/differentiable.md#L434-L453):

| Loss Function Type | Formula | Applications | Implementation Notes |
|-------------------|---------|--------------|----------------------|
| Distance-based | |L d_final - d_ref|| | Simple geometry targets | Easy to implement, bounded |
| RMSD | √(Σᵢ‖rᵢ - rᵢ,ref‖² / N) | Structure matching | Requires alignment for rigid body motion |
| dRMSD | √(Σᵢⱼ‖dᵢⱼ - dᵢⱼ,ref‖²) | Distance map matching | Invariant to rotation/translation |
| RDF | ∫‖g(r) - g_ref(r)‖ dr | Phase matching | Currently requires manual implementation |
| Radius of Gyration | ‖R_g - R_g,ref‖ | Compactness control | Computationally inexpensive |
| Energy | |E - E_ref|| | Thermodynamic targets | Can be noisy in short simulations |

For intermediate observables throughout the simulation (e.g., RDF, time-dependent distances), the recommended strategy is to split the simulation into chunks and accumulate the loss across chunks. This provides temporal resolution while avoiding the memory overhead of storing full trajectories for gradient computation.

```julia
function loss_with_trajectory(..., n_chunks)
    loss_total = 0.0
    coords_current = copy(coords)
    velocities_current = copy(velocities)
    
    for chunk in 1:n_chunks
        # Run short simulation chunk
        simulate!(sys, simulator, chunk_steps)
        
        # Compute chunk contribution to loss
        chunk_loss = compute_loss(sys.coords)
        loss_total += chunk_loss
        
        # Prepare for next chunk
        coords_current = sys.coords
        velocities_current = sys.velocities
        sys = System(...)  # Reset with current state
    end
    
    return loss_total / n_chunks
end
```

> [!TIP]
> The magnitude of gradients may be less important than their sign in stochastic systems. Consider sampling gradients across different sources of stochasticity such as initial velocities, conformations, and thermostat random seeds. Use gradient averaging or sign-based optimization strategies to improve convergence robustness.

## Practical Considerations

Several practical considerations affect the success of gradient-based force field optimization. Exploding gradients can be problematic with velocity Verlet integration in the NVE ensemble; thermostatted simulations (Langevin dynamics, Berendsen thermostat, or velocity rescaling) typically provide more stable gradient behavior [differentiable.md](docs/src/differentiable.md#L453-L458).

Performance considerations include:
- Reverse mode AD: Memory intensive, suitable for many parameters
- Forward mode AD: Memory efficient (constant in step count), suitable for few parameters
- GPU execution: Often 5-10× slower than non-differentiable simulations
- CPU execution: Typically within one order of magnitude of primal performance

The table below summarizes the trade-offs between AD modes and computational backends:

| Configuration | Memory | Speed | Parameter Count Suitability | Use Cases |
|--------------|--------|-------|----------------------------|-----------|
| CPU Reverse | High | Good | 10-1000+ | Standard force field refinement |
| CPU Forward | Low | Moderate | 1-10 | Few-parameter optimization |
| GPU Reverse | Very High | Slow | 10-1000+ | Large system, many parameters (requires checkpointing) |
| GPU Forward | Moderate | Moderate | 1-10 | Few parameters on GPU |

Gradient checkpointing with Checkpointing.jl can mitigate memory constraints for large simulations by storing only selected intermediate states and recomputing others during the backward pass, trading compute for memory.

The golden rule of differentiable programming applies: always validate gradients against finite differencing when possible, especially for custom interactions or loss functions. Molly.jl's test suite provides examples of gradient verification patterns that can be adapted to custom optimization workflows.

## Next Steps

For comprehensive understanding of differentiable simulation, explore [Automatic Differentiation with Enzyme](24-automatic-differentiation-with-enzyme) to learn about Enzyme integration and advanced AD techniques. To understand the stochastic aspects of differentiable simulations, see [Differentiable Stochastic Simulations](26-differentiable-stochastic-simulations). For practical optimization strategies and performance considerations, refer to [Performance Considerations for AD](27-performance-considerations-for-ad).

To implement custom loss functions or understand interaction parameter injection, examine [Custom Interaction Implementation](12-custom-interaction-implementation). For integration with machine learning workflows, explore [Neural Network Potentials](25-gradient-based-force-field-optimization#neural-network-potentials) and related packages like Flux.jl.