

Umbrella sampling is a powerful enhanced sampling technique that enables calculation of free energy profiles along collective variables (CVs) by applying biasing potentials (umbrella windows) and using statistical reweighting methods like MBAR to recover the unbiased free energy landscape. This page covers the implementation of custom restraints, setup of umbrella sampling simulations, and free energy analysis using the Multistate Bennett Acceptance Ratio (MBAR) method in Molly.

## Custom Restraint Interactions

Molly provides flexible mechanisms for implementing custom restraint interactions that bias simulations along desired collective variables. While the library includes a built-in `HarmonicPositionRestraint` for positional restraints, advanced users can implement specialized restraints for dihedral angles, distances, or other CVs.

### Harmonic Position Restraints

The `HarmonicPositionRestraint` constrains an atom to a reference position using a harmonic potential, defined mathematically as:

```math
V(\boldsymbol{x}) = \frac{1}{2} k |\boldsymbol{x} - \boldsymbol{x}_0|^2
```

The force is calculated as `F = -k(r - r₀)` in the direction from the atom to the restraint center, with special handling for zero-distance cases to avoid numerical instability. This implementation provides a straightforward way to restrain atoms near specific coordinates during simulations.

Sources: [harmonic_position_restraint.jl](src/interactions/harmonic_position_restraint.jl#L1-L33)

### Custom Dihedral Restraints

For torsional CVs, custom dihedral restraints can be implemented as specific interactions. The documentation demonstrates a `DihedralRestraint` that uses a flat-bottom harmonic potential, enabling umbrella sampling along dihedral angles without excessive force at potential well boundaries:

```julia
struct DihedralRestraint{A, K}
    ϕ0::A  # Reference dihedral angle (radians)
    Δϕ::A  # Width of flat bottom region (radians)
    k::K   # Force constant
end
```

The potential energy function applies the harmonic restraint only when the dihedral angle deviates beyond the flat-bottom region, with the force calculation requiring careful geometric derivatives using vector cross products and gradient decomposition. This approach prevents large forces that could destabilize simulations when the CV is within acceptable bounds.

Sources: [free_energy.md](docs/src/free_energy.md#L21-L193)

## Thermodynamic State Management

The MBAR implementation uses `ThermoState` objects to encapsulate the thermodynamic parameters (inverse temperature β and pressure p) along with the `System` object used for energy evaluation. Each umbrella window corresponds to a different thermodynamic state, potentially with different restraint parameters or simulation conditions.

The `ThermoState` constructor performs unit consistency checks to ensure β and pressure are provided in compatible dimensions with the system's energy units. This validation prevents subtle unit conversion errors that could corrupt free energy calculations. The inverse temperature β = 1/k_BT is the critical parameter that converts potential energies to reduced potentials used in MBAR equations.

Sources: [mbar.jl](src/free_energy/mbar.jl#L9-L64)

## Umbrella Sampling Workflow

### Pulling Simulations

Before running umbrella sampling production simulations, pulling simulations gradually move the system along the CV to generate equilibrated starting structures for each umbrella window. This sequential approach ensures each window begins from a configuration already near its target CV value, avoiding large initial forces and improving equilibration efficiency.

The workflow involves:
1. Loading an equilibrated starting structure
2. Calculating the initial CV value
3. Iteratively running biased simulations with restraint centers shifted by increments
4. Using final coordinates from each window as starting coordinates for the next

For dihedral umbrella sampling spanning 360° with 60 windows, each window advances the restraint center by 6°, using a flat-bottom potential with 3° width to provide overlap between adjacent windows.

Sources: [free_energy.md](docs/src/free_energy.md#L193-L328)

### Umbrella Window Production

Once pulling simulations complete, each umbrella window runs independently as a production simulation. These simulations use identical force fields, thermostats, barostats, and integrators as the pulling simulations, but run for much longer durations to adequately sample the CV distribution within each window.

The independent nature of umbrella windows enables parallel execution, significantly reducing wall-clock time for the total sampling effort. Each window writes trajectory files at regular intervals (e.g., every picosecond) for subsequent analysis.

Sources: [free_energy.md](docs/src/free_energy.md#L328-L420)

## Statistical Subsampling

Before MBAR analysis, trajectory data must be subsampled to remove correlated frames and provide statistically independent samples. The `statistical_inefficiency` function estimates the integrated autocorrelation time using the Initial Positive Sequence (IPS) method with finite-sample tapering.

The statistical inefficiency g determines the subsampling stride as `ceil(Int, g)`, with effective sample size calculated as `floor(N / stride)`. This preprocessing step ensures MBAR's assumption of independent samples is approximately satisfied and prevents overestimation of statistical precision.

Sources: [stats.jl](src/free_energy/stats.jl#L1-L109)

## MBAR Data Assembly

The `assemble_mbar_inputs` function prepares the reduced potential matrix u of size N×K, where N is the total number of frames across all windows and K is the number of thermodynamic states (umbrella windows). For each frame n and state k, the reduced potential is:

```math
u_{n,k} = \beta_k U_k(\mathbf{x}_n) + \beta_k p_k V_n
```

where U_k is the potential energy under state k's Hamiltonian (including restraints), V_n is the system volume, and the second term is included only for constant-pressure simulations.

The function also optionally computes reduced potentials for a target (unbiased) state and can apply per-frame energy shifts to improve numerical stability by preventing exponent overflow in the MBAR equations.

Sources: [mbar.jl](src/free_energy/mbar.jl#L66-L199)

## Solving MBAR Equations

The core MBAR algorithm solves the self-consistent equations:

```math
f_k = -\ln \sum_{n=1}^N \frac{N_k \exp(f_k - u_{n,k})}{\sum_{j=1}^K N_j \exp(f_j - u_{n,j})}}
```

where f_k are the dimensionless free energies for each state. The `iterate_mbar` function uses an initial guess based on mean reduced potentials and iterates until the relative change falls below a specified tolerance or the maximum iteration count is reached.

The iteration uses numerically stable log-sum-exp operations to avoid overflow, with gauge-fixing setting f[1] = 0 to handle the indeterminate additive constant in free energies.

Sources: [mbar.jl](src/free_energy/mbar.jl#L201-L297)

## Potential of Mean Force Calculation

The `pmf_with_uncertainty` function computes a 1D free energy profile along a collective variable using MBAR-reweighted probabilities. The workflow involves:

1. Computing MBAR weights W[n,k] and target-state weights w_target[n]
2. Binning the collective variable values into histogram bins
3. Estimating probability density p(r) in each bin using weighted samples
4. Converting to PMF as F(r) = -k_BT ln p(r)

The implementation also computes asymptotic uncertainties using MBAR's covariance propagation framework, providing error estimates for each point on the PMF.

Sources: [mbar.jl](src/free_energy/mbar.jl#L503-L841)

## High-Level MBAR Interface

Molly provides a convenient high-level function that orchestrates the complete MBAR analysis workflow. This function accepts coordinate trajectories, boundaries, thermodynamic states, target state definition, and collective variable values, then automatically:

1. Assembles MBAR input matrices
2. Solves for state free energies
3. Computes reweighting weights
4. Bins the collective variable
5. Returns the PMF with uncertainties

The returned `PMF` struct contains bin centers, widths, edges, PMF values in both k_BT and energy units, uncertainties, and underlying probability densities.

Sources: [mbar.jl](src/free_energy/mbar.jl#L800-L841)

> [!TIP]
> When implementing custom restraints, always include numerical tolerances that scale with system dimensions to handle edge cases like collinear geometries in dihedral calculations. The dihedral restraint example demonstrates this pattern with `tol_b` and `tol_xy` scaled by characteristic length Ls.

> [!TIP]
> Use the `shift=true` option in `assemble_mbar_inputs` when energy differences between states are large, as this subtracts per-frame row minima from the reduced potential matrix to prevent numerical overflow in the MBAR exponent calculations.

## Next Steps

For comprehensive free energy analysis, continue to [Multistate Bennett Acceptance Ratio (MBAR)](28-multistate-bennett-acceptance-ratio-mbar) for deeper theoretical understanding, or proceed to [Potential of Mean Force Calculations](30-potential-of-mean-force-calculations) for practical applications. If you need to optimize force field parameters using free energy data, explore [Gradient-Based Force Field Optimization](25-gradient-based-force-field-optimization).

For understanding the foundational integration methods used in umbrella sampling simulations, review [Integration Methods: Verlet, Velocity Verlet, and Langevin](15-integration-methods-verlet-velocity-verlet-and-langevin).