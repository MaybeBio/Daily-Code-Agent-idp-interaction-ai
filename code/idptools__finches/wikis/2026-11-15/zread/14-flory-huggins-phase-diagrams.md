

Finches predicts homotypic liquid-liquid phase separation (LLPS) phase diagrams by bridging **mean-field interaction parameters** (ε) with the **analytical solution to Flory-Huggins theory**. The pipeline converts a protein sequence into a temperature–concentration (T vs. φ) phase diagram with both binodal and spinodal curves, enabling rapid, sequence-resolved predictions of phase behavior without simulation.

## Theoretical Foundation

The phase diagram engine is built on the closed-form solution to the two-component Flory-Huggins model derived by Qian, Michaels & Knowles (2022), which provides **analytical expressions for the binodal and spinodal** curves that previously required numerical root-finding. In Flory-Huggins theory, the central interaction parameter is χ (chi), defined as:

> **χ = ε / (k_B × T)**

where ε is the site-to-site contact energy, k_B is Boltzmann's constant, and T is temperature. Larger positive χ values correspond to stronger attractive interactions and thus a greater driving force for phase separation. The **critical point** — the thermodynamic boundary below which no phase separation occurs regardless of concentration — is given analytically as:

| Quantity | Expression |
|---|---|
| Critical volume fraction (φ_c) | 1 / (1 + √n) |
| Critical chi (χ_c) | 0.5 × (1 + 1/√n)² |

where n is the polymer chain length. For n = 1 (symmetric solute-solvent), this recovers the classic result χ_c = 2, φ_c = 0.5.

Sources: [backend.py](finches/analytical_fh/backend.py#L31-L48), [floryhuggins.py](finches/analytical_fh/floryhuggins.py#L1-L30)

## Architecture: From Sequence to Phase Diagram

The conversion from a protein sequence to a phase diagram proceeds through a well-defined computational pipeline. Finches first computes the mean-field interaction parameter ε from the sequence using a forcefield model, then converts ε to a per-residue contact energy (δ_ε), builds the χ vs. φ diagram using the analytical FH solution, and finally recasts χ into temperature via T = δ_ε / χ.

```mermaid
flowchart LR
    A["Protein Sequence"] --> B["ε Calculation<br/>(forcefield model)"]
    B --> C["δ_ε = -ε / L<br/>(per-residue contact energy)"]
    C --> D["χ vs. φ Diagram<br/>(analytical FH binodal/spinodal)"]
    D --> E["T vs. φ Diagram<br/>(T = δ_ε / χ)"]
    
    style A fill:#e8f4fd,stroke:#2980b9
    style B fill:#fdf2e8,stroke:#e67e22
    style C fill:#fdf2e8,stroke:#e67e22
    style D fill:#e8f8e8,stroke:#27ae60
    style E fill:#e8f8e8,stroke:#27ae60
```

The **sign convention** is critical: finches defines negative ε as attractive (promoting phase separation). When converting to Flory-Huggins χ, the sign is reversed because χ is positive for attractive interactions. The per-residue division (`-ε / len(seq)`) is necessary because ε is a chain-level quantity, but the FH solution already accounts for polymer length through n, so only the average per-residue contribution should feed into χ.

Sources: [epsilon_to_FHtheory.py](finches/epsilon_to_FHtheory.py#L399-L448), [floryhuggins.py](finches/analytical_fh/floryhuggins.py#L72-L130)

## Three Binodal Computation Modes

The analytical FH backend offers three modes for computing the binodal, all sharing the same input signature `(chi, n)`:

| Mode | Function | Description | When to Use |
|---|---|---|---|
| `analytic_binodal` | `backend.analytic_binodal()` | Closed-form analytic solution from the Qian et al. paper | **Default and recommended** for all cases |
| `binodal` | `backend.binodal()` | Self-consistent iterative solution (5 iterations by default) with optional improved Newton map | Validation or comparison against the analytic result |
| `GL_binodal` | `backend.GL_binodal()` | Ginzburg-Landau expansion near the critical point | Only valid near χ_c; provided for completeness |

The `analytic_binodal` mode uses distinct functional forms for n = 1 (symmetric mixture, using `tanh`-based expressions) and n > 1 (asymmetric, using `coth`-based expressions with the dimensionless distance parameter D = (χ - χ_c) / χ_c). For n > 1, the solution involves auxiliary variables `a = n^(1/4)`, hyperbolic cotangent terms, and exponential mapping to recover the dense and dilute volume fractions.

> [!TIP]
> If ε is positive (net repulsive), finches overrides it to -0.01 — a negligibly attractive value — because a repulsive ε would invert the χ-to-T conversion, yielding unphysical "inverse" phase diagrams. This is a pragmatic workaround; sequences with repulsive ε simply do not phase-separate.

Sources: [backend.py](finches/analytical_fh/backend.py#L200-L316), [floryhuggins.py](finches/analytical_fh/floryhuggins.py#L55-L97)

## API Reference

### Low-Level: `analytical_fh` Module

These functions operate directly in χ–φ space and require you to specify χ values explicitly.

**`floryhuggins.calculate_binodal(L, mode, chi_min, chi_max, n_points)`** — Sweeps χ from `chi_min` to `chi_max` with `n_points` evenly-spaced values, computing the binodal at each. Returns a 5-element tuple: `(chi_values, dilute_φ, dense_φ, critical_φ, critical_χ)`.

**`floryhuggins.calculate_spinodal(L, chi_min, chi_max, n_points)`** — Same sweep structure for the spinodal curve. Returns the same 5-element tuple format.

**`backend.critical(n)`** — Returns `[φ_c, χ_c]` for polymer length n.

**`backend.spinodal(chi, n)`** / **`backend.binodal(chi, n)`** / **`backend.analytic_binodal(chi, n)`** / **`backend.GL_binodal(chi, n)`** — Compute a single point on the respective curve. Each returns `[dense_φ, dilute_φ]` (or raises `ValueError` if χ < χ_c).

Sources: [floryhuggins.py](finches/analytical_fh/floryhuggins.py#L8-L197), [backend.py](finches/analytical_fh/backend.py#L31-L316)

### Mid-Level: `epsilon_to_FHtheory` Module

These functions bridge sequence-level ε calculations with the FH backend, producing T–φ diagrams.

**`return_phase_diagram(seq, X_class)`** — Computes ε from the sequence using the provided `InteractionMatrixConstructor`, then delegates to `epsilon_to_phase_diagram`. This is the most common entry point for single-sequence diagrams.

**`epsilon_to_phase_diagram(seq, epsilon)`** — Given a pre-computed ε value and sequence (used only for length), builds the full T–φ phase diagram. Returns an 8-element list:

| Index | Content |
|---|---|
| `[0]` | Dilute phase φ (binodal) |
| `[1]` | Dense phase φ (binodal) |
| `[2]` | `[critical_φ, critical_T]` (binodal) |
| `[3]` | Temperature array (binodal) |
| `[4]` | Dilute phase φ (spinodal) |
| `[5]` | Dense phase φ (spinodal) |
| `[6]` | `[critical_φ, critical_T]` (spinodal) |
| `[7]` | Temperature array (spinodal) |

Sources: [epsilon_to_FHtheory.py](finches/epsilon_to_FHtheory.py#L241-L399)

### High-Level: Frontend Methods

The frontend objects (`Mpipi_frontend`, `CALVADOS_frontend`) provide the most convenient API.

**`frontend.build_phase_diagram(seq)`** — Computes ε internally via `self.epsilon(seq, seq)` and calls `epsilon_to_phase_diagram`. Returns the same 8-element list as above.

**`frontend.plot_phase_diagram(seq, ...)`** — Builds and renders the phase diagram in a single call. Accepts `line_color`, `line_width`, `xlim`, `ylim`, `xlog`, `width`, `height`, and `filename` for customization. Returns `[phase_data, fig, ax]`.

**`frontend.plot_multiple_phase_diagrams(seqdict, tc_ref, ...)`** — Overlays multiple phase diagrams on one plot. `seqdict` maps sequence names to `[sequence_string, color]`, and `tc_ref` specifies which sequence's critical temperature is used to normalize the y-axis (T/T_c).

```python
from finches.frontend.mpipi_frontend import Mpipi_frontend

mf = Mpipi_frontend()
seq = 'MASNDYTQQATQSYGAYPTQPGQGYSQQSSQPYGQQSYSGYSQSTDTSGYGQSSYSSYGQSQNTGYGTQSTPQGYGSTGGYGSSQSSQSSYGQQSSYPGYGQQPAPSSTSGSYGSSSQSSSYGQPQSGSYSQQPSYGGQQQSYGQQQSYNPPQGYGQQNQYNS'

# Build and plot in one call
result = mf.plot_phase_diagram(seq, ylim=[0.7, 1.05])

# Or build data separately for custom plotting
B = mf.build_phase_diagram(seq)
plt.plot(B[0], B[3], 'blue', label='Binodal (dilute arm)')
plt.plot(B[1], B[3], 'blue', label='Binodal (dense arm)')
plt.plot(B[4], B[7], 'blue', ls='--', label='Spinodal (dilute)')
plt.plot(B[5], B[7], 'blue', ls='--', label='Spinodal (dense)')
plt.plot(B[2][0], B[2][1], 'ko', ms=4, label='Critical point')
```

Sources: [frontend_base.py](finches/frontend/frontend_base.py#L1000-L1199), [epsilon_to_FHtheory.py](finches/epsilon_to_FHtheory.py#L241-L330)

## Condition-Dependent Phase Diagrams

Finches provides three builder functions that iterate over environmental conditions, re-initializing the forcefield at each step to capture how phase boundaries shift with solution chemistry. Each returns `[condition_list, diagrams_dict, epsilons_dict]`.

| Function | Condition Parameter | Effect on Phase Behavior |
|---|---|---|
| `build_SALT_dependent_phase_diagrams(seq, X_class, condition_list)` | `X_class.parameters.salt` | Higher salt screens electrostatic repulsion, generally expanding the two-phase region |
| `build_PH_dependent_phase_diagrams(seq, X_class, condition_list, ...)` | `X_class.parameters.pH` | Alters residue charge states (especially His, Asp, Glu), reshaping electrostatic contributions |
| `build_DIELECTRIC_dependent_phase_diagrams(seq, X_class, condition_list, ...)` | `X_class.parameters.dielectric` | Changes the solvent dielectric constant, modulating all electrostatic terms |

Each function validates that the requested condition parameter exists on the forcefield model, iterates through the condition list updating the lookup dictionary at each step, computes both the phase diagram and the ε value, and **restores the original parameters** after completion. For salt-dependent diagrams, `build_SALT_dependent_phase_diagrams` also reads and preserves the `charge_prefactor` and `null_interaction_baseline` from the `X_class`.

```python
from finches import epsilon_to_FHtheory
from finches.frontend.mpipi_frontend import Mpipi_frontend

mf = Mpipi_frontend()
salts = [0.01, 0.05, 0.1, 0.15, 0.2, 0.3]

result = epsilon_to_FHtheory.build_SALT_dependent_phase_diagrams(
    seq, mf.IMC_object, salts
)

# result[0] = salt values
# result[1] = dict mapping salt -> phase diagram data
# result[2] = dict mapping salt -> epsilon value
```

For a deeper treatment of how phase diagrams respond to environmental conditions, see [Condition-Dependent Phase Behavior](15-condition-dependent-phase-behavior).

Sources: [epsilon_to_FHtheory.py](finches/epsilon_to_FHtheory.py#L44-L237)

## Important Caveats

**Repulsive ε handling**: When ε > 0 (net repulsive), the code sets ε = -0.01 before computing the diagram. This avoids a divide-by-zero error in the χ-to-T conversion and prevents unphysical inverted phase diagrams. The resulting diagram will show an extremely narrow two-phase region that effectively signals "no phase separation." If you need to distinguish truly marginal phase-separating sequences from non-separating ones, inspect the raw ε value before calling the diagram builder.

**χ range defaults**: The `epsilon_to_phase_diagram` function uses `chi_max=0.8` with `n_points=50000` for the binodal and `n_points=10000` for the spinodal. For very long polymers (n > 1000), the critical χ approaches 0.5 from above, so the default `chi_min=0.5` is appropriate. However, for short polymers or extremely strong interactions, you may need to increase `chi_max` by calling `floryhuggins.calculate_binodal` directly.

**Mean-field approximation**: These diagrams are derived from a mean-field lattice model (Flory-Huggins) and do not capture concentration fluctuations, finite-size effects, or multi-component interactions. They are best interpreted as qualitative predictions of the *tendency* and *approximate boundaries* of phase separation.

Sources: [epsilon_to_FHtheory.py](finches/epsilon_to_FHtheory.py#L399-L448)

## Suggested Reading Path

- For the ε calculation that feeds into phase diagrams: [Epsilon Calculation and Weighting](9-epsilon-calculation-and-weighting)
- For how environmental conditions reshape phase boundaries: [Condition-Dependent Phase Behavior](15-condition-dependent-phase-behavior)
- For the frontend API that wraps these functions: [Mpipi and CALVADOS Frontends](5-mpipi-and-calvados-frontends)
- For forcefield parameter details: [Mpipi Forcefield Parameters](11-mpipi-forcefield-parameters)