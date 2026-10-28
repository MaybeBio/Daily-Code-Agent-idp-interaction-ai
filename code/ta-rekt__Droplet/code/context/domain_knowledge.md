# Domain Knowledge: Biomolecular Condensate Simulations

**Last Updated**: January 23, 2026  
**Purpose**: Physical and biological context for experimental design

---

## 1. Scientific Background

### 1.1 Research Question

**Core Hypothesis**: Biomolecular condensates (phase-separated droplets) can enhance cellular sensing capabilities beyond what diffusion alone permits.

**Specific Context**: 
- Transcription factors (TF) bind to enhancer DNA sequences
- This recruits RNA polymerase (RNAP) to promoter sequences
- TF concentration is a signal to be "sensed"
- Condensation may amplify weak signals or sharpen responses

### 1.2 Physical System

**Lattice Gas Model**:
- 3D cubic lattice (side length L = 16-32 sites)
- Particles occupy single sites (hard-core exclusion)
- Periodic boundary conditions (no walls)
- Metropolis Monte Carlo dynamics (detailed balance)

**Energy Function**:
```
H = Σ_ij ε_ij n_i n_j
```
- `ε_ij`: Pairwise interaction energy between species i and j
- `n_i`: Occupancy (0 or 1) at site i
- Negative `ε`: Attractive interactions
- Positive `ε`: Repulsive interactions (rare)

**Temperature**:
- Set to 1 (in units of kB·T)
- Can vary via interaction strengths: `ε_effective = ε_0 / T`

---

## 2. Key Physical Concepts

### 2.1 Phase Separation / Condensation

**What it is**: Spontaneous demixing into dense and dilute phases due to attractive self-interactions.

**When it occurs**:
- Self-interaction `ε < 0` (attractive)
- Below critical temperature: `|ε| > ε_critical ≈ 0.5` (rule of thumb)
- Sufficient particle density (typically > 10% filling fraction)

**Observable signatures**:
- Bimodal density distribution (droplet vs bulk)
- Large spatial density fluctuations
- Cluster formation (detected via connected components)

**Example**:
```julia
tf = FastBulk(length=2048)  # ~50% filling in 16³ lattice
intPot!(tf, -1.0)           # Attractive self-interaction
# Expect: TF forms droplet(s) after equilibration
```

### 2.2 Prewetting

**What it is**: Preferential condensation near an attractive surface (1st-order transition precursor).

**Physical setup**:
- Bulk species (TF) with self-attraction
- Linear surface (DNA, enhancer) with stronger TF binding
- As TF concentration increases, forms layer on surface before bulk droplets appear

**Key parameters**:
- `ε_bulk = ε_TF-TF`: Controls bulk phase separation
- `ε_surface = ε_TF-enhancer`: Controls surface affinity
- Typically: `|ε_surface| > |ε_bulk|` (e.g., -2.0 vs -1.0)

**Experimental observable**:
- TF density profile perpendicular to surface
- Enhanced binding near surface relative to bulk

**Example**:
```julia
enhancer = FastBulk(length=6, initShape="line")
tf = FastBulk(length=4096)
intPot!(tf, -1.0)            # Bulk condensation
intPot!(tf, enhancer, -2.0)  # Surface binding
# Expect: TF layer around enhancer, thickness depends on concentration
```

### 2.3 Cooperative Binding

**What it is**: Binding affinity modulated by presence of other bound species.

**Mechanisms**:
- Direct TF-RNAP interaction: `ε_TF-RNAP < 0` stabilizes co-localization
- Condensate-mediated: TF droplet creates high local concentration → enhances RNAP recruitment

**Observable signature**:
- Non-linear dose-response curves (Hill coefficient > 1)
- Binding thresholds (switch-like behavior)

**Example**:
```julia
intPot!(tf, enhancer, -2.0)      # TF binds enhancer
intPot!(rnap, promoter, -2.0)    # RNAP binds promoter
intPot!(tf, rnap, -2.0)          # Cooperative interaction
# Expect: RNAP binding increases super-linearly with TF concentration
```

### 2.4 On-Site Coupling Interactions (Two-Way and Three-Way)

**What it is**: On-site coupling interactions where energy depends on multiple species co-located at the same lattice site. Two-way coupling rewards/penalizes two species sharing a site; three-way coupling rewards/penalizes three species sharing a site. These are distinct from nearest-neighbor `intMat` interactions.

**Physical formula**: For species `self` at position `r`, the coupling contribution is:

```
# Two-way (j == k): on-site co-location of self and j
H_coupling_2way(self, r) = Σ_j coupling[self, j, j] × cplpq(j, j, r)

# Three-way (j < k): all three co-located at the same site
H_coupling_3way(self, r) = Σ_{j<k} coupling[self, j, k] × cplpq(j, k, r)
```

where `cplpq(j, k, pos)` returns 1 if both species `j` and `k` occupy site `pos`, 0 otherwise. The coupling matrix has two symmetries:
- **j,k symmetry**: `coupling[i,j,k] = coupling[i,k,j]`
- **First-two-index symmetry**: All BulkMonomer participants share the same coupling value. E.g., `couplingPot!(a, b, c, -0.5)` sets `coupling[a,b,c] = coupling[b,a,c] = coupling[c,a,b] = -0.5`.

**Key design principle**: Coupling energy is "owned by self" — only the species that owns the coupling (has a `couplingMat` field) has its Hamiltonian updated. All BulkMonomer participants own the coupling due to first-two-index symmetry. Partners' `deltaHP` does NOT include coupling contributions.

**API**:
```julia
# Three-way coupling: all three species co-located at the same site
couplingPot!(a, b, c, energy)

# Two-way coupling: two species co-located at the same site (on-site, distinct from intMat)
couplingPot!(a, b, energy)

# a, b, c: any AbstractMonomer (but at least one must be BulkMonomer to own it)
# energy: Float64 — coupling strength (negative = attractive/cooperative)
```

**Restrictions**:
- `FastBulkMonomer` does NOT support couplings (throws `ErrorException`)
- Only `BulkMonomer` can own couplings (has `couplingMat`, `couplingIdx` fields)
- Partners can be any species (BulkMonomer, PolyMonomer, etc.)
- Two-way `couplingMat[i,j,j]` is distinct from nearest-neighbor `intMat[i,j]`

**Example: Cooperative binding with on-site coupling**:
```julia
tf1 = BulkMonomer(name="TF1", length=100)
tf2 = BulkMonomer(name="TF2", length=50)
dna = poly.monomers[1]  # PolyMonomer from a Polymer

# Nearest-neighbor interactions (standard intMat)
intPot!(tf1, -1.0)               # TF1 self-interaction
intPot!(tf2, -1.0)               # TF2 self-interaction
intPot!(tf1, dna, -0.5)          # TF1-DNA nearest-neighbor binding
intPot!(tf2, dna, -0.5)          # TF2-DNA nearest-neighbor binding

# Two-way on-site coupling: TF1 and TF2 at the same site
couplingPot!(tf1, tf2, -0.2)     # On-site co-localization bonus

# Three-way on-site coupling: TF1, TF2, and DNA all at the same site
couplingPot!(tf1, tf2, dna, -0.3)  # Extra stabilization when all three co-locate
```

**Observable signature**:
- Enhanced co-localization at sites where multiple coupled species overlap
- Energy distributions shifted more negative at multi-species co-location sites
- Two-way coupling shifts the on-site penalty/reward independently of nearest-neighbor `intMat`

---

## 3. Key Metrics & Interpretations

### 3.1 Hamiltonian (Energy)

**Definition**: `H = Σ_ij ε_ij n_i n_j + Σ_{i,j<k} coupling_i[j,k] × nnpq_terms` (sum over all occupied neighbor pairs, plus three-way coupling terms)

**What it tells you**:
- More negative = more favorable interactions
- Time series shows equilibration (plateau) vs drift
- Large fluctuations indicate phase transition proximity

**Typical values**:
- Small system (512 particles): H ~ -500 to -5000
- Sign of convergence: H variance decreases after burn-in

### 3.2 Acceptance Rate

**Definition**: Fraction of Monte Carlo moves accepted (Metropolis criterion)

**What it tells you**:
- Too high (> 60%): Weak interactions, fast dynamics, may under-sample
- Too low (< 5%): System jammed, slow equilibration
- Optimal: 20-40% for efficient sampling

**Troubleshooting**:
- Low acceptance + drift: Increase runtime or decrease particle density
- High acceptance: System may be in gas phase (check interaction strengths)

### 3.3 Binding Events

**Definition**: `nspeciesBound(species1, species2_list)` counts species1 particles with ≥1 neighbor from species2

**What it tells you**:
- Time series: Binding dynamics (on/off kinetics)
- Mean value: Equilibrium binding occupancy
- Variance: Binding fluctuations (relates to information capacity)

**Example interpretation**:
```julia
binding = nspeciesBound(tf, [enhancer])
# binding ∈ [0, min(tf.length, enhancer.length × 26)]
# (26 = max neighbors per site in 3D)
# High mean + low variance = stable binding
# Low mean + high variance = transient binding
```

### 3.4 Spatial Density Fields

**Definition**: 3D histogram of particle positions (accumulated over trajectory)

**What it tells you**:
- Droplet shapes and sizes
- Spatial organization (e.g., TF clustered near enhancers)
- Symmetry breaking

**Analysis**:
- Threshold to identify high-density regions (droplets)
- Radial density profiles: `ρ(r)` around reference point
- Structure factors (Fourier transform) for periodicity

**Typical output**:
```julia
tf_density = zeros(Float64, (L, L, L))
# After normalization: tf_density[i,j,k] ∈ [0, 1]
# Represents probability of finding TF at site (i,j,k)
```

### 3.5 Fisher Information Rate

**Definition** (conceptual): Sensitivity of binding statistics to parameter changes

**Calculation** (post-processing):
```python
# From time series B(t) of binding events at different concentrations c
# Fisher info ∝ (∂⟨B⟩/∂c)² / Var(B)
```

**Interpretation**:
- Higher FI = better sensing of concentration
- Condensates may enhance FI via variance reduction (ultra-sensitivity)

---

## 4. Dimensional Analysis & Parameter Choices

### 4.1 Lattice Size

**Trade-offs**:
- Larger L → less finite-size effects, higher computational cost
- Smaller L → faster runs, but boundary effects significant

**Guidelines**:
- L = 16: Quick tests, single droplets
- L = 20-24: Standard production runs
- L = 32: High resolution, multiple droplets
- L > 40: Rare (very expensive, use only if necessary)

**Rule of thumb**: Droplet diameter should be < L/2 to avoid self-interaction through periodic boundaries

### 4.2 Particle Counts

**Filling Fraction**:
```
ϕ = N_total / L³
```
- ϕ < 0.1: Dilute (gas phase likely)
- ϕ = 0.1-0.4: Moderate (prewetting, partial condensation)
- ϕ = 0.4-0.5: Dense (full condensation, slow dynamics)
- ϕ > 0.5: Overcrowded (jamming, avoid unless specific goal)

**Typical setups**:
- Bulk studies: 512-2048 particles in 16³ lattice (ϕ ~ 0.1-0.5)
- Prewetting: 1000-4000 TF in 32³ lattice (ϕ ~ 0.03-0.12)
- Polymers: Fewer particles due to chain connectivity

### 4.3 Interaction Strengths

**Energy scale reference**: All in units of k_B T (Boltzmann constant × Temperature)

**Typical ranges**:
- Nonspecific interactions: `ε = -0.1 to -0.5`
  - Example: TF-DNA backbone (weak, non-selective)
- Specific binding: `ε = -1.0 to -2.0`
  - Example: TF-enhancer (moderate, selective)
- Strong binding: `ε = -2.0 to -4.0`
  - Example: Highly specific protein-protein (rare in bulk systems)

**Physical meaning**:
- `ε = -1.0`: Binding probability ~ e^(-1) ≈ 37% (comparable to thermal fluctuations)
- `ε = -2.0`: Binding probability ~ e^(-2) ≈ 13.5% (strongly favored)
- `ε = -4.0`: Binding probability ~ e^(-4) ≈ 1.8% (essentially locked)

**Warning**: `|ε| > 5.0` often causes numerical issues (overflow, extreme jamming)

**Three-way coupling strengths**:
- Cooperative enhancement: `coupling = -0.1 to -0.5` (modest context-dependent stabilization)
- Strong cooperative coupling: `coupling = -0.5 to -1.0`
- Typically `|coupling| < |ε_two-way|` (coupling is a correction to two-way interactions)
- Sign convention: negative = cooperative (enhanced by context), positive = competitive (suppressed)

### 4.4 Time Scales

**Monte Carlo Time**:
- 1 MC sweep = N attempted moves (N = total particles)
- Real time mapping: Not directly physical (model-dependent)

**Simulation runtime**:
- `t_max = 5 min`: Quick tests, may not equilibrate
- `t_max = 10 min`: Standard runs, most systems equilibrate
- `t_max = 20+ min`: Slow dynamics (polymers, dense systems)

**Burn-in fraction**:
- `t_burn_in = 0.5 × t_max`: Conservative (definitely equilibrated)
- `t_burn_in = 0.8 × t_max`: Aggressive (maximizes statistics)
- Monitor Hamiltonian: Equilibrated when H variance becomes stationary

**Sampling interval**:
- `step = 100`: Fine temporal resolution (large data)
- `step = 1000`: Standard (balance statistics vs file size)
- `step = 5000`: Coarse (only equilibrium averages)

---

## 5. Experimental Design Principles

### 5.0 Start Small, Scale Up

**CRITICAL PRINCIPLE**: Always begin with minimal resources for rapid iteration

**Initial proposal requirements**:
- **Runtime**: 1 minute (not 5, not 10)
- **CPU count**: 1 (not 16)
- **Reason**: Fast feedback loop for validating proposal format, spec file syntax, and basic physics

**Workflow**:
1. **Draft proposal** → 1 min, 1 CPU test run
2. **Iterate** based on errors, visualization, preliminary results
3. **Scale up** ONLY after validation: increase runtime (5-10 min), then CPUs (parallel)
4. **Production** runs: Long time, many replicates

**Why this matters**:
- Catch errors early (syntax, initialization, metric bugs)
- Test visualization pipeline before investing compute time
- Verify physical behavior (equilibration, acceptance rate) quickly
- Iterate on proposal format without wasting hours

**Template override**: When copying from existing spec files, ALWAYS replace their runtime/CPU settings with small initial values, regardless of what the template used.

### 5.1 Hypothesis-Driven Sweeps

**Good sweep design**:
1. **Single variable**: Change one parameter at a time
2. **Physical range**: Cover regime from weak to strong effect
3. **Resolution**: Enough points to see smooth trends (typically 10-20 points)

**Example**:
```julia
# Hypothesis: TF condensation enhances RNAP binding
# Sweep: TF concentration from dilute to dense
vals = range(1300, stop=4000, step=180)  # 16 points
# Measure: RNAP-promoter binding vs TF concentration
```

**Bad example**:
```julia
# Sweep both concentration AND interaction strength
alphas = [1000, 2000, 3000]
betas = [-0.5, -1.0, -2.0]
# Problem: Can't isolate which effect causes observed changes
```

### 5.2 Control Experiments

**Always include**:
- **Null condition**: No interactions (ε = 0) → baseline from diffusion alone
- **Decoupled condition**: Remove cooperative interactions → isolate condensation effect
- **Saturation**: High enough parameter that effect plateaus

**Example**:
```julia
# To prove condensation matters, compare:
# 1. Full model: ε_TF-TF = -1.0, ε_TF-enhancer = -2.0
# 2. No condensation: ε_TF-TF = 0.0, ε_TF-enhancer = -2.0
# 3. No binding: ε_TF-TF = -1.0, ε_TF-enhancer = 0.0
```

### 5.3 Replication & Statistics

**When needed**:
- Phase transitions (large fluctuations)
- Stochastic binding (few events per run)
- Small signals (need to average out noise)

**When not needed**:
- Smooth dose-response curves (deterministic limit)
- Large systems (self-averaging)
- Exploratory parameter scans

**Typical approach**:
```julia
nrepeats = 3  # For publication-quality data
nrepeats = 1  # For parameter exploration
```

---

## 6. Common Experimental Patterns

### 6.1 Phase Diagram

**Goal**: Map out condensation regimes in parameter space

**Setup**:
- Sweep temperature (via `ε = ε_0 / T`) or density
- Measure: Density variance, cluster sizes
- Output: Phase boundary separating gas/liquid/solid

**Application**: Understanding when condensates form

### 6.2 Dose-Response Curve

**Goal**: Quantify how output (binding, transcription) depends on input (TF concentration)

**Setup**:
- Sweep TF particle count
- Measure: Binding occupancy, event frequency
- Output: Sigmoidal curve → Hill coefficient, EC50

**Application**: Characterizing sensitivity and dynamic range

### 6.3 Geometry Scan

**Goal**: Test effect of spatial organization (distances, orientations)

**Setup**:
- Sweep distance between enhancer and promoter
- Measure: RNAP binding, TF density profiles
- Output: Optimal spacing for communication

**Application**: Architectural constraints on gene regulation

### 6.4 Perturbation Response

**Goal**: Test robustness or adaptation

**Setup**:
- Run system to equilibrium
- Change parameter (TF concentration)
- Measure: Relaxation time to new equilibrium
- Output: Timescales for sensing

**Application**: Kinetic aspects of information processing

---

## 7. Biological Interpretation

### 7.1 Transcription Factor (TF)

**Role**: Regulatory proteins that bind DNA to control gene expression

**Model representation**:
- `FastBulk` species (free monomers)
- Self-attractive (`ε < 0`) → can form condensates
- Binds enhancers with moderate affinity

**Biological parameter range**:
- Copy number: 10-10,000 molecules per cell
- Concentration: nM to μM range
- Binding affinity: Kd ~ 1-100 nM (ε ~ -1 to -3)

### 7.2 RNA Polymerase (RNAP)

**Role**: Enzyme that transcribes DNA to RNA

**Model representation**:
- `FastBulk` species
- Lower copy number than TF (typically 10-100)
- Binds promoters strongly when recruited

**Biological context**:
- Often recruited by TF-enhancer complexes
- Cooperative binding with cofactors

### 7.3 Enhancers and Promoters

**Role**: DNA regulatory sequences

**Model representation**:
- `FastBulk` with `initShape="line"` (rigid, immobile)
- Enhancers: Bind TF
- Promoters: Bind RNAP

**Biological context**:
- Enhancers can be distant from promoters (kb-Mb in genome)
- 3D proximity enables regulation (chromosome looping)
- Model simplification: Place nearby on lattice

### 7.4 Condensates in Biology

**Known examples**:
- Transcriptional condensates (Mediator, Pol II clusters)
- Nucleoli (ribosome biogenesis)
- P granules (C. elegans germline)

**Hypothesized functions**:
- Concentrating factors (mass action)
- Excluding inhibitors (partitioning)
- Sensing concentration changes (Fisher information)

**Model goals**:
- Test whether condensation enhances sensing
- Identify parameter regimes where it matters
- Generate predictions for experiments

---

## 8. Validation & Sanity Checks

### 8.1 Physical Realism

**Energy conservation**:
- Hamiltonian should be extensive (scale with system size)
- No perpetual motion (H should fluctuate around mean, not monotonic drift)

**Detailed balance**:
- Metropolis algorithm ensures this (built into `move!()`)
- Acceptance rate should be symmetric for +ΔE and -ΔE moves

**Ergodicity**:
- System must explore configuration space
- Check: Do different initial conditions converge to same observables?

### 8.2 Numerical Stability

**Overflow check**:
- If `|ε| × N_neighbors > 20`, exponentials in Metropolis may overflow
- Manifests as NaN in Hamiltonian or acceptance rate

**Underflow check**:
- Very weak interactions (|ε| < 0.01) may be rounded to zero
- System behaves as non-interacting gas

### 8.3 Convergence Tests

**Equilibration**:
- Plot H(t): Should plateau after burn-in
- Plot ⟨observable⟩ in time windows: Should be constant

**Finite-size effects**:
- Run same experiment with L = 20 and L = 32
- Results should be qualitatively similar (quantitative deviations OK)

**Sampling sufficiency**:
- Autocorrelation time τ: Measure from time series
- Effective samples: N_eff ~ T_total / τ
- Rule: N_eff > 100 for reliable statistics

---

**End of Domain Knowledge**

**Key Takeaway**: Good experimental design requires balancing physical realism (appropriate ε, ϕ, L), computational feasibility (t_max, step), and scientific goals (testable hypotheses, informative metrics).
