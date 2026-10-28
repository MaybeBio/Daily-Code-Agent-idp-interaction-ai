# Experiment Type Catalog

**Last Updated**: January 23, 2026  
**Purpose**: Reference for existing experiment patterns and when to use them

---

## Overview

This document catalogs the experiment types currently in `spec/`, organized by scientific goal. Use this to identify which pattern fits a new experimental idea.

---

## 1. Bulk Phase Studies

### 1.1 bulk_phases.jl

**Scientific Goal**: Study phase separation of a single species as function of interaction strength

**Physical Setup**:
- Lattice: 16³ sites, 1 species
- Species: Bulk monomer (512 particles, ~12% filling)
- No spatial structure (uniform initial distribution)

**Parameter Sweep**:
- `alpha`: Interaction potential (-1.0 to +1.0 in 0.25 steps)
- `beta`: Chemical potential (fixed at 1.0)

**Key Metrics**:
- Hamiltonian (phase transition indicator)
- Acceptance rate (dynamics measure)
- Spatial density distribution

**When to Use**:
- Testing basic condensation physics
- Validating interaction potential implementations
- Quick sanity checks

**Runtime**: Fast (~5 min on 1 CPU)

---

### 1.2 flory_huggins.jl

**Scientific Goal**: Temperature-dependent phase separation (Flory-Huggins theory)

**Physical Setup**:
- Lattice: 16³ sites, 1 species
- Species: FastBulk (2048 particles, ~50% filling)
- Uniform initial distribution

**Parameter Sweep**:
- `alpha`: Effective temperature (0.1 to 1.6 in 0.1 steps)
- Interaction: `intPot = -1.0 / temp` (temperature scaling)

**Key Metrics**:
- Hamiltonian vs temperature
- Droplet size distribution (minsize = 50)
- Density fluctuations

**When to Use**:
- Mapping phase diagrams (gas-liquid transitions)
- Understanding temperature effects
- Comparing to mean-field theory

**Runtime**: Moderate (10 min × 16 CPUs)

---

### 1.3 poly_phases.jl

**Scientific Goal**: Phase behavior of polymer chains

**Physical Setup**:
- Lattice: Polymer species with connectivity constraints
- Interactions: Bead-bead attractions

**Parameter Sweep**:
- Typically interaction strength or temperature

**Key Metrics**:
- Chain conformations (radius of gyration)
- Inter-chain clustering
- Phase separation kinetics

**When to Use**:
- Studying polymeric condensates
- Chromatin models
- Systems where chain connectivity matters

**Runtime**: Slow (polymers have restricted moves)

---

## 2. Prewetting Studies

**Concept**: Surface-induced condensation before bulk phase separation

### 2.1 prewetting_depletion.jl

**Scientific Goal**: Depletion layer dynamics near surfaces

**Physical Setup**:
- Lattice: Multi-species with rigid surface elements
- Bulk species with surface interactions

**Parameter Sweep**:
- Surface-bulk interaction strength

**Key Metrics**:
- Density profile perpendicular to surface
- Layer thickness vs parameters

**When to Use**:
- Understanding surface effects
- DNA-protein interface studies

---

### 2.2 filament_prewetting.jl

**Scientific Goal**: Condensation along linear DNA-like structures

**Physical Setup**:
- Linear "filaments" (FastBulk with initShape="line")
- Bulk species that binds filaments

**Parameter Sweep**:
- Filament-bulk interaction strength or bulk density

**Key Metrics**:
- Binding occupancy along filament
- Cluster sizes near filament

**When to Use**:
- Modeling DNA-protein interactions
- Chromatin fiber condensation

---

## 3. Regulatory Systems (Gene Regulation Models)

### 3.1 prewet_transcription.jl ⭐ (Complex Example)

**Scientific Goal**: Test if TF condensation enhances RNAP recruitment to promoters

**Physical Setup**:
- Lattice: 32³ sites, 8 species
- DNA architecture:
  - Rod 1: `dna1` - `enhancer` - `dna2` (horizontal)
  - Rod 2: `dna1_l` - `promoter` - `dna2_l` (horizontal, offset by 3 sites)
- Mobile species:
  - TF (1300-4000 particles, sweep parameter)
  - RNAP (32 particles, fixed)

**Interaction Network**:
- TF self-attraction: -1.0 (condensation)
- TF-enhancer: -2.0 (strong binding)
- RNAP-promoter: -2.0 (strong binding)
- TF-RNAP: -2.0 (cooperative)
- Weak nonspecific: -0.2 (TF/RNAP to DNA)

**Parameter Sweep**:
- `alpha`: TF particle count (1300 → 4000, 16 points)
- Hypothesis: More TF → more RNAP binding

**Key Metrics**:
- `binding`: TF-enhancer occupancy over time
- `transcription`: RNAP-promoter occupancy over time
- `tf_density_lattice`: 3D spatial distribution
- `rnap_density_lattice`: 3D spatial distribution

**Expected Results**:
- Sigmoidal dose-response: RNAP binding vs TF concentration
- TF droplet formation at enhancer at high concentration
- Enhanced RNAP recruitment when TF condensed

**When to Use**:
- Testing condensate-mediated gene regulation
- Studying cooperative TF-RNAP interactions
- Dose-response curves for sensing

**Runtime**: Moderate-long (10 min × 16 CPUs ≈ 30 min wall-clock)

**Data Size**: ~50 MB per run (density lattices are large)

---

### 3.2 prewet_enhancers.jl

**Scientific Goal**: Test effect of enhancer-enhancer spacing on TF condensation

**Physical Setup**:
- Lattice: 32³ sites, 7 species
- Two parallel DNA rods with enhancers
- TF (4096 particles, high concentration)

**Parameter Sweep**:
- `alpha`: Vertical separation between rods (1 → 16 lattice units)
- Tests: Can TF droplet bridge two enhancers?

**Interaction Network**:
- TF self-attraction: -1.0
- TF-enhancer: -2.0
- Weak TF-DNA: -0.5

**Key Metrics**:
- TF density distribution
- Droplet morphology (single vs separate droplets)

**Expected Results**:
- Close spacing (< 5 units): Single bridging droplet
- Intermediate spacing: Elongated droplet
- Far spacing (> 10 units): Separate droplets

**When to Use**:
- Studying 3D genome architecture effects
- Enhancer-enhancer communication
- Geometric constraints on condensation

**Runtime**: Moderate (10 min × 16 CPUs)

---

## 4. Pattern Selection Guide

### Decision Tree

**Q1: Single species or multiple?**
- Single → Bulk phase studies (bulk_phases, flory_huggins)
- Multiple → Go to Q2

**Q2: Any polymers (chains)?**
- Yes → poly_phases pattern
- No → Go to Q3

**Q3: Rigid structures (DNA, surfaces)?**
- No → Bulk phase studies
- Yes → Go to Q4

**Q4: What's the main question?**
- Surface layering → prewetting_depletion, filament_prewetting
- Regulatory function → prewet_transcription, prewet_enhancers
- Geometric effects → prewet_enhancers
- Dose-response → prewet_transcription

---

## 5. Common Modifications

### 5.1 Changing Sweep Parameter

**Original** (e.g., prewet_transcription):
```julia
tf.length = alpha  # Sweep TF concentration
```

**Modified** (sweep interaction strength instead):
```julia
tf.length = 2000  # Fixed concentration
intPot!(tf, enhancer, alpha)  # Sweep binding affinity
```

### 5.2 Adding Species

**Steps**:
1. Increase `nspecies` in Lattice
2. Define new species (FastBulk, etc.)
3. Add to `slist`
4. Define interactions with `intPot!()`
5. Update `movingsps` if mobile
6. Add relevant metrics in loop

**Example**: Adding a repressor protein
```julia
# Physical setup
repressor = FastBulk(length=100, initShape="uniform")

# Interactions
intPot!(repressor, promoter, -2.0)  # Competes with RNAP
intPot!(repressor, tf, +0.5)        # Repulsive to TF

# slist
slist = [..., repressor]

# movingsps
movingsps = [(tf, tf.intIdx), (rnap, rnap.intIdx), (repressor, repressor.intIdx)]

# Metrics
repressor_binding = []
push!(repressor_binding, nspeciesBound(repressor, [promoter]))
```

### 5.3 Multi-Dimensional Sweeps

**Avoid unless necessary** (confounds interpretation)

**If required**:
```julia
# 2D sweep: TF concentration × interaction strength
vals_alpha = range(1000, stop=4000, step=500)  # 7 points
vals_beta = range(-0.5, stop=-2.0, step=-0.5)  # 4 points

alphas = []
betas = []
for a in vals_alpha
    for b in vals_beta
        push!(alphas, a)
        push!(betas, b)
    end
end

nparams = length(alphas)  # = 28 total conditions
```

### 5.4 Replicates

**Change**:
```julia
nrepeats = 3  # Instead of 1
```

**Effect**: Each parameter value run 3 times with different random seeds

**When**: Stochastic systems, phase transitions, rare events

---

## 6. Experimental Design Checklist

Before creating a new experiment, verify:

- [ ] **Hypothesis**: Is there a clear, testable prediction?
- [ ] **Metrics**: Are the measurements sufficient to test hypothesis?
- [ ] **Controls**: Can you isolate the effect (e.g., no condensation baseline)?
- [ ] **Parameter range**: Does sweep cover relevant physical regime?
- [ ] **Lattice size**: Is L large enough to avoid finite-size artifacts?
- [ ] **Particle density**: Is filling fraction reasonable (< 0.5)?
- [ ] **Interaction strengths**: Are |ε| values in physical range (0.1 - 2.0)?
- [ ] **Runtime**: Is t_max sufficient for equilibration?
- [ ] **Sampling**: Is step size appropriate for dynamics?
- [ ] **Pattern match**: Does it resemble an existing experiment? Use that as template!

---

## 7. Quick Reference Table

| Experiment | Lattice | # Species | Main Sweep | Key Metric | Runtime |
|------------|---------|-----------|------------|------------|---------|
| bulk_phases | 16³ | 1 | Interaction | Hamiltonian | Fast |
| flory_huggins | 16³ | 1 | Temperature | Droplet size | Medium |
| poly_phases | var | 1 | Interaction | R_g, clusters | Slow |
| prewet_depletion | var | 2+ | Surface binding | Density profile | Medium |
| filament_prewet | var | 2+ | Bulk density | Binding occupancy | Medium |
| prewet_transcription | 32³ | 8 | TF count | RNAP binding | Medium-Long |
| prewet_enhancers | 32³ | 7 | Geometry | Droplet shape | Medium |

---

## 8. Future Experiment Ideas

**Potential extensions** (not yet implemented):

1. **Kinetic proofreading**: Multi-step binding with intermediate states
2. **Chromatin loops**: Polymer chains with specific binding sites
3. **Phase separation kinetics**: Nucleation and growth dynamics
4. **Multi-condensate systems**: Competition between different droplets
5. **Dynamic enhancers**: Mobile DNA elements (loop extrusion)
6. **Crowding agents**: Inert volume excluders
7. **Active processes**: Energy-driven non-equilibrium dynamics

---

**End of Catalog**

**Usage**: When planning a new experiment, first check if it fits an existing pattern. If so, copy that spec file and modify. If not, use the checklist and template guide to build from scratch.
