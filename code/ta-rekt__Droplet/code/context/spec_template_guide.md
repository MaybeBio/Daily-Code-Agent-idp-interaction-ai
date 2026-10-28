# Specification File Template Guide

**Last Updated**: January 24, 2026  
**Purpose**: Reference for generating valid Julia specification files for Monte Carlo simulations

---

## 0. Key Principles

### 0.1 Minimize Changes

**When creating new spec files**: Always start from the most similar existing spec file and modify as little as possible. Copy-paste liberally from existing specs rather than writing from scratch.

**Workflow**:
1. Identify most similar existing spec (use `experiment_types.md`)
2. Copy that file as starting point
3. Modify only what's necessary for new experiment
4. Keep same structure, naming conventions, patterns

### 0.2 Start Small, Scale Up

**CRITICAL**: First proposals always use minimal resources for rapid iteration.

**Initial test configuration**:
- **t_max = 1**: 1 minute runtime (NOT 5, NOT 10)
- **ncpus = 1**: Single CPU (NOT parallel)
- **Reason**: Fast validation of proposal, spec syntax, physics, and visualization

**Override template values**: When copying from existing spec file, ALWAYS replace runtime/CPU with small initial values, regardless of what the original used.

**Scaling sequence**:
1. **Test**: 1 min, 1 CPU → validate everything works
2. **Iterate**: Fix errors, refine visualization
3. **Medium**: 5-10 min, 1-4 CPUs → check equilibration
4. **Production**: 10+ min, 16+ CPUs → final data collection

### 0.3 Development Workflow (from researcher)

1. **dev.jl phase**: Develop new features, test with quick runs
2. **viz.ipynb iteration**: Visualize results, refine parameters
3. **spec file creation**: Copy from dev.jl or similar spec → create spec file
4. **Test runs**: Small scale, iterate on errors (5-50 iterations typical)
5. **Scale up**: More CPUs, larger parameter ranges, longer runtime
6. **Production**: Deploy to cluster if needed

**Common failure modes**:
- Hamiltonian not converging (insufficient steps)
- JSON read/write errors (spec file → viz.ipynb)
- Syntax errors in spec file
- Incorrect parameter ranges (doesn't show desired behavior)

---

## 1. File Structure Overview

Every spec file has three main sections:

```julia
### Physical Parameters
# Define lattice and species

### Simulation Parameters  
# Runtime, sweep ranges, parallelization

### Functions
# target!() function with Monte Carlo loop
```

---

## 2. Physical Parameters Section

### 2.1 Lattice Definition

**Always required as first element**:

```julia
l = Lattice(;
    boxlength = 32,      # Linear dimension (total volume = boxlength³)
    nspecies = 8         # Number of distinct species types
)
```

**Guidelines**:
- `boxlength`: Typically 16, 20, 24, 32 (power of 2 or multiple of 4)
- `nspecies`: Count all species in slist (including DNA segments, TF, RNAP, etc.)

### 2.2 Species Types

#### FastBulk (Most Common)

**Default choice** - use this unless you need polymers. (Bulk is deprecated legacy type.)

For non-polymer particles (free monomers, transcription factors, etc.):

```julia
tf = FastBulk(;
    length = 4096,              # Number of particles
    initShape = "uniform"       # Initial distribution
)
```

**Init Shapes**:
- `"uniform"`: Random placement across lattice
- `"line"`: Linear arrangement (requires `initPos`, `initParams`)
- `"box"`: Cubic cluster (requires `initPos`)

**With positioning**:
```julia
dna = FastBulk(;
    length = 13,
    initShape = "line",
    initParams = 1,              # Direction parameter
    initPos = [1, 10, 16]       # [x, y, z] starting position
)
```

#### Bulk (Legacy)

Similar to FastBulk but slower. Use FastBulk unless specific reason.

```julia
b = Bulk(;
    length = 512,
    intPot = -1.0,      # Self-interaction potential
    chemPot = 0.0       # Chemical potential
)
```

#### Polymer

For chain molecules with connectivity constraints:

```julia
poly = Polymer(;
    length = 100,               # Number of beads
    maxlength = 100,
    nMonomerTypes = 2,          # Distinct bead types
    initShape = "line",
    initPos = [8, 8, 8]
)
```

**Note**: Polymers require additional initialization via `initMonomers!()` - see Section 4.

### 2.3 Species List

Always define at end of section:

```julia
slist = [dna1, dna2, enhancer, promoter, tf, rnap]
```

**Order matters**: This list determines species indices used in interactions.

---

## 3. Simulation Parameters Section

### 3.1 Runtime Configuration

```julia
t_max = 10              # Total runtime in minutes
step = 1000             # Sampling interval (Monte Carlo steps)
t_burn_in = 0.8 * t_max # Equilibration time (discard initial data)
```

**Guidelines**:
- `t_max`: Typically 5-20 minutes for production runs
- `step`: 100-5000 depending on system dynamics
- `t_burn_in`: Usually 0.5-0.8 × t_max

### 3.2 Parameter Sweep Setup

**Pattern for parallelized sweeps**:

```julia
vals = range(1300, stop=4000, step=180)  # Parameter values to sweep
alphas = zeros(16)                        # Pre-allocate array
nrepeats = 1                              # Replicates per parameter value

for (i, v) in enumerate(vals)
    alphas[(i-1)*nrepeats+1:i*nrepeats] = ones(nrepeats) .* v
end

betas = [0.0]           # Secondary parameter (often unused)
nparams = length(alphas) # Total number of sweep points
```

**For single parameter runs**:
```julia
alphas = [4096]
betas = [0.0]
nparams = 1
```

### 3.3 Analysis Parameters

```julia
minsize = 5     # Minimum droplet size for cluster detection
```

---

## 4. Target Function Section

### 4.1 Function Signature

**Always use this exact signature**:

```julia
function target!(res, args...)
    t_start = time_ns()  # Start monotonic timer
    t_current = 0

    # Unpack arguments
    slist, alpha, beta = args
    
    # Extract individual species from slist
    dna1, dna2, enhancer, tf, rnap = slist
```

### 4.2 Parameter Assignment

Map sweep parameters to physical quantities:

```julia
    # Example: alpha controls TF particle count
    tf.length = alpha
    
    # Example: alpha controls geometric parameter
    enhancer.initPos[2] = enhancer.initPos[2] + alpha
```

### 4.3 Interaction Potentials

**Self-interactions** (attractive for condensation):
```julia
intPot!(tf, -1.0)       # Negative = attractive
intPot!(rnap, -1.0)
```

**Multi-species interactions**:
```julia
intPot!(tf, enhancer, -2.0)     # TF binds strongly to enhancer
intPot!(tf, dna, -0.2)          # TF weakly binds DNA
intPot!(rnap, promoter, -2.0)   # RNAP binds strongly to promoter
intPot!(tf, rnap, -2.0)         # Cooperative TF-RNAP interaction
```

**Interaction strength guidelines**:
- `-0.1 to -0.5`: Weak/nonspecific interactions
- `-1.0 to -2.0`: Moderate/specific binding
- `-2.0 to -4.0`: Strong/highly specific
- Positive values: Repulsive (rarely used)

### 4.4 Initialization

**Standard pattern**:
```julia
    slist = [dna1, dna2, enhancer, tf, rnap]  # Update with actual species
    movingsps = [(tf, tf.intIdx), (rnap, rnap.intIdx)]  # Species that move
    
    init!(l, slist)  # Initialize lattice and all species
```

**For polymers**, add monomer initialization:
```julia
    initMonomers!(poly, monomer_list)
    init!(l, slist)
```

### 4.5 Data Collection Setup

**Initialize accumulators**:
```julia
    nmoves = [0, 0]  # Move counters (one per moving species)
    frame = 1
    deltaH = zeros(Float64, l.nspecies)
    
    # Metrics
    hams = []           # Hamiltonian time series
    accs = []           # Acceptance rate time series
    binding = []        # Binding event counts
    transcription = []  # Transcription event counts
    
    # Density fields
    tf_density = zeros(Float64, (l.boxlength, l.boxlength, l.boxlength))
    tf_total = 0.0
    niter = 0
```

**Initial Hamiltonian calculation**:
```julia
    for (sp, idx) in tf.intSpecies
        push!(hams, sp.ham)
        ham!(sp)
        sp.deltaH = 0.0
    end
```

### 4.6 Monte Carlo Loop

**Standard time-based loop**:
```julia
    while t_current <= t_max
        t_current = (time_ns() - t_start) / 6e10  # Convert to minutes
        niter += 1
        
        for (s, i) in movingsps
            if move!(s)  # Attempt move, returns true if accepted
                nmoves[i] += 1
                
                # Sample data at specified intervals
                if mod(nmoves[i], step) == 0
                    # Update energy changes
                    deltaH[i] += s.deltaH
                    
                    # Update partner species energies
                    for (p, j) in s.intSpecies
                        if j != s.intIdx
                            deltaH[j] += p.deltaHP
                        end
                    end
                    
                    # Record metrics
                    ham!(s)
                    push!(hams, s.ham)
                    push!(accs, niter == 0 ? 0 : nmoves[i] / niter)
                    
                    deltaH[i] = 0
                    frame += 1
                    
                    # Application-specific metrics
                    push!(binding, nspeciesBound(tf, [enhancer]))
                    push!(transcription, nspeciesBound(rnap, [promoter]))
                    
                    # Spatial density (post burn-in only)
                    if t_current >= t_burn_in
                        for (_, val) in s.pos
                            tf_density[val[1], val[2], val[3]] += 1.0
                            tf_total += 1.0
                        end
                    end
                end
            end
        end
    end
```

**Alternative: Step-based loop** (for faster systems):
```julia
    nsteps_total = 1_000_000
    for step in 1:nsteps_total
        for (s, i) in movingsps
            if move!(s)
                nmoves[i] += 1
                # ... sampling logic ...
            end
        end
    end
```

### 4.7 Results Export

**Normalize accumulated data**:
```julia
    tf_density = tf_density ./ tf_total  # Convert counts to probability
```

**Build results dictionary**:
```julia
    res["species_name"] = tf.name
    res["nframes"] = frame
    res["step"] = step
    res["temperature"] = 1.0
    res["beta"] = beta
    res["interaction_matrix"] = tf.intMat
    res["boxlength"] = l.boxlength
    res["length"] = tf.length
    
    # Experimental metrics
    res["binding"] = binding
    res["transcription"] = transcription
    res["tf_density_lattice"] = tf_density
    res["hamiltonian"] = hams
    res["acceptance_rate"] = accs
    
    # Final state
    coords_tf = Matrix{Int32}(undef, tf.length, 4)
    res["tf_final_coordinates"] = extractCoords!(tf, coords_tf)
    
    # Metadata
    res["runtime_min"] = (time_ns() - t_start) / 6e10
    
    return res
end
```

**Required result fields**:
- `species_name`: Identifier string
- `nframes`: Number of samples collected
- `step`: Sampling interval
- `boxlength`: Lattice size
- `runtime_min`: Actual wall-clock time

**Common optional fields**:
- `binding`, `transcription`: Event time series
- `*_density_lattice`: 3D spatial distributions
- `hamiltonian`: Energy time series
- `acceptance_rate`: Move acceptance fraction
- `*_final_coordinates`: End-state particle positions

---

## 5. Common Patterns

### 5.1 Bulk Phase Separation

**Use case**: Study condensation of single species

```julia
# Simple bulk system
l = Lattice(boxlength=16, nspecies=1)
b = FastBulk(length=2048, initShape="uniform")
slist = [b]

# Sweep over temperature (via interaction strength)
vals = range(0.1, stop=1.6, step=0.1)
# In target!(): intPot!(b, -1.0 / temp)
```

### 5.2 Prewetting on Surfaces

**Use case**: Condensate formation near attractive surfaces (DNA, enhancers)

```julia
# Linear DNA elements + bulk TF
dna = FastBulk(length=13, initShape="line", initPos=[1,10,16])
enhancer = FastBulk(length=6, initShape="line", initPos=[14,10,16])
tf = FastBulk(length=4096, initShape="uniform")

# Strong TF-enhancer binding
intPot!(tf, -1.0)            # Self-condensation
intPot!(tf, enhancer, -2.0)  # Surface binding
intPot!(tf, dna, -0.2)       # Nonspecific background
```

### 5.3 Multi-Component Binding

**Use case**: Competition between species (TF vs RNAP for binding sites)

```julia
# Two mobile species
tf = FastBulk(length=1300)
rnap = FastBulk(length=32)

# Both compete for promoter
intPot!(tf, promoter, -2.0)
intPot!(rnap, promoter, -2.0)
intPot!(tf, rnap, -2.0)  # May cooperate or compete

movingsps = [(tf, tf.intIdx), (rnap, rnap.intIdx)]
```

### 5.4 Geometric Sweeps

**Use case**: Vary distances between elements

```julia
# Two parallel filaments
enhancer1 = FastBulk(length=6, initPos=[13, 10, 16])
enhancer2 = FastBulk(length=6, initPos=[13, 10+alpha, 16])

# Sweep alpha from 1 to 16 (separation distance)
vals = range(1, stop=16, step=1)
```

---

## 6. Debugging Tips

### 6.1 Common Syntax Errors

**Missing semicolon in keyword arguments**:
```julia
# Wrong:
l = Lattice(boxlength = 32, nspecies = 8)

# Correct:
l = Lattice(; boxlength = 32, nspecies = 8)
```

**Incorrect argument unpacking**:
```julia
# Wrong:
slist, alpha = args

# Correct:
slist, alpha, beta = args  # Always include beta
```

### 6.2 Runtime Issues

**Acceptance rate too low** (< 10%):
- System may be jammed (too many particles)
- Interaction potentials too strong (try -1.0 instead of -3.0)

**No sampling data**:
- Check `t_max` is large enough (> 2 minutes)
- Verify `step` divides `nmoves` reasonably
- Ensure species are in `movingsps` list

**NaN in Hamiltonian**:
- Overflow from excessive interaction strength
- Bug in custom interaction logic

### 6.3 Physical Validity

**Overcrowding check**:
```julia
total_particles = sum([s.length for s in slist])
available_sites = l.boxlength^3
filling_fraction = total_particles / available_sites

# Should be < 0.5 for equilibration
```

**Energy scale check**:
```julia
# Typical: |intPot| ~ 0.1 to 2.0 (in units of kT)
# Avoid |intPot| > 5.0 unless very specific reason
```

---

## 7. Quick Reference

### Available Functions (from `src/`)

**Initialization**:
- `init!(lattice, species_list)`: Set up system
- `initMonomers!(polymer, monomer_list)`: For polymers

**Interaction Setup**:
- `intPot!(species, potential)`: Self-interaction
- `intPot!(species1, species2, potential)`: Pairwise

**Monte Carlo**:
- `move!(species)`: Attempt random move, returns Bool
- `ham!(species)`: Calculate Hamiltonian

**Analysis**:
- `nspeciesBound(species1, species2_list)`: Count bound particles
- `extractCoords!(species, array)`: Export positions

### Typical File Sizes

- Simple bulk: ~50 lines
- Prewetting: ~150 lines  
- Complex multi-species: ~200 lines

### Naming Conventions

- Spec files: `category_descriptor.jl` (e.g., `prewet_transcription.jl`)
- Species: Lowercase, descriptive (`tf`, `rnap`, `enhancer`, `dna1`)
- Global params: Lowercase (`l`, `slist`, `alphas`, `t_max`)
- Functions: CamelCase with `!` for mutation (`target!`, `init!`)

---

**End of Template Guide**
