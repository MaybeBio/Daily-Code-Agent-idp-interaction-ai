# Validation Checklist for Experiment Proposals

**Last Updated**: January 23, 2026  
**Purpose**: Ensure experiment proposals are syntactically valid, physically sensible, and scientifically sound

---

## 1. Syntax Validation (Automated)

### 1.1 Julia Compilation Check

**Command**:
```bash
julia --check-bounds=yes --project=@. -e 'include("spec/proposal.jl")'
```

**Expected output**: No errors, only potential warnings

**Common errors**:
- Missing semicolon in keyword arguments: `Lattice(; boxlength=...)`
- Undefined variables: Check all species defined before use
- Type mismatches: Ensure Int32 for indices, Float64 for energies
- Missing closing brackets/parentheses

**Pass criteria**: ✅ Clean compilation

---

### 1.2 Required Structure Check

**Section 1 - Physical Parameters**:
- [ ] Lattice defined first with `boxlength` and `nspecies`
- [ ] All species instantiated (FastBulk, Bulk, or Polymer)
- [ ] `slist` array created with all species

**Section 2 - Simulation Parameters**:
- [ ] `t_max` or `nsteps_total` defined
- [ ] `step` (sampling interval) defined
- [ ] `alphas` array defined
- [ ] `betas` array defined
- [ ] `nparams = length(alphas)` set correctly

**Section 3 - Target Function**:
- [ ] Function signature: `function target!(res, args...)`
- [ ] Argument unpacking: `slist, alpha, beta = args`
- [ ] Interaction potentials set via `intPot!()`
- [ ] `init!(l, slist)` called
- [ ] Monte Carlo loop present
- [ ] Results dictionary populated
- [ ] Function returns `res`

**Pass criteria**: ✅ All sections present with required elements

---

## 2. Physical Validity (Automated)

### 2.1 Particle Density Check

**Calculation**:
```julia
total_particles = sum([s.length for s in slist])
available_sites = boxlength^3
filling_fraction = total_particles / available_sites
```

**Thresholds**:
- ✅ `filling_fraction < 0.5`: Safe (good mobility)
- ⚠️ `0.5 ≤ filling_fraction < 0.7`: Dense (slow equilibration expected)
- ❌ `filling_fraction ≥ 0.7`: Overcrowded (likely jamming)

**Action if warning/fail**:
- Reduce particle counts
- Increase `boxlength`
- Extend `t_max` and `t_burn_in`

---

### 2.2 Interaction Strength Check

**Scan all `intPot!()` calls**:
```julia
# Extract all ε values
```

**Thresholds**:
- ✅ `|ε| ≤ 2.0`: Typical physical range
- ⚠️ `2.0 < |ε| ≤ 4.0`: Strong, verify intentional
- ❌ `|ε| > 5.0`: Extreme (potential overflow/jamming)

**Special cases**:
- Positive ε (repulsive): Rare, ensure intentional
- Zero ε: Verify not a bug (no interaction)

**Action if warning/fail**:
- Justify strong interactions (e.g., "irreversible binding")
- Check for typos (e.g., `-20.0` instead of `-2.0`)

---

### 2.3 Time Scale Check

**Runtime sufficiency**:
- ✅ `t_max ≥ 5 min`: Likely sufficient for most systems
- ⚠️ `2 min ≤ t_max < 5 min`: Short, monitor equilibration
- ❌ `t_max < 2 min`: Too short (except for quick tests)

**Burn-in fraction**:
- ✅ `0.5 ≤ t_burn_in/t_max ≤ 0.8`: Standard
- ⚠️ `t_burn_in/t_max < 0.5`: May include transient data
- ⚠️ `t_burn_in/t_max > 0.8`: Few samples for statistics

**Sampling interval**:
- ✅ `100 ≤ step ≤ 5000`: Reasonable
- ⚠️ `step < 100`: Very fine (large data files)
- ⚠️ `step > 5000`: Coarse (may miss dynamics)

**Action if warning**:
- For slow systems (polymers, dense): Increase `t_max`
- For fast systems (dilute bulk): Can decrease `t_max`

---

### 2.4 Sweep Configuration Check

**Parameter range**:
- ✅ Minimum 5 points: Enough to see trends
- ✅ Maximum 30 points: Keeps runtime reasonable
- ⚠️ < 5 points: May miss features
- ⚠️ > 30 points: Long runtime, verify necessary

**Replicates**:
- ✅ `nrepeats = 1`: Standard for parameter exploration
- ✅ `nrepeats = 3-5`: Good for publication
- ⚠️ `nrepeats > 5`: Likely overkill (unless high variance expected)

**Computational cost estimate**:
```julia
ncpus = [user specified]
nparams = length(alphas)
wall_clock_time = (t_max × nparams) / ncpus  # minutes
```

**Thresholds**:
- ✅ `wall_clock_time < 60 min`: Quick
- ⚠️ `60 min ≤ wall_clock_time < 240 min`: Moderate (few hours)
- ❌ `wall_clock_time > 240 min`: Long (verify necessary)

---

## 3. Quick Functionality Test (Automated)

### 3.1 Single-Parameter Test Run

**Modify proposal temporarily**:
```julia
alphas = [alphas[1]]  # Only first parameter
betas = [betas[1]]
t_max = 0.5  # 30 seconds
nparams = 1
```

**Run**:
```bash
julia --project=@. sim/sweep.jl spec/proposal.jl 1 1 /tmp/test_output
```

**Expected**:
- Execution completes within 60 seconds
- Output file created: `/tmp/test_output/output_cpu1.json`
- JSON is valid (can be parsed)

**Check JSON structure**:
```bash
python3 -c "import json; data=json.load(open('/tmp/test_output/output_cpu1.json')); print(list(data.keys()))"
```

**Expected keys**: `["1"]` (parameter index)

**Nested structure**:
```bash
python3 -c "import json; d=json.load(open('/tmp/test_output/output_cpu1.json')); print(list(d['1']['1'].keys()))"
```

**Expected fields**: At minimum `["species_name", "nframes", "runtime_min", "hamiltonian"]`

**Pass criteria**: ✅ Runs without error, produces valid JSON with expected structure

---

### 3.2 Equilibration Smoke Test

**Parse Hamiltonian time series**:
```python
import json
data = json.load(open('/tmp/test_output/output_cpu1.json'))
hams = data['1']['1']['hamiltonian']
```

**Visual check**:
```python
import matplotlib.pyplot as plt
plt.plot(hams)
plt.xlabel('Frame')
plt.ylabel('Hamiltonian')
plt.title('Quick equilibration check')
plt.savefig('/tmp/ham_test.png')
```

**Expected**:
- Initial transient (first 20% of frames)
- Plateau or fluctuations around mean (remaining frames)

**Red flags**:
- Monotonic trend (not equilibrated)
- NaN or Inf values (numerical overflow)
- Constant value (system frozen)

**Note**: 30-second run won't fully equilibrate, but should show reasonable behavior

---

## 4. Scientific Validity (Researcher Review)

### 4.1 Hypothesis Clarity

**Question**: Is there a single, testable hypothesis?

**Good example**: 
> "TF condensation at enhancers increases RNAP binding to promoters, quantified by a positive correlation between TF concentration and RNAP-promoter occupancy."

**Bad example**:
> "See what happens with different parameters."

**Checklist**:
- [ ] Hypothesis stated clearly in proposal
- [ ] Specific prediction (increases, decreases, non-monotonic, etc.)
- [ ] Metrics chosen directly test hypothesis

---

### 4.2 Parameter Range Justification

**Question**: Does the sweep cover the relevant physical regime?

**Good justification**:
> "TF concentration from 1300-4000 spans dilute (ϕ=0.04) to dense (ϕ=0.12), crossing expected condensation threshold at ϕ≈0.08."

**Bad justification**:
> "Trying a bunch of values."

**Checklist**:
- [ ] Range includes expected transition or effect
- [ ] Range covers control regime (e.g., low concentration baseline)
- [ ] Range includes saturation regime (effect plateaus)

---

### 4.3 Control Experiments

**Question**: Can you isolate the mechanism?

**Example** (testing if condensation matters):
- **Full model**: `intPot!(tf, -1.0)` (condensation ON)
- **Control**: `intPot!(tf, 0.0)` (condensation OFF)
- **Comparison**: Does RNAP binding differ between conditions?

**Checklist**:
- [ ] Null hypothesis condition included (no effect expected)
- [ ] Alternative mechanism tested (e.g., diffusion-only baseline)
- [ ] Positive control (known outcome for validation)

---

### 4.4 Metric Appropriateness

**Question**: Do the metrics answer the question?

**Examples**:
- Hypothesis about binding → Measure `nspeciesBound()`
- Hypothesis about droplet size → Measure density distribution + clustering
- Hypothesis about kinetics → Measure time series, autocorrelation

**Checklist**:
- [ ] Metrics directly relate to hypothesis
- [ ] Sufficient data collected (time series, spatial, final state)
- [ ] Metrics are unambiguous (clear definition)

---

### 4.5 Expected Results

**Question**: What would confirm/reject the hypothesis?

**Good expectation**:
> "Expect sigmoidal dose-response: RNAP binding low at TF<2000, increases sharply at 2000-3000, plateaus at TF>3000. Hill coefficient > 1 indicates cooperativity."

**Bad expectation**:
> "Something interesting should happen."

**Checklist**:
- [ ] Qualitative prediction stated
- [ ] Quantitative threshold specified (if applicable)
- [ ] Failure modes considered (what if hypothesis is wrong?)

---

## 5. Computational Feasibility

### 5.1 Resource Check

**Memory estimate**:
```
Per CPU ≈ (boxlength³ × nspecies × 4 bytes) + (nframes × nmetrics × 8 bytes)
        ≈ few hundred MB per CPU (typically < 1 GB)
```

**Disk space estimate**:
```
Per output file ≈ nframes × nmetrics × 8 bytes + density_lattice (if saved)
                ≈ 1-50 MB per CPU (larger if density lattices saved)
```

**Checklist**:
- [ ] Total memory (ncpus × per-CPU memory) < available RAM
- [ ] Total disk space < available storage
- [ ] Acceptable wall-clock time (see section 2.4)

---

### 5.2 Parallelization Sanity

**Check**:
```julia
nparams % ncpus == 0  # Ideal: evenly divisible
```

**Example**:
- `nparams = 16`, `ncpus = 16` → 1 param per CPU ✅
- `nparams = 16`, `ncpus = 10` → Uneven (6 CPUs do 2 params, 4 CPUs do 1) ⚠️
- `nparams = 16`, `ncpus = 32` → Wasted CPUs (only 16 will be used) ❌

**Action if warning**:
- Adjust `nparams` to be multiple of `ncpus`
- Or adjust `ncpus` to divide `nparams` evenly

---

## 6. Pre-Approval Checklist

**Before submitting proposal to researcher**:

### Automated Checks
- [ ] Julia compilation succeeds
- [ ] Required structure present
- [ ] Filling fraction < 0.7
- [ ] Interaction strengths |ε| < 5.0
- [ ] Runtime t_max ≥ 2 min
- [ ] Quick test run completes
- [ ] Valid JSON output produced

### Researcher Review Items
- [ ] Clear hypothesis stated
- [ ] Parameter range justified
- [ ] Control experiments included
- [ ] Metrics appropriate
- [ ] Expected results specified
- [ ] Computational cost reasonable

---

## 7. Post-Execution Validation

**After simulation completes** (Phase 2/3):

### 7.1 Data Integrity

- [ ] All `output_cpu{N}.json` files present (N = 1 to ncpus)
- [ ] All files parse as valid JSON
- [ ] All parameter indices present (check for crashes)
- [ ] File sizes reasonable (not empty, not absurdly large)

### 7.2 Equilibration Check

**For each parameter**:
- [ ] Hamiltonian reaches plateau
- [ ] Acceptance rate in reasonable range (5-60%)
- [ ] Runtime close to `t_max` (no early termination)

**Red flags**:
- Monotonic H(t): Increase runtime
- Acceptance < 5%: System jammed, reduce density or |ε|
- Acceptance > 80%: Interactions too weak, not interesting regime

### 7.3 Scientific Validation

- [ ] Results qualitatively match expectations (or interesting deviation)
- [ ] Dose-response curves smooth (not jagged from poor statistics)
- [ ] No unphysical values (NaN, Inf, negative densities)

### 7.4 Reproducibility

**If replicates run** (`nrepeats > 1`):
- [ ] Variance within expected range
- [ ] No outliers suggesting bugs

---

## 8. Validation Failure Actions

### Syntax Errors
→ **Fix code**, rerun compilation check

### Physical Invalidity (jamming, overflow)
→ **Adjust parameters** (reduce density or |ε|), regenerate proposal

### Quick Test Fails
→ **Debug Julia code** (print statements, check logic), fix and retest

### Hypothesis Unclear
→ **Refine with researcher**, iterate proposal

### Computational Infeasibility
→ **Reduce scope** (smaller lattice, fewer parameters, shorter runtime), regenerate

### Post-Execution Failures
→ **Diagnose** (poor equilibration, bugs, wrong parameters), propose fix, rerun if necessary

---

## 9. Agent Validation Workflow

**Step-by-step process**:

1. **Generate spec file** from proposal
2. **Run syntax check** (Julia compilation)
   - If fails: Report errors, suggest fixes, iterate
3. **Run physical validity checks** (density, interactions, timescales)
   - If warnings: Report to researcher, get approval to proceed or adjust
4. **Run quick test** (30-second single-parameter run)
   - If fails: Debug code, report issue
   - If succeeds: Include output summary in proposal
5. **Compile validation report**:
   ```
   ✅ Syntax: Clean compilation
   ✅ Physical: Filling fraction = 0.12, all |ε| < 2.0
   ✅ Quick test: Ran successfully, H converged
   ⚠️ Note: 16 parameters × 10 min = 2.7 hours wall-clock
   ```
6. **Present to researcher** with proposal for approval

---

## 10. Validation Report Template

```markdown
# Validation Report: [Experiment Name]

**Date**: [YYYYMMDD]
**Spec File**: `spec/[name].jl`

## Automated Checks

### Syntax ✅ | ⚠️ | ❌
- Julia compilation: [PASS/FAIL]
- Required structure: [PASS/FAIL]
- Details: [Any notes]

### Physical Validity ✅ | ⚠️ | ❌
- Filling fraction: [value] ([status])
- Interaction strengths: [min-max range] ([status])
- Runtime: [t_max] min ([status])
- Details: [Any warnings or notes]

### Quick Test ✅ | ⚠️ | ❌
- Execution: [PASS/FAIL]
- JSON output: [VALID/INVALID]
- Hamiltonian: [Reasonable behavior? Y/N]
- Details: [Any issues]

## Computational Cost

- Parameters: [nparams]
- CPUs: [ncpus]
- Wall-clock time: [estimate]
- Status: [Reasonable/Long/Too Long]

## Researcher Review Required

**Scientific Checks**:
- Hypothesis clarity: [Your assessment]
- Parameter justification: [Your assessment]
- Control experiments: [Present? Y/N]
- Metric appropriateness: [Your assessment]

**Recommendations**:
- [Any suggestions for improvement]

## Overall Status

[READY FOR APPROVAL | NEEDS REVISION | INFEASIBLE]

[Summary statement]
```

---

**End of Validation Checklist**

**Usage**: Agents should run through sections 1-6 automatically before presenting proposals. Researcher reviews section 4 (scientific validity) during approval process. Sections 7-8 guide post-execution validation.
