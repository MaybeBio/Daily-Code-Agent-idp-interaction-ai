# Experiment Proposal: heterogeneous_polymer_looping

**Date**: 20260128  
**Status**: Draft

---

## 1. Scientific Goal

```yaml
hypothesis: "Bulk species interactions with patterned polymer regions can promote looping (pattern proximity) even at weak inter-pattern interaction strengths"

research_question: "Can differential bulk-polymer interactions compensate for weak direct pattern-pattern attraction to induce polymer looping?"

expected_outcome: |
  At strong inter-pattern attraction (ε = -2.0), patterns should contact frequently regardless of bulk.
  At weak or zero inter-pattern attraction (ε ≥ -0.5), bulk-mediated effects should significantly 
  increase contact frequency and reduce mean pattern distance compared to no-bulk control.
  Time to first contact should decrease with both stronger inter-pattern attraction AND bulk presence.
```

**Scientific Context**:  
In biological systems, polymer looping (e.g., DNA enhancer-promoter contacts) is often mediated by proteins that interact differentially with specific sequence regions. This experiment tests whether a similar mechanism can occur through phase separation: bulk species that preferentially interact with patterned polymer regions may create local concentration gradients that bring distant patterns together, even when direct pattern-pattern interaction is weak or unfavorable.

**Comparison to Control**:  
By comparing systems with and without bulk interactions (control achievable by setting bulk-polymer interaction to 0), we can isolate the contribution of bulk-mediated effects from direct inter-pattern attraction.

---

## 2. Physical Setup

```yaml
lattice:
  boxlength: 42
  nspecies: 4  # base polymer, pattern1, pattern2, bulk

species:
  - name: poly
    type: Polymer
    length: 32  # initial length
    maxlength: 128  # final length
    nMonomerTypes: 3  # base (type 1), pattern1 (type 2), pattern2 (type 3)
    initShape: line
    initPos: [5, 5, 5]  # 5 sites from box edge
    mobile: true
    pattern:
      architecture: |
        Bead positions and types:
        - Beads 1-28: type 1 (base)
        - Beads 29-36: type 2 (pattern 1, length 8)
        - Beads 37-92: type 1 (base, separation length 56)
        - Beads 93-100: type 3 (pattern 2, length 8)
        - Beads 101-128: type 1 (base)
      
      summary: "Two patterned regions (8 beads each) separated by 56 beads, each 28 beads from polymer ends"
    
  - name: bulk
    type: FastBulk
    length: 4096
    initShape: uniform
    mobile: true

interactions:
  interaction_matrix: |
    4x4 matrix (base, pattern1, pattern2, bulk):
         [  0.5    0.5    0.5   -0.1 ]  # base (type 1)
         [  0.5    0.5   SWEEP  -1.0 ]  # pattern1 (type 2)
         [  0.5   SWEEP  0.5   -1.0 ]  # pattern2 (type 3)
         [ -0.1   -1.0   -1.0   -1.5 ]  # bulk
  
  notes:
    - "SWEEP = inter-pattern interaction (pattern1-pattern2), varies from -2.0 to 0.0"
    - "Base self-interaction: 0.5 (slightly repulsive, prevents collapse)"
    - "Pattern self-interactions: 0.5 (same as base)"
    - "Bulk-base: -0.1 (weak attraction)"
    - "Bulk-pattern: -1.0 (strong attraction, 10x stronger than bulk-base)"
    - "Bulk self-interaction: -1.5 (promotes bulk condensation)"
```

**Pattern Architecture Details**:
- Total polymer: 128 beads
- Pattern 1 (type 2): beads 29-36 (length 8)
- Pattern 2 (type 3): beads 93-100 (length 8)
- Separation: 56 base beads between patterns
- Distance from ends: 28 beads from each terminus (symmetric placement)

**Interaction Design Rationale**:
- Repulsive polymer self-interaction (0.5) prevents global collapse, isolating looping to pattern-pattern attraction
- Strong bulk-pattern attraction (-1.0) vs weak bulk-base (-0.1): creates preferential localization of bulk near patterns
- Strong bulk self-interaction (-1.5): ensures bulk phase separates, creating concentrated regions
- Sweep inter-pattern interaction from attractive (-2.0) to neutral (0.0) to test threshold for bulk-mediated looping

---

## 3. Simulation Parameters

```yaml
runtime:
  t_max: 1  # minute (initial test)
  t_burn_in: 0.2  # 0.2 minutes = 12 seconds
  step: 10000  # MC steps between samples

sweep:
  parameter: "inter_pattern_interaction (pattern1-pattern2)"
  range:
    start: -2.0
    stop: 0.0
    step: 0.25
  replicates: 1
  total_conditions: 9  # -2.0, -1.75, -1.5, -1.25, -1.0, -0.75, -0.5, -0.25, 0.0
```

---

## 4. Execution Configuration

```yaml
execution:
  cpu_count: 1
  estimated_runtime: "1 minute"
  output_directory: "data/20260128/heterogeneous_polymer_looping_42_n2-0"

output_data:
  format: "JSON (output_cpu1.json)"
  structure: |
    results[param_idx]["1"] = {
      "species_name": "poly",
      "nframes": int,
      "step": 10000,
      "temperature": 1.0,
      "beta": 0.0,
      "interaction_matrix": [[...]],  # 4x4 matrix
      "boxlength": 42,
      "length_poly": int,  # polymer length time series
      "length_bulk": int,  # bulk length (constant)
      "radius_of_gyration": [...],  # polymer Rg time series
      "pattern_distances": [...],  # euclidean distance between pattern centers
      "contact_events": [...],  # binary: 1 if patterns in contact, 0 otherwise
      "first_contact_frame": int,  # frame number of first contact (-1 if never)
      "hamiltonian_poly": [...],
      "hamiltonian_bulk": [...],
      "acceptance_rate": [...],
      "final_coordinates_poly": [...],
      "final_coordinates_bulk": [...],
      "runtime_min": float
    }
  
  key_metrics:
    - radius_of_gyration: "Polymer Rg time series"
    - pattern_distances: "Euclidean distance between pattern1 and pattern2 centers"
    - contact_events: "Binary time series (1 = in contact, 0 = not in contact)"
    - first_contact_frame: "Frame index when patterns first contact (-1 if never)"
    - hamiltonian_poly: "Polymer energy time series"
    - hamiltonian_bulk: "Bulk energy time series"
```

**Contact Definition**:  
Two patterns are "in contact" if ANY bead from pattern 1 is adjacent (26 neighbors in 3D) to ANY bead from pattern 2.

---

## 5. Analysis & Visualization

```yaml
derived_metrics:
  - name: mean_pattern_distance
    calculation: "np.mean(pattern_distances[burn_in_idx:]) for each parameter"
    purpose: "Average separation between patterns at equilibrium"
    
  - name: contact_frequency
    calculation: "np.mean(contact_events[burn_in_idx:]) for each parameter"
    purpose: "Fraction of time patterns are in contact (0 to 1)"
    
  - name: time_to_first_contact
    calculation: "first_contact_frame * step (convert frame to MC steps)"
    purpose: "Kinetic measure: how quickly do patterns find each other?"
    note: "Value is -1 if no contact occurs during simulation"
    
  - name: mean_hamiltonian_poly
    calculation: "np.mean(hamiltonian_poly[burn_in_idx:])"
    purpose: "Check polymer equilibration"
    
  - name: mean_hamiltonian_bulk
    calculation: "np.mean(hamiltonian_bulk[burn_in_idx:])"
    purpose: "Check bulk equilibration"

plots:
  - name: mean_distance_vs_interaction
    type: scatter
    description: |
      Mean pattern distance vs inter-pattern interaction strength.
      Compare WITH bulk (experiment) vs WITHOUT bulk (control).
      Expect: distance decreases with stronger attraction; bulk presence reduces distance.
    note: "Deferred to later visualization phase"
    
  - name: contact_frequency_vs_interaction
    type: scatter
    description: |
      Contact frequency vs inter-pattern interaction strength.
      Compare WITH bulk vs WITHOUT bulk.
      Expect: frequency increases with stronger attraction; bulk enhances frequency at weak attraction.
    note: "Deferred to later visualization phase"
    
  - name: time_to_contact_vs_interaction
    type: scatter
    description: |
      Time to first contact vs inter-pattern interaction.
      Compare WITH bulk vs WITHOUT bulk.
      Expect: faster contact with stronger attraction; bulk accelerates contact.
    note: "Deferred to later visualization phase"

dashboard:
  filename: "TBD - deferred to post-experiment phase"
  note: "User will handle visualization separately"
```

---

## 6. Comparison to Similar Experiments

```yaml
most_similar: poly_bulk_phases.jl

key_differences:
  - "poly_bulk_phases: homogeneous polymer (single monomer type)"
  - "heterogeneous_polymer_looping: heterogeneous polymer (3 monomer types with spatial pattern)"
  - "poly_bulk_phases: sweeps polymer self-interaction"
  - "heterogeneous_polymer_looping: sweeps inter-pattern interaction (off-diagonal interaction matrix element)"
  - "poly_bulk_phases: measures radius of gyration only"
  - "heterogeneous_polymer_looping: measures pattern distances, contact events, first contact time"
  - "poly_bulk_phases: 16³ lattice"
  - "heterogeneous_polymer_looping: 42³ lattice (larger to accommodate longer polymer)"

changes_from_template:
  added:
    - "Heterogeneous polymer with 3 monomer types (nMonomerTypes = 3)"
    - "Custom pattern specification (two 8-bead regions separated by 56 beads)"
    - "4x4 interaction matrix (3 polymer types + bulk)"
    - "Pattern distance tracking (euclidean distance between pattern centers)"
    - "Contact event detection (binary time series)"
    - "First contact frame metric (kinetic measure)"
  
  modified:
    - "Lattice size: 16 → 42 (larger for longer polymer)"
    - "Polymer initial length: 16 → 32"
    - "Polymer max length: 64 → 128"
    - "Bulk particle count: 1024 → 4096 (maintain similar density)"
    - "Sweep parameter: polymer self-interaction → inter-pattern interaction"
    - "Interaction matrix: 2x2 → 4x4 (heterogeneous polymer)"
    - "Runtime: 20 min → 1 min (initial test)"
  
  kept_same:
    - "Species types: Polymer + FastBulk"
    - "CPU count: 1 (single CPU for initial test)"
    - "Output format: JSON with time series"
    - "Burn-in strategy: fraction of runtime"
```

---

## 7. Technical Notes

```yaml
notes:
  - "INITIAL TEST: 1 minute runtime, 1 CPU - validate before scaling"
  - "Pattern definition uses :custom pattern type with explicit bead type assignments"
  - "Polymer initialization: initMonomers! must be called before init!"
  - "Contact detection: check if any pattern1 bead has pattern2 bead in 26-neighbor shell"
  - "First contact frame: track minimum frame index where contact occurs, -1 if never"
  - "Burn-in: 0.2 min = 12 seconds (conservative for quick test)"
  - "Sampling: step=10000 (coarse, suitable for slow polymer dynamics)"
  - "Lattice size 42³ chosen to give polymer space (max extent ~32 lattice sites for extended chain)"
  - "Bulk density: 4096 / 42³ ≈ 0.055 (moderate filling fraction)"

pattern_setup_code: |
  # Example pattern specification for heterogeneous polymer
  pattern = zeros(Int, 128)
  pattern[1:28] .= 1      # base type
  pattern[29:36] .= 2     # pattern 1
  pattern[37:92] .= 1     # base type
  pattern[93:100] .= 3    # pattern 2
  pattern[101:128] .= 1   # base type

control_experiment:
  description: "Compare with and without bulk-polymer interactions"
  setup: "Set bulk-polymer interactions (matrix[4,1:3]) to 0.0 instead of [-0.1, -1.0, -1.0]"
  expected: |
    Without bulk: looping only occurs at strong inter-pattern attraction (ε < -1.0)
    With bulk: looping occurs even at weak inter-pattern attraction (ε ≈ -0.5)
  comparison: "Demonstrates bulk-mediated looping enhancement"
```

---

**End of Proposal**
