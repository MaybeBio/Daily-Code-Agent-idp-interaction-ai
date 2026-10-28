# Experiment Proposal: poly_bulk_phases

**Date**: 20260124  
**Status**: Draft

---

## 1. Scientific Goal

```yaml
hypothesis: "Bulk-polymer interactions can induce polymer collapse at positive self-interaction values that would not collapse in isolation"

research_question: "Can a polymer with repulsive self-interaction (ε_poly > 0) still collapse via attractive interactions with a condensing bulk species?"

expected_outcome: |
  Polymer with ε_poly = 0.5 (repulsive) should NOT collapse alone, but SHOULD collapse 
  when interacting with bulk species (ε_bulk-poly = -1.0, ε_bulk-bulk = -1.0).
  Radius of gyration should show phase transition despite positive self-interaction.
```

---

## 2. Physical Setup

```yaml
lattice:
  boxlength: 16
  nspecies: 2

species:
  - name: poly
    type: Polymer
    length: 16
    maxlength: 64
    initShape: line
    initPos: [1, 8, 8]
    mobile: true
    
  - name: bulk
    type: FastBulk
    length: 1024
    initShape: uniform
    mobile: true

interactions:
  self_interactions:
    - species: poly
      potential: SWEEP  # -1.0 to 1.0
      rationale: "Sweep parameter - testing collapse with repulsive self-interaction"
      
    - species: bulk
      potential: -1.0
      rationale: "Attractive for bulk condensation"
  
  pairwise_interactions:
    - species_pair: [bulk, poly]
      potential: -1.0
      rationale: "Attractive bulk-polymer interaction to induce collapse"
```

---

## 3. Simulation Parameters

```yaml
runtime:
  t_max: 1  # minute
  t_burn_in: 0.5  # fraction
  step: 10000  # MC steps between samples

sweep:
  parameter: "poly self-interaction potential"
  range:
    start: -1.0
    stop: 1.0
    step: 0.25
  replicates: 1
  total_conditions: 9  # -1.0, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5, 0.75, 1.0
```

---

## 4. Execution Configuration

```yaml
execution:
  cpu_count: 1
  estimated_runtime: "1 minute"
  output_directory: "data/20260124/poly_bulk_phases_16_n1-1"

output_data:
  format: "JSON (output_cpu1.json)"
  structure: |
    results[param_idx]["1"] = {
      "species_name": "poly",
      "nframes": int,
      "step": 10000,
      "temperature": 1.0,
      "beta": 0.0,
      "interaction_matrix": [...],
      "boxlength": 16,
      "length": poly length,
      "radius_of_gyration": [...],  # time series
      "hamiltonian_poly": [...],     # polymer energy time series
      "hamiltonian_bulk": [...],     # bulk energy time series
      "acceptance_rate": [...],
      "final_coordinates_poly": [...],  # polymer bead positions
      "final_coordinates_bulk": [...],  # bulk particle positions
      "runtime_min": float
    }
  
  key_metrics:
    - radius_of_gyration: "Polymer Rg time series"
    - hamiltonian_poly: "Polymer energy time series"
    - hamiltonian_bulk: "Bulk energy time series"
    - final_coordinates_poly: "Final polymer conformation"
    - final_coordinates_bulk: "Final bulk positions"
```

---

## 5. Analysis & Visualization

```yaml
derived_metrics:
  - name: mean_radius_of_gyration
    calculation: "np.mean(radius_of_gyration[burn_in_idx:]) for each parameter"
    purpose: "Equilibrium polymer size as function of self-interaction"
    
  - name: std_radius_of_gyration
    calculation: "np.std(radius_of_gyration[burn_in_idx:])"
    purpose: "Fluctuations in polymer size (indicates phase transition)"
    
  - name: mean_hamiltonian_poly
    calculation: "np.mean(hamiltonian_poly[burn_in_idx:])"
    purpose: "Check equilibration of polymer"
    
  - name: mean_hamiltonian_bulk
    calculation: "np.mean(hamiltonian_bulk[burn_in_idx:])"
    purpose: "Check equilibration of bulk"

plots:
  - name: rg_vs_self_interaction
    type: scatter
    x_axis:
      metric: polymer_self_interaction
      label: "Polymer Self-Interaction (ε_poly)"
      scale: linear
      range: [-1.0, 1.0]
    y_axis:
      metric: mean_radius_of_gyration
      calculation: "mean of Rg after burn-in"
      label: "Radius of Gyration ⟨Rg⟩"
      scale: linear
      range: auto
    description: |
      Should show phase transition: collapsed (low Rg) at negative ε_poly,
      and possibly also at positive ε_poly due to bulk interactions.
      Key test: Does Rg remain low at ε_poly = 0.5?
      
  - name: rg_fluctuations
    type: scatter
    x_axis:
      metric: polymer_self_interaction
      label: "Polymer Self-Interaction (ε_poly)"
      scale: linear
      range: [-1.0, 1.0]
    y_axis:
      metric: std_radius_of_gyration
      calculation: "std of Rg after burn-in"
      label: "Rg Fluctuations (σ_Rg)"
      scale: linear
      range: auto
    description: "Large fluctuations indicate proximity to phase transition"
    
  - name: hamiltonian_equilibration
    type: line
    description: |
      Show hamiltonian time series for polymer and bulk at selected parameters
      (e.g., ε_poly = -1.0, 0.0, 0.5, 1.0) to verify equilibration
    
  - name: final_conformations
    type: 3d_scatter
    description: |
      Visualize final polymer and bulk configurations for key parameter values
      (e.g., ε_poly = 0.5) to confirm visual collapse despite repulsive self-interaction

dashboard:
  filename: "dashboard_poly_bulk_phases.py"
  function_name: "dashboard_poly_bulk_phases"
  arguments:
    - data_dir: "Path to output directory"
    - param_range: "list(np.arange(-1.0, 1.25, 0.25))"
  returns: "matplotlib figure with 4 subplots"
  template: "Create new (closest is bulk_phases, but needs polymer-specific plots)"

notebook_entry:
  date_header: "## 20260124"
  experiment_header: "### Poly-Bulk Phase Transition"
  code: |
    from dashboard_poly_bulk_phases import dashboard_poly_bulk_phases
    import numpy as np
    
    data_dir = "data/20260124/poly_bulk_phases_16_n1-1"
    param_range = list(np.arange(-1.0, 1.25, 0.25))
    
    fig = dashboard_poly_bulk_phases(data_dir, param_range)
    plt.show()
```

---

## 6. Comparison to Similar Experiments

```yaml
most_similar: bulk_phases.jl

key_differences:
  - "bulk_phases: single Bulk species, sweeps interaction potential"
  - "poly_bulk_phases: Polymer + FastBulk, sweeps polymer self-interaction"
  - "bulk_phases: deprecated (uses legacy Bulk type)"
  - "poly_bulk_phases: uses Polymer + FastBulk"
  - "bulk_phases: measures only bulk density"
  - "poly_bulk_phases: measures radius of gyration + both species' energies"

changes_from_template:
  added:
    - "Polymer species with connectivity constraints"
    - "FastBulk species (bulk is now mobile second species)"
    - "Bulk-polymer interaction potential"
    - "radius_of_gyration metric"
    - "Separate hamiltonian tracking for each species"
    - "Final coordinates for both species"
  
  modified:
    - "Sweep parameter: now polymer self-interaction (not bulk)"
    - "nspecies: 1 → 2"
    - "Updated from deprecated Bulk to FastBulk"
  
  kept_same:
    - "Lattice size: 16³"
    - "Runtime: 1 min (reduced from template's 10 min for initial testing)"
    - "Acceptance rate tracking"
```

---

## 7. Technical Notes

```yaml
notes:
  - "bulk_phases.jl is deprecated - this updates to modern FastBulk + Polymer types"
  - "Control experiment: Also run WITHOUT bulk species (bulk.length = 0) at ε_poly = 0.5 to verify no collapse"
  - "Key test case: ε_poly = 0.5 should show collapsed Rg WITH bulk, expanded Rg WITHOUT bulk"
  - "Polymer initialization requires initMonomers! call before init!"
  - "Track both species' hamiltonians separately for equilibration verification"
  - "Burn-in fraction = 0.5 (conservative, polymers equilibrate slowly)"
  - "Sampling interval = 10000 steps (coarse, suitable for slow polymer dynamics)"
  - "INITIAL TEST: 1 minute runtime, 1 CPU - iterate on proposal before scaling up"

control_run:
  description: "Sanity check - polymer alone at ε_poly = 0.5"
  setup:
    - "Same polymer (length=64, ε_poly=0.5)"
    - "NO bulk species (remove from slist)"
  expected: "Rg should be LARGE (expanded/swollen polymer)"
  comparison: "Rg_control >> Rg_with_bulk demonstrates bulk-induced collapse"
```

---

# End of Proposal
