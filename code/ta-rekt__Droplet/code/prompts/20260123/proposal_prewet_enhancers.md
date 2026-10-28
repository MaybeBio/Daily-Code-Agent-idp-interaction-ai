# Experiment Proposal: prewet_enhancers

**Date**: 20251218 (original experiment)  
**Status**: Reference example (reverse-engineered from existing spec)

---

## 1. Scientific Goal

```yaml
hypothesis: "Close spacing between parallel DNA rods allows single bridging TF droplet; far spacing produces separate droplets"
research_question: "How does geometric distance between enhancer sites affect TF condensate morphology?"
expected_outcome: "Transition from single bridging droplet (low distance) to separate droplets (high distance)"
```

---

## 2. Physical Setup

```yaml
lattice:
  boxlength: 32
  nspecies: 7

species:
  - name: tf
    type: FastBulk
    length: 4096
    initShape: uniform
    mobile: true
    
  - name: dna1
    type: FastBulk
    length: 13
    initShape: line
    initPos: [1, 10, 16]
    initParams: 1
    mobile: false
    
  - name: enhancer
    type: FastBulk
    length: 6
    initShape: line
    initPos: [13, 10, 16]
    initParams: 1
    mobile: false
    
  - name: dna2
    type: FastBulk
    length: 13
    initShape: line
    initPos: [19, 10, 16]
    initParams: 1
    mobile: false
    
  - name: dna1_l
    type: FastBulk
    length: 13
    initShape: line
    initPos: [1, "10+alpha", 16]  # Second rod, y-position varies
    initParams: 1
    mobile: false
    
  - name: enhancer_l
    type: FastBulk
    length: 6
    initShape: line
    initPos: [13, "10+alpha", 16]
    initParams: 1
    mobile: false
    
  - name: dna2_l
    type: FastBulk
    length: 13
    initShape: line
    initPos: [19, "10+alpha", 16]
    initParams: 1
    mobile: false

interactions:
  self_interactions:
    - species: tf
      potential: -1.0
      rationale: "Attractive for condensation"
  
  pairwise_interactions:
    - species_pair: [tf, enhancer]
      potential: -2.0
      rationale: "Strong specific binding to enhancer sites"
      
    - species_pair: [tf, enhancer_l]
      potential: -2.0
      rationale: "Strong specific binding to enhancer sites"
      
    - species_pair: [tf, dna1]
      potential: -0.5
      rationale: "Weak nonspecific binding to DNA"
      
    - species_pair: [tf, dna2]
      potential: -0.5
      rationale: "Weak nonspecific binding to DNA"
      
    - species_pair: [tf, dna1_l]
      potential: -0.5
      rationale: "Weak nonspecific binding to DNA"
      
    - species_pair: [tf, dna2_l]
      potential: -0.5
      rationale: "Weak nonspecific binding to DNA"
```

---

## 3. Simulation Parameters

```yaml
runtime:
  t_max: 10  # minutes
  t_burn_in: 0.8  # fraction of t_max
  step: 1000  # MC steps between samples

sweep:
  parameter: "rod_separation (alpha)"
  range:
    start: 1
    stop: 16
    step: 1
  replicates: 1
  total_conditions: 16
  
  note: "alpha modifies y-position of second rod: dna1_l.initPos[2] = 10 + alpha"
```

---

## 4. Execution Configuration

```yaml
execution:
  cpu_count: 16
  estimated_runtime: "~10 minutes"
  output_directory: "data/20251218/prewet_enhancers_32_1-16"

output_data:
  format: "JSON (output_cpu1.json through output_cpu16.json)"
  structure: |
    results[param_idx]["1"] = {
      "species_name": "tf",
      "nframes": int,
      "step": 1000,
      "temperature": 1.0,
      "beta": 0.0,
      "boxlength": 32,
      "length": 4096,
      "density_lattice": [32][32][32] array,
      "hamiltonian": [...],
      "acceptance_rate": [...],
      "final_coordinates": [...],
      "runtime_min": float
    }
  
  key_metrics:
    - density_lattice: "3D TF probability density (post burn-in)"
    - hamiltonian: "Energy time series"
    - acceptance_rate: "Move acceptance time series"
    - final_coordinates: "TF positions at end of simulation"
```

---

## 5. Analysis & Visualization

```yaml
derived_metrics:
  - name: mean_hamiltonian
    calculation: "np.mean(hamiltonian) for each parameter value"
    purpose: "Average energy as function of rod separation"
    
  - name: mean_force
    calculation: "np.gradient(mean_hamiltonian, rod_separation)"
    purpose: "Mean force between rods: F = -dH/dr"

plots:
  - name: mean_force_vs_distance
    type: line
    x_axis:
      metric: rod_separation
      label: "Rod Distance (lattice units)"
      scale: linear
      range: [1, 16]
    y_axis:
      metric: mean_force
      calculation: "gradient of mean hamiltonian"
      label: "Mean Force (dH/dr)"
      scale: linear
      range: auto
    description: "Effective force between rods as function of separation"
    
  - name: density_cross_section_yz
    type: heatmap
    description: "YZ cross-section of TF density at x=16 for various rod separations"
    note: "Show multiple separations (e.g., alpha=1, 5, 10, 16) in subplots"
    
  - name: density_cross_section_xz
    type: heatmap
    description: "XZ cross-section of TF density at y=10+alpha/2 (midplane between rods)"
    note: "Show droplet bridging vs separation"
    
  - name: hamiltonian_time_series
    type: line
    x_axis:
      metric: frame
      label: "Frame Number"
      scale: linear
    y_axis:
      metric: hamiltonian
      label: "Hamiltonian"
      scale: linear
    description: "Check equilibration for selected parameter values"
    note: "Show alpha=1 and alpha=16 to verify convergence"

dashboard:
  filename: "dashboard_prewet_enhancers.py"
  function_name: "dashboard_prewet_enhancers"
  arguments:
    - data_dir: "Path to output directory"
    - param_range: "list(range(1, 17))"
  returns: "matplotlib figure with 4 subplots"
  template: "Create new (no similar dashboard exists yet)"

notebook_entry:
  date_header: "## 20251218"
  experiment_header: "### Prewet Enhancers"
  code: |
    from dashboard_prewet_enhancers import dashboard_prewet_enhancers
    
    data_dir = "data/20251218/prewet_enhancers_32_1-16"
    param_range = list(range(1, 17))
    
    fig = dashboard_prewet_enhancers(data_dir, param_range)
    plt.show()
```

---

## 6. Comparison to Similar Experiments

```yaml
most_similar: prewet_transcription.jl

key_differences:
  - "prewet_transcription has 8 species (TF + RNAP + 2 rods)"
  - "prewet_enhancers has 7 species (TF only + 2 rods)"
  - "prewet_transcription sweeps TF particle count"
  - "prewet_enhancers sweeps geometric parameter (rod separation)"
  - "prewet_transcription measures binding/transcription events"
  - "prewet_enhancers measures only spatial density"

changes_from_template:
  removed:
    - "RNAP species"
    - "promoter DNA segments"
    - "binding and transcription metrics (nspeciesBound calls)"
  
  modified:
    - "Sweep parameter: alpha now modifies geometry (rod y-position) instead of particle count"
    - "Reduced to single mobile species (tf only)"
    - "nspecies: 8 → 7"
  
  kept_same:
    - "Two parallel rod architecture"
    - "TF self-interaction: -1.0"
    - "TF-enhancer interaction: -2.0"
    - "Runtime: 10 min, burn-in: 0.8"
    - "Lattice: 32³"
```

---

## 7. Technical Notes

```yaml
notes:
  - "Rod separation controlled via alpha: dna1.initPos[2] = dna1_l.initPos[2] + alpha"
  - "Both rods have identical structure: dna1 (13) + enhancer (6) + dna2 (13)"
  - "High TF count (4096) ensures significant condensation"
  - "density_lattice is normalized: values represent probability density"
```

---

**End of Proposal**
