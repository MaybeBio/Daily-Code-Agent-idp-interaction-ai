# Experiment Proposal: [EXPERIMENT_NAME]

**Date**: [YYYYMMDD]  
**Status**: [Draft | Testing | Complete]

---

## 1. Scientific Goal

```yaml
hypothesis: "[What physical/biological principle are you testing?]"

research_question: "[What specific question does this experiment answer?]"

expected_outcome: "[What do you predict will happen and why?]"
```

---

## 2. Physical Setup

```yaml
lattice:
  boxlength: [int - lattice size, typically 16 or 32]
  nspecies: [int - number of distinct species]

species:
  - name: [species_name]
    type: [FastBulk | Polymer]
    length: [int - number of particles/monomers]
    maxlength: [int - POLYMER ONLY: maximum chain length]
    initShape: [uniform | line]
    initPos: [x, y, z]  # optional for uniform
    initParams: [int]   # optional: direction for line (1=x, 2=y, 3=z)
    mobile: [true | false]
  
  # Add more species as needed...

interactions:
  self_interactions:
    - species: [species_name]
      potential: [float or SWEEP]
      rationale: "[Why this value?]"
  
  pairwise_interactions:
    - species_pair: [species1, species2]
      potential: [float]
      rationale: "[Why this value?]"
  
  # Add more interactions as needed...
```

---

## 3. Simulation Parameters

```yaml
runtime:
  t_max: [int - minutes, START SMALL: 1 minute for initial testing]
  t_burn_in: [float - fraction of t_max, typically 0.5-0.8]
  step: [int - MC steps between samples, typically 1000-10000]

sweep:
  parameter: "[Which parameter are you sweeping?]"
  range:
    start: [float]
    stop: [float]
    step: [float]
  replicates: [int - typically 1 for initial testing]
  total_conditions: [int - number of parameter values]
```

---

## 4. Execution Configuration

```yaml
execution:
  cpu_count: [int - START WITH 1 FOR INITIAL TESTING]
  estimated_runtime: "[e.g., '1 minute' - keep short for iteration]"
  output_directory: "data/[YYYYMMDD]/[experiment_name]_[boxlength]_[param_info]"

output_data:
  format: "JSON (output_cpu1.json through output_cpu[N].json)"
  structure: |
    results[param_idx]["1"] = {
      "species_name": "[primary species]",
      "nframes": int,
      "step": int,
      "temperature": 1.0,
      "beta": 0.0,
      "boxlength": int,
      "length": int,
      [... metric time series ...],
      "hamiltonian": [...],
      "acceptance_rate": [...],
      "final_coordinates": [...],
      "runtime_min": float
    }
  
  key_metrics:
    - [metric_name]: "[Description of what this measures]"
    # List all output metrics from spec file
```

---

## 5. Analysis & Visualization

```yaml
derived_metrics:
  - name: [metric_name]
    calculation: "[How to compute from raw data, e.g., np.mean(...)]"
    purpose: "[What does this tell you?]"
  
  # Add more derived metrics...

plots:
  - name: [plot_name]
    type: [scatter | line | heatmap | 3d_scatter]
    x_axis:
      metric: [parameter or derived metric]
      label: "[Axis label with units]"
      scale: [linear | log]
      range: [auto | [min, max]]
    y_axis:
      metric: [derived metric]
      calculation: "[How to compute]"
      label: "[Axis label with units]"
      scale: [linear | log]
      range: [auto | [min, max]]
    description: "[What should this plot show? What patterns indicate success?]"
  
  # Add more plots...

dashboard:
  filename: "dashboard_[experiment_name].py"
  function_name: "dashboard_[experiment_name]"
  arguments:
    - data_dir: "Path to output directory"
    - param_range: "[e.g., list(range(1, 17)) or np.arange(-1.0, 1.25, 0.25)]"
  returns: "matplotlib figure with [N] subplots"
  template: "[Which existing dashboard to copy from, or 'Create new']"

notebook_entry:
  date_header: "## [YYYYMMDD]"
  experiment_header: "### [Descriptive Title]"
  code: |
    from dashboard_[experiment_name] import dashboard_[experiment_name]
    
    data_dir = "data/[YYYYMMDD]/[output_directory]"
    param_range = [...]
    
    fig = dashboard_[experiment_name](data_dir, param_range)
    plt.show()
```

---

## 6. Comparison to Similar Experiments

```yaml
most_similar: [filename.jl - which existing spec is closest?]

key_differences:
  - "[What's different scientifically?]"
  - "[What's different technically?]"
  # List major differences...

changes_from_template:
  added:
    - "[New species, metrics, or features]"
    # List additions...
  
  modified:
    - "[What changed and how]"
    # List modifications...
  
  kept_same:
    - "[What was preserved from template]"
    # List unchanged elements...
```

---

## 7. Technical Notes

```yaml
notes:
  - "ALWAYS START SMALL: 1 minute runtime, 1 CPU for initial testing"
  - "Scale up ONLY after proposal is validated and results look reasonable"
  - "[Any special initialization requirements?]"
  - "[Any technical gotchas or requirements?]"
  - "[Dependencies between species or parameters?]"
  # Add experiment-specific notes...

control_run:  # Optional, if applicable
  description: "[What is the control experiment?]"
  setup:
    - "[How to set up the control]"
  expected: "[What should the control show?]"
  comparison: "[How does control validate the main experiment?]"
```

---

**End of Proposal**
