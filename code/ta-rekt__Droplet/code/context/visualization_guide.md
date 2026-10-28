# Visualization & Analysis Guide

**Last Updated**: January 24, 2026  
**Purpose**: Guidelines for creating dashboard functions and analysis workflows

---

## 1. Overview

**Pipeline**: Spec file (raw output) → Dashboard function (refined analysis) → viz.ipynb (plots)

Each experiment produces:
1. **Spec file** (`spec/*.jl`): Generates raw JSON output data
2. **Dashboard function** (separate Python file): Loads JSON, computes derived metrics, creates plots
3. **viz.ipynb entry**: Markdown header + code cell calling dashboard function

---

## 2. Dashboard Functions

### 2.1 Key Principle: Reuse and Minimize Changes

**When creating dashboards**: Start from most similar existing dashboard and modify minimally.

**Keep it simple**: Dashboard should do only what's needed - no unnecessary complexity.

### 2.2 File Organization

**Location**: Create separate Python file for each dashboard (allows reuse)

**Naming**: `dashboard_[experiment_name].py`

**Example**: `dashboard_prewet_enhancers.py`

### 2.3 Typical Dashboard Structure

```python
import json
import numpy as np
import matplotlib.pyplot as plt

def dashboard_prewet_enhancers(data_dir, param_range):
    """
    Dashboard for prewet_enhancers experiment.
    
    Args:
        data_dir: Path to data directory (e.g., "data/20251218/prewet_enhancers_...")
        param_range: List of parameter values swept
    
    Returns:
        fig: Matplotlib figure with plots
    """
    # 1. Load data from JSON files
    results = {}
    for cpu in range(1, 17):  # Adjust based on CPU count
        filepath = f"{data_dir}/output_cpu{cpu}.json"
        with open(filepath, 'r') as f:
            data = json.load(f)
            results.update(data)
    
    # 2. Extract raw metrics
    hamiltonians = []
    density_lattices = []
    for param_idx in sorted(results.keys()):
        hamiltonians.append(results[param_idx]['1']['hamiltonian'])
        density_lattices.append(results[param_idx]['1']['density_lattice'])
    
    # 3. Compute derived metrics
    mean_hams = [np.mean(h) for h in hamiltonians]
    mean_force = np.gradient(mean_hams, param_range)  # dH/dr
    
    # 4. Create plots
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    # Plot 1: Mean force vs distance
    axes[0].plot(param_range, mean_force, 'o-')
    axes[0].set_xlabel('Rod Distance (lattice units)')
    axes[0].set_ylabel('Mean Force (dH/dr)')
    axes[0].set_title('Force between rods')
    
    # Plot 2: Example density cross-section
    # ... (add more plots as needed)
    
    plt.tight_layout()
    return fig
```

### 2.4 Common Patterns

**Loading data**:
```python
# Pattern 1: Load from multiple CPU files
results = {}
for cpu in range(1, ncpus+1):
    with open(f"{data_dir}/output_cpu{cpu}.json") as f:
        results.update(json.load(f))

# Pattern 2: Sort by parameter index
sorted_keys = sorted(results.keys(), key=int)
```

**Extracting metrics**:
```python
# Get time series
hamiltonians = [results[k]['1']['hamiltonian'] for k in sorted_keys]

# Get spatial data
density = [results[k]['1']['density_lattice'] for k in sorted_keys]

# Get final states
final_coords = [results[k]['1']['final_coordinates'] for k in sorted_keys]
```

**Derived metrics** (general pattern):
- **Means**: `mean_ham = [np.mean(h) for h in hamiltonians]`
- **Gradients/Forces**: `force = np.gradient(mean_values, param_range)`
- **Correlations**: `corr = np.corrcoef(metric1, metric2)`
- **Statistical measures**: Will be added as experiments evolve

---

## 3. viz.ipynb Structure

### 3.1 Organization

**Hierarchical structure**:
```
# YYYYMMDD (date header, level 1)
## Experiment Name (level 2)
[Code cell with dashboard call]
```

### 3.2 Example Entry

```python
# In viz.ipynb, add:

## 20251218
### Prewet Enhancers

# Import dashboard
from dashboard_prewet_enhancers import dashboard_prewet_enhancers

# Set parameters
data_dir = "data/20251218/prewet_enhancers_32_1-16"
param_range = list(range(1, 17))

# Run dashboard
fig = dashboard_prewet_enhancers(data_dir, param_range)
plt.show()
```

### 3.3 Iterative Refinement

**During development**:
1. Run experiment with small parameter range
2. Create initial dashboard in viz.ipynb
3. Iterate: adjust plots, add/remove metrics
4. Once stable: move dashboard to separate file
5. Scale up experiment, refine dashboard

**Typical iteration**: 5-50 rounds of refinement (same as spec files)

---

## 4. Common Plot Types

### 4.1 Parametric Plots (Most Common)

**Pattern**: Quantity vs swept parameter

```python
# Example: Density vs temperature
plt.plot(temperatures, densities, 'o-')
plt.xlabel('Temperature')
plt.ylabel('Density')
```

### 4.2 Cross-Sectional Density Plots

**Pattern**: 2D slice of 3D density lattice

```python
# YZ cross-section at x=16
density_yz = density_lattice[16, :, :]
plt.imshow(density_yz, origin='lower', cmap='viridis')
plt.colorbar(label='TF Density')
plt.xlabel('Z')
plt.ylabel('Y')
```

### 4.3 Time Series (Future)

**Pattern**: Metric evolution over simulation time

```python
# Example: Hamiltonian vs frame
plt.plot(hamiltonians)
plt.xlabel('Frame')
plt.ylabel('Hamiltonian')
```

---

## 5. Analysis Workflow

### 5.1 Raw Output → Refined Metrics

**Spec file produces (examples, not exclusive)**:
- `hamiltonian`: Energy time series
- `acceptance_rate`: Move acceptance fraction
- `density_lattice`: 3D spatial distribution
- `binding`, `transcription`: Event counts (if applicable)
- `final_coordinates`: End-state positions

**Dashboard computes**:
- **Means**: Time-averaged quantities
- **Derivatives**: Forces, gradients
- **Spatial features**: Cross-sections, profiles
- **Statistical measures**: Variances, correlations (future)

### 5.2 Typical Analysis Chain

```
1. Load JSON data
2. Extract raw metrics (hamiltonian, density, etc.)
3. Compute derived quantities (means, forces, etc.)
4. Create multiple plots in single figure
5. Return figure for display in notebook
```

---

## 6. Best Practices

### 6.1 Code Reuse

- **Copy existing dashboards**: Start from similar experiment
- **Modify minimally**: Change only parameter names, calculations, plot labels
- **Keep structure**: Maintain same function signature, return types

### 6.2 Simplicity

- **One dashboard per experiment**: Don't try to make universal functions
- **Focused purpose**: Only compute/plot what's needed for this experiment
- **Clear variable names**: `mean_force` not `mf`, `rod_separation` not `alpha`

### 6.3 Documentation

**Docstring should include**:
- Brief description of experiment
- Args with types and meanings
- Returns description
- Example usage (optional)

---

## 7. Quick Reference

### Dashboard Template

```python
import json
import numpy as np
import matplotlib.pyplot as plt

def dashboard_[experiment_name](data_dir, param_range, **kwargs):
    """Brief description of experiment."""
    
    # Load data
    results = {}
    for cpu in range(1, ncpus+1):
        with open(f"{data_dir}/output_cpu{cpu}.json") as f:
            results.update(json.load(f))
    
    # Extract + compute
    sorted_keys = sorted(results.keys(), key=int)
    raw_metric = [results[k]['1']['metric_name'] for k in sorted_keys]
    derived_metric = compute_something(raw_metric)
    
    # Plot
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(param_range, derived_metric, 'o-')
    ax.set_xlabel('Parameter')
    ax.set_ylabel('Metric')
    ax.set_title('Title')
    
    return fig
```

### viz.ipynb Entry Template

```markdown
## YYYYMMDD
### Experiment Name
```

```python
from dashboard_experiment_name import dashboard_experiment_name

data_dir = "data/YYYYMMDD/experiment_name_params"
param_range = list(range(start, stop, step))

fig = dashboard_experiment_name(data_dir, param_range)
plt.show()
```

---

**End of Guide**
