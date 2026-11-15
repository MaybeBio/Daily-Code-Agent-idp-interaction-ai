# Interaction Maps

<details>
<summary>Relevant source files</summary>

The following files were used as context for generating this wiki page:

- [demo/nucleic_acid/protein_nucleic_acid_binding.ipynb](demo/nucleic_acid/protein_nucleic_acid_binding.ipynb)
- [demo/phase_diagrams/phase_diagram_demo.ipynb](demo/phase_diagrams/phase_diagram_demo.ipynb)
- [demo/protein_matrix/interaction_matrix_demo.ipynb](demo/protein_matrix/interaction_matrix_demo.ipynb)
- [docs/getting_started.rst](docs/getting_started.rst)
- [docs/idr_idr.rst](docs/idr_idr.rst)
- [docs/phase_diagrams.rst](docs/phase_diagrams.rst)
- [finches/interaction_vector.py](finches/interaction_vector.py)

</details>



This section covers the generation and interpretation of 2D interaction maps (intermaps) between protein sequences, particularly intrinsically disordered regions (IDRs). Interaction maps provide a visual representation of predicted intermolecular interactions between two sequences, showing how residues from each sequence are likely to interact with one another.

For information about calculating single epsilon values between sequences, see [Epsilon Calculations](#3.2). For analysis of IDR interactions with folded domains, see [IDR-Folded Domain Analysis](#3.4).

## What are Interaction Maps

Interaction maps (intermaps) are 2D visualizations where each pixel represents the predicted interaction strength between residues from two different sequences. The x-axis corresponds to one sequence, the y-axis to the other sequence, and the color intensity indicates the strength of attraction (typically purple) or repulsion (typically green) between residue pairs.

```mermaid
graph TB
    subgraph "Input"
        S1["Sequence 1<br/>(x-axis)"] 
        S2["Sequence 2<br/>(y-axis)"]
    end
    
    subgraph "Processing"
        SW["Sliding Window<br/>Analysis"]
        EP["Epsilon Parameter<br/>Calculation"]
        IM["Interaction Matrix<br/>Generation"]
    end
    
    subgraph "Output"
        MAP["2D Interaction Map<br/>(intermap)"]
        MAT["Raw Matrix Data"]
    end
    
    S1 --> SW
    S2 --> SW
    SW --> EP
    EP --> IM
    IM --> MAP
    IM --> MAT
```

*Sources: [docs/idr_idr.rst:5-8](), [demo/protein_matrix/interaction_matrix_demo.ipynb:99-101]()*

## Basic Usage

### Generating Interaction Maps

The primary interface for generating interaction maps is through the frontend classes using the `interaction_figure()` method:

```python
from finches import Mpipi_frontend, CALVADOS_frontend

# Initialize frontend objects
mf = Mpipi_frontend()
cf = CALVADOS_frontend()

# Generate interaction map
mf.interaction_figure(sequence1, sequence2)
```

The `interaction_figure()` method creates a visual plot showing the predicted interactions between two sequences.

```mermaid
graph LR
    subgraph "Frontend Interface"
        MF["Mpipi_frontend"]
        CF["CALVADOS_frontend"]
    end
    
    subgraph "Core Methods"
        IF["interaction_figure()"]
        IIM["intermolecular_idr_matrix()"]
    end
    
    subgraph "Backend Processing"
        IMC["InteractionMatrixConstructor"]
        ES["epsilon_stateless functions"]
    end
    
    MF --> IF
    CF --> IF
    MF --> IIM
    CF --> IIM
    IF --> IMC
    IIM --> IMC
    IMC --> ES
```

*Sources: [docs/idr_idr.rst:13-32](), [finches/frontend/](), [demo/protein_matrix/interaction_matrix_demo.ipynb:54-62]()*

### Accessing Raw Matrix Data

For numerical analysis, use `intermolecular_idr_matrix()` to obtain the underlying data:

```python
# Get raw matrix data
(M, d1, d2) = mf.intermolecular_idr_matrix(sequence1, sequence2)
```

Where:
- `M`: Interaction matrix tuple containing sliding epsilon values and position indices
- `d1`: Disorder profile for sequence 1
- `d2`: Disorder profile for sequence 2

*Sources: [docs/idr_idr.rst:40-51]()*

## Customization Parameters

### Window Size

The `window_size` parameter controls the sliding window used for calculation:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `window_size` | 31 | Number of amino acids in sliding window |
| `tic_frequency` | varies | Spacing of axis tick marks |
| `vmin`/`vmax` | model-dependent | Color scale limits |

```python
# Smaller window for shorter sequences
mf.interaction_figure(seq1, seq2, window_size=13, tic_frequency=15)
```

*Sources: [docs/idr_idr.rst:57-62]()*

### Color Scale Adjustment

Interaction maps use clipped color scales with model-specific ranges. Values exceeding the min/max thresholds appear at the extreme colors, which can mask relative differences between strongly interacting regions.

```python
# Custom color scale
mf.interaction_figure(seq1, seq2, vmin=-2, vmax=2)
```

*Sources: [docs/idr_idr.rst:54-56]()*

## Implementation Architecture

```mermaid
graph TB
    subgraph "User Interface"
        UI["interaction_figure()<br/>intermolecular_idr_matrix()"]
    end
    
    subgraph "Matrix Construction"
        IMC["InteractionMatrixConstructor"]
        SW["Sliding Window<br/>Processing"]
        EP["Epsilon Calculation"]
    end
    
    subgraph "Visualization"
        PLT["matplotlib<br/>Plotting"]
        CM["Colormap<br/>Application"]
    end
    
    subgraph "Data Output"
        MAT["Matrix Tuple<br/>(M, d1, d2)"]
        FIG["Figure Object"]
    end
    
    UI --> IMC
    IMC --> SW
    SW --> EP
    EP --> PLT
    EP --> MAT
    PLT --> CM
    CM --> FIG
```

*Sources: [finches/epsilon_calculation.py](), [finches/interaction_vector.py](), High-level system architecture diagram*

## Advanced Features

### Disorder Filtering

Frontend objects can apply disorder filtering to focus analysis on truly disordered regions:

```python
# Enable disorder filtering
mf = Mpipi_frontend(disorder_1=True, disorder_2=True)
result = mf.intermolecular_idr_matrix(seq1, seq2)
```

When enabled, the disorder profiles `d1` and `d2` contain predicted disorder scores rather than all 1s.

### Interaction Vectors

For one-dimensional analysis, use functions from the `interaction_vector` module:

```mermaid
graph LR
    subgraph "Vector Functions"
        SSI["show_sequence_interaction_vector()"]
        SFD["show_folded_domain_interaction_on_sequence()"]
        MIV["make_interaction_vector_plot()"]
    end
    
    subgraph "Core Calculations" 
        GSE["get_sequence_epsilon_vectors()"]
        GIE["get_interdomain_epsilon_vectors()"]
    end
    
    SSI --> GSE
    SFD --> GIE
    GSE --> MIV
    GIE --> MIV
```

*Sources: [finches/interaction_vector.py:100-144](), [finches/interaction_vector.py:19-95]()*

## Output Interpretation

### Matrix Structure

The interaction matrix output `M` is a tuple containing:
1. Matrix of sliding epsilon values
2. Sequence 1 position indices mapping to matrix coordinates  
3. Sequence 2 position indices mapping to matrix coordinates

### Color Coding

- **Purple/Attractive**: Negative epsilon values indicating favorable interactions
- **Green/Repulsive**: Positive epsilon values indicating unfavorable interactions
- **Intensity**: Magnitude of interaction strength within the color scale limits

### Spatial Resolution

Each pixel represents interactions centered on specific residue pairs, smoothed over the sliding window. The effective resolution depends on the window size parameter.

*Sources: [docs/idr_idr.rst:47-51](), [docs/idr_idr.rst:5-8]()*

---