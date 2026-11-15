# Frontend Classes

<details>
<summary>Relevant source files</summary>

The following files were used as context for generating this wiki page:

- [docs/api.rst](docs/api.rst)
- [docs/conf.py](docs/conf.py)
- [finches/data/fingerprints.py](finches/data/fingerprints.py)
- [finches/epsilon_stateless.py](finches/epsilon_stateless.py)
- [finches/frontend/calvados_frontend.py](finches/frontend/calvados_frontend.py)
- [finches/frontend/frontend_base.py](finches/frontend/frontend_base.py)
- [finches/frontend/mpipi_frontend.py](finches/frontend/mpipi_frontend.py)

</details>



This page provides detailed API documentation for the frontend classes in FINCHES. These classes serve as the primary user interface for calculating epsilon values, generating interaction matrices, and creating visualization plots. The frontend classes abstract away the complexity of the underlying calculation engine and provide a consistent, easy-to-use interface regardless of which forcefield model is being used.

For information about the underlying calculation engine, see [Calculation Engine](#5.2). For examples of using these classes, see [Basic Usage Examples](#6.1).

## Class Hierarchy and Architecture

The frontend system is built around an inheritance hierarchy with a common base class that ensures consistent interfaces across different forcefield implementations.

### Frontend Class Inheritance

```mermaid
classDiagram
    class FinchesFrontend {
        <<abstract>>
        +IMC_object: InteractionMatrixConstructor
        +intermolecular_idr_matrix()
        +epsilon()
        +epsilon_vectors()
        +interaction_figure()
        +per_residue_attractive_vector()
        +per_residue_repulsive_vector()
        +protein_nucleic_vector()
        +plot_protein_nucleic_vector()
    }
    
    class Mpipi_frontend {
        +model: Mpipi_model
        +__init__(salt, dielectric)
        +intermolecular_idr_matrix()
        +interaction_figure()
    }
    
    class CALVADOS_frontend {
        +model: calvados_model
        +__init__(salt, pH, temp)
        +intermolecular_idr_matrix()
        +epsilon()
        +interaction_figure()
        +protein_nucleic_vector()
    }
    
    FinchesFrontend <|-- Mpipi_frontend
    FinchesFrontend <|-- CALVADOS_frontend
    
    note for CALVADOS_frontend "Uses @RNA_check decorator\nto prevent RNA handling"
    note for Mpipi_frontend "Supports RNA sequences\n(poly-U)"
```

Sources: [finches/frontend/frontend_base.py:20-41](), [finches/frontend/mpipi_frontend.py:12-21](), [finches/frontend/calvados_frontend.py:28-43]()

### Component Integration

```mermaid
graph TB
    subgraph "Frontend Layer"
        MF["Mpipi_frontend"]
        CF["CALVADOS_frontend"]
        FB["FinchesFrontend"]
    end
    
    subgraph "Model Layer"
        MM["Mpipi_model"]
        CM["calvados_model"]
    end
    
    subgraph "Calculation Engine"
        IMC["InteractionMatrixConstructor"]
        ES["epsilon_stateless functions"]
    end
    
    subgraph "External Dependencies"
        META["metapredict"]
        MPL["matplotlib"]
        NP["numpy"]
    end
    
    FB --> MF
    FB --> CF
    MF --> MM
    CF --> CM
    MF --> IMC
    CF --> IMC
    IMC --> ES
    
    MF --> META
    CF --> META
    FB --> MPL
    FB --> NP
    
    note1["Handles RNA sequences\n(sequences with 'U')"]
    note2["RNA sequences blocked\nvia @RNA_check decorator"]
    
    MF -.-> note1
    CF -.-> note2
```

Sources: [finches/frontend/frontend_base.py:1-11](), [finches/frontend/mpipi_frontend.py:4-21](), [finches/frontend/calvados_frontend.py:1-43]()

## Base Class: FinchesFrontend

The `FinchesFrontend` class serves as an abstract base class that defines the common interface for all frontend implementations. It cannot be instantiated directly and enforces implementation of required attributes in derived classes.

### Core Attributes

| Attribute | Type | Description |
|-----------|------|-------------|
| `IMC_object` | `InteractionMatrixConstructor` | Core calculation engine, must be defined in subclasses |

### Key Methods Overview

```mermaid
graph LR
    subgraph "Calculation Methods"
        EM["epsilon()"]
        EV["epsilon_vectors()"]
        IIDM["intermolecular_idr_matrix()"]
    end
    
    subgraph "Visualization Methods"
        IF["interaction_figure()"]
        PPNV["plot_protein_nucleic_vector()"]
    end
    
    subgraph "Analysis Methods"
        PRAV["per_residue_attractive_vector()"]
        PRRV["per_residue_repulsive_vector()"]
        PNV["protein_nucleic_vector()"]
    end
    
    FB["FinchesFrontend"] --> EM
    FB --> EV
    FB --> IIDM
    FB --> IF
    FB --> PPNV
    FB --> PRAV
    FB --> PRRV
    FB --> PNV
```

Sources: [finches/frontend/frontend_base.py:206-282](), [finches/frontend/frontend_base.py:288-565](), [finches/frontend/frontend_base.py:570-828]()

### Core Calculation Methods

#### epsilon()
Calculates the overall epsilon value between two sequences.

**Parameters:**
- `seq1` (str): First input sequence
- `seq2` (str): Second input sequence  
- `use_aliphatic_weighting` (bool): Enable aliphatic residue weighting (default: True)
- `use_charge_weighting` (bool): Enable charge-based weighting (default: True)

**Returns:** Float epsilon value

#### intermolecular_idr_matrix()
Generates a 2D interaction matrix using sliding window analysis.

**Key Parameters:**
- `window_size` (int): Sliding window size (default: 31)
- `use_cython` (bool): Use optimized Cython implementation (default: True)
- `disorder_1`, `disorder_2` (bool): Generate disorder profiles (default: True)
- `null_shuffle` (bool/int): Apply sequence shuffling for baseline correction

**Returns:** Tuple containing interaction matrix, disorder profile 1, disorder profile 2

#### epsilon_vectors()
Returns attractive and repulsive interaction vectors between sequences.

**Returns:** Tuple of (attractive_vector, repulsive_vector)

Sources: [finches/frontend/frontend_base.py:206-242](), [finches/frontend/frontend_base.py:46-199](), [finches/frontend/frontend_base.py:246-282]()

### Visualization Methods

#### interaction_figure()
Creates comprehensive interaction matrix plots with disorder profiles.

**Key Parameters:**
- `vmin`, `vmax` (float): Color scale limits (default: -3, 3)
- `cmap` (str): Colormap for visualization (default: 'PRGn')
- `seq1_domains`, `seq2_domains` (list): Domain boundaries for highlighting
- `fname` (str): Output filename (None for display)
- `zero_folded` (bool): Zero out folded regions (default: True)

**Returns:** Tuple of matplotlib figure and axes objects

Sources: [finches/frontend/frontend_base.py:288-565]()

## Mpipi_frontend

The `Mpipi_frontend` class implements the Mpipi forcefield model for protein-protein and protein-RNA interactions.

### Constructor Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `salt` | float | 0.150 | Salt concentration in M |
| `dielectric` | float | 80.0 | Dielectric constant |

### RNA Sequence Support

Mpipi supports RNA sequences containing 'U' nucleotides. The implementation automatically detects RNA sequences and adjusts disorder profile generation accordingly.

```mermaid
flowchart TD
    INPUT["Input Sequences"]
    CHECK{"Contains 'U'?"}
    RNA_PATH["RNA Handling"]
    PROTEIN_PATH["Protein Handling"]
    DISORDER_OFF["disorder_profile = False"]
    DISORDER_ON["disorder_profile = True"]
    
    INPUT --> CHECK
    CHECK -->|Yes| RNA_PATH
    CHECK -->|No| PROTEIN_PATH
    RNA_PATH --> DISORDER_OFF
    PROTEIN_PATH --> DISORDER_ON
```

Sources: [finches/frontend/mpipi_frontend.py:13-21](), [finches/frontend/mpipi_frontend.py:105-126]()

### Method Overrides

The `Mpipi_frontend` overrides specific methods to handle RNA sequences:

- `intermolecular_idr_matrix()`: Automatically disables disorder profiles for sequences containing 'U'
- `interaction_figure()`: Adjusts disorder profile generation for RNA sequences

Sources: [finches/frontend/mpipi_frontend.py:27-126](), [finches/frontend/mpipi_frontend.py:128-291]()

## CALVADOS_frontend

The `CALVADOS_frontend` class implements the CALVADOS2 forcefield model for protein-protein interactions only.

### Constructor Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `salt` | float | 0.150 | Salt concentration in M |
| `pH` | float | 7.4 | Solution pH |
| `temp` | float | 288 | Temperature in K |

### RNA Sequence Restrictions

CALVADOS does not support RNA sequences. The implementation uses the `@RNA_check` decorator to prevent RNA handling.

```mermaid
graph TD
    INPUT["Method Call"]
    DECORATOR["@RNA_check"]
    CHECK{"Contains 'U'?"}
    ERROR["ValueError:\nCALVADOS2 cannot handle RNA"]
    PROCEED["Continue Execution"]
    
    INPUT --> DECORATOR
    DECORATOR --> CHECK
    CHECK -->|Yes| ERROR
    CHECK -->|No| PROCEED
```

**RNA Check Implementation:**
```python
def RNA_check(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        if 'U' in args[1] or 'U' in args[2]:
            raise ValueError("CALVADOS2 cannot handle RNA ('U')")
        return func(*args, **kwargs)
    return wrapper
```

Sources: [finches/frontend/calvados_frontend.py:16-26](), [finches/frontend/calvados_frontend.py:34-43]()

### Decorated Methods

The following methods are protected by the `@RNA_check` decorator:

- `intermolecular_idr_matrix()`
- `epsilon()`
- `interaction_figure()`

The `protein_nucleic_vector()` method is stubbed and raises an exception since CALVADOS cannot handle RNA.

Sources: [finches/frontend/calvados_frontend.py:51-139](), [finches/frontend/calvados_frontend.py:145-184](), [finches/frontend/calvados_frontend.py:344-349]()

## Method Parameter Patterns

### Common Parameters Across Methods

| Parameter | Type | Default | Purpose |
|-----------|------|---------|---------|
| `seq1`, `seq2` | str | Required | Input amino acid sequences |
| `use_aliphatic_weighting` | bool | True | Weight aliphatic residue clusters |
| `use_charge_weighting` | bool | True | Weight charged residue interactions |
| `use_cython` | bool | True | Use optimized Cython implementation |
| `window_size` | int | 31 | Sliding window size for matrix calculations |

### Visualization-Specific Parameters

| Parameter | Type | Default | Purpose |
|-----------|------|---------|---------|
| `vmin`, `vmax` | float | Model-dependent | Color scale limits |
| `cmap` | str | 'PRGn' | Matplotlib colormap |
| `tic_frequency` | int | 100 | Axis tick spacing |
| `fname` | str | None | Output filename |
| `zero_folded` | bool | True | Zero out predicted folded regions |

Sources: [finches/frontend/frontend_base.py:46-120](), [finches/frontend/frontend_base.py:317-396]()

## Usage Patterns

### Basic Epsilon Calculation
```python
# Mpipi model
mf = Mpipi_frontend(salt=0.150, dielectric=80.0)
epsilon_value = mf.epsilon(seq1, seq2)

# CALVADOS model  
cf = CALVADOS_frontend(salt=0.150, pH=7.4, temp=288)
epsilon_value = cf.epsilon(seq1, seq2)
```

### Interaction Matrix Generation
```python
# Generate matrix with disorder profiles
matrix_data = mf.intermolecular_idr_matrix(seq1, seq2, window_size=31)
interaction_matrix, disorder1, disorder2 = matrix_data

# Create visualization
fig_data = mf.interaction_figure(seq1, seq2, fname="output.png")
```

### Model-Specific Considerations
```python
# Mpipi supports RNA
rna_seq = "UUUUUUUUUUUUUUUUUUUU"
epsilon_rna = mf.epsilon(protein_seq, rna_seq)  # Works

# CALVADOS blocks RNA
try:
    epsilon_rna = cf.epsilon(protein_seq, rna_seq)  # Raises ValueError
except ValueError as e:
    print("CALVADOS cannot handle RNA sequences")
```

Sources: [finches/frontend/mpipi_frontend.py:1-292](), [finches/frontend/calvados_frontend.py:1-350]()

---