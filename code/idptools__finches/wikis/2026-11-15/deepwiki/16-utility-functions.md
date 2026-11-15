# Utility Functions

<details>
<summary>Relevant source files</summary>

The following files were used as context for generating this wiki page:

- [finches/interaction_vector.py](finches/interaction_vector.py)
- [finches/parsing_aminoacid_sequences.py](finches/parsing_aminoacid_sequences.py)
- [finches/sequence_tools.py](finches/sequence_tools.py)

</details>



This page documents the supporting utilities and helper functions that provide foundational functionality for sequence analysis, weighting calculations, and visualization across the FINCHES package. These utilities are primarily used internally by the core calculation engine and frontend interfaces.

For complete API documentation of the main calculation classes, see [Calculation Engine](#5.2). For structure analysis utilities, see [Structure Analysis](#5.3).

## Overview

The utility functions in FINCHES are organized into three main categories:

1. **Sequence Analysis Utilities** - Basic sequence property calculations and manipulations
2. **Weighting and Masking Functions** - Advanced weighting schemes for charge and aliphatic clustering effects  
3. **Visualization and Plotting Utilities** - Tools for creating interaction vector plots and sequence visualizations

## Sequence Analysis Utilities

### Core Sequence Properties

The fundamental sequence analysis functions calculate basic biophysical properties of protein sequences.

```mermaid
graph TB
    subgraph "sequence_tools.py"
        NCPR["calculate_NCPR()"]
        FCR["calculate_FCR()"] 
        BOTH["calculate_FCR_and_NCPR()"]
        MASK["mask_sequence()"]
        WINDOW["get_neighbors_window_of3()"]
    end
    
    subgraph "Input"
        SEQ["Amino Acid Sequence"]
        TARGET["Target Residues"]
    end
    
    subgraph "Output"
        NET_CHARGE["Net Charge per Residue"]
        FRAC_CHARGE["Fraction Charged Residues"]
        BINARY_MASK["Binary Mask [0,1,0,1...]"]
        LOCAL_FRAG["Local Fragment"]
    end
    
    SEQ --> NCPR
    SEQ --> FCR
    SEQ --> BOTH
    SEQ --> MASK
    SEQ --> WINDOW
    TARGET --> MASK
    
    NCPR --> NET_CHARGE
    FCR --> FRAC_CHARGE
    BOTH --> NET_CHARGE
    BOTH --> FRAC_CHARGE
    MASK --> BINARY_MASK
    WINDOW --> LOCAL_FRAG
```

**Sources:** [finches/sequence_tools.py:12-114]()

| Function | Purpose | Return Type |
|----------|---------|-------------|
| `calculate_NCPR()` | Net charge per residue (positive - negative)/length | float |
| `calculate_FCR()` | Fraction of charged residues | float |
| `calculate_FCR_and_NCPR()` | Both calculations in one call | list[float, float] |
| `mask_sequence()` | Binary mask for target residues | list[int] |
| `get_neighbors_window_of3()` | Extract 3-residue window around position | str |

### Advanced Sequence Masking

More sophisticated masking operations support clustering analysis and fragment extraction.

```mermaid
graph TB
    subgraph "Fragment Analysis Pipeline"
        MASK["Binary Mask"]
        EXTRACT["extract_fragments()"]
        NEIGHBORS["MASK_n_closest_nearest_neighbors()"]
        GROUPS["Aliphatic Groups"]
    end
    
    subgraph "Parameters"
        MAX_SEP["max_separation=1"]
        MAX_DIST["max_distance=4"]
    end
    
    subgraph "Processing Steps"
        SPLIT["Split by separator"]
        COUNT["Count neighbors"]
        WEIGHT["Apply weighting"]
    end
    
    MASK --> EXTRACT
    EXTRACT --> SPLIT
    SPLIT --> NEIGHBORS
    MAX_SEP --> EXTRACT
    MAX_DIST --> NEIGHBORS
    NEIGHBORS --> COUNT
    COUNT --> WEIGHT
    WEIGHT --> GROUPS
```

**Sources:** [finches/sequence_tools.py:201-317]()

The `extract_fragments()` function splits binary masks into fragments based on gap tolerance, while `MASK_n_closest_nearest_neighbors()` counts local clustering density for each position.

## Weighting and Masking Functions

### Charge Weighting System

The charge weighting system implements context-dependent charge effects by analyzing local charge environments around interacting residue pairs.

```mermaid
graph TB
    subgraph "get_charge_weighted_mask() Pipeline"
        SEQ1["sequence1"]
        SEQ2["sequence2"]
        CHARGED["Check if both charged"]
        FRAGMENT["Extract 6-residue fragment"]
        CALC["Calculate |NCPR/FCR|"]
        WEIGHT["Charge Weight"]
    end
    
    subgraph "Fragment Construction"
        WIN1["get_neighbors_window_of3(i, seq1)"]
        WIN2["get_neighbors_window_of3(j, seq2)"]
        CONCAT["Concatenate fragments"]
    end
    
    subgraph "Matrices"
        ATTR["Attractive Matrix"]
        REPUL["Repulsive Matrix"]
    end
    
    SEQ1 --> CHARGED
    SEQ2 --> CHARGED
    CHARGED --> WIN1
    CHARGED --> WIN2
    WIN1 --> CONCAT
    WIN2 --> CONCAT
    CONCAT --> FRAGMENT
    FRAGMENT --> CALC
    CALC --> WEIGHT
    WEIGHT --> REPUL
    WEIGHT --> ATTR
```

**Sources:** [finches/parsing_aminoacid_sequences.py:25-194]()

The charge weighting algorithm:
1. Identifies charged residue pairs between sequences
2. Extracts 3-residue windows around each charged residue  
3. Concatenates fragments and calculates `|NCPR/FCR|`
4. Uses this as a weighting factor for like-charge cluster effects

### Aliphatic Clustering Weighting

Aliphatic weighting accounts for cooperative hydrophobic effects when aliphatic residues cluster together.

```mermaid
graph TB
    subgraph "Aliphatic Weighting System"
        ALI_GROUPS["get_aliphatic_groups()"]
        ALI_MASK["get_aliphatic_weighted_mask()"]
        MULTIPLIER["multiplier_weighting"]
    end
    
    subgraph "Clustering Algorithm"
        BINARY["mask_sequence(['A','V','I','L','M'])"]
        NEIGHBORS["MASK_n_closest_nearest_neighbors()"]
        GROUP_SIZE["Group sizes 1,2,3"]
    end
    
    subgraph "Weight Matrix"
        SINGLE["1_1: 1.0"]
        MEDIUM["2_2: 1.5"]
        CLUSTER["3_3: 3.0"]
    end
    
    ALI_GROUPS --> BINARY
    BINARY --> NEIGHBORS
    NEIGHBORS --> GROUP_SIZE
    GROUP_SIZE --> ALI_MASK
    ALI_MASK --> MULTIPLIER
    
    MULTIPLIER --> SINGLE
    MULTIPLIER --> MEDIUM
    MULTIPLIER --> CLUSTER
```

**Sources:** [finches/parsing_aminoacid_sequences.py:252-341]()

The aliphatic weighting scheme assigns higher interaction strengths to clusters of aliphatic residues, with weights ranging from 1.0 (isolated) to 3.0 (highly clustered).

### Folded Domain Specific Weighting

The `get_charge_weighted_FD_mask()` function provides specialized charge weighting for folded domain-IDR interactions.

**Sources:** [finches/parsing_aminoacid_sequences.py:200-246]()

Key difference from standard charge weighting:
- Treats folded domain residues individually (no clustering with neighbors)
- Combines single FD residue with 3-residue IDR window
- Used specifically for surface-accessible folded domain (SAFD) calculations

## Visualization and Plotting Utilities

### Interaction Vector Plotting

The interaction vector module provides comprehensive plotting capabilities for visualizing epsilon-based interactions.

```mermaid
graph TB
    subgraph "High-Level Plot Functions"
        FD_PLOT["show_folded_domain_interaction_on_sequence()"]
        SEQ_PLOT["show_sequence_interaction_vector()"]
    end
    
    subgraph "Core Plotting Engine"
        MAKE_PLOT["make_interaction_vector_plot()"]
        HTML_SEQ["show_sequence_HTML()"]
    end
    
    subgraph "Data Processing"
        SAFD["pdb_to_SDFDresidues_and_xyzs()"]
        EPSILON["get_interdomain_epsilon_vectors()"]
        SEQ_EPS["get_sequence_epsilon_vectors()"]
    end
    
    subgraph "Vector Data"
        ATTR["attractive_vector"]
        REPUL["repulsive_vector"]
    end
    
    FD_PLOT --> SAFD
    FD_PLOT --> EPSILON
    SEQ_PLOT --> SEQ_EPS
    
    EPSILON --> ATTR
    EPSILON --> REPUL
    SEQ_EPS --> ATTR
    SEQ_EPS --> REPUL
    
    ATTR --> MAKE_PLOT
    REPUL --> MAKE_PLOT
    MAKE_PLOT --> HTML_SEQ
```

**Sources:** [finches/interaction_vector.py:18-314]()

### Sequence HTML Visualization

The `show_sequence_HTML()` function creates rich HTML representations of protein sequences with customizable coloring and formatting.

**Sources:** [finches/sequence_tools.py:321-455]()

Key features:
- Amino acid specific color coding using `AA_COLOR` dictionary
- Customizable block sizes and line wrapping
- Bold highlighting for specific positions or residue types
- Opacity effects for de-emphasized regions

| Parameter | Purpose | Default |
|-----------|---------|---------|
| `blocksize` | Residues per block | 10 |
| `newline` | Residues per line | 50 |
| `fontsize` | Font size in pixels | 14 |
| `bold_positions` | Positions to bold | [] |
| `colors` | Custom color overrides | {} |

### Advanced Vector Analysis

The `make_interaction_vector_for_folded_domain()` function generates full-length domain interaction vectors by mapping surface-accessible residue data back to complete domain sequences.

**Sources:** [finches/interaction_vector.py:258-315]()

This function:
1. Extracts SAFD sequence and coordinates from PDB
2. Calculates interaction vectors for SAFD-IDR pairs  
3. Maps vectors back to full folded domain using `map_SAFD_vector_to_full_folded_domain()`
4. Returns vectors with length equal to complete folded domain

## Integration with Core Systems

These utility functions integrate with the main FINCHES calculation pipeline:

```mermaid
graph TB
    subgraph "Frontend Layer"
        MPIPI["Mpipi_frontend"]
        CALV["CALVADOS_frontend"]
    end
    
    subgraph "Calculation Engine"
        IMC["InteractionMatrixConstructor"]
        EPSILON_CALC["epsilon_calculation module"]
    end
    
    subgraph "Utility Functions"
        CHARGE_WEIGHT["get_charge_weighted_mask()"]
        ALI_WEIGHT["get_aliphatic_weighted_mask()"]
        SEQ_PROPS["calculate_NCPR/FCR()"]
        PLOTS["interaction_vector plotting"]
    end
    
    subgraph "Data Flow"
        MATRICES["Weighted Matrices"]
        VECTORS["Interaction Vectors"]
        PLOTS_OUT["Visualizations"]
    end
    
    MPIPI --> IMC
    CALV --> IMC
    IMC --> EPSILON_CALC
    EPSILON_CALC --> CHARGE_WEIGHT
    EPSILON_CALC --> ALI_WEIGHT
    EPSILON_CALC --> SEQ_PROPS
    
    CHARGE_WEIGHT --> MATRICES
    ALI_WEIGHT --> MATRICES
    MATRICES --> VECTORS
    VECTORS --> PLOTS
    PLOTS --> PLOTS_OUT
```

**Sources:** [finches/parsing_aminoacid_sequences.py:1-388](), [finches/sequence_tools.py:1-456](), [finches/interaction_vector.py:1-318]()

The utilities serve as the computational foundation for advanced weighting schemes and visualization capabilities that distinguish FINCHES from simpler mean-field approaches.

---