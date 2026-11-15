# Overview

<details>
<summary>Relevant source files</summary>

The following files were used as context for generating this wiki page:

- [README.md](README.md)
- [changelog.md](changelog.md)
- [docs/index.rst](docs/index.rst)
- [finches/__init__.py](finches/__init__.py)

</details>



## Purpose and Scope

This document provides a high-level introduction to the FINCHES software system, its core architecture, and primary capabilities. FINCHES (**F**irst-principle **I**nteractions via **CHE**mical **S**pecificity) is a Python package for predicting intermolecular interactions driven by intrinsically disordered regions (IDRs) in proteins using chemical physics principles extracted from coarse-grained force fields.

For detailed installation and setup instructions, see [Getting Started](#2). For comprehensive usage examples, see [User Guide](#3). For complete API documentation, see [API Reference](#5).

## What is FINCHES

FINCHES enables computational prediction of IDR-mediated interactions through sequence-based analysis. The system transforms amino acid sequences into quantitative interaction predictions using established forcefield parameters, supporting multiple analysis modes including pairwise interaction maps, mean-field interaction parameters (epsilon values), and phase behavior predictions.

**Sources:** [README.md:7-8](), [docs/index.rst:9-23]()

## System Architecture

The FINCHES system follows a layered architecture with distinct separation between user interfaces, computational engines, and underlying forcefield models:

```mermaid
graph TB
    subgraph "User Interface Layer"
        MF[Mpipi_frontend]
        CF[CALVADOS_frontend] 
        FB[FinchesFrontend]
        FB --> MF
        FB --> CF
    end
    
    subgraph "Computational Engine"
        IMC[InteractionMatrixConstructor]
        ES[epsilon_stateless]
        MM[matrix_manipulation.pyx]
        IMC --> ES
        IMC --> MM
    end
    
    subgraph "Forcefield Models"
        MP[mpipi_model]
        CAL[calvados_model]
        CUST[custom_model]
    end
    
    subgraph "Analysis Tools"
        FD[FoldedDomain]
        PST[PDB_structure_tools]
        SP[parsing_aminoacid_sequences]
        ST[sequence_tools]
    end
    
    MF --> IMC
    CF --> IMC
    IMC --> MP
    IMC --> CAL
    IMC --> CUST
    IMC --> FD
    IMC --> SP
    
    FD --> PST
    SP --> ST
```

**Sources:** [finches/__init__.py:16-17](), [README.md:25-34]()

## User Interface Layer

FINCHES provides two primary frontend classes that implement identical interfaces for different forcefield models:

```mermaid
graph LR
    subgraph "Frontend Classes"
        MF["Mpipi_frontend"]
        CF["CALVADOS_frontend"]
        FB["FinchesFrontend (base)"]
        
        FB --> MF
        FB --> CF
    end
    
    subgraph "Core Methods"
        EP["epsilon()"]
        IF["interaction_figure()"]
        IM["interaction_map()"]
        PD["phase_diagram()"]
        RV["per_residue_attractive_vector()"]
    end
    
    subgraph "Analysis Types"
        IDR_IDR["IDR-IDR interactions"]
        IDR_FD["IDR-folded domain"]
        PHASE["Phase behavior"]
        VECTORS["Interaction vectors"]
    end
    
    MF --> EP
    CF --> EP
    EP --> IDR_IDR
    
    MF --> IF
    CF --> IF
    IF --> IDR_IDR
    IF --> IDR_FD
    
    MF --> IM
    CF --> IM
    IM --> IDR_IDR
    IM --> IDR_FD
    
    MF --> PD
    CF --> PD
    PD --> PHASE
    
    MF --> RV
    CF --> RV
    RV --> VECTORS
```

**Sources:** [finches/__init__.py:16-17](), [README.md:36-43]()

## Core Capabilities

FINCHES provides four primary analysis capabilities, each accessible through the frontend interfaces:

| Capability | Output | Use Case |
|------------|--------|----------|
| **Epsilon Calculation** | Scalar interaction parameter | Quantifying overall interaction strength between sequence pairs |
| **Interaction Maps** | 2D heatmap matrices | Visualizing residue-specific interaction patterns |
| **IDR-Folded Domain Analysis** | Surface interaction profiles | Analyzing IDR interactions with structured protein domains |
| **Phase Diagrams** | Concentration-temperature plots | Predicting liquid-liquid phase separation behavior |

**Sources:** [docs/index.rst:13-22]()

## Computational Models

FINCHES implements multiple forcefield models through a unified interface, enabling comparative analysis using different physical frameworks:

```mermaid
graph TB
    subgraph "Forcefield Interface"
        FFI["compute_interaction_parameter()"]
    end
    
    subgraph "Mpipi Implementation"
        MPIPI["mpipi_model"]
        WF["Wang-Frenkel potential"]
        COUL["Coulombic interactions"]
        MPIPI --> WF
        MPIPI --> COUL
    end
    
    subgraph "CALVADOS Implementation" 
        CALV["calvados_model"]
        AH["Ashbaugh-Hatch potential"]
        YUK["Yukawa electrostatics"]
        CALV --> AH
        CALV --> YUK
    end
    
    subgraph "Custom Implementation"
        CUSTOM["custom_model"]
        USER["User-defined parameters"]
        CUSTOM --> USER
    end
    
    FFI --> MPIPI
    FFI --> CALV  
    FFI --> CUSTOM
    
    subgraph "Parameter Data"
        MPIPI_DATA["Mpipi parameters"]
        CALV_DATA["CALVADOS parameters"]
        CUSTOM_DATA["Custom interaction tables"]
    end
    
    MPIPI --> MPIPI_DATA
    CALV --> CALV_DATA
    CUSTOM --> CUSTOM_DATA
```

**Sources:** [README.md:62-71](), [changelog.md:26-30]()

## Data Processing Pipeline

The system processes biological sequence and structural data through a standardized pipeline that transforms input data into quantitative interaction predictions:

```mermaid
flowchart TD
    subgraph "Input Data"
        SEQ1["Protein sequence 1"]
        SEQ2["Protein sequence 2"]
        PDB["PDB structure files"]
        PARAMS["Model parameters"]
    end
    
    subgraph "Processing Functions"
        PARSE["parsing_aminoacid_sequences"]
        STRUCT["PDB_structure_tools"]
        VALID["sequence validation"]
    end
    
    subgraph "Core Calculations"
        MATRIX["Pairwise matrix generation"]
        EPSILON["epsilon_stateless functions"]
        SLIDING["Sliding window analysis"]
        BASELINE["Null baseline correction"]
    end
    
    subgraph "Output Generation"
        INTERMAP["2D interaction maps"]
        EPS_VAL["Epsilon scalar values"]
        PHASE_PRED["Phase diagram predictions"]
        VECTORS["Per-residue interaction vectors"]
    end
    
    SEQ1 --> PARSE
    SEQ2 --> PARSE
    PDB --> STRUCT
    PARAMS --> VALID
    
    PARSE --> MATRIX
    STRUCT --> MATRIX
    VALID --> MATRIX
    
    MATRIX --> EPSILON
    MATRIX --> SLIDING
    EPSILON --> BASELINE
    
    SLIDING --> INTERMAP
    BASELINE --> EPS_VAL
    EPS_VAL --> PHASE_PRED
    EPSILON --> VECTORS
```

**Sources:** [docs/index.rst:14-21](), [README.md:36-43]()

## Performance Architecture

FINCHES employs strategic performance optimizations through Cython compilation for computationally intensive matrix operations:

```mermaid
graph LR
    subgraph "Python Layer"
        API["Frontend APIs"]
        WRAP["Method wrappers"]
        PY_CALC["Python calculations"]
    end
    
    subgraph "Optimized Core"
        CY_OPT["matrix_manipulation.pyx"]
        LOOKUP["Parameter lookups"]
    end
    
    subgraph "External Dependencies"
        NP["NumPy arrays"]
        SP["SciPy functions"]
        MD["MDTraj structures"]
        MP["MetaPredict disorder"]
    end
    
    API --> WRAP
    WRAP --> PY_CALC
    WRAP --> CY_OPT
    
    PY_CALC --> LOOKUP
    CY_OPT --> LOOKUP
    
    PY_CALC --> NP
    CY_OPT --> NP
    PY_CALC --> SP
    
    WRAP --> MD
    WRAP --> MP
```

**Sources:** [README.md:94-108](), [changelog.md:145]()

## Current Status and Availability

FINCHES is available in multiple deployment formats to accommodate different user preferences and computational environments:

| Deployment | Access Method | Target Users |
|------------|---------------|--------------|
| **Python Package** | `pip install` from GitHub | Developers, power users |
| **Google Colab** | Pre-configured notebooks | Educational, quick analysis |
| **Web Interface** | https://finches-online.com | General users, demonstrations |

The software is currently in public beta status, with active development focused on performance improvements and feature additions.

**Sources:** [README.md:12-17](), [README.md:16-17](), [changelog.md:4-6]()

---