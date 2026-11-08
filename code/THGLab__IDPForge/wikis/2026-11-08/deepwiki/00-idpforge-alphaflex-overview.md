# IDPForge & AlphaFlex Overview

<details>
<summary>Relevant source files</summary>

The following files were used as context for generating this wiki page:

- [AlphaFlex/Data_Inputs/knot_screening.json](AlphaFlex/Data_Inputs/knot_screening.json)
- [AlphaFlex/README.md](AlphaFlex/README.md)
- [LICENSE](LICENSE)
- [README.md](README.md)
- [configs/sample.yml](configs/sample.yml)
- [data/sic1_pre_exp.txt](data/sic1_pre_exp.txt)
- [environment.yml](environment.yml)
- [esm/__init__.py](esm/__init__.py)
- [setup.py](setup.py)

</details>



IDPForge is a generative transformer-based protein language diffusion model designed to create all-atom ensembles for Intrinsically Disordered Proteins (IDPs) and proteins containing Intrinsically Disordered Regions (IDRs). While traditional structural biology tools often struggle with the conformational heterogeneity of disordered proteins, IDPForge leverages ESM2-based embeddings and diffusion processes to sample valid structural ensembles that can be optionally guided by experimental data.

AlphaFlex (AFX-IDPForge) is an extension pipeline built upon IDPForge to handle large-scale ensemble generation for proteins with complex architectures, such as those containing multiple folded domains connected by disordered linkers.

### System Relation & Purpose

| System | Primary Purpose | Key Entities |
| :--- | :--- | :--- |
| **IDPForge** | Core generative model and sampling engine. | `idpforge/model.py`, `idpforge/diffuse.py` |
| **AlphaFlex** | Orchestration pipeline for multi-domain/IDR proteins. | `Step_1_case_label.py` through `Step_4_ldr_stitch.py` |

Sources: [README.md:1-4](), [AlphaFlex/README.md:1-5](), [setup.py:5-8]()

---

## Core Capabilities

IDPForge and AlphaFlex provide a comprehensive toolkit for IDP/IDR modeling:

*   **Generative Diffusion:** Uses a diffusion framework to sample protein conformations in both Euclidean (backbone coordinates) and torsional (dihedral angles) space.
*   **Experimental Guidance:** Incorporates gradient-based potentials during reverse diffusion to bias ensembles toward experimental observables such as Paramagnetic Relaxation Enhancement (PRE) or smFRET.
*   **Folded Context Retention:** Specifically designed to model disordered regions while maintaining the structural integrity of pre-existing folded domains.
*   **Topological Validation:** Includes a structural validation suite to ensure generated ensembles are free of knots (unless native) and physical violations.
*   **Large-Scale Orchestration:** The AlphaFlex pipeline automates the classification, template preparation, sampling, and kinematic stitching of complex proteins.

Sources: [README.md:1-4](), [configs/sample.yml:34-40](), [AlphaFlex/README.md:49-115]()

---

## High-Level Architecture

The system transitions from protein sequence and secondary structure definitions into a 3D coordinate space through a series of neural modules and diffusion steps.

### Conceptual Workflow: Natural Language to Code Entity Space

The following diagram maps the high-level scientific concepts to the specific code entities that implement them.

**Diagram: IDPForge Generative Flow**
```mermaid
graph TD
    subgraph "Input Space"
        SEQ["Sequence & SS Encoding"]
        EXP["Experimental Potentials"]
    end

    subgraph "Generative Core (idpforge/)"
        TRUNK["FoldingTrunk (model.py)"]
        DIFF["Diffuser (diffuse.py)"]
        POT["Potential (utils/potential.py)"]
    end

    subgraph "Output Space"
        ENS["Protein Ensemble (.pdb)"]
        VAL["Validation (utils/structure_validation.py)"]
    end

    SEQ --> TRUNK
    TRUNK --> DIFF
    EXP --> POT
    POT -- "Gradient Guidance" --> DIFF
    DIFF --> ENS
    ENS --> VAL
```
Sources: [idpforge/model.py:1-100](), [idpforge/diffuse.py:1-50](), [idpforge/utils/potential.py:1-50](), [idpforge/utils/structure_validation.py:1-50]()

### AlphaFlex Pipeline Architecture

The AlphaFlex pipeline manages the lifecycle of generating ensembles for proteins with folded domains (templates) and disordered regions (IDRs).

**Diagram: AlphaFlex Multi-Step Pipeline**
```mermaid
graph LR
    subgraph "Data & Labeling"
        S1["Step_1_case_label.py"]
        DB["AlphaFlex_database_Nov2025.json"]
    end

    subgraph "Template Preparation"
        S2["Step_2_mk_ldr_template.py"]
        MK_LDR["mk_ldr_template.py"]
        MK_FLEX["mk_flex_template.py"]
    end

    subgraph "Sampling & Assembly"
        S3["Step_3_sample_conformer.py"]
        S4["Step_4_ldr_stitch.py"]
        STITCH["assemble_kinematic_chain"]
    end

    DB --> S1
    S1 --> S2
    S2 --> MK_LDR
    S2 --> MK_FLEX
    MK_LDR --> S3
    MK_FLEX --> S3
    S3 --> S4
    S4 --> STITCH
```
Sources: [AlphaFlex/README.md:49-115](), [AlphaFlex/config.py:1-50]()

---

## Detailed Documentation Sections

For in-depth technical information, please refer to the following child pages:

### [1.1 Getting Started: Installation & Configuration](#)
Covers the environment setup using `environment.yml`, the integration of `OpenFold` resources, and the primary configuration files (`sample.yml`, `train.yml`) used to control model hyperparameters and diffusion schedules.
*   **Key Files:** `environment.yml`, `setup.py`, `configs/sample.yml`
*   **Details:** [README.md:15-97]()

### [1.2 Data Files & Example Data](#)
Documents the required data structures for training and inference, including secondary structure encodings, IGSO3 diffusion caches, and experimental data formats (e.g., `sic1_pre_exp.txt`).
*   **Key Files:** `data/example_data.pkl`, `data/sic1_pre_exp.txt`, `AlphaFlex/Data_Inputs/knot_screening.json`
*   **Details:** [README.md:200-205](), [AlphaFlex/README.md:7-16]()

---
*Note: This is a high-level parent page. For implementation details, loss functions, and specific model layers, see the [IDPForge Core Model Architecture](#) section.*

---