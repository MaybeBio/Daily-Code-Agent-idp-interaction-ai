# AlphaFlex Pipeline

<details>
<summary>Relevant source files</summary>

The following files were used as context for generating this wiki page:

- [AlphaFlex/Data_Inputs/knot_screening.json](AlphaFlex/Data_Inputs/knot_screening.json)
- [AlphaFlex/README.md](AlphaFlex/README.md)
- [AlphaFlex/config.py](AlphaFlex/config.py)

</details>



The **AlphaFlex** pipeline is a four-step automated workflow designed for the large-scale generation of disordered protein ensembles. It leverages AlphaFold2 (AF2) structures as templates and uses the IDPForge diffusion model to sample conformational ensembles of Intrinsically Disordered Regions (IDRs) while maintaining the structural integrity of folded domains.

The pipeline is designed to handle proteins ranging from fully disordered (Category 0) to complex multi-domain proteins with disordered tails, linkers, and loops (Categories 1-3) [AlphaFlex/README.md:51-58]().

### Pipeline Architecture

The workflow progresses from database labeling to final ensemble assembly. Each step is orchestrated by a dedicated script that consumes the outputs of the previous stage.

```mermaid
graph TD
    subgraph "Step 1: Case Labeling"
        S1["Step_1_case_label.py"] --> S1B["Step_1B_subset_label.py"]
    end

    subgraph "Step 2: Template Prep"
        S2["Step_2_mk_ldr_template.py"]
    end

    subgraph "Step 3: Sampling"
        S3["Step_3_sample_conformer.py"]
    end

    subgraph "Step 4: Assembly"
        S4["Step_4_ldr_stitch.py"]
    end

    DB[("AlphaFlex_database.json")] --> S1
    S1B -->|"Filtered ID List"| S2
    S2 -->|"NPZ Templates"| S3
    S3 -->|"Raw PDB Pools"| S4
    S4 -->|"Final Ensembles"| Output[/Protein Ensembles/]

    style S1 fill:none,stroke-dasharray: 5 5
    style S1B fill:none,stroke-dasharray: 5 5
    style S2 fill:none,stroke-dasharray: 5 5
    style S3 fill:none,stroke-dasharray: 5 5
    style S4 fill:none,stroke-dasharray: 5 5
```
**Sources:** [AlphaFlex/README.md:3-16](), [AlphaFlex/config.py:27-101]()

---

### Step 1: Case Labeling & Subset Filtering
This stage classifies proteins based on their IDR topology using the `AlphaFlex_database_Nov2025.json` [AlphaFlex/README.md:9-10](). 

- **Case Labeling**: `Step_1_case_label.py` augments the database with IDR types: **Tails** (N- or C-terminal), **Linkers** (between domains), and **Loops** (internal to a domain) [AlphaFlex/README.md:51-58]().
- **Subset Filtering**: `Step_1B_subset_label.py` generates specific lists of UniProt IDs for processing, allowing users to filter by total protein length, IDR count, or specific IDR lengths [AlphaFlex/config.py:34-51]().

For details on classification logic and schema, see [Step 1: Case Labeling & Subset Filtering](#5.1).

**Sources:** [AlphaFlex/config.py:27-51](), [AlphaFlex/README.md:49-80]()

---

### Step 2: Template Generation
`Step_2_mk_ldr_template.py` prepares the structural starting points for the diffusion model. It handles two primary template types:

1. **Static Templates**: Created by `mk_ldr_template.py` for Tails and Loops, where the folded environment is held rigid [AlphaFlex/README.md:86-87]().
2. **Flexible Templates**: Created by `mk_flex_template.py` for Linkers, where adjacent folded domains are randomly translated and rotated to sample the inter-domain space [AlphaFlex/README.md:88-89]().

For details on template `.npz` formats and size-capping logic, see [Step 2: Template Generation](#5.2).

**Sources:** [AlphaFlex/config.py:53-74](), [AlphaFlex/README.md:82-94]()

---

### Step 3: Conformer Sampling
`Step_3_sample_conformer.py` manages the high-throughput execution of the IDPForge model (`sample_ldr.py` or `sample_idp.py`).

- **Parallelization**: Supports `split_index` for distributing workloads across HPC nodes [AlphaFlex/config.py:76-96]().
- **Quality Gates**: Implements curvature gates (`SAMPLE_JUNCTION_KAPPA`) and knot screening using `knot_screening.json` to ensure generated disordered regions do not create unphysical topologies with the folded domains [AlphaFlex/config.py:87-89](), [AlphaFlex/Data_Inputs/knot_screening.json:1-21]().

For details on the sampling loop and state persistence, see [Step 3: Conformer Sampling](#5.3).

**Sources:** [AlphaFlex/config.py:76-96](), [AlphaFlex/README.md:96-115]()

---

### Step 4: Kinematic Stitching & Ensemble Assembly
The final stage, `Step_4_ldr_stitch.py`, assembles individual IDR conformers back into the full-length protein.

- **Kinematic Assembly**: Uses Monte Carlo stitching to join disordered segments to folded domains while minimizing clashes [AlphaFlex/config.py:120-122]().
- **Relaxation**: Performs AMBER energy minimization with harmonic restraints (`RELAX_STIFFNESS`) on the folded domains to ensure chemical validity at the junctions [AlphaFlex/config.py:109-114]().

For details on the stitching algorithm and structural validation, see [Step 4: Kinematic Stitching & Ensemble Assembly](#5.4).

**Sources:** [AlphaFlex/config.py:99-122](), [AlphaFlex/README.md:117-133]()

---

### Configuration and Inputs
The pipeline is controlled by a centralized configuration file and several reference datasets.

| Entity | Code Reference | Role |
| :--- | :--- | :--- |
| **Pipeline Config** | `AlphaFlex/config.py` | Global paths, batch sizes, and physics thresholds. |
| **Knot Screening** | `knot_screening.json` | Pre-calculated native topology for AF2 templates. |
| **Master Database** | `AlphaFlex_database_Nov2025.json` | UniProt IDR boundaries and PAE-based interactions. |
| **Weights** | `weights/mdl.ckpt` | The trained IDPForge diffusion model weights. |

For details on schema and configuration parameters, see [AlphaFlex Configuration & Data Inputs](#5.5).

**Sources:** [AlphaFlex/config.py:1-25](), [AlphaFlex/Data_Inputs/knot_screening.json:1-21](), [AlphaFlex/README.md:7-16]()

---