# Advanced Analysis Examples

<details>
<summary>Relevant source files</summary>

The following files were used as context for generating this wiki page:

- [demo/nucleic_acid/protein_nucleic_acid_binding.ipynb](demo/nucleic_acid/protein_nucleic_acid_binding.ipynb)
- [demo/phase_diagrams/phase_diagram_demo.ipynb](demo/phase_diagrams/phase_diagram_demo.ipynb)
- [demo/protein_matrix/interaction_matrix_demo.ipynb](demo/protein_matrix/interaction_matrix_demo.ipynb)
- [finches/epsilon_to_FHtheory.py](finches/epsilon_to_FHtheory.py)

</details>



This page demonstrates complex analytical workflows using FINCHES, including phase diagram generation with Flory-Huggins theory, protein-nucleic acid interaction analysis, and multi-condition parameter sweeps. For basic usage patterns such as simple epsilon calculations and interaction maps, see [Basic Usage Examples](#6.1).

## Phase Diagram Analysis with Flory-Huggins Theory

FINCHES integrates analytical solutions to the Flory-Huggins model for generating phase diagrams from computed epsilon values. This functionality bridges molecular-level interaction parameters with thermodynamic phase behavior.

### Core Phase Diagram Workflow

```mermaid
flowchart TD
    SEQ["Protein Sequence"] --> EPSILON["epsilon_stateless.get_sequence_epsilon_value()"]
    EPSILON --> CONVERT["epsilon_to_FHtheory.epsilon_to_phase_diagram()"]
    CONVERT --> BINODAL["floryhuggins.calculate_binodal()"]
    CONVERT --> SPINODAL["floryhuggins.calculate_spinodal()"]
    
    BINODAL --> PHASE_DATA["Phase Separation Boundaries"]
    SPINODAL --> PHASE_DATA
    
    PHASE_DATA --> PLOT["Phase Diagram Visualization"]
    
    subgraph "Conversion Logic"
        EPSILON_CALC["Negative epsilon = attractive"]
        LENGTH_NORM["epsilon/sequence_length"]
        CHI_CONV["chi = delta_eps / (kB * T)"]
    end
    
    CONVERT --> EPSILON_CALC
    EPSILON_CALC --> LENGTH_NORM
    LENGTH_NORM --> CHI_CONV
```

The phase diagram calculation involves several key transformations. When epsilon values are positive (repulsive), they are set to a small attractive value (-0.01) to prevent mathematical issues in the Flory-Huggins calculations. The epsilon value is normalized by sequence length since it represents chain-level interactions while the Flory-Huggins model requires per-residue energies.

**Sources:** [finches/epsilon_to_FHtheory.py:356-448]()

### Multi-Condition Phase Analysis

FINCHES provides specialized functions for analyzing phase behavior across different environmental conditions:

```mermaid
graph TB
    INPUT["Sequence + Frontend"] --> SALT_SWEEP["build_SALT_dependent_phase_diagrams()"]
    INPUT --> PH_SWEEP["build_PH_dependent_phase_diagrams()"]
    INPUT --> DIELEC_SWEEP["build_DIELECTRIC_dependent_phase_diagrams()"]
    
    SALT_SWEEP --> SALT_LOOP["Loop: condition_list values"]
    PH_SWEEP --> PH_LOOP["Loop: condition_list values"] 
    DIELEC_SWEEP --> DIELEC_LOOP["Loop: condition_list values"]
    
    SALT_LOOP --> UPDATE_PARAMS["X_class.parameters.salt = value"]
    PH_LOOP --> UPDATE_PARAMS2["X_class.parameters.pH = value"]
    DIELEC_LOOP --> UPDATE_PARAMS3["X_class.parameters.dielectric = value"]
    
    UPDATE_PARAMS --> LOOKUP_UPDATE["X_class._update_lookup_dict()"]
    UPDATE_PARAMS2 --> LOOKUP_UPDATE
    UPDATE_PARAMS3 --> LOOKUP_UPDATE
    
    LOOKUP_UPDATE --> PHASE_CALC["return_phase_diagram()"]
    PHASE_CALC --> RESULTS["Condition-indexed results"]
```

Each condition sweep function follows the same pattern: iterating through condition values, updating the forcefield parameters, refreshing the lookup tables, and computing phase diagrams for each condition.

**Sources:** [finches/epsilon_to_FHtheory.py:33-271]()

## Protein-Nucleic Acid Interaction Analysis

FINCHES supports heterotypic interactions between proteins and nucleic acids, extending beyond traditional protein-protein interactions. This capability is particularly relevant for analyzing RNA-binding proteins and DNA-binding domains.

### Nucleic Acid Parameter Integration

```mermaid
flowchart LR
    PROTEIN_SEQ["Protein Sequence"] --> MPIPI["Mpipi_frontend"]
    RNA_SEQ["RNA/DNA Sequence"] --> CUSTOM_PARAMS["Custom Nucleic Acid Parameters"]
    
    MPIPI --> MATRIX["InteractionMatrixConstructor"]
    CUSTOM_PARAMS --> MATRIX
    
    MATRIX --> EPSILON["Cross-species epsilon calculation"]
    MATRIX --> INTERMAP["Heterotypic interaction maps"]
    
    subgraph "Parameter Sources"
        MPIPI_PARAMS["Mpipi amino acid parameters"]
        NA_PARAMS["Nucleic acid interaction parameters"]
    end
    
    MPIPI_PARAMS --> MATRIX
    NA_PARAMS --> MATRIX
```

The analysis involves extending the standard amino acid parameter sets to include nucleic acid residues, enabling calculation of protein-nucleic acid interaction energies using the same computational framework.

**Sources:** [demo/nucleic_acid/protein_nucleic_acid_binding.ipynb:1-200]()

### Heterotypic Interaction Matrices

When analyzing protein-nucleic acid systems, interaction matrices reveal binding specificity patterns:

```mermaid
graph TD
    PROT_RES["Protein Residues (rows)"] --> MATRIX["Interaction Matrix"]
    NA_RES["Nucleic Acid Residues (cols)"] --> MATRIX
    
    MATRIX --> HOTSPOTS["Binding Hotspots"]
    MATRIX --> SPECIFICITY["Sequence Specificity"]
    MATRIX --> AVIDITY["Multivalent Interactions"]
    
    subgraph "Analysis Types"
        SPATIAL["Spatial interaction patterns"]
        ENERGETIC["Energy landscape mapping"]
        COMPETITIVE["Competitive binding analysis"]
    end
    
    HOTSPOTS --> SPATIAL
    SPECIFICITY --> ENERGETIC
    AVIDITY --> COMPETITIVE
```

**Sources:** [demo/nucleic_acid/protein_nucleic_acid_binding.ipynb:1-200]()

## Complex Interaction Matrix Analysis

Advanced interaction matrix analysis goes beyond simple pairwise epsilon calculations to examine spatially-resolved interaction patterns and multi-domain systems.

### IDR-IDR Interaction Networks

```mermaid
flowchart TB
    IDR1["med1_idr sequence"] --> FRONTEND["Mpipi_frontend/CALVADOS_frontend"]
    IDR2["med14_idr sequence"] --> FRONTEND
    IDR3["med15_idr sequence"] --> FRONTEND
    
    FRONTEND --> MATRIX_CALC["generate_interaction_matrix()"]
    
    MATRIX_CALC --> HOMO1["med1-med1 homotypic"]
    MATRIX_CALC --> HOMO2["med14-med14 homotypic"] 
    MATRIX_CALC --> HOMO3["med15-med15 homotypic"]
    
    MATRIX_CALC --> HETERO12["med1-med14 heterotypic"]
    MATRIX_CALC --> HETERO13["med1-med15 heterotypic"]
    MATRIX_CALC --> HETERO23["med14-med15 heterotypic"]
    
    subgraph "Analysis Outputs"
        HEATMAPS["Interaction heatmaps"]
        VECTORS["Per-residue interaction vectors"]
        WINDOWS["Sliding window analysis"]
    end
    
    HOMO1 --> HEATMAPS
    HETERO12 --> HEATMAPS
    HOMO1 --> VECTORS
    HETERO12 --> VECTORS
    HOMO1 --> WINDOWS
    HETERO12 --> WINDOWS
```

The interaction matrix analysis reveals position-specific interaction preferences across different IDR pairs, enabling identification of binding motifs and cooperative interaction regions.

**Sources:** [demo/protein_matrix/interaction_matrix_demo.ipynb:80-150]()

### Multi-Scale Analysis Pipeline

```mermaid
graph TB
    SEQUENCE["Input Sequences"] --> GLOBAL["Global epsilon calculation"]
    SEQUENCE --> LOCAL["Local interaction matrices"]
    SEQUENCE --> DOMAINS["Domain-specific analysis"]
    
    GLOBAL --> PHASE["Phase behavior prediction"]
    LOCAL --> BINDING["Binding site identification"]
    DOMAINS --> MODULARITY["Modular interaction analysis"]
    
    subgraph "Resolution Levels"
        CHAIN["Chain-level (epsilon)"]
        RESIDUE["Residue-level (matrices)"]
        SEGMENT["Segment-level (domains)"]
    end
    
    PHASE --> CHAIN
    BINDING --> RESIDUE
    MODULARITY --> SEGMENT
    
    subgraph "Integration"
        COMBINED["Multi-scale interaction model"]
    end
    
    CHAIN --> COMBINED
    RESIDUE --> COMBINED
    SEGMENT --> COMBINED
```

**Sources:** [demo/protein_matrix/interaction_matrix_demo.ipynb:1-200](), [finches/epsilon_to_FHtheory.py:275-351]()

## Custom Forcefield Integration

Advanced users can integrate custom interaction parameters for specialized systems or novel amino acid modifications.

### Custom Parameter Workflow

```mermaid
flowchart LR
    CUSTOM_DICT["Custom interaction dictionary"] --> CUSTOM_MODEL["custom_model.CustomModel"]
    LITERATURE["Literature parameters"] --> CUSTOM_DICT
    EXPERIMENT["Experimental data"] --> CUSTOM_DICT
    
    CUSTOM_MODEL --> CONSTRUCTOR["InteractionMatrixConstructor"]
    STANDARD_FF["Mpipi/CALVADOS"] --> CONSTRUCTOR
    
    CONSTRUCTOR --> HYBRID["Hybrid parameter calculations"]
    HYBRID --> RESULTS["Analysis with custom parameters"]
    
    subgraph "Parameter Types"
        SIGMA["Size parameters (sigma)"]
        EPSILON_PARAM["Energy parameters (epsilon)"]
        CHARGE["Charge modifications"]
        CUSTOM_PROPS["Custom properties"]
    end
    
    CUSTOM_DICT --> SIGMA
    CUSTOM_DICT --> EPSILON_PARAM
    CUSTOM_DICT --> CHARGE
    CUSTOM_DICT --> CUSTOM_PROPS
```

Custom forcefields enable analysis of non-standard systems including modified amino acids, synthetic polymers, or hybrid protein-polymer systems.

**Sources:** [finches/forcefields/custom_model.py:1-100]()

## Integration with External Tools

FINCHES analysis workflows can be integrated with structure prediction and experimental validation tools.

### Analysis Integration Pipeline

```mermaid
graph TB
    METAPREDICT["metapredict disorder prediction"] --> SEQUENCES["IDR identification"]
    SPARROW["sparrow structure analysis"] --> SEQUENCES
    
    SEQUENCES --> FINCHES["FINCHES interaction analysis"]
    
    FINCHES --> PHASE_PRED["Phase separation prediction"]
    FINCHES --> BINDING_PRED["Binding specificity prediction"]
    
    PHASE_PRED --> VALIDATION["Experimental validation"]
    BINDING_PRED --> VALIDATION
    
    subgraph "External Integration"
        STRUCT_PRED["Structure prediction tools"]
        EXP_DATA["Experimental datasets"]
        ML_MODELS["Machine learning models"]
    end
    
    VALIDATION --> STRUCT_PRED
    VALIDATION --> EXP_DATA
    VALIDATION --> ML_MODELS
```

**Sources:** [demo/protein_matrix/interaction_matrix_demo.ipynb:55-65](), [demo/nucleic_acid/protein_nucleic_acid_binding.ipynb:55-65]()

---