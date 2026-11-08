# JudgementDay Analysis Pipeline

<details>
<summary>Relevant source files</summary>

The following files were used as context for generating this wiki page:

- [analysis/analysis.py](analysis/analysis.py)
- [analysis/fig_4_s5_chimera.py](analysis/fig_4_s5_chimera.py)
- [analysis/get_other_method_preds.py](analysis/get_other_method_preds.py)
- [analysis/side_analysis.py](analysis/side_analysis.py)

</details>



## Purpose and Scope

The JudgementDay Analysis Pipeline is a comprehensive evaluation framework for assessing the performance of Disobind predictions against AlphaFold2, AlphaFold3, and other methods on out-of-distribution (OOD) test sets. This pipeline orchestrates the comparison of multiple prediction methods, generates specialized binary masks for disorder and amino acid analysis, calculates performance metrics across various categories, and produces visualization outputs.

For information about generating Disobind predictions for analysis, see [Generating OOD Predictions](#5.2). For details on processing AlphaFold predictions, see [Processing AlphaFold Predictions](#5.3). For information about specific evaluation metrics, see [Performance Metrics and Evaluation](#5.4).

## System Overview

The JudgementDay pipeline comprises four major components that work together to perform comprehensive model evaluation:

```mermaid
graph TB
    subgraph "Input Generation"
        PredClass["Prediction class<br/>get_disobind_predictions.py"]
        AF2MClass["AF2MPredictions class<br/>get_af_prediction.py"]
        OtherClass["Othermethods class<br/>get_other_method_preds.py"]
    end
    
    subgraph "Core Orchestration"
        JDClass["JudgementDay class<br/>analysis.py"]
    end
    
    subgraph "Prediction Sources"
        DisoNPY["Disobind_Predictions.npy<br/>6 tasks × CG levels"]
        AF2NPY["Predictions_af2m_results.npy<br/>Structure-based contacts"]
        AF3NPY["Predictions_af3_results.npy<br/>Structure-based contacts"]
        OtherNPY["other_methods.npy<br/>AIUPred/DeepDisoBind/MORFchibi"]
        TargetH5["Target_bcmap_test.h5<br/>Ground truth contacts"]
    end
    
    subgraph "Analysis Outputs"
        ResultsCSV["Results_OOD_set.csv<br/>Metrics by category"]
        ConfPlot["AF_confidence_plot.png<br/>ipTM comparison"]
        SparsityPlot["Sparsity_F1_plot.png<br/>Dataset characteristics"]
        CaseCSV["Case_sp_analysis.csv<br/>Sample-wise metrics"]
    end
    
    PredClass --> DisoNPY
    AF2MClass --> AF2NPY
    AF2MClass --> AF3NPY
    OtherClass --> OtherNPY
    
    DisoNPY --> JDClass
    AF2NPY --> JDClass
    AF3NPY --> JDClass
    OtherNPY --> JDClass
    TargetH5 --> JDClass
    
    JDClass --> ResultsCSV
    JDClass --> ConfPlot
    JDClass --> SparsityPlot
    JDClass --> CaseCSV
```

**Sources:** [analysis/analysis.py:1-87](), [analysis/get_disobind_predictions.py:1-50](), [analysis/get_af_prediction.py:1-77](), [analysis/get_other_method_preds.py:1-35]()

## Core Classes and Components

### JudgementDay Class

The `JudgementDay` class is the primary orchestrator for the analysis pipeline. It manages configuration, loads predictions from all sources, and coordinates evaluation across multiple categories.

```mermaid
graph LR
    subgraph "JudgementDay Configuration"
        Init["__init__<br/>lines 29-97"]
        Config["Configuration Fields"]
    end
    
    subgraph "Core Methods"
        Forward["forward()<br/>lines 107-123<br/>Main entry point"]
        GetTasks["get_tasks()<br/>lines 125-131<br/>Returns task list"]
        CreateTensors["create_ood_set_tensors()<br/>lines 192-316<br/>Generator function"]
        EvalPerf["eval_performance()<br/>lines 438-508<br/>Calculate metrics"]
    end
    
    subgraph "Analysis Methods"
        BaseModel["get_base_model_preds()<br/>lines 319-331"]
        Combine["combine_diso_af_preds()<br/>lines 334-351"]
        DisorderPreds["get_preds_for_disorder_order_residues()<br/>lines 354-382"]
        IntTypePreds["get_preds_for_interaction_types()<br/>lines 385-412"]
        OtherMethods["get_other_method_preds()<br/>lines 416-434"]
    end
    
    Init --> Config
    Forward --> CreateTensors
    CreateTensors --> EvalPerf
    EvalPerf --> BaseModel
    EvalPerf --> Combine
    EvalPerf --> DisorderPreds
    EvalPerf --> IntTypePreds
    EvalPerf --> OtherMethods
```

**Key Configuration Parameters:**

| Parameter | Default | Purpose |
|-----------|---------|---------|
| `data_version` | 19 | Dataset version identifier |
| `mode` | "ood" | Dataset type (ood/misc) |
| `max_len` | 100 (ood) / 200 (misc) | Maximum protein length |
| `contact_threshold` | 0.5 | Cutoff for binary contact prediction |
| `iptm_cutoff` | 0.0 | Minimum ipTM score for confident AF predictions |
| `device` | "cuda" | Computation device |

**Sources:** [analysis/analysis.py:29-97]()

### Othermethods Class

The `Othermethods` class is responsible for parsing and standardizing predictions from external interface predictors: AIUPred, DeepDisoBind, and MORFchibi. It ensures these predictions are aligned with the Disobind output format for direct comparison.

```mermaid
graph TB
    subgraph "Othermethods Workflow"
        OMInit["__init__<br/>lines 20-35"]
        OMForward["forward()<br/>lines 38-47"]
        
        AIU["get_aiupred_predictions()<br/>lines 65-97"]
        DeepDiso["get_deepdisobind_predictions()<br/>lines 100-122"]
        MORF["get_morfchibi_predictions()<br/>lines 184-210"]
        Assemble["assemble_interfaces_for_ood_entries()<br/>lines 213-264"]
    end
    
    subgraph "MORF Processing"
        MORFParse["parse_morfchibi_output()<br/>lines 125-140"]
        MORFValid["valid_morf()<br/>lines 143-157"]
        MORFProc["process_morfchibi_pred()<br/>lines 160-181"]
    end
    
    OMForward --> AIU
    OMForward --> DeepDiso
    OMForward --> MORF
    MORF --> MORFParse
    MORFParse --> MORFProc
    MORFProc --> MORFValid
    OMForward --> Assemble
```

**External Method Specifications:**

| Method | Implementation / Parsing Logic | Cutoff / Processing |
|-----------|--------------------------------|---------------------|
| **AIUPred** | Uses `aiupred_lib.predict_binding` [analysis/get_other_method_preds.py:86-89]() | Threshold: 0.5 [analysis/get_other_method_preds.py:91]() |
| **DeepDisoBind** | Parses FASTA result files [analysis/get_other_method_preds.py:110-118]() | Binary string extraction |
| **MORFchibi** | Suggests 0.775 cutoff [analysis/get_other_method_preds.py:168]() | Sections < 4 residues filtered [analysis/get_other_method_preds.py:172-178]() |

**Sources:** [analysis/get_other_method_preds.py:14-264]()

### AF2MPredictions Class

The `AF2MPredictions` class processes AlphaFold2 and AlphaFold3 predictions, extracting contact maps from structures with confidence filtering based on pLDDT and PAE scores.

**Key Processing Steps:**

1. **Model Selection:** Ranks AF models by ipTM+pTM scores, selects best model [analysis/get_af_prediction.py:96-132]()
2. **Structure Parsing:** Extracts C-alpha coordinates and pLDDT values [analysis/get_af_prediction.py:163-195]()
3. **Contact Map Creation:** Computes inter-residue distances, applies 8Å threshold [analysis/get_af_prediction.py:198-214]()
4. **Confidence Filtering:** Creates binary masks for pLDDT≥70 and PAE≤5 [analysis/get_af_prediction.py:315-361]()
5. **Coarse-Graining:** Applies MaxPool2d for CG levels 5 and 10 [analysis/get_af_prediction.py:267-312]()

**Sources:** [analysis/get_af_prediction.py:28-498]()

## Analysis Workflow

### Data Loading and Preparation

The pipeline begins by loading predictions from all sources and preparing them for evaluation:

```mermaid
graph LR
    subgraph "Load Predictions"
        LoadAF2["Load AF2 predictions<br/>af2m_preds = np.load()"]
        LoadAF3["Load AF3 predictions<br/>af3_preds = np.load()"]
        LoadDiso["Load Disobind predictions<br/>disobind_preds = np.load()"]
        LoadOther["Load other methods<br/>other_methods = np.load()"]
        LoadTarget["Load ground truth<br/>target_cmap = h5py.File()"]
    end
    
    subgraph "Tensor Creation"
        CreateOOD["create_ood_set_tensors()<br/>Generator yields task data"]
    end
    
    subgraph "Task Loop"
        Tasks["For each task:<br/>interaction_1/5/10<br/>interface_1/5/10"]
    end
    
    LoadAF2 --> CreateOOD
    LoadAF3 --> CreateOOD
    LoadDiso --> CreateOOD
    LoadOther --> CreateOOD
    LoadTarget --> CreateOOD
    
    CreateOOD --> Tasks
```

The `create_ood_set_tensors()` method is a generator that yields data for each task, assembling all predictions, masks, and targets into a unified dictionary structure.

**Sources:** [analysis/analysis.py:107-123](), [analysis/analysis.py:192-316]()

### Side Analysis Utilities

The `side_analysis.py` module provides utility functions for segregating the dataset into specific biological categories for targeted evaluation.

| Function | Logic | File Reference |
|----------|-------|----------------|
| `get_dor_ddr_complexes` | Segregates Dynamic (DDR) vs. Static (DOR) complexes based on contact map variance across conformers. | [analysis/side_analysis.py:55-107]() |
| `get_frac_disordered` | Calculates the fraction of residues in a sequence fragment that are annotated as disordered. | [analysis/side_analysis.py:125-138]() |
| `get_full_idr_complexes` | Identifies complexes where Protein 1 is 100% disordered. | [analysis/side_analysis.py:141-183]() |

**Sources:** [analysis/side_analysis.py:55-183]()

### Evaluation Categories

The pipeline evaluates predictions across five main categories:

1. **Base Model Evaluation** [analysis/analysis.py:319-331]()
2. **Combined Predictions** [analysis/analysis.py:334-351]()
   - Formula: `combined = max(disobind, af_pred)`
3. **Disorder-Specific Analysis** [analysis/analysis.py:354-382]()
   - Evaluates performance on IDR-IDR vs ordered interactions.
4. **Interaction Type Analysis** [analysis/analysis.py:385-412]()
   - Applies amino acid property masks (aromatic, hydrophobic, etc.).
5. **Other Methods Comparison** [analysis/analysis.py:416-434]()
   - Comparison with AIUPred, DeepDisoBind, and MORFchibi.

**Sources:** [analysis/analysis.py:438-508]()

## Visualization and Visualization Utilities

### Structural Visualization (ChimeraX)

The `fig_4_s5_chimera.py` script provides automation for ChimeraX to visualize predicted vs. target interfaces, specifically for NMR ensembles.

```python
def ensemble_in_background(pdb_file, chain1, chain2, predicted_interface1, predicted_interface2, out_path, name):
    # Opens mmCIF, sets preset 1, and handles transparency for NMR models
    # Colors predicted interface residues (Blue for Chain 1, Red for Chain 2)
```

**Key Visualization Features:**
- **Transparency:** Sets 95% transparency for non-representative NMR models to show structural ensembles [analysis/fig_4_s5_chimera.py:37-38]().
- **Interface Mapping:** Maps binary interface predictions back to PDB residue positions for 3D visualization [analysis/fig_4_s5_chimera.py:69-80]().

**Sources:** [analysis/fig_4_s5_chimera.py:9-130]()

### Sparsity vs F1 Plot

Visualizes the relationship between dataset sparsity and Disobind F1-score performance. The pipeline loads contact density from `fraction_positives.json` and correlates it with performance metrics.

**Sources:** [analysis/analysis.py:670-705]()

## Configuration and Customization

### Mode Selection

The pipeline supports two operational modes:

**OOD Mode (Out-of-Distribution):**
- Dataset: `v_{version}/Target_bcmap_test_v_{version}.h5`
- Max length: 100 residues
- Includes comparisons with external methods.

**Misc Mode (Miscellaneous Dataset):**
- Dataset: `Misc/misc_test_target.h5`
- Max length: 200 residues
- Focuses on sample-wise analysis and top prediction extraction.

**Sources:** [analysis/analysis.py:42-67]()

### Exclusion of Redundant Entries

The pipeline automatically excludes entries with sequence redundancy to PDB70 at 20% identity to ensure the integrity of the OOD test set.

```python
# Excluded entries (sequence redundant with PDB70)
excluded = ["P0DTC9--P0DTD1_2", "Q96PU5--Q96PU5_0", "P0AG11--P0AG11_4",
            "Q9IK92--Q9IK91_0", "Q16236--O15525_0", "P12023--P12023_0",
            "O85041--O85043_0", "P25024--P10145_0"]
```

**Sources:** [analysis/analysis.py:223-239](), [analysis/side_analysis.py:75-79]()

---