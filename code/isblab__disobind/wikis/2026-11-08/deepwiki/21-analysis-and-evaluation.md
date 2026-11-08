# Analysis and Evaluation

<details>
<summary>Relevant source files</summary>

The following files were used as context for generating this wiki page:

- [analysis/README.md](analysis/README.md)
- [analysis/analysis.py](analysis/analysis.py)
- [analysis/get_other_method_preds.py](analysis/get_other_method_preds.py)
- [dataset/README.md](dataset/README.md)

</details>



## Purpose and Scope

This document describes the comprehensive analysis and evaluation system for Disobind predictions. The analysis pipeline evaluates Disobind's performance on out-of-distribution (OOD) test sets and compares it against AlphaFold2/3 predictions and competing methods (AIUPred, DeepDisoBind, MORFchibi). The system generates predictions with specialized masks for disorder regions, amino acid types, and Linear Interacting Peptides (LIPs), then calculates performance metrics across multiple evaluation categories.

For information about running Disobind predictions, see [Running Predictions](#2). For details on model training, see [Model Architecture](#4). The subsections cover specific components: [JudgementDay Analysis Pipeline](#5.1), [Generating OOD Predictions](#5.2), [Processing AlphaFold Predictions](#5.3), [Performance Metrics and Evaluation](#5.4), [Disorder-Specific Analysis](#5.5), [Comparing with Other Methods](#5.6), and [IDPPI Dataset Evaluation](#5.7).

## System Overview

The analysis system operates in three stages:

1. **Prediction Generation**: Generate Disobind predictions with disorder and amino acid masks for the OOD test set using `get_disobind_predictions.py`.
2. **AlphaFold Processing**: Extract contact maps from AF2/3 structures with confidence filtering using `get_af_prediction.py`.
3. **Comparative Evaluation**: Calculate metrics across multiple dimensions and generate analysis results using `analysis.py`.

All scripts are located in the `analysis/` directory and process predictions from the dataset version 19 by default [analysis/analysis.py:31-32]().

**Diagram: Analysis Pipeline Architecture**

```mermaid
graph TB
    subgraph Input["Input Data Sources"]
        TestSet["OOD Test Set<br/>../database/v_19/"]
        DisorderDBs["Disorder Databases<br/>DisProt/IDEAL/MobiDB"]
        AF2Files["AlphaFold2 Structures<br/>AF2_preds/"]
        AF3Files["AlphaFold3 Structures<br/>AF3_preds/"]
        OtherInputs["Other Method Inputs<br/>../database/other_methods/"]
    end
    
    subgraph Generation["Prediction Generation"]
        PredClass["Prediction class<br/>get_disobind_predictions.py"]
        LoadModels["load_model()<br/>params.py"]
        CreateMasks["create_lip_masks()<br/>create_aa_masks()"]
        GetDisorder["get_disorder_matrix()"]
        PredForward["predict()<br/>Run all 6 tasks"]
        
        PredClass --> LoadModels
        PredClass --> CreateMasks
        PredClass --> GetDisorder
        LoadModels --> PredForward
        CreateMasks --> PredForward
        GetDisorder --> PredForward
    end
    
    subgraph AFProcessing["AlphaFold Processing"]
        AF2MClass["AF2MPredictions class<br/>get_af_prediction.py"]
        GetBestModel["get_best_model()<br/>ranking_debug.json"]
        GetCoords["get_coordinates()<br/>PDBParser/MMCIFParser"]
        CreateCMap["create_contact_map()<br/>8Å threshold"]
        GetPAE["get_PAE_matrix()"]
        ApplyConf["get_plddt_pae_mat()<br/>pLDDT≥70, PAE≤5"]
        CoarseGrainAF["coarse_grain()<br/>CG 1/5/10"]
        
        AF2MClass --> GetBestModel
        GetBestModel --> GetCoords
        GetBestModel --> GetPAE
        GetCoords --> CreateCMap
        CreateCMap --> ApplyConf
        GetPAE --> ApplyConf
        ApplyConf --> CoarseGrainAF
    end
    
    subgraph OtherMethods["Other Methods Processing"]
        OtherClass["Othermethods class<br/>get_other_method_preds.py"]
        AIUPred["get_aiupred_predictions()<br/>aiupred_lib"]
        DeepDiso["get_deepdisobind_predictions()<br/>Parse FASTA"]
        MorfChibi["get_morfchibi_predictions()<br/>Parse TXT"]
        AssembleInt["assemble_interfaces_for_ood_entries()"]
        
        OtherClass --> AIUPred
        OtherClass --> DeepDiso
        OtherClass --> MorfChibi
        AIUPred --> AssembleInt
        DeepDiso --> AssembleInt
        MorfChibi --> AssembleInt
    end
    
    subgraph Evaluation["Evaluation & Analysis"]
        JDClass["JudgementDay class<br/>analysis.py"]
        CreateTensors["create_ood_set_tensors()<br/>Generator function"]
        GetBasePreds["get_base_model_preds()"]
        CombinePreds["combine_diso_af_preds()<br/>Max operation"]
        GetDisorderPreds["get_preds_for_disorder_order_residues()"]
        GetAATypePreds["get_preds_for_interaction_types()"]
        CalcMetrics["calculate_metrics()<br/>torch_metrics()"]
        
        JDClass --> CreateTensors
        CreateTensors --> GetBasePreds
        CreateTensors --> CombinePreds
        CreateTensors --> GetDisorderPreds
        CreateTensors --> GetAATypePreds
        GetBasePreds --> CalcMetrics
        CombinePreds --> CalcMetrics
        GetDisorderPreds --> CalcMetrics
        GetAATypePreds --> CalcMetrics
    end
    
    subgraph Outputs["Output Files"]
        DisoNPY["Disobind_Predictions.npy<br/>6 tasks × masks"]
        AF2NPY["Predictions_af2m_results_*.npy"]
        AF3NPY["Predictions_af3_results_*.npy"]
        OtherNPY["other_methods.npy"]
        ResultsCSV["Results_OOD_set_*.csv<br/>All metrics"]
        ConfPlot["AF_confidence_plot_*.png"]
        SparsityPlot["Sparsity_F1_plot_*.png"]
    end
    
    TestSet --> PredClass
    DisorderDBs --> PredClass
    AF2Files --> AF2MClass
    AF3Files --> AF2MClass
    OtherInputs --> OtherClass
    
    PredForward --> DisoNPY
    CoarseGrainAF --> AF2NPY
    CoarseGrainAF --> AF3NPY
    AssembleInt --> OtherNPY
    
    DisoNPY --> JDClass
    AF2NPY --> JDClass
    AF3NPY --> JDClass
    OtherNPY --> JDClass
    
    CalcMetrics --> ResultsCSV
    CalcMetrics --> ConfPlot
    CalcMetrics --> SparsityPlot
```

Sources: [analysis/analysis.py:25-97](), [analysis/README.md:1-68](), [analysis/get_other_method_preds.py:14-48]()

## Core Classes and Components

The analysis system is built around four primary classes, each responsible for a distinct phase of the evaluation pipeline.

**Table: Core Analysis Classes**

| Class | File | Primary Responsibility | Key Methods |
|-------|------|----------------------|-------------|
| `Prediction` | `get_disobind_predictions.py` | Generate Disobind predictions with masks | `predict()`, `create_lip_masks()`, `create_aa_masks()`, `get_disorder_matrix()` |
| `AF2MPredictions` | `get_af_prediction.py` | Process AlphaFold2/3 structures | `get_af2m_prediction()`, `get_af3_prediction()`, `get_plddt_pae_mat()`, `coarse_grain()` |
| `Othermethods` | `get_other_method_preds.py` | Parse competing method predictions | `get_aiupred_predictions()`, `get_deepdisobind_predictions()`, `get_morfchibi_predictions()` |
| `JudgementDay` | `analysis.py` | Orchestrate evaluation and metrics | `create_ood_set_tensors()`, `eval_performance()`, `calculate_metrics()` |

Sources: [analysis/analysis.py:25-97](), [analysis/get_other_method_preds.py:14-36]()

### JudgementDay Pipeline Details

The `JudgementDay` class is the central orchestrator for evaluating results across multiple protein-protein interaction (PPI) tasks. It supports two primary modes: `"ood"` for the standard out-of-distribution test set and `"misc"` for case-specific datasets [analysis/analysis.py:42]().

**Diagram: JudgementDay Evaluation Logic**

```mermaid
graph TB
    subgraph Logic["JudgementDay::eval_performance()"]
        Tasks["get_tasks()<br/>interaction_1/5/10<br/>interface_1/5/10"]
        Generator["create_ood_set_tensors()"]
        
        Tasks --> Generator
        Generator --> Loop["For each task batch"]
        
        subgraph Models["Model Processing"]
            Base["get_base_model_preds()<br/>Disobind, AF2, AF3"]
            Combined["combine_diso_af_preds()<br/>Disobind+AF2, Disobind+AF3"]
            DisoOrder["get_preds_for_disorder_order_residues()<br/>IDR-IDR, IDR-any, order"]
            AAInter["get_preds_for_interaction_types()<br/>Aromatic, Hydrophobic, etc."]
            Others["get_other_method_preds()<br/>AIUPred, DeepDisoBind, MORFchibi"]
        end
        
        Loop --> Models
        Models --> Metrics["calculate_metrics()<br/>Recall, Precision, F1, AUROC, MCC"]
        Metrics --> Results["Results_OOD_set_19.csv"]
    end
```

Sources: [analysis/analysis.py:125-131](), [analysis/analysis.py:192-316](), [analysis/analysis.py:438-509]()

## Task and Evaluation Categories

The analysis evaluates predictions across 6 tasks (interaction vs interface × 3 coarse-graining levels) and multiple evaluation categories.

**Table: Evaluation Categories Detail**

| Field Key | Data Type | Description |
|-----------|-----------|-------------|
| `AF2_pLDDT_PAE` | np.array | AF2 predictions with confidence filtering [analysis/analysis.py:221]() |
| `Disobind` | np.array | Disobind model predictions [analysis/analysis.py:213]() |
| `IDR-IDR` | np.array | Binary mask for IDR-IDR interactions [analysis/analysis.py:254]() |
| `order` | np.array | Binary mask for ordered region interactions [analysis/analysis.py:255]() |
| `lips` | np.array | Binary mask for LIP interactions [analysis/analysis.py:261]() |
| `Aiupred` | np.array | AIUPred interface predictions (interface_1 only) [analysis/analysis.py:302]() |

Sources: [analysis/analysis.py:192-316]()

## Competing Method Comparison

The `Othermethods` class provides a unified interface for comparing Disobind against established disordered binding predictors. It handles the specific nuances of each tool, such as MORFchibi's minimum length requirement (26 residues) and thresholding heuristics [analysis/get_other_method_preds.py:188-190]().

**Table: Competing Methods Integration**

| Method | Source | Thresholding Logic |
|--------|--------|-------------------|
| AIUPred | `aiupred_lib` | 0.5 propensity cutoff [analysis/get_other_method_preds.py:91]() |
| DeepDisoBind | FASTA result parsing | `protein_binary` field [analysis/get_other_method_preds.py:102-105]() |
| MORFchibi | TSV result parsing | 0.775 cutoff + 4-residue window filter [analysis/get_other_method_preds.py:168-181]() |

Sources: [analysis/get_other_method_preds.py:65-218]()

## Dataset-Specific Evaluations

The pipeline includes specialized logic for different test sets:

*   **OOD Set**: The standard non-redundant test set derived from PDB and disorder databases [analysis/analysis.py:51-54]().
*   **Misc Dataset**: Case-specific entries (e.g., specific IDR complexes) where `JudgementDay` generates additional JSON summaries for visualization [analysis/analysis.py:63-66]().
*   **IDPPI**: Evaluation on the IDPPI protein-protein interaction dataset is facilitated by `prep_idppi_input2.py` and evaluated using scripts like `idppi_preds.py` [dataset/README.md:43-47]().

Sources: [analysis/analysis.py:56-66](), [dataset/README.md:43-56]()

## Performance Metrics

Evaluation is performed using a standard suite of metrics implemented in `torch_metrics`.

**Table: Calculated Metrics**

| Metric | Code Symbol | Description |
|--------|-------------|-------------|
| Recall | `recall` | Sensitivity to true contacts [analysis/analysis.py:726]() |
| Precision | `precision` | Accuracy of predicted contacts [analysis/analysis.py:727]() |
| F1-score | `f1` | Harmonic mean of Precision and Recall [analysis/analysis.py:728]() |
| AUROC | `auroc` | Area under the Receiver Operating Characteristic curve [analysis/analysis.py:731]() |
| MCC | `mcc` | Matthews Correlation Coefficient [analysis/analysis.py:730]() |

Sources: [analysis/analysis.py:709-736]()

---