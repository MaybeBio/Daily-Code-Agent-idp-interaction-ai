# Disorder-Specific Analysis

<details>
<summary>Relevant source files</summary>

The following files were used as context for generating this wiki page:

- [analysis/analysis.py](analysis/analysis.py)
- [analysis/get_af_prediction.py](analysis/get_af_prediction.py)
- [analysis/get_disobind_predictions.py](analysis/get_disobind_predictions.py)
- [analysis/get_other_method_preds.py](analysis/get_other_method_preds.py)
- [analysis/params.py](analysis/params.py)

</details>



## Purpose and Scope

This page explains how Disobind evaluates predictions specifically for different disorder contexts: interactions between intrinsically disordered regions (IDR-IDR) versus interactions between ordered regions. This analysis reveals how well Disobind and competing methods perform on disordered versus structured protein interactions.

For the overall analysis pipeline, see [JudgementDay Analysis Pipeline](#5.1). For generating the base predictions that are filtered by disorder type, see [Generating OOD Predictions](#5.2). For comparison with partner-independent interface predictors, see [Comparing with Other Methods](#5.6).

**Key Points:**
- Disorder-specific analysis is performed only at **CG=1** (residue-level resolution) [analysis/analysis.py:464-465]().
- Evaluates performance on **IDR-IDR** and **ordered** interaction categories [analysis/analysis.py:354-360]().
- Uses specialized masks for **LIPs**, **aromatic**, **hydrophobic**, **polar**, and **disorder-promoting** amino acids [analysis/get_disobind_predictions.py:6-9]().
- Uses disorder annotations from DisProt, IDEAL, and MobiDB databases [analysis/get_disobind_predictions.py:78-81]().
- Compares Disobind, AlphaFold2, AlphaFold3, and combined predictions [analysis/analysis.py:361-362]().

Sources: [analysis/analysis.py:1-8](), [analysis/get_disobind_predictions.py:1-9](), [analysis/analysis.py:464-465](), [analysis/analysis.py:354-360]()

## Overview

Disorder-specific analysis addresses a fundamental question: does Disobind perform better on disordered versus ordered protein interactions? By creating binary masks that identify disordered and ordered residues, the system filters predictions into distinct categories and calculates separate performance metrics for each.

The analysis evaluates:
1. **IDR-IDR interactions**: Both proteins have disordered regions interacting [analysis/get_disobind_predictions.py:538-540]().
2. **Ordered interactions**: Both proteins have ordered regions interacting [analysis/get_disobind_predictions.py:560-562]().
3. **Residue-type specific interactions**: Performance on aromatic, hydrophobic, polar, and disorder-promoting residues [analysis/get_disobind_predictions.py:6-7]().
4. **LIPs (Linear Interaction Peptides)**: Specific analysis for motifs that mediate IDR interactions [analysis/get_disobind_predictions.py:8]().

```mermaid
graph TB
    subgraph "Disorder Database Sources"
        DisProt["DisProt Database<br/>(disprot_path)"]
        IDEAL["IDEAL Database<br/>(ideal_path)"]
        MobiDB["MobiDB Database<br/>(mobidb_path)"]
    end
    
    subgraph "Disorder Mask Creation"
        LoadDBs["load_disorder_dbs()<br/>Load all databases"]
        FindRegions["find_disorder_regions()<br/>Extract disorder annotations"]
        GetPositions["get_disordered_positions()<br/>Convert ranges to positions"]
        CreateMatrix["get_disorder_matrix()<br/>Create binary masks"]
    end
    
    subgraph "Binary Masks"
        IDRIDR["IDR-IDR Matrix<br/>p1_pos * p2_pos.T"]
        IDRAny["IDR-any Matrix<br/>p1_pos + p2_pos.T"]
        OrderMat["Order Matrix<br/>p1_ordered * p2_ordered.T"]
    end
    
    subgraph "Prediction Filtering"
        DisoPreds["Disobind Predictions"]
        AF2Preds["AF2 Predictions"]
        AF3Preds["AF3 Predictions"]
        CombinedPreds["Combined Predictions"]
        
        FilterIDR["Filter by IDR-IDR mask<br/>preds * IDR-IDR matrix"]
        FilterOrder["Filter by Order mask<br/>preds * Order matrix"]
    end
    
    subgraph "Evaluation Results"
        IDRMetrics["IDR-IDR Metrics<br/>Recall/Precision/F1"]
        OrderMetrics["Ordered Metrics<br/>Recall/Precision/F1"]
        Comparison["Compare Performance<br/>Disordered vs Ordered"]
    end
    
    DisProt --> LoadDBs
    IDEAL --> LoadDBs
    MobiDB --> LoadDBs
    LoadDBs --> FindRegions
    FindRegions --> GetPositions
    GetPositions --> CreateMatrix
    
    CreateMatrix --> IDRIDR
    CreateMatrix --> IDRAny
    CreateMatrix --> OrderMat
    
    IDRIDR --> FilterIDR
    OrderMat --> FilterOrder
    
    DisoPreds --> FilterIDR
    DisoPreds --> FilterOrder
    AF2Preds --> FilterIDR
    AF2Preds --> FilterOrder
    AF3Preds --> FilterIDR
    AF3Preds --> FilterOrder
    CombinedPreds --> FilterIDR
    CombinedPreds --> FilterOrder
    
    FilterIDR --> IDRMetrics
    FilterOrder --> OrderMetrics
    
    IDRMetrics --> Comparison
    OrderMetrics --> Comparison
```

**Diagram: Disorder-Specific Analysis Workflow**

Sources: [analysis/get_disobind_predictions.py:134-179](), [analysis/analysis.py:354-412](), [analysis/get_disobind_predictions.py:538-562]()

## Disorder Database Integration

The disorder-specific analysis relies on three disorder databases to identify intrinsically disordered regions:

| Database | Purpose | File Path |
|----------|---------|-----------|
| **DisProt** | Experimentally validated disorder annotations | `../database/input_files/DisProt.csv` [analysis/get_disobind_predictions.py:79]() |
| **IDEAL** | Intrinsically disordered proteins with function | `../database/input_files/IDEAL.csv` [analysis/get_disobind_predictions.py:80]() |
| **MobiDB** | Computational disorder predictions and LIPs | `../database/input_files/MobiDB.csv` [analysis/get_disobind_predictions.py:81]() |

The `Prediction` class loads these databases during initialization:

```python
# From analysis/get_disobind_predictions.py
self.disprot, self.ideal, self.mobidb = load_disorder_dbs(
    disprot_path = self.disprot_path, 
    ideal_path = self.ideal_path, 
    mobidb_path = self.mobidb_path 
)
```
[analysis/get_disobind_predictions.py:142-144]()

The disorder annotations are consolidated across databases using `find_disorder_regions()` from `dataset.utility`, which merges overlapping disorder annotations and returns unified disorder ranges for each UniProt ID [analysis/get_disobind_predictions.py:423-430]().

Sources: [analysis/get_disobind_predictions.py:78-144](), [analysis/get_disobind_predictions.py:422-450]()

## Creating Specialized Masks

The `Prediction` class creates binary masks for different amino acid properties and interaction motifs to enable fine-grained analysis.

### LIP and Amino Acid Type Masks
The pipeline creates masks for specific biological features:
1. **LIPs (Linear Interaction Peptides)**: Identified using `create_lip_masks()` [analysis/get_disobind_predictions.py:164]().
2. **Amino Acid Types**: Identified using `create_aa_masks()` [analysis/get_disobind_predictions.py:165]().
   - **Disorder-promoting**: A, R, G, Q, S, P, E, K.
   - **Aromatic**: F, W, Y, H.
   - **Hydrophobic**: V, I, L, M, F, W, Y.
   - **Polar**: S, T, N, Q.

### Interaction Matrices
For each binary complex, binary vectors are created where `1` indicates the presence of a feature and `0` indicates its absence. These are then converted into interaction matrices:

```python
# IDR-IDR interactions (both proteins disordered)
disorder_mat1 = p1_pos * p2_pos.T [analysis/get_disobind_predictions.py:538-540]()

# Ordered interactions (both proteins ordered)
order_mat = p1_ordered * p2_ordered.T [analysis/get_disobind_predictions.py:560-562]()
```

For **interface** tasks, the vectors are concatenated:
```python
# target: ( N, L1, L2 ) --> ( N, L1+L2 )
target = torch.cat( ( p1_target, p2_target ), axis = 0 ) [analysis/analysis.py:176]()
```

Sources: [analysis/get_disobind_predictions.py:538-588](), [analysis/analysis.py:176](), [analysis/get_disobind_predictions.py:164-165]()

## Interaction Type Categories

The disorder-specific analysis evaluates two primary categories:

### IDR-IDR Interactions
**Definition**: Interactions where **both** proteins have disordered residues at the interacting positions [analysis/get_disobind_predictions.py:538-540]().
**Biological Significance**: These represent interactions between intrinsically disordered regions, which are challenging for structure-based methods like AlphaFold since disorder implies lack of stable structure.

### Ordered Interactions
**Definition**: Interactions where **both** proteins have ordered (structured) residues at the interacting positions [analysis/get_disobind_predictions.py:560-562]().
**Biological Significance**: These represent traditional protein-protein interactions between structured domains, which AlphaFold is designed to predict accurately.

Sources: [analysis/get_disobind_predictions.py:538-563](), [analysis/analysis.py:363-364]()

## Filtering Predictions by Disorder Type

The `JudgementDay` class filters predictions by disorder type using the `get_preds_for_disorder_order_residues()` method [analysis/analysis.py:354]().

```mermaid
graph LR
    subgraph "Input Predictions"
        DisoPred["Disobind Predictions<br/>ood_dict['Disobind']"]
        AF2Pred["AF2 Predictions<br/>ood_dict['AF2_pLDDT_PAE']"]
        AF3Pred["AF3 Predictions<br/>ood_dict['AF3_pLDDT_PAE']"]
        CombPred["Combined Predictions<br/>preds_dict['Disobind_AF2']"]
    end
    
    subgraph "Disorder Masks"
        IDRMask["IDR-IDR Mask<br/>ood_dict['IDR-IDR']"]
        OrderMask["Order Mask<br/>ood_dict['order']"]
        TargetCmap["Target Contact Map<br/>ood_dict['targets']"]
    end
    
    subgraph "Filtering Operation"
        ApplyIDR["Element-wise Multiply<br/>pred * IDR-IDR mask"]
        ApplyOrder["Element-wise Multiply<br/>pred * Order mask"]
        FilterTarget["Filter Target<br/>target * mask"]
    end
    
    subgraph "Filtered Results"
        IDRResults["Disobind_IDR-IDR<br/>AF2_pLDDT_PAE_IDR-IDR<br/>AF3_pLDDT_PAE_IDR-IDR<br/>Disobind_AF2_IDR-IDR"]
        OrderResults["Disobind_order<br/>AF2_pLDDT_PAE_order<br/>AF3_pLDDT_PAE_order<br/>Disobind_AF2_order"]
    end
    
    DisoPred --> ApplyIDR
    AF2Pred --> ApplyIDR
    AF3Pred --> ApplyIDR
    CombPred --> ApplyIDR
    
    DisoPred --> ApplyOrder
    AF2Pred --> ApplyOrder
    AF3Pred --> ApplyOrder
    CombPred --> ApplyOrder
    
    IDRMask --> ApplyIDR
    OrderMask --> ApplyOrder
    TargetCmap --> FilterTarget
    
    ApplyIDR --> IDRResults
    ApplyOrder --> OrderResults
    FilterTarget --> IDRResults
    FilterTarget --> OrderResults
```

**Diagram: Prediction Filtering by Disorder Type**

The filtering process applies the binary masks to both the prediction and the target contact map to isolate specific interaction types [analysis/analysis.py:375-376]().

Sources: [analysis/analysis.py:354-382]()

## Model Comparison Across Disorder Categories

The disorder-specific analysis evaluates four model types on both IDR-IDR and ordered interactions [analysis/analysis.py:361-362]():

| Model Type | Description | Key |
|------------|-------------|-----|
| **Disobind** | Sequence-based Disobind predictions | `"Disobind"` |
| **AF2** | AlphaFold2 predictions (with confidence filtering) | `"AF2_pLDDT_PAE"` |
| **AF3** | AlphaFold3 predictions (with confidence filtering) | `"AF3_pLDDT_PAE"` |
| **Disobind+AF2** | Combined Disobind and AF2 predictions | `"Disobind_AF2"` |

### Evaluation Scope
**Performed only at CG=1**: Disorder-specific analysis is restricted to residue-level predictions [analysis/analysis.py:464-465]().

**Metrics Calculated**: For each model-disorder combination, the following metrics are computed using `torch_metrics()` [analysis/analysis.py:476]():
- Recall, Precision, F1-score
- Average Precision, MCC, AUROC, Accuracy [src/metrics.py:20-21]() (referenced in [analysis/analysis.py:20]())

Sources: [analysis/analysis.py:438-509](), [analysis/analysis.py:464-470](), [analysis/analysis.py:361-362]()

## Implementation Code Flow

The complete disorder-specific analysis involves coordination between two classes:

```mermaid
graph TB
    subgraph "Prediction Class<br/>(get_disobind_predictions.py)"
        InitPred["__init__()<br/>Set paths to disorder DBs"]
        LoadDBs["forward()<br/>Load DisProt, IDEAL, MobiDB"]
        PredLoop["predict() loop<br/>For each OOD entry"]
        GetDisOrder["get_disorder_matrix()<br/>For each task (obj, cg)"]
        GetDisPos["get_disordered_positions()<br/>Query disorder databases"]
        CreateMasks["Create binary masks<br/>IDR-IDR, IDR-any, Order"]
        StorePreds["Store in predictions dict<br/>with disorder masks"]
    end
    
    subgraph "JudgementDay Class<br/>(analysis.py)"
        InitJudge["__init__()<br/>Set paths to prediction files"]
        LoadPreds["forward()<br/>Load all predictions"]
        CreateOOD["create_ood_set_tensors()<br/>Generator for each task"]
        EvalPerf["eval_performance()<br/>Main evaluation loop"]
        CheckCG["Check if task is CG=1<br/>interaction_1 or interface_1"]
        FilterDis["get_preds_for_disorder_order_residues()<br/>Apply masks to predictions"]
        CalcMetrics["calculate_metrics()<br/>For each model × disorder type"]
        SaveResults["Save to CSV<br/>full_results_file"]
    end
    
    InitPred --> LoadDBs
    LoadDBs --> PredLoop
    PredLoop --> GetDisOrder
    GetDisOrder --> GetDisPos
    GetDisPos --> CreateMasks
    CreateMasks --> StorePreds
    
    InitJudge --> LoadPreds
    LoadPreds --> CreateOOD
    CreateOOD --> EvalPerf
    EvalPerf --> CheckCG
    CheckCG -->|"Yes (CG=1)"| FilterDis
    FilterDis --> CalcMetrics
    CalcMetrics --> SaveResults
    
    StorePreds -.->|"predictions.npy<br/>with disorder masks"| LoadPreds
```

**Diagram: Code Flow for Disorder-Specific Analysis**

### Key Methods

**In `Prediction` class** ([analysis/get_disobind_predictions.py]()):

| Method | Lines | Purpose |
|--------|-------|---------|
| `get_disordered_positions()` | 422-450 | Query disorder databases for a UniProt ID |
| `get_disorder_matrix()` | 454-592 | Create IDR-IDR, IDR-any, and Order matrices |
| `predict()` | 810-896 | Main prediction loop that calls disorder matrix creation |

**In `JudgementDay` class** ([analysis/analysis.py]()):

| Method | Lines | Purpose |
|--------|-------|---------|
| `create_ood_set_tensors()` | 192-316 | Load predictions and masks into unified format |
| `get_preds_for_disorder_order_residues()` | 354-382 | Filter predictions by disorder type |
| `eval_performance()` | 438-509 | Calculate metrics for all model-disorder combinations |

Sources: [analysis/get_disobind_predictions.py:134-896](), [analysis/analysis.py:25-509]()

## Output Files

Disorder-specific analysis results are written to:

**Main Results**: `Analysis_OOD_{version}_{iptm}/Results_OOD_set_{version}.csv` [analysis/analysis.py:84-85]()
- Contains metrics for all models including disorder-specific variants.

**Logs**: `Predictions_ood_v_{version}/Logs.json` [analysis/get_disobind_predictions.py:129]()
- Contains disorder statistics: `"IDR-IDR_interactions"`, `"IDR-any_interactions"`, `"Ordered_interactions"` [analysis/get_disobind_predictions.py:575-577]().

**Disorder Dictionary**: `Predictions_ood_v_{version}/disorder_dict.json` [analysis/get_disobind_predictions.py:125]()
- Cached disorder positions for each UniProt ID to avoid redundant database queries [analysis/get_disobind_predictions.py:153-155]().

Sources: [analysis/get_disobind_predictions.py:84-129](), [analysis/analysis.py:81-96](), [analysis/get_disobind_predictions.py:575-577]()

---