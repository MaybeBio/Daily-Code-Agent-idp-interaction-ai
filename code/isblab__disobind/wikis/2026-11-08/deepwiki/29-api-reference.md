# API Reference

<details>
<summary>Relevant source files</summary>

The following files were used as context for generating this wiki page:

- [analysis/get_af_prediction.py](analysis/get_af_prediction.py)
- [analysis/get_disobind_predictions.py](analysis/get_disobind_predictions.py)
- [analysis/params.py](analysis/params.py)
- [example/test.fasta](example/test.fasta)
- [run_disobind.py](run_disobind.py)

</details>



This page provides comprehensive API documentation for the key classes, functions, and modules in the Disobind codebase. The API is organized by functional area, covering prediction, analysis, dataset creation, and utility functions.

For practical usage examples, see [Examples and Tutorials](#7). For conceptual understanding of the system architecture, see [Overview](#1). For details on model architecture, see [Model Architecture](#4).

---

## Core System Architecture: Code Entity Mapping

The following diagram maps the major system components to their code entities:

```mermaid
graph TB
    subgraph "Prediction System"
        RunDiso["run_disobind.py"]
        DisoClass["Disobind class"]
        AfClass["AfPrediction class"]
        
        RunDiso --> DisoClass
        RunDiso --> AfClass
    end
    
    subgraph "Analysis System"
        AnalysisPy["analysis/analysis.py"]
        JudgeClass["JudgementDay class"]
        PredPy["get_disobind_predictions.py"]
        PredClass["Prediction class"]
        AfPredPy["get_af_prediction.py"]
        AF2MClass["AF2MPredictions class"]
        OtherPy["get_other_method_preds.py"]
        OtherClass["Othermethods class"]
        
        AnalysisPy --> JudgeClass
        PredPy --> PredClass
        AfPredPy --> AF2MClass
        OtherPy --> OtherClass
    end
    
    subgraph "Dataset System"
        EmbPy["create_input_embeddings.py"]
        EmbClass["Embeddings class"]
        UtilPy["dataset/utility.py"]
        UtilFuncs["Utility Functions"]
        
        EmbPy --> EmbClass
        UtilPy --> UtilFuncs
    end
    
    subgraph "Configuration"
        ParamsPy["analysis/params.py"]
        ParamFunc["parameter_files()"]
        
        ParamsPy --> ParamFunc
    end
    
    DisoClass --> EmbClass
    PredClass --> EmbClass
    DisoClass --> AfClass
    JudgeClass --> PredClass
    JudgeClass --> AF2MClass
    JudgeClass --> OtherClass
```

**Sources:** [run_disobind.py:44-45](), [analysis/analysis.py:25-26](), [analysis/get_disobind_predictions.py:38-39](), [analysis/get_af_prediction.py:28-29](), [analysis/get_other_method_preds.py:14-15](), [dataset/create_input_embeddings.py:18-19](), [analysis/params.py:6-7]()

---

## Prediction Classes

### Disobind Class

The main entry point for running Disobind predictions on protein pairs. Handles input processing, embedding generation, model loading, and prediction output.

**Location:** [run_disobind.py:44-826]()

**Key Attributes:**

| Attribute | Type | Description |
|-----------|------|-------------|
| `input_file` | str | Path to CSV/FASTA file containing protein pairs [run_disobind.py:51-51]() |
| `cores` | int | Number of CPU cores for parallel processing [run_disobind.py:55-55]() |
| `predict_cmap` | bool | Whether to predict contact maps (default: False) [run_disobind.py:57-57]() |
| `required_cg` | int | Coarse-graining level (0, 1, 5, or 10) [run_disobind.py:59-59]() |
| `output_dir` | str | Directory for output files [run_disobind.py:61-61]() |
| `device` | str | Device for inference ("cpu" or "cuda") [run_disobind.py:68-68]() |
| `threshold` | float | Contact probability threshold (default: 0.5) [run_disobind.py:70-70]() |
| `predictions` | dict | Dictionary storing all predictions [run_disobind.py:84-84]() |

**Key Methods:**

```mermaid
graph LR
    forward["forward()"] --> process["process_input_file()"]
    forward --> pairs["process_input_pairs()"]
    forward --> getpred["get_predictions()"]
    
    getpred --> create["create_embeddings()"]
    getpred --> predict["predict()"]
    
    predict --> load["load_model()"]
    predict --> apply["apply_settings()"]
    predict --> getinput["get_input_tensors()"]
```

**Main Methods:**

- **`forward()`** [run_disobind.py:111-126]()
  - Main execution pipeline orchestrating input processing and prediction.
  - Saves final predictions to `.npy` file [run_disobind.py:126-126]().

- **`process_input_file()`** [run_disobind.py:212-258]()
  - Supports CSV (multiple jobs) and FASTA (single pair) formats.
  - Returns: `(entry_ids, af_dict)` where `af_dict` contains AlphaFold metadata.

- **`create_embeddings(batch)`** [run_disobind.py:375-395]()
  - Generates T5 embeddings using the `Embeddings` class.
  - Stores embeddings in HDF5 format [run_disobind.py:106-106]().

- **`predict(required_tasks, af_dict)`** [run_disobind.py:532-661]()
  - Generates predictions for all tasks.
  - Combines Disobind and AlphaFold predictions if structural files are provided.

For details, see [Core Prediction Classes](#6.1).

**Sources:** [run_disobind.py:44-826]()

---

### AfPrediction Class

Processes AlphaFold2/3 structure predictions and extracts confident contacts based on pLDDT and PAE metrics.

**Location:** [run_disobind.py:831-1172]()

**Key Methods:**

- **`get_confident_interactions(prot1_res, prot2_res)`** [run_disobind.py:986-1128]()
  - Extracts contacts from AF2/3 structures.
  - Applies filters: pLDDT ≥ 70 [run_disobind.py:74-74]() and PAE ≤ 5 [run_disobind.py:76-76]().

**Sources:** [run_disobind.py:831-1172]()

---

## Analysis Classes

### JudgementDay Class

Orchestrates comprehensive evaluation of Disobind, AlphaFold2/3, and other methods on out-of-distribution (OOD) test sets.

**Location:** [analysis/analysis.py:25-742]()

**Key Methods:**

- **`forward()`** [analysis/analysis.py:107-123]()
  - Main execution pipeline for loading predictions and calculating performance metrics.

- **`calculate_metrics(pred, target, ...)`** [analysis/analysis.py:709-736]()
  - Computes Recall, Precision, F1, AvgPrecision, MCC, AUROC, and Accuracy.

For details, see [Analysis Classes](#6.2).

**Sources:** [analysis/analysis.py:25-742]()

---

### Prediction Class

Generates Disobind predictions on OOD/test datasets with specialized masks for disorder regions, LIPs, and amino acid types.

**Location:** [analysis/get_disobind_predictions.py:38-1031]()

**Key Methods:**

- **`create_lip_masks(headers)`** [analysis/get_disobind_predictions.py:639-684]()
  - Creates binary masks for Linear Interacting Peptides (LIPs).

- **`create_aa_masks(headers)`** [analysis/get_disobind_predictions.py:741-805]()
  - Creates masks for aromatic, hydrophobic, polar, and disorder-promoting residues.

**Sources:** [analysis/get_disobind_predictions.py:38-1031]()

---

### AF2MPredictions Class

Processes pre-existing AlphaFold2/3 predictions to extract contact maps at multiple coarse-graining levels (1, 5, 10).

**Location:** [analysis/get_af_prediction.py:28-498]()

**Key Methods:**

- **`get_best_model(model_path, header)`** [analysis/get_af_prediction.py:96-132]()
  - Selects the top-ranked model from AF2/3 outputs based on confidence scores.

**Sources:** [analysis/get_af_prediction.py:28-498]()

---

## Dataset Creation Classes

### Embeddings Class

Generates T5 protein embeddings and manages dataset splits (train/dev/test) with necessary padding.

**Location:** [dataset/create_input_embeddings.py:18-713]()

**Key Methods:**

- **`forward()`** [dataset/create_input_embeddings.py:126-159]()
  - High-level workflow for creating dataset partitions and saving to disk.

- **`apply_padding(prot1, prot2, cmap, key)`** [dataset/create_input_embeddings.py:371-408]()
  - Pads embeddings and contact maps to `max_len` [dataset/create_input_embeddings.py:46-46]().

For details, see [Dataset Creation Classes](#6.3).

**Sources:** [dataset/create_input_embeddings.py:18-713]()

---

## Utility Functions

Disobind relies on a suite of utility functions for structural analysis and sequence processing.

### Structural Utilities

- **`load_PDB(pdb, pdb_path)`** [dataset/utility.py:26-59]()
  - Loads structure files in PDB or mmCIF format.

- **`get_contact_map(coords1, coords2, contact_threshold)`** [dataset/utility.py:108-124]()
  - Computes binary contact maps from CA coordinates.

### Sequence & Disorder Utilities

- **`load_disorder_dbs(...)`** [dataset/utility.py:850-866]()
  - Loads DisProt, IDEAL, and MobiDB databases for region identification.

- **`consolidate_regions(positions, min_len)`** [dataset/utility.py:810-831]()
  - Merges fragmented residue lists into (start, end) tuples.

For details, see [Utility Functions](#6.4).

**Sources:** [dataset/utility.py:1-1200+](), [run_disobind.py:36-37]()

---

## Parameter Configuration

### parameter_files Function

Provides model checkpoint paths and metadata for specific model versions (e.g., version 19).

**Location:** [analysis/params.py:6-60]()

**Mapping Example:**
- `interaction_1` maps to model `Epsilon_3_6.2` [analysis/params.py:19-22]().
- `interface_1` maps to model `Epsilon_3_16` [analysis/params.py:43-46]().

**Sources:** [analysis/params.py:1-60]()

---