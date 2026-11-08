# Examples and Tutorials

<details>
<summary>Relevant source files</summary>

The following files were used as context for generating this wiki page:

- [README.md](README.md)
- [example/output.tar.gz](example/output.tar.gz)
- [example/pae_model_4_multimer_v3_pred_4.json](example/pae_model_4_multimer_v3_pred_4.json)
- [example/test.csv](example/test.csv)
- [example/unrelaxed_model_4_multimer_v3_pred_4.pdb](example/unrelaxed_model_4_multimer_v3_pred_4.pdb)

</details>



This page provides practical examples and tutorials demonstrating common Disobind workflows. Each example includes concrete input formats, command-line invocations, and expected outputs with references to actual code paths and data files.

**Scope**: This page covers high-level summaries and quick-start examples. For detailed step-by-step tutorials, see the following child pages:
- [Example: Basic Prediction](#7.1) — Complete walkthrough of running predictions on a simple protein pair from input to output using the `example/test.csv` and `example/test.fasta` files.
- [Example: AlphaFold-Enhanced Predictions](#7.2) — Tutorial showing how to combine Disobind with AlphaFold predictions for improved accuracy using the example PDB and PAE JSON files.
- [Example: Custom Dataset Creation](#7.3) — Step-by-step guide for creating a custom dataset from PDB structures for training, walking through the numbered pipeline scripts 1-4.
- [Example: Performance Analysis](#7.4) — Walkthrough of running comprehensive performance analysis on predictions using `JudgementDay`, from generating OOD predictions to interpreting results CSV.

For installation and initial setup, see [Installation and Setup](#1.1). For detailed API documentation, see [API Reference](#6).

---

## Overview of Example Workflows

The Disobind repository includes complete example data in the `example/` directory demonstrating two primary use cases: standalone sequence-based prediction and structure-enhanced prediction.

### Natural Language to Code Entity Space: Prediction Pipeline
The following diagram maps the logical steps of the prediction workflow to the specific code entities responsible for them.

```mermaid
flowchart TD
    subgraph "Natural Language Space"
        Input["User provides protein pairs<br/>(UniProt IDs or FASTA)"]
        Proc["System retrieves sequences<br/>and generates embeddings"]
        Inference["Model predicts contacts<br/>and interface residues"]
        Output["User receives CSV/NPY<br/>prediction files"]
    end

    subgraph "Code Entity Space"
        CSV["example/test.csv"]
        Main["run_disobind.py"]
        UniProt["from_APIs_with_love.get_uniprot_seq()"]
        Emb["Embeddings.initialize()"]
        E3["Epsilon_3.forward()"]
        Results["Predictions.npy"]
    end

    Input -.-> CSV
    CSV --> Main
    Main --> UniProt
    UniProt --> Emb
    Emb --> E3
    E3 --> Results
    Results -.-> Output
```
**Sources**: [example/test.csv:1-2](), [run_disobind.py:43-122](), [run_disobind.py:831-1177](), [dataset/from_APIs_with_love.py]()

---

## Quick Start: Basic Disobind Prediction

The simplest prediction requires only UniProt IDs and residue ranges in CSV format. Each row represents a protein pair where Protein 1 must be an IDR and Protein 2 can be ordered or disordered.

### Command Execution
```bash
python run_disobind.py -i csv -f ./example/test.csv
```

This command triggers the `Disobind` class in `run_disobind.py`, which orchestrates sequence retrieval, embedding generation via the `Embeddings` class, and inference using the `Epsilon_3` model architecture. For a full walkthrough of this scenario, see [Example: Basic Prediction](#7.1).

**Sources**: [README.md:79-80](), [run_disobind.py:108-122](), [run_disobind.py:299-327]()

---

## AlphaFold-Enhanced Prediction

Disobind can combine its sequence-based predictions with structural confidence from AlphaFold2 or AlphaFold3. This process uses the `AfPrediction` class to parse PDB files and PAE (Predicted Aligned Error) JSON files to filter for high-confidence structural interactions.

### Code Entity Mapping: AF2 Integration
This diagram illustrates how external AlphaFold data is ingested by the `AfPrediction` system.

```mermaid
flowchart LR
    subgraph "External AF Data"
        PDB["example/unrelaxed_model_4...pdb"]
        PAE["example/pae_model_4...json"]
    end

    subgraph "AfPrediction Class"
        Parser["Bio.PDB.PDBParser"]
        Conf["get_confident_interactions()"]
        Threshold["dist_threshold=8, plddt_cutoff=70, pae_cutoff=5"]
    end

    subgraph "Integration"
        MaxOp["np.max(axis=1)"]
        Diso["Disobind.predict()"]
    end

    PDB --> Parser
    PAE --> Conf
    Parser --> Conf
    Threshold --> Conf
    Conf --> MaxOp
    Diso --> MaxOp
```

For details on configuring offsets and chain IDs for this workflow, see [Example: AlphaFold-Enhanced Predictions](#7.2).

**Sources**: [example/test.csv:2](), [run_disobind.py:831-845](), [run_disobind.py:1118-1128](), [README.md:57-62]()

---

## Dataset Creation and Analysis

Beyond running predictions, Disobind provides a pipeline for users to create their own datasets and evaluate model performance.

### Custom Dataset Creation
Users can build training datasets from PDB structures using a four-step pipeline located in the `dataset/` directory. This involves:
1. Collecting PDB IDs from disorder databases.
2. Generating binary complexes and mapping sequences using SIFTS.
3. Reducing redundancy using `MMSeqs2`.
4. Generating T5 embeddings.

For details, see [Example: Custom Dataset Creation](#7.3).

**Sources**: [README.md:117-118](), [dataset/README.md]()

### Performance Analysis
The `JudgementDay` class in the `analysis/` directory allows for comprehensive evaluation of predictions against ground truth, including the calculation of AUROC, MCC, and F1 scores. It also supports comparing Disobind against other methods like `AIUPred` or `DeepDisoBind`.

For details, see [Example: Performance Analysis](#7.4).

**Sources**: [README.md:123-124](), [analysis/README.md]()

---

## Command-Line Options Summary

The `run_disobind.py` script is the primary entry point for examples.

| Flag | Description | Default |
|------|-------------|---------|
| `-i` | Input file type (`csv` or `fasta`) | `csv` |
| `-f` | Path to the input file | Required |
| `-cm`| Predict inter-protein contact maps | `False` |
| `-cg`| Coarse-grained resolution (0, 1, 5, 10) | `1` |
| `-d` | Device (`cpu` or `cuda`) | `cpu` |

**Sources**: [README.md:85-94](), [run_disobind.py:44-82]()

---