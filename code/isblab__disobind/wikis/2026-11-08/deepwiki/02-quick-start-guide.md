# Quick Start Guide

<details>
<summary>Relevant source files</summary>

The following files were used as context for generating this wiki page:

- [README.md](README.md)
- [example/output.tar.gz](example/output.tar.gz)
- [example/pae_model_4_multimer_v3_pred_4.json](example/pae_model_4_multimer_v3_pred_4.json)
- [example/test.csv](example/test.csv)
- [example/unrelaxed_model_4_multimer_v3_pred_4.pdb](example/unrelaxed_model_4_multimer_v3_pred_4.pdb)

</details>



## Purpose and Scope

This guide provides a practical walkthrough for running your first Disobind predictions. It covers the essential steps from preparing input files to interpreting output predictions for binary protein complexes containing at least one intrinsically disordered region (IDR).

For detailed information about installation and environment setup, see [Installation and Setup](1.1). For comprehensive documentation on prediction options and parameters, see [Running Predictions](2). For output file formats and interpretation details, see [Output Files and Interpretation](2.4).

**Sources:** [README.md:1-129]()

---

## Prerequisites Check

Before running Disobind, verify your installation:

```bash
# Activate the conda environment
conda activate disobind

# Verify Python packages
python -c "import torch; import Bio; print('Dependencies OK')"

# Check model files exist
ls models/Epsilon_3_6/Version_0/
ls models/Epsilon_3_16/Version_2/
```

The Disobind models should be located in the `models/` directory with the following structure:

| Model Type | CG Level | Directory | Purpose |
|------------|----------|-----------|---------|
| Epsilon_3_6 | 10 | `models/Epsilon_3_6/Version_0/` | Interaction prediction (contact maps) |
| Epsilon_3_6 | 5 | `models/Epsilon_3_6/Version_1/` | Interaction prediction (contact maps) |
| Epsilon_3_6 | 1 | `models/Epsilon_3_6/Version_2/` | Interaction prediction (contact maps) |
| Epsilon_3_16 | 10 | `models/Epsilon_3_16/Version_0/` | Interface residue prediction |
| Epsilon_3_16 | 5 | `models/Epsilon_3_16/Version_1/` | Interface residue prediction |
| Epsilon_3_16 | 1 | `models/Epsilon_3_16/Version_2/` | Interface residue prediction |

**Sources:** [run_disobind.py:86-87](), [README.md:14-38]()

---

## Basic Prediction Workflow

The following diagram illustrates the complete workflow from input preparation to final predictions:

**Basic Disobind Prediction Flow**

```mermaid
flowchart TD
    Start["User Prepares<br/>Input CSV"]
    ParseCSV["run_disobind.py<br/>Parse CSV Input"]
    CheckCache{"Embeddings<br/>Cached?"}
    DownloadSeq["Download UniProt<br/>Sequences via API"]
    LoadCache["Load Cached<br/>Embeddings from HDF5"]
    CreateEmbed["Create T5 Embeddings<br/>Embeddings class"]
    LoadModels["Load Epsilon_3 Models<br/>get_model()"]
    RunPred["Generate Predictions<br/>Disobind.predict()"]
    SaveOutput["Save Results<br/>CSV + NPY files"]
    End["Predictions Complete"]
    
    Start --> ParseCSV
    ParseCSV --> CheckCache
    CheckCache -->|No| DownloadSeq
    CheckCache -->|Yes| LoadCache
    DownloadSeq --> CreateEmbed
    CreateEmbed --> LoadModels
    LoadCache --> LoadModels
    LoadModels --> RunPred
    RunPred --> SaveOutput
    SaveOutput --> End
```

**Sources:** [run_disobind.py:43-123](), [run_disobind.py:164-203]()

---

## Understanding Input Format

Disobind accepts protein pairs in CSV or FASTA format. Each entry represents one binary complex to predict.

### Input Format for Disobind Only (CSV)

For predictions using only Disobind (sequence-based):

```
UniProt_ID1,start1,end1,UniProt_ID2,start2,end2
```

**Example:**
```csv
P04273,95,193,P04273,95,192
```

This specifies a homodimer where:
- Protein 1: UniProt P04273, residues 95-193
- Protein 2: UniProt P04273, residues 95-192

### Input Format for Disobind+AF2 (CSV)

To combine Disobind predictions with AlphaFold2 structural information:

```
UniProt_ID1,start1,end1,UniProt_ID2,start2,end2,AF2_struct_file,AF2_json_file,chain1,chain2,offset1,offset2
```

**Example:**
```csv
P04273,95,193,P04273,95,193,./example/unrelaxed_model_4_multimer_v3_pred_4.pdb,./example/pae_model_4_multimer_v3_pred_4.json,B,C,0,0
```

### FASTA format

If proteins do not have a UniProt accession, use FASTA format:
```
>UniProt_ID1, start1, end1, AF2_struct_file_path, AF2_pae_file_path, chain1, offset1
Protein 1 Sequence
>UniProt_ID2, start2, end2, AF2_struct_file_path, AF2_pae_file_path, chain2, offset2
Protein 2 Sequence
```

**Important Constraints:**
- Protein 1 must be an IDR.
- Protein 2 may or may not be an IDR.
- Only binary complexes (AB) are supported. For non-binary complexes (ABC), convert them into binary pairs (AB, BC, AC).

**Sources:** [run_disobind.py:208-255](), [README.md:42-73](), [example/test.csv:1-2]()

---

## Running Your First Prediction

### Example 1: Basic Disobind Prediction

Use the provided example file `example/test.csv`. Run Disobind with default settings (interface prediction only, CG=1):

```bash
python run_disobind.py -i csv -f ./example/test.csv 
```

**Key Command-Line Options:**

| Flag | Default | Description |
|------|---------|-------------|
| `-i` | `csv` | Input file type (`csv` or `fasta`) |
| `-f` | Required | Path to the input file |
| `-o` | `output` | Output directory name |
| `-d` | `cpu` | Device: `cpu` or `cuda` |
| `-c` | `2` | Number of cores for downloading UniProt sequences |
| `-cm` | `False` | Predict inter-protein contact maps |
| `-cg` | `1` | Coarse-graining level: `0` (all), `1`, `5`, or `10` |

**Sources:** [run_disobind.py:126-161](), [README.md:76-93]()

---

## Understanding the Output

### Output Files

After running Disobind, the following files are generated in the output directory:

1. **CSV output file**: Contains residue-level interaction predictions.
2. **`Predictions.npy`**: Contains predictions for all input sequence fragment pairs in a nested dictionary.

### CSV Output Interpretation

| Protein1 | Residue1 | Protein2 | Residue2 |
|----------|----------|----------|----------|
| X1 | 10 | X2 | 40 |

- **Contact Map Prediction (`-cm`)**: Residue 10 in protein X1 interacts with residue 40 in protein X2.
- **Interface Prediction (Default)**: Residues in Protein1 (left) may interact with one or more residues in Protein2 (right).

**Sources:** [README.md:95-114]()

---

## Code Components Reference

The following diagram maps the main code entities involved in the prediction workflow:

**Disobind Class Structure and Workflow**

```mermaid
graph TB
    subgraph "Main Entry Point"
        RunScript["run_disobind.py<br/>Main script"]
        DisobindClass["Disobind class<br/>__init__: line 44<br/>forward: line 108"]
    end
    
    subgraph "Input Processing"
        ReadCSV["read_csv_input()<br/>line 208-255"]
        GetUniIDs["get_unique_uni_ids()<br/>line 280-295"]
        DownloadSeqs["download_uniprot_seq()<br/>line 299-327"]
        ProcessPairs["process_input_pairs()<br/>line 333-368"]
    end
    
    subgraph "Embedding Generation"
        EmbClass["Embeddings class<br/>create_input_embeddings.py"]
        CreateEmb["create_embeddings()<br/>line 375-395"]
    end
    
    subgraph "Model Loading"
        LoadModel["load_model()<br/>line 399-429"]
        GetModel["get_model()<br/>src/models/get_model.py"]
    end
    
    subgraph "Prediction"
        Predict["predict()<br/>line 532-661"]
        GetInputTensors["get_input_tensors()<br/>line 475-529"]
    end
    
    subgraph "AlphaFold Integration"
        AfPredClass["AfPrediction class<br/>line 831-1166"]
        GetConfident["get_confident_interactions()<br/>line 1150-1166"]
    end
    
    RunScript --> DisobindClass
    DisobindClass --> ReadCSV
    ReadCSV --> ProcessPairs
    ProcessPairs --> GetUniIDs
    GetUniIDs --> DownloadSeqs
    ProcessPairs --> CreateEmb
    
    DisobindClass --> LoadModel
    LoadModel --> GetModel
    
    DisobindClass --> Predict
    Predict --> GetInputTensors
    
    Predict --> AfPredClass
    AfPredClass --> GetConfident
```

**Key Classes and Methods:**

| Class/Function | Location | Purpose |
|----------------|----------|---------|
| `Disobind` | [run_disobind.py:43-826]() | Main orchestrator class |
| `Disobind.forward()` | [run_disobind.py:108-122]() | Entry point for prediction pipeline |
| `Disobind.predict()` | [run_disobind.py:532-661]() | Runs models and generates predictions |
| `AfPrediction` | [run_disobind.py:831-1166]() | Processes AlphaFold2 structures |

**Sources:** [run_disobind.py:43-1166]()

---

## AlphaFold2 Integration Strategy

Disobind can combine sequence-based predictions with AlphaFold structural contacts to improve accuracy.

**Combination Strategy:**
Disobind and AF2 predictions are combined using an element-wise maximum operation in `Disobind.predict()` [run_disobind.py:619-643]().

**Confidence Filtering:**
The `AfPrediction` class filters structural contacts based on:
- **pLDDT**: Per-residue confidence [run_disobind.py:1118-1128]().
- **PAE**: Predicted Aligned Error [run_disobind.py:1130-1148]().

**Sources:** [run_disobind.py:831-1166](), [run_disobind.py:619-643]()

---