# Generating OOD Predictions

<details>
<summary>Relevant source files</summary>

The following files were used as context for generating this wiki page:

- [analysis/get_af_prediction.py](analysis/get_af_prediction.py)
- [analysis/get_disobind_predictions.py](analysis/get_disobind_predictions.py)
- [analysis/params.py](analysis/params.py)

</details>



## Purpose and Scope

This page describes how to use the `Prediction` class to generate Disobind predictions on out-of-distribution (OOD) datasets. The `Prediction` class generates predictions across all tasks (interaction and interface) and coarse-graining levels (1, 5, 10), while simultaneously creating binary masks for disorder-specific, amino acid type-specific, and LIP-specific analysis.

For information about analyzing these predictions after generation, see [JudgementDay Analysis Pipeline](#5.1). For processing AlphaFold predictions to compare with Disobind, see [Processing AlphaFold Predictions](#5.3).

**Sources:** [analysis/get_disobind_predictions.py:1-9]()

## Overview

The `Prediction` class in [analysis/get_disobind_predictions.py:38-132]() serves as the prediction engine for evaluating Disobind models on OOD test sets. It performs three main functions:

1. **Model Inference**: Runs trained Disobind models across all 6 task configurations (Interaction/Interface at CG 1, 5, 10).
2. **Mask Generation**: Creates binary masks for specialized analysis including disorder regions, amino acid types, and Linear Interacting Peptides (LIPs).
3. **Output Packaging**: Assembles predictions and masks into a structured dictionary format saved as a `.npy` file.

The class supports two operational modes:
- **`ood` mode**: For the standard OOD test set evaluation [analysis/get_disobind_predictions.py:85-94]().
- **`misc` mode**: For case-specific datasets with specialized padding requirements (e.g., max length 200) [analysis/get_disobind_predictions.py:96-105]().

**Sources:** [analysis/get_disobind_predictions.py:38-132]()

## Workflow Diagram

```mermaid
flowchart TD
    Start([Initialize Prediction Class]) --> LoadDisorder[Load Disorder Databases<br/>load_disorder_dbs]
    
    LoadDisorder --> ReadInput[Read Input CSV<br/>read_input]
    ReadInput --> LoadSeq[Load UniProt Sequences<br/>from JSON]
    
    LoadSeq --> CheckMotifs{Motifs Dict<br/>Exists?}
    CheckMotifs -->|No| GetMotifs[get_motifs<br/>Extract LIP Motifs from MobiDB]
    CheckMotifs -->|Yes| LoadMotifs[Load Cached Motifs]
    GetMotifs --> CreateLIPMasks
    LoadMotifs --> CreateLIPMasks[create_lip_masks<br/>Binary LIP masks]
    
    CreateLIPMasks --> CreateAAMasks[create_aa_masks<br/>AA Type Masks]
    
    CreateAAMasks --> CreateEmb[create_embeddings<br/>via Embeddings class]
    
    CreateEmb --> LoopTasks{For Each Task<br/>6 total}
    
    LoopTasks --> ApplySettings[apply_settings<br/>Load model for obj+CG]
    ApplySettings --> PrepInputs[get_input_tensors<br/>Pad, CG, to Device]
    PrepInputs --> RunModel[Model Forward Pass<br/>with torch.no_grad]
    
    RunModel --> GetDisorder[get_disorder_matrix<br/>IDR-IDR, IDR-any, order]
    GetDisorder --> AssemblePred[Assemble Prediction Dict<br/>Disobind output + masks]
    
    AssemblePred --> MoreTasks{More<br/>Tasks?}
    MoreTasks -->|Yes| LoopTasks
    MoreTasks -->|No| SaveOutput[Save Predictions NPY<br/>Disobind_Predictions.npy]
    
    SaveOutput --> SaveLogs[Save Logs JSON<br/>counts per category]
    SaveLogs --> End([Complete])
```

**Sources:** [analysis/get_disobind_predictions.py:134-179](), [analysis/get_disobind_predictions.py:810-896]()

## Input Preparation

### Input File Format

The `Prediction` class accepts input in two formats specified by `self.input_type` [analysis/get_disobind_predictions.py:181-209]():

| Format | Description | Use Case |
|--------|-------------|----------|
| **CSV** | Comma-separated entry IDs | Multiple binary complexes (recommended) |
| **FASTA** | FASTA format with headers | Single binary complex |

**CSV Format [analysis/get_disobind_predictions.py:211-231]():**
`UniID1:start1:end1--UniID2:start2:end2_num,UniID1:start1:end1--UniID2:start2:end2_num,...`

**Sources:** [analysis/get_disobind_predictions.py:181-250]()

### Required Input Files

The class requires several pre-existing files depending on the mode:

**For `ood` mode [analysis/get_disobind_predictions.py:85-94]():**
* `prot_1-2_test_v_19.csv`: Input file containing prot1/2 headers.
* `Uniprot_seq.json`: UniProt sequences for all IDs.
* `Target_bcmap_test_v_19.h5`: Target contact maps.

**For `misc` mode [analysis/get_disobind_predictions.py:96-105]():**
* `misc_test_input.csv`: Input file for misc headers.
* `Uniprot_seq_misc.json`: Sequences for misc entries.
* `misc_test_target.h5`: Target maps for misc entries.

**Sources:** [analysis/get_disobind_predictions.py:85-107]()

## Prediction Configuration

### Key Parameters

The `Prediction` class constructor defines critical configuration parameters:

```mermaid
graph LR
    Config[Prediction Configuration] --> Data[Data Parameters]
    Config --> Model[Model Parameters]
    Config --> Proc[Processing Parameters]
    
    Data --> DataVer["data_version = 19"]
    Data --> EmbType["embedding_type = 'T5'"]
    Data --> Scope["scope = 'global'"]
    
    Model --> ModelVer["model_version = 19"]
    Model --> Device["device = 'cuda'"]
    Model --> Thresh["threshold = 0.5"]
    
    Proc --> Mode["mode = 'ood' | 'misc'"]
    Proc --> MaxLen["max_len = 100 (ood)<br/>200 (misc)"]
    Proc --> MultidimAvg["multidim_avg = 'global'"]
```

**Sources:** [analysis/get_disobind_predictions.py:39-65]()

### Model Parameter Files

Model parameters are loaded via the `parameter_files()` function from [analysis/params.py:6-11](). The function returns specific checkpoints and hyperparameter strings for each task and coarse-graining (CG) level.

| Task | CG Level | Model Version | File Prefix |
|------|----------|---------------|-------------|
| Interaction | 1 | 19 | `Epsilon_3_6.2` |
| Interaction | 5 | 19 | `Epsilon_3_6.1` |
| Interaction | 10 | 19 | `Epsilon_3_6` |
| Interface | 1 | 19 | `Epsilon_3_16` |
| Interface | 5 | 19 | `Epsilon_3_16.1` |
| Interface | 10 | 19 | `Epsilon_3_16.2` |

**Sources:** [analysis/params.py:14-59]()

## Running Predictions

### Execution

To generate OOD predictions, run the script from the `analysis` directory:

```bash
python get_disobind_predictions.py
```

The script initializes the `Prediction` class and calls the `forward()` method [analysis/get_disobind_predictions.py:134-179]().

### Per-Entry Processing

For each entry in the OOD set, the `predict()` method [analysis/get_disobind_predictions.py:810-896]() executes a nested loop over objectives (Interaction, Interface) and CG levels (1, 5, 10).

```mermaid
flowchart TD
    LoopEntry[For each entry_id] --> GetEmb[Load prot1_emb, prot2_emb]
    GetEmb --> GetTarget[Load target_cmap]
    
    GetTarget --> LoopObj{For each objective}
    
    LoopObj --> LoopCG{For each CG}
    
    LoopCG --> ApplySettings[apply_settings<br/>Load specific model]
    ApplySettings --> GetInputs[get_input_tensors<br/>Padding & Coarse Graining]
    
    GetInputs --> ModelForward[Model Forward Pass]
    ModelForward --> ApplyMask[Apply target_mask]
    
    ApplyMask --> ExtractOutput[extract_model_output<br/>Reshape to task format]
    ExtractOutput --> GetDisorderMat[get_disorder_matrix<br/>IDR-IDR, IDR-any, order]
    
    GetDisorderMat --> AssembleDict[Assemble Task Dict]
    AssembleDict --> NextTask{Next Task?}
```

**Sources:** [analysis/get_disobind_predictions.py:810-895]()

## Output Structure

### Prediction Dictionary Format

The `Prediction` class generates a nested dictionary saved as `Disobind_Predictions.npy` [analysis/get_disobind_predictions.py:175]():

```mermaid
graph TB
    Root[predictions dict] --> Entry1[entry_id_1]
    
    Entry1 --> Task1[interaction_1]
    Entry1 --> Task2[interface_5]
    
    Task1 --> Diso[Disobind<br/>np.array]
    Task1 --> TMask[target_mask]
    Task1 --> IDRIDR[IDR-IDR mask]
    Task1 --> P1AA[prot1_aa_mask<br/>dict]
    Task1 --> P1LIP[prot1_lips_mask]
    
    P1AA --> DPAA1[disorder_promoting_aa]
    P1AA --> Arom1[aromatic_aa]
```

**Sources:** [analysis/get_disobind_predictions.py:875-886]()

## Binary Mask Generation

### Disorder Masks

The `get_disorder_matrix()` method [analysis/get_disobind_predictions.py:454-591]() creates three types of disorder masks:
1. **IDR-IDR**: Residue pairs where both residues are disordered.
2. **IDR-any**: Residue pairs where at least one residue is disordered.
3. **order**: Residue pairs where both residues are ordered.

These are generated by querying `find_disorder_regions` [analysis/get_disobind_predictions.py:437-444]() and then applying `MaxPool` if the CG level is greater than 1 [analysis/get_disobind_predictions.py:532-555]().

**Sources:** [analysis/get_disobind_predictions.py:454-591]()

### LIP and Amino Acid Masks

* **LIP Masks**: Extracted from MobiDB annotations using `get_motifs()` [analysis/get_disobind_predictions.py:594-651]() and `create_lip_masks()` [analysis/get_disobind_predictions.py:653-685]().
* **AA Type Masks**: Categorized in `create_aa_masks()` [analysis/get_disobind_predictions.py:687-807]() into:
    * **Disorder-promoting**: R, P, Q, E, G, S, A, K [analysis/get_disobind_predictions.py:722-723]().
    * **Aromatic**: F, Y, W [analysis/get_disobind_predictions.py:724-725]().
    * **Hydrophobic**: A, V, L, I, P, M, F, W [analysis/get_disobind_predictions.py:726-727]().
    * **Polar**: S, T, C, N, Q, Y, D, E, K, R, H [analysis/get_disobind_predictions.py:728-729]().

**Sources:** [analysis/get_disobind_predictions.py:594-807]()

## Model Loading and Settings

### apply_settings()

The `apply_settings()` method [analysis/get_disobind_predictions.py:306-344]() configures the model based on the prediction task:
* Sets `objective0` to `"interaction"` or `"interface"`.
* Sets `objective1` to the CG level (1, 5, or 10).
* Sets `objective2` (binning) to `True` if CG > 1.
* Loads the model weights using `get_model` and `load_state_dict` [analysis/get_disobind_predictions.py:284-304]().

**Sources:** [analysis/get_disobind_predictions.py:284-344]()

### extract_model_output()

After inference, `extract_model_output()` [analysis/get_disobind_predictions.py:898-924]() reshapes the raw model output:
* For **interaction** tasks: Reshapes to `(max_len/CG, max_len/CG)`.
* For **interface** tasks: Reshapes to `(2 * max_len/CG, 1)`.

**Sources:** [analysis/get_disobind_predictions.py:898-924]()

---