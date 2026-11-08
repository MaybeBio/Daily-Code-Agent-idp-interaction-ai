# Training Your Own Models

<details>
<summary>Relevant source files</summary>

The following files were used as context for generating this wiki page:

- [src/build_model.py](src/build_model.py)
- [src/hparams_search.py](src/hparams_search.py)
- [src/model_versions.py](src/model_versions.py)
- [src/models/Epsilon_3.py](src/models/Epsilon_3.py)
- [src/models/get_model.py](src/models/get_model.py)
- [src/utils.py](src/utils.py)

</details>



This page provides a step-by-step guide for training custom Disobind models on your own datasets. This includes preparing training data, configuring model hyperparameters, running the training pipeline, and evaluating results.

**Scope:** This page focuses on the practical workflow for training models. For details about the Epsilon_3 architecture, see [Epsilon_3 Network Design](#4.1). For configuration file specifications, see [Model Configuration Files](#4.2). For internal training mechanics, see [Training Pipeline](#4.3). For dataset preparation from PDB files, see [Dataset Creation Pipeline](#3).

---

## Training Workflow Overview

```mermaid
flowchart TD
    Start([Start Training Process]) --> CheckData{Dataset<br/>Ready?}
    
    CheckData -->|No| PrepData[Prepare Dataset<br/>See Section 3]
    PrepData --> CreateEmb[Generate Embeddings<br/>create_input_embeddings.py]
    CreateEmb --> SplitData[Split Train/Dev/Test<br/>90%/5%/5%]
    
    CheckData -->|Yes| LoadData[Load Existing Dataset]
    SplitData --> LoadData
    
    LoadData --> ConfigModel[Configure Model<br/>model_versions.py]
    ConfigModel --> SetObjective{Select Task}
    
    SetObjective --> Interaction[interaction_bin<br/>Contact Map Prediction]
    SetObjective --> Interface[interface_bin<br/>Interface Residue Prediction]
    
    Interaction --> SetCG{Coarse-Graining<br/>Level?}
    Interface --> SetCG
    
    SetCG --> CG1[CG = 1<br/>Residue-level]
    SetCG --> CG5[CG = 5<br/>5-residue bins]
    SetCG --> CG10[CG = 10<br/>10-residue bins]
    
    CG1 & CG5 & CG10 --> CreateYAML[Generate Config YAML<br/>version_*.yml]
    
    CreateYAML --> RunTrain[Run Training<br/>hparams_search.py]
    
    RunTrain --> TrainLoop[Training Loop<br/>25 epochs default]
    TrainLoop --> ValidStep[Validation Step<br/>Calculate Metrics]
    ValidStep --> CheckEpoch{More<br/>Epochs?}
    CheckEpoch -->|Yes| TrainLoop
    CheckEpoch -->|No| Calibrate[Calibrate Model<br/>Platt/Beta/Isotonic]
    
    Calibrate --> TestEval[Test Set Evaluation]
    TestEval --> SaveModel[Save Model Checkpoint<br/>.pth file]
    SaveModel --> SaveMetrics[Save Metrics<br/>CSV + Summary]
    
    SaveMetrics --> End([Training Complete])
```

**Purpose:** This workflow illustrates the complete process from dataset preparation through model training to checkpoint generation. The key decision points are task selection (interaction vs interface) and coarse-graining level (1, 5, or 10).

**Sources:** [src/model_versions.py:1-152](), [src/build_model.py:25-48]()

---

## Prerequisites

Before training models, ensure you have:

| Requirement | Description | Location |
|------------|-------------|----------|
| **Training Dataset** | Non-redundant protein complex dataset with contact maps | `../database/v_21/` |
| **Embeddings** | T5/ESM2/BERT embeddings for all sequences | `../database/v_21/T5/global-None/` |
| **Split Files** | Pre-split train/dev/test sets (.npy format) | Same as embeddings directory |
| **Configuration Template** | Python script to generate YAML config | `src/model_versions.py` |
| **Training Script** | Main training execution script | `src/hparams_search.py` |
| **GPU Resources** | CUDA-enabled GPU recommended | Set in `model_versions.py` |

**Important:** The dataset must be prepared using the pipeline described in [Dataset Creation Pipeline](#3). Training requires pre-computed embeddings generated via [Generating Embeddings](#3.4).

**Sources:** [src/model_versions.py:17-101](), [src/hparams_search.py:81-91]()

---

## Step 1: Prepare Training Dataset

### Dataset File Structure

Your training data should follow this structure:

```
database/
└── v_21/                    # Dataset version
    └── T5/                  # Embedding type
        └── global-None/     # Embedding configuration
            ├── Train_set_global_v_21.npy
            ├── Dev_set_global_v_21.npy
            └── Test_set_global_v_21.npy
```

### Dataset Format

Each `.npy` file contains a NumPy array where each row represents one protein complex. The `Trainer.get_inputs` method handles slicing these arrays.

| Component | Logic | Description |
|-----------|-------|-------------|
| `prot1_embedding` | `batch[:,:,:1024]` | T5 embeddings for protein 1 |
| `prot2_embedding` | `batch[:,:,1024:-2*max_len]` | T5 embeddings for protein 2 |
| `target_contact_map` | `target[:,:,:max_len]` | Binary contact matrix |
| `target_mask` | `target[:,:,max_len:]` | Mask for padding positions |

**Note:** All sequences are typically padded to `max_len=200` residues as defined in `model_versions.py`.

**Sources:** [src/model_versions.py:19-22](), [src/build_model.py:103-153]()

---

## Step 2: Configure Model Parameters

### Basic Configuration

Edit `model_versions.py` to set training parameters:

```mermaid
graph LR
    subgraph "Configuration Parameters"
        Dataset["dataset_version = 21"]
        MaxLen["max_len = 200"]
        Model["model = Epsilon_3"]
        Emb["emb = T5"]
    end
    
    subgraph "Task Selection"
        TaskType{Task Type}
        TaskType --> CM["interaction_bin"]
        TaskType --> IR["interface_bin"]
    end
    
    subgraph "Model Architecture"
        ProjDim["projection_dim = 128"]
        HidLayers["num_hid_layers = [[0,0,0,0,'vanilla','']]"]
        Activation["activation1 = [['elu', None]]"]
    end
    
    Dataset & MaxLen & Model & Emb --> Generate["OmegaConf.save()"]
    TaskType --> Generate
    ProjDim & HidLayers & Activation --> Generate
```

**Sources:** [src/model_versions.py:17-88](), [src/model_versions.py:151]()

### Key Configuration Parameters

#### Task Objective Configuration

The `objective` parameter in `model_versions.py` controls what the model predicts:

```python
# Format: [task_type, coarse_grain_level, pool_type, bin_post_proj, bin_input, single_output]
objective = ["interface_bin", 10, "avg", False, True, False]
```

**Parameter Explanation:**

| Index | Parameter | Options | Description |
|-------|-----------|---------|-------------|
| 0 | `task_type` | `"interaction_bin"`, `"interface_bin"` | Contact map vs interface prediction |
| 1 | `bin_size` | `1`, `5`, `10` | Coarse-graining level |
| 2 | `pool_type` | `"avg"` | Pooling method |
| 4 | `bin_input` | `True`/`False` | If True, averages input embeddings |

**Sources:** [src/model_versions.py:23-41](), [src/utils.py:92-120]()

#### Model Architecture Parameters

The `model_params` dictionary defines the `Epsilon_3` network structure.

```python
"model_params": {
    "emb_size": 1024,                              
    "projection_layer": [[128, "ln2", True, 1, ""]], 
    "input_layer": ["op-od", "vanilla", "avg2d"],  # avg2d is used for interface tasks
    "num_hid_layers": [[0, 0, 0, 0, "vanilla", ""]], 
    "dropouts": [[0.2, 0, 0, 0, 0]],               
    "activation1": [["elu", None]],                
    "norm": [True, "LN"],                          
}
```

**Key Logic:**
- **`projection_layer`**: Used in `create_projection_layers` to reduce embedding dimensionality. [src/models/Epsilon_3.py:64-80]()
- **`input_layer`**: The `avg2d` parameter in the third index is specifically for interface models to aggregate the 2D interaction tensor into 1D interface profiles. [src/models/Epsilon_3.py:83-90]()

**Sources:** [src/model_versions.py:51-88](), [src/models/Epsilon_3.py:15-128]()

#### Training Parameters

The `train_params` dictionary defines the optimization strategy.

```python
"train_params": {
    "loss": "se_loss",                    
    "optimizer": "AdamW",                 
    "learning_rate": [2e-4],              
    "max_epochs": 25,                     
    "scheduler": {
        "apply": True,
        "name": "exp",
        "gamma": [0.9],
    }
}
```

**Sources:** [src/model_versions.py:102-147](), [src/build_model.py:50-100]()

---

## Step 3: Run Training

### Command Line Execution

```bash
cd src
python hparams_search.py -f version_14.2.yml -m manual
```

**Arguments:**
- `-f`: Path to configuration YAML file. [src/hparams_search.py:40-42]()
- `-m`: Mode (`manual` for single configuration). [src/hparams_search.py:43-45]()

### Training Execution Flow

```mermaid
sequenceDiagram
    participant Script as hparams_search.py
    participant Loader as DatasetLoader
    participant Model as Epsilon_3
    participant Trainer as Trainer
    
    Script->>Loader: load_dataset()
    Loader-->>Script: Train/Dev/Test sets
    Script->>Model: get_model(config)
    Script->>Trainer: Trainer(config, device)
    
    loop max_epochs
        Script->>Trainer: training_step(model, train_set)
        Trainer->>Model: forward(p1, p2, mask)
        Model-->>Trainer: preds
        Trainer->>Trainer: calculate_loss_n_metrics()
        Trainer->>Trainer: optimizer.step()
    end
    
    Script->>Trainer: test_step(model, test_set)
    Trainer->>Trainer: calibrate_model()
    Trainer-->>Script: Final Metrics
```

**Sources:** [src/hparams_search.py:173-205](), [src/build_model.py:219-294]()

---

## Step 4: Monitor and Evaluate

### Calibration and Reliability

The training pipeline includes an optional calibration step using `calibration_curve` from `sklearn`.

```python
def plot_reliabity_diagram( uncal_preds, cal_preds, target, file_name ):
    # Calculates fraction of positives vs mean predicted probability
    target_prob1, cal_pred_prob = calibration_curve( target, cal_preds, n_bins = 10 )
```

**Sources:** [src/utils.py:21-63]()

### Output Files

Upon completion, the `Trainer` and `HparamSearch` classes save results:
- **Checkpoints**: Saved via `torch.save` if `save_model` is True. [src/model_versions.py:144-145]()
- **Metrics**: Dumped into CSV and text summaries. [src/hparams_search.py:202-205]()
- **Plots**: Reliability diagrams and performance curves. [src/utils.py:21-63]()

---

## Advanced: Coarse-Graining Implementation

The `prepare_input` function in `utils.py` handles the data flow for custom coarse-grained models:

```mermaid
flowchart LR
    subgraph "Input Processing"
        P1["Protein 1 Emb"]
        P2["Protein 2 Emb"]
        T["Target Map"]
    end

    subgraph "Coarse Graining"
        P1 -->|AvgPool1d| P1_bin["Binned Emb"]
        P2 -->|AvgPool1d| P2_bin["Binned Emb"]
        T -->|MaxPool2d| T_bin["Binned Target"]
    end

    subgraph "Task Projection"
        T_bin -->|interaction| Flat["Flattened Map"]
        T_bin -->|interface| Cat["Concat P1+P2 Interface"]
    end
```

**Key Code Entity:**
`utils.prepare_input` implements the logic for `bin_size` (1, 5, 10) by applying `nn.AvgPool1d` to embeddings and `nn.MaxPool2d` to target contact maps.

**Sources:** [src/utils.py:92-156](), [src/utils.py:169-194]()

---