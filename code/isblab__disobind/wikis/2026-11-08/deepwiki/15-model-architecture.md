# Model Architecture

<details>
<summary>Relevant source files</summary>

The following files were used as context for generating this wiki page:

- [src/README.md](src/README.md)
- [src/model_versions.py](src/model_versions.py)
- [src/models/Epsilon_3.py](src/models/Epsilon_3.py)
- [src/models/get_model.py](src/models/get_model.py)

</details>



This document provides a comprehensive description of the **Epsilon_3** neural network architecture used in Disobind for predicting protein-protein interactions and interface residues. It covers the model's design, configuration system, and integration with the training pipeline.

**Scope**: This page focuses on the architecture design, configuration files, and model instantiation. For details on training procedures, see [Training Pipeline](#4.3). For information on running predictions with trained models, see [Running Predictions](#2).

---

## Architecture Overview

Disobind uses a neural network called **Epsilon_3** that processes protein embeddings to predict either contact maps (interaction task) or interface residues (interface task) at three coarse-graining levels (CG 1, 5, 10). The same base architecture is used for all six prediction tasks, with configuration parameters specifying task-specific behavior.

### Core Design Philosophy

The Epsilon_3 architecture follows a three-stage processing pipeline:

1. **Projection Block**: Projects high-dimensional protein embeddings (1024-dim T5 embeddings) to lower-dimensional representations.
2. **Interaction Tensor Construction**: Combines projected embeddings from two proteins to capture pairwise residue relationships.
3. **Output Block**: Processes the interaction tensor through optional hidden layers and produces predictions.

```mermaid
graph TB
    Input1["Protein 1 Embedding<br/>(N, L1, 1024)"]
    Input2["Protein 2 Embedding<br/>(N, L2, 1024)"]
    
    subgraph ProjectionBlock["Projection Block"]
        Proj1["projection_layer1<br/>Linear(1024 → 128/256)"]
        Proj2["projection_layer2<br/>Optional separate layer"]
        Drop1["dropout1<br/>p=0.2"]
        Z1["z1 (N, L1, 128/256)"]
        Z2["z2 (N, L2, 128/256)"]
    end
    
    subgraph InteractionBlock["Interaction Tensor Construction"]
        Capture["capture_interaction()<br/>op-od operations"]
        IntTensor["interaction_tensor<br/>(N, L1, L2, 2C)"]
    end
    
    subgraph ProcessingBlock["Processing Block"]
        IntTask{"Task Type?"}
        InterfaceProj["interface_block()<br/>avg2d/avg1d/lin"]
        ContactReshape["Reshape to<br/>(N, L1×L2, 2C)"]
        IntermediateTensor["o (N, Seq, C)"]
    end
    
    subgraph HiddenBlock["Hidden Layers (Optional)"]
        Downsample["downsampling_layers<br/>3 layers, scale=2"]
        Upsample["upsampling_layers<br/>LayerNorm + ELU"]
    end
    
    subgraph OutputBlock["Output Block"]
        OutputLinear["contact<br/>Linear(C → 1)"]
        Sigmoid["sigmoid activation"]
        FinalOutput["Predictions"]
    end
    
    Input1 --> Proj1
    Input2 --> Proj1
    Proj1 --> Drop1
    Drop1 --> Z1
    Drop1 --> Z2
    
    Z1 --> Capture
    Z2 --> Capture
    Capture --> IntTensor
    
    IntTensor --> IntTask
    IntTask -->|"interface"| InterfaceProj
    IntTask -->|"interaction"| ContactReshape
    InterfaceProj --> IntermediateTensor
    ContactReshape --> IntermediateTensor
    
    IntermediateTensor --> Downsample
    Downsample --> Upsample
    Upsample --> OutputLinear
    OutputLinear --> Sigmoid
    Sigmoid --> FinalOutput
```

**Sources**: [src/models/Epsilon_3.py:15-131](), [src/models/Epsilon_3.py:132-354]()

---

## Model Components

### Input Processing and Projection

The model begins by projecting T5 embeddings from 1024 dimensions to a lower-dimensional space (typically 128 or 256 dimensions). The projection layer is created dynamically based on configuration.

| Component | Purpose | Configuration Key | Code Location |
|-----------|---------|-------------------|---------------|
| `projection_layer1` | Projects protein 1 embeddings | `projection_layer` | [src/models/Epsilon_3.py:64-70]() |
| `projection_layer2` | Optional separate projection for protein 2 | `projection_layer[4]` | [src/models/Epsilon_3.py:71-80]() |
| `dropout1` | Regularization after projection | `dropouts[0]` | [src/models/Epsilon_3.py:34]() |

The projection layer type is often specified as `"ln2"`, which creates a Linear → LayerNorm → Activation sequence.

```mermaid
graph LR
    E1["e1<br/>(N, L1, 1024)"]
    E2["e2<br/>(N, L2, 1024)"]
    
    PL1["projection_layer1"]
    PL2["projection_layer2 or<br/>projection_layer1"]
    
    subgraph "Projection Layer (ln2 type)"
        Lin["Linear<br/>(1024 → proj_dim)"]
        LN["LayerNorm<br/>(proj_dim)"]
        Act["activation1<br/>(ELU)"]
    end
    
    D1["dropout1"]
    Z1["z1<br/>(N, L1, proj_dim)"]
    Z2["z2<br/>(N, L2, proj_dim)"]
    
    E1 --> PL1
    E2 --> PL2
    PL1 --> Lin
    PL2 --> Lin
    Lin --> LN
    LN --> Act
    Act --> D1
    D1 --> Z1
    D1 --> Z2
```

**Sources**: [src/models/Epsilon_3.py:64-81](), [src/model_versions.py:59-59]()

---

### Interaction Tensor Construction

The core innovation of Epsilon_3 is how it combines embeddings from two proteins to create an interaction tensor. The model supports multiple aggregation operations specified via the `input_layer` configuration parameter.

#### Aggregation Operations

The `input_layer` parameter format is `["op1-op2", "mode", "interface_method"]`:

| Operation | Mathematical Definition | Shape Transformation |
|-----------|-------------------------|---------------------|
| `op` (outer product) | $z_1 \otimes z_2$ | (N,L1,C) × (N,L2,C) → (N,L1,L2,C) |
| `od` (outer difference) | $|z_1 - z_2|$ | (N,L1,C) × (N,L2,C) → (N,L1,L2,C) |
| `os` (outer sum) | $z_1 + z_2$ | (N,L1,C) × (N,L2,C) → (N,L1,L2,C) |

The standard configuration `"op-od"` concatenates outer product and outer difference operations.

```mermaid
graph TB
    Z1["z1 (N, L1, C)"]
    Z2["z2 (N, L2, C)"]
    
    subgraph CaptureInteraction["capture_interaction()"]
        Unsq1["z1.unsqueeze(2)<br/>(N, L1, 1, C)"]
        Unsq2["z2.unsqueeze(1)<br/>(N, 1, L2, C)"]
        
        OP["z1 * z2<br/>outer product"]
        OD["abs(z1 - z2)<br/>outer difference"]
    end
    
    Concat["torch.cat([z_i1, z_i2], dim=-1)"]
    IntTensor["interaction_tensor<br/>(N, L1, L2, 2C)"]
    
    Z1 --> Unsq1
    Z2 --> Unsq2
    Unsq1 --> OP
    Unsq1 --> OD
    Unsq2 --> OP
    Unsq2 --> OD
    
    OP --> Concat
    OD --> Concat
    Concat --> IntTensor
```

**Sources**: [src/models/Epsilon_3.py:228-296]()

---

### Task-Specific Processing

After constructing the interaction tensor, processing diverges based on the prediction objective:

#### Interaction (Contact Map) Task
For contact map prediction, the interaction tensor (N, L1, L2, 2C) is reshaped to (N, L1×L2, 2C) to treat each residue pair independently.

#### Interface (Residue Label) Task
For interface prediction, the model aggregates over one dimension to identify which residues are involved in interactions. The `interface_block()` method supports three aggregation strategies:

| Strategy | Configuration Value | Operation | Output Shape |
|----------|-------------------|-----------|--------------|
| Average 2D | `"avg2d"` | `avg2d()` with mask handling | (N, L1+L2, 2C) |
| Average 1D | `"avg1d"` | `torch.mean(axis=2)` | (N, L1+L2, 2C) |
| Linear projection | `"lin"` | Learnable linear layer | (N, L1+L2, 2C) |

The `avg2d()` method (default for interface tasks) computes masked averages to ignore padding residues.

**Sources**: [src/models/Epsilon_3.py:181-225]()

---

### Hidden Layers

The model optionally processes the intermediate tensor through downsampling and upsampling layers. These are configured via `num_hid_layers = [#US, #DS, #blocks, scale_factor, block_type, residual_type]`.

Standard configurations use:
- 0 upsampling layers
- 3 downsampling layers  
- Scale factor of 2
- Vanilla block (no residual connections)

**Sources**: [src/models/Epsilon_3.py:92-108](), [src/model_versions.py:73-73]()

---

## Configuration System

### YAML Configuration Files

Model configurations are defined in `src/model_versions.py` and saved as YAML files. These files define all architectural and training parameters.

```mermaid
graph TB
    subgraph ConfigFile["version_X.yml"]
        Version["Version: '14.2'"]
        ModelParams["conf.model_params"]
        DatasetParams["conf.dataset"]
        TrainParams["conf.train_params"]
    end
    
    subgraph ModelParamsDetail["model_params"]
        EmbSize["emb_size: 1024"]
        ProjLayer["projection_layer: [dim, type, bias, mult, sep]"]
        Objective["objective: [task, CG, pool, ...]"]
        InputLayer["input_layer: [op-od, vanilla, avg2d]"]
        Activation["activation1: [elu, null]"]
    end
    
    ConfigFile --> ModelParams
    ConfigFile --> TrainParams
    ModelParams --> ModelParamsDetail
```

**Sources**: [src/model_versions.py:43-151]()

### Model Version Mapping

Each configuration file corresponds to a specific task. For details, see [Model Configuration Files](#4.2).

### Configuration to Model Instantiation

The configuration flow from YAML to model instance:

```mermaid
graph LR
    YAML["version_X.yml<br/>YAML file"]
    OmegaConf["OmegaConf.save()<br/>in model_versions.py"]
    Config["config.conf.model_params<br/>DictConfig object"]
    GetModel["get_model(config)<br/>src/models/get_model.py"]
    Epsilon3["Epsilon_3(...)<br/>Model instance"]
    
    YAML --> OmegaConf
    OmegaConf --> Config
    Config --> GetModel
    GetModel --> Epsilon3
```

The `get_model()` function in [src/models/get_model.py:4-26]() reads the configuration and instantiates `Epsilon_3` with the specified parameters.

---

## Training Pipeline Integration

### Trainer Class Architecture

The `Trainer` class (detailed in [Training Pipeline](#4.3)) orchestrates the complete training loop, integrating the Epsilon_3 model with optimization and loss computation.

### Input Preparation

The training pipeline prepares input tensors by:
1. Extracting embeddings and targets from the batch.
2. Applying coarse-graining (CG) if specified in the `objective`.
3. Reformatting based on task objective.

For details on the training loop and loss functions, see [Training Pipeline](#4.3).

**Sources**: [src/model_versions.py:102-147]()

---

## Model Activation and Normalization

### Activation Functions
The architecture uses two activation functions specified in configuration:
- `activation1`: Hidden layers, typically `["elu", null]`.
- `activation2`: Output layer, typically `["sigmoid", true]`.

**Sources**: [src/models/Epsilon_3.py:39-43](), [src/model_versions.py:83-84]()

### Normalization Strategy
Layer Normalization (LN) is applied in the projection layers and downsampling layers when `norm` is set to `[True, "LN"]`.

**Sources**: [src/models/Epsilon_3.py:57-63](), [src/model_versions.py:78-78]()

---

## Summary of Key Code Entities

| Component | Class/Function | File | Purpose |
|-----------|----------------|------|---------|
| Neural Network | `Epsilon_3` | [src/models/Epsilon_3.py:15]() | Main model architecture |
| Model Factory | `get_model()` | [src/models/get_model.py:4]() | Instantiate model from config |
| Interaction Ops | `capture_interaction()` | [src/models/Epsilon_3.py:228]() | Combine protein embeddings |
| Interface Aggregation | `avg2d()` | [src/models/Epsilon_3.py:181]() | Masked averaging for interface |
| Layer Creation | `create_projection_layers` | [src/models/Epsilon_3.py:10]() | Utility for projection blocks |

**Sources**: [src/models/Epsilon_3.py:1-25](), [src/models/get_model.py:1-27]()

---