# Training Pipeline

<details>
<summary>Relevant source files</summary>

The following files were used as context for generating this wiki page:

- [idpforge/loader.py](idpforge/loader.py)
- [idpforge/loss.py](idpforge/loss.py)
- [idpforge/utils/validation_metrics.py](idpforge/utils/validation_metrics.py)
- [idpforge/wrapper.py](idpforge/wrapper.py)
- [train.py](train.py)

</details>



The training pipeline for IDPForge is built on **PyTorch Lightning**, providing a structured framework for training the diffusion model on intrinsically disordered protein (IDP) ensembles. The pipeline handles the orchestration of rotational and Euclidean diffusion schedules, complex loss calculations including structural violations and secondary structure-aware distance losses, and model weight management via Exponential Moving Average (EMA).

## Training Architecture Overview

The training process is centralized in `train.py`, which initializes the diffusion parameters, data modules, and the model wrapper.

### Training Flow
The following diagram illustrates the relationship between the data loading, the diffusion process, and the training loop.

```mermaid
graph TD
    subgraph "Data Preparation"
        A["Pickle Files (.pkl)"] --> B["DiffDataset"]
        B --> C["BatchCollator"]
    end

    subgraph "Training Loop (IDPForgeWrapper)"
        C --> D["training_step()"]
        D --> E["IDPForge Model"]
        E --> F["calc_loss()"]
        F --> G["Optimizer/Scheduler"]
        G --> H["EMA Update"]
    end

    subgraph "Entities"
        B -- "uses" --> I["Diffuser.diffuse_pose()"]
        E -- "calls" --> J["FoldingTrunk / StructureModule"]
        F -- "computes" --> K["FAPE / Torsion / Dist Loss"]
    end
```
**Sources:** [train.py:21-39](), [idpforge/wrapper.py:56-83](), [idpforge/loader.py:18-45]()

---

## Training Entrypoint & Data Loading

The `train.py` script serves as the primary entrypoint. It configures the `Diffuser` and `Denoiser` based on the `n_tsteps` defined in the YAML configuration and passes them to the `IDPloader`.

*   **Data Loading:** `IDPloader` (a `LightningDataModule`) manages `DiffDataset` instances that load pre-processed protein ensembles.
*   **Diffusion in Dataset:** The `DiffDataset` applies noise to ground-truth coordinates on-the-fly using `diffuser.diffuse_pose()`, generating the $x_t$ (Euclidean) and $\alpha_t$ (torsional) inputs required for the denoising task.
*   **Batching:** `BatchCollator` utilizes `input_process()` to handle sequence encoding and secondary structure mapping for variable-length proteins in a single batch.

For details on data formats and the loading mechanism, see [Training Entrypoint & Data Loading](#3.1).

**Sources:** [train.py:25-35](), [idpforge/loader.py:62-99](), [idpforge/loader.py:102-115]()

---

## Loss Functions & Training Wrapper

Training logic is encapsulated in `IDPForgeWrapper`, which inherits from `pl.LightningModule`. It manages the forward pass, loss calculation, and validation metrics.

### Key Components
*   **IDPForgeWrapper:** Orchestrates the `training_step` and `validation_step`. It implements **self-conditioning**, where the model occasionally uses its own previous prediction as an input to improve refinement.
*   **Loss Function (`calc_loss`):** A composite loss function including:
    *   **FAPE Loss:** Frame Aligned Point Error for backbone and sidechain atoms.
    *   **Torsion Loss:** Measures the error in predicted dihedral angles.
    *   **CB Distance Loss:** A custom distance loss that enforces secondary structure constraints and CA-CA connectivity.
*   **EMA & Scheduling:** Uses `ExponentialMovingAverage` to maintain a smoothed version of model weights for validation and inference. The learning rate is managed by `AlphaFoldLRScheduler`.

For details on specific loss terms and the PyTorch Lightning implementation, see [Loss Functions & Training Wrapper](#3.2).

**Sources:** [idpforge/wrapper.py:17-36](), [idpforge/wrapper.py:65-76](), [idpforge/loss.py:42-89](), [idpforge/wrapper.py:133-135]()

---

## Training Logic to Code Mapping

This diagram maps high-level training concepts to the specific classes and functions in the IDPForge codebase.

```mermaid
graph LR
    subgraph "Natural Language Space"
        Direction["Diffusion Training"]
        Data["Data Augmentation"]
        Optimization["Loss Optimization"]
        Metric["Validation Metric"]
    end

    subgraph "Code Entity Space"
        Direction --- IDPForge["IDPForge.forward()"]
        Data --- DiffusePose["Diffuser.diffuse_pose()"]
        Optimization --- CalcLoss["idpforge/loss.py:calc_loss()"]
        Metric --- RgDist["rg_dist_per_group()"]
    end

    IDPForge -- "calls" --> Trunk["FoldingTrunk"]
    CalcLoss -- "includes" --> FAPE["compute_fape()"]
    CalcLoss -- "includes" --> CBDist["cb_dist_loss()"]
```
**Sources:** [idpforge/model.py:6-7](), [idpforge/loss.py:42-52](), [idpforge/utils/diff_utils.py:11-12](), [idpforge/utils/validation_metrics.py:32-38]()

## Summary Table of Training Configuration

| Component | Code Reference | Description |
| :--- | :--- | :--- |
| **Model Wrapper** | `IDPForgeWrapper` | PyTorch Lightning module for training/validation. |
| **Loss Calculation** | `calc_loss()` | Aggregates FAPE, torsion, distance, and violation losses. |
| **Data Module** | `IDPloader` | Handles training/validation datasets and multi-worker loading. |
| **LR Scheduler** | `AlphaFoldLRScheduler` | Custom scheduler with warmup and decay steps. |
| **EMA** | `ExponentialMovingAverage` | Maintains stable weights for evaluation. |

**Sources:** [idpforge/wrapper.py:17](), [idpforge/loss.py:42](), [idpforge/loader.py:102](), [idpforge/wrapper.py:152](), [idpforge/wrapper.py:27]()

---