# Loss Functions & Training Wrapper

<details>
<summary>Relevant source files</summary>

The following files were used as context for generating this wiki page:

- [idpforge/loss.py](idpforge/loss.py)
- [idpforge/utils/validation_metrics.py](idpforge/utils/validation_metrics.py)
- [idpforge/wrapper.py](idpforge/wrapper.py)

</details>



This page details the training infrastructure of IDPForge, centered around the `IDPForgeWrapper` PyTorch Lightning module. It covers the loss computation logic, optimization strategies, and validation metrics used to train the model for intrinsically disordered protein (IDP) ensemble generation.

## Training Wrapper Architecture

The `IDPForgeWrapper` encapsulates the `IDPForge` model, the diffusion process, and the optimization logic. It manages training and validation steps, weight averaging via EMA, and learning rate scheduling.

### Entity Mapping: Training Logic to Code

The following diagram maps high-level training concepts to their corresponding classes and functions in the codebase.

**Training Wrapper Overview**
```mermaid
graph TD
    subgraph "Lightning Module: IDPForgeWrapper"
        A["IDPForgeWrapper (idpforge/wrapper.py)"] --> B["training_step()"]
        A --> C["validation_step()"]
        A --> D["configure_optimizers()"]
    end

    subgraph "Optimization & Weights"
        D --> E["AlphaFoldLRScheduler (openfold/utils/lr_schedulers.py)"]
        A --> F["ExponentialMovingAverage (openfold/utils/exponential_moving_average.py)"]
        B --> G["on_before_zero_grad() -> ema.update()"]
    end

    subgraph "Loss Computation"
        B --> H["calc_loss() (idpforge/loss.py)"]
        H --> I["backbone_loss (FAPE)"]
        H --> J["torsion_loss"]
        H --> K["cb_dist_loss"]
        H --> L["viol_loss"]
    end
```
**Sources:** [idpforge/wrapper.py:17-36](), [idpforge/wrapper.py:133-135](), [idpforge/loss.py:42-90]()

### Key Components

*   **Self-Conditioning:** During training, the model optionally uses its own predictions from a previous diffusion step as input for the current step. This is implemented in `training_step` with a 50% probability [idpforge/wrapper.py:64-71]().
*   **EMA Weight Management:** The wrapper maintains an Exponential Moving Average (EMA) of the model weights. During validation, it caches the current training weights, loads the EMA weights, and restores the training weights once the validation epoch ends [idpforge/wrapper.py:86-91](), [idpforge/wrapper.py:128-131]().
*   **Learning Rate Scheduler:** It utilizes the `AlphaFoldLRScheduler`, which implements a linear warmup followed by a decay phase, specifically tuned for protein folding models [idpforge/wrapper.py:152-153]().

---

## Loss Functions

The `calc_loss` function aggregates several structural and geometric loss components. The total loss is a weighted sum of these components as defined in the training configuration [idpforge/loss.py:88-89]().

### Loss Component Breakdown

| Loss Component | Function / Source | Description |
| :--- | :--- | :--- |
| **FAPE** | `backbone_loss` & `sidechain_fape` | Frame Aligned Point Error for backbone and sidechain frames [idpforge/loss.py:54-68](). |
| **Torsion** | `torsion_loss` | Mean absolute error of predicted torsion angles (backbone and $\chi$) [idpforge/loss.py:44-52](), [idpforge/loss.py:159-166](). |
| **CB Distance** | `cb_dist_loss` | Penalizes deviations in $C\beta-C\beta$ distances, with specific logic for secondary structure types [idpforge/loss.py:124-146](). |
| **Violation** | `viol_loss` | Penalizes steric clashes and bond length/angle violations using OpenFold utilities [idpforge/loss.py:82-86](). |
| **CA Connectivity**| `ca_connectivity_loss`| Ensures physical continuity of the peptide chain by penalizing $C\alpha-C\alpha$ distance deviations [idpforge/loss.py:135-136](). |

### Data Flow: Loss Computation
The diagram below illustrates how predicted coordinates and frames are processed through various loss functions.

**Loss Data Flow**
```mermaid
graph LR
    subgraph "Model Output"
        P1["pred['positions']"]
        P2["pred['frames']"]
        P3["pred['angles']"]
    end

    subgraph "Loss Components"
        L1["FAPE Loss"]
        L2["Torsion Loss"]
        L3["CB Distance Loss"]
        L4["Violation Loss"]
    end

    P2 --> L1
    P3 --> L2
    P1 --> L3
    P1 --> L4

    subgraph "Ground Truth"
        T1["true_coords"]
        T2["true_rigids"]
        T3["true_torsions"]
    end

    T2 --> L1
    T3 --> L2
    T1 --> L3
```
**Sources:** [idpforge/loss.py:42-87](), [idpforge/loss.py:124-146]()

### Specialized IDP Losses

*   **Secondary Structure Aware Distance Loss:** In `cb_dist_loss`, the model applies different penalties based on the secondary structure type (`sstype`). For residues predicted as coils, it clamps the error to allow for the inherent flexibility of IDPs, while applying stricter penalties to regions assigned to structured states [idpforge/loss.py:138-145]().
*   **CA Connectivity:** A specific term is added to the distance loss to reinforce the 3.8Å $C\alpha-C\alpha$ distance, preventing "breaks" in the generated chain [idpforge/loss.py:135-136]().

---

## Validation & Metrics

Validation is performed using the `recon` method of the model, which performs a reconstruction from noise [idpforge/wrapper.py:93]().

### Validation Metrics

1.  **dRMSD:** Calculates the distance root-mean-square deviation of $C\alpha$ atoms [idpforge/wrapper.py:102-103]().
2.  **Structural Violations:** Monitors the physical plausibility of the generated ensembles [idpforge/wrapper.py:104-105]().
3.  **Rg Distribution (`rg_error`):** A critical metric for IDPs. It computes the divergence between the predicted and true Radius of Gyration ($R_g$) distributions.
    *   **Logic:** The `rg_dist_per_group` function groups conformers by sequence identity and compares the mean $R_g$ of the predicted ensemble against the ground truth [idpforge/utils/validation_metrics.py:32-53]().

### Metric Implementations

| Metric | Code Reference |
| :--- | :--- |
| **Radius of Gyration** | `calc_rg_with_mask` [idpforge/utils/validation_metrics.py:4-29]() |
| **Grouped Rg Error** | `rg_dist_per_group` [idpforge/utils/validation_metrics.py:32-53]() |
| **Violation Metric** | `viol_loss(..., return_metric=True)` [idpforge/wrapper.py:104]() |

**Sources:** [idpforge/wrapper.py:96-123](), [idpforge/utils/validation_metrics.py:1-53]()

---