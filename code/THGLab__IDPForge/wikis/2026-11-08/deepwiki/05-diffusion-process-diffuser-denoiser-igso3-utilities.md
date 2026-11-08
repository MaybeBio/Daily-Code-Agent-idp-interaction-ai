# Diffusion Process: Diffuser, Denoiser & IGSO3 Utilities

<details>
<summary>Relevant source files</summary>

The following files were used as context for generating this wiki page:

- [idpforge/utils/diff_utils.py](idpforge/utils/diff_utils.py)
- [idpforge/utils/igso3_utils.py](idpforge/utils/igso3_utils.py)
- [idpforge/utils/pre_relax.py](idpforge/utils/pre_relax.py)

</details>



This page details the diffusion framework within IDPForge, covering the mathematical schedules for Euclidean and torsional noise, the implementation of SO(3) rotational diffusion for protein backbones, and the mechanics of the reverse diffusion sampling loop.

## Overview of the Diffusion Framework

IDPForge utilizes a diffusion probabilistic model to generate protein ensembles. The process involves two primary components:
1.  **Forward Diffusion**: Gradually adding noise to protein structures (coordinates and rotations) according to predefined schedules.
2.  **Reverse Diffusion (Denosing)**: Using the neural network to predict the clean structure ($x_0$) from a noisy structure ($x_t$), then stepping back to $x_{t-1}$.

### System Entity Map

The following diagram maps the conceptual diffusion steps to the specific classes and functions implemented in the codebase.

"Entity Mapping: Diffusion Logic to Code"
```mermaid
graph TD
    subgraph "Noise Schedules"
        A["Beta Schedule"] --> B["linear_beta_schedule"]
        C["Noise Scale"] --> D["get_noise_schedule"]
    end

    subgraph "Diffusion Engine"
        E["Diffuser Class"] --> F["idpforge/model.py"]
        G["SO3 Diffuser"] --> H["SO3Diffuser Class"]
    end

    subgraph "Sampling Operations"
        I["Reverse Step"] --> J["get_next_ca"]
        I --> K["get_next_frames"]
        I --> L["get_next_chi_angles"]
        M["Alignment"] --> N["align_coords (Kabsch)"]
    end

    B --- [idpforge/utils/diff_utils.py:56-61]
    D --- [idpforge/utils/diff_utils.py:64-76]
    J --- [idpforge/utils/diff_utils.py:98-125]
    K --- [idpforge/utils/diff_utils.py:126-161]
    L --- [idpforge/utils/diff_utils.py:162-185]
    N --- [idpforge/utils/diff_utils.py:25-54]
```
**Sources:** [idpforge/utils/diff_utils.py:1-185](), [idpforge/model.py:1-100]()

---

## Noise Schedules and Beta Calculations

IDPForge supports linear schedules for the variance of the noise added at each timestep.

*   **Linear Beta Schedule**: Computes a sequence of $\beta$ values from $\beta_0$ to $\beta_T$. The values are scaled by the number of steps to maintain consistency across different sampling lengths [idpforge/utils/diff_utils.py:56-61]().
*   **Alphabar ($\bar{\alpha}$)**: Represents the cumulative product of $(1 - \beta_t)$, determining the total amount of noise at time $t$.
*   **Coordinate Scaling**: Coordinates are typically scaled by a factor (e.g., `crd_scale=0.25`) during diffusion to normalize the variance [idpforge/utils/diff_utils.py:187]().

**Sources:** [idpforge/utils/diff_utils.py:56-96]()

---

## IGSO3 Rotational Diffusion

For backbone orientations, IDPForge uses the **Isotropic Gaussian on SO(3)** (IGSO3) distribution. This allows for principled diffusion on the manifold of rotation matrices rather than in Euclidean space.

### Key Utilities
*   **Exp/Log Maps**: Functions to map between the Lie algebra $\mathfrak{so}(3)$ (rotation vectors) and the Lie group $SO(3)$ (rotation matrices) [idpforge/utils/igso3_utils.py:27-40]().
*   **Numerical Approximation**: Since the IGSO3 density involves an infinite sum, `f_igso3` provides a truncated power series approximation [idpforge/utils/igso3_utils.py:42-64]().
*   **Pre-computation**: `calculate_igso3` discretizes the $\sigma$ and $\omega$ (angle) spaces to pre-compute CDFs and score norms, which are stored in `diff_igso3.pkl` for fast lookup during sampling [idpforge/utils/igso3_utils.py:88-132]().

"Data Flow: SO(3) Sampling"
```mermaid
flowchart LR
    subgraph "Pre-computation"
        A["calculate_igso3"] --> B["diff_igso3.pkl"]
    end
    
    subgraph "Inference"
        B --> C["SO3Diffuser.sample_vec"]
        C --> D["Rotation Matrix (R_t)"]
        D --> E["get_next_frames"]
    end
```
**Sources:** [idpforge/utils/igso3_utils.py:1-132](), [idpforge/utils/diff_utils.py:126-161]()

---

## Reverse Diffusion Sampling Loop

The sampling process iteratively refines a structure from pure noise ($t=T$) to a folded ensemble ($t=0$).

### 1. Initialization
The process starts at `init_sample`, which places residues at the origin and applies random translations and IGSO3-sampled rotations [idpforge/utils/diff_utils.py:187-196]().

### 2. Denoising Steps
At each step $t$, the model predicts $x_0$. The following functions then calculate $x_{t-1}$:
*   **`get_next_ca`**: Updates $C\alpha$ positions using the Gaussian posterior mean $\mu$ and variance $\sigma$ [idpforge/utils/diff_utils.py:98-125]().
*   **`get_next_frames`**: Updates residue orientations. It converts coordinates to rigid frames using `rigid_from_3_points_np`, calculates the rotation transition via the `SO3Diffuser`, and applies it to the backbone atoms [idpforge/utils/diff_utils.py:126-161]().
*   **`get_next_chi_angles`**: Diffuses sidechain torsion angles ($\chi$) while ensuring they stay within the $[-\pi, \pi]$ range using `wrap_rad` [idpforge/utils/diff_utils.py:162-185]().

### 3. Alignment (Kabsch)
To prevent the protein from drifting during diffusion, `align_coords` uses the Kabsch algorithm to align the predicted $C\alpha$ coordinates to a reference (usually the current noisy state), ensuring the center of mass and orientation are consistent [idpforge/utils/diff_utils.py:25-54]().

### 4. Continuity Screening
Before AMBER relaxation, the system performs a "pre-relax" check. `check_backbone_continuity` calculates $C\alpha-C\alpha$ distances. If the distance exceeds thresholds (e.g., $9.12\text{\AA}$ for backbone or $6.46\text{\AA}$ for junctions), the sample may be flagged as physically unrealistic [idpforge/utils/pre_relax.py:10-29]().

| Parameter | Threshold ($\text{\AA}$) | Source |
| :--- | :--- | :--- |
| `_JUNCTION_CA_THRESHOLD` | 6.46 | [idpforge/utils/pre_relax.py:6]() |
| `_BACKBONE_CA_THRESHOLD` | 9.12 | [idpforge/utils/pre_relax.py:7]() |

**Sources:** [idpforge/utils/diff_utils.py:25-196](), [idpforge/utils/pre_relax.py:1-29]()

---