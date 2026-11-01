# Sampling Algorithms

<details>
<summary>Relevant source files</summary>

The following files were used as context for generating this wiki page:

- [starling/frontend/__init__.py](starling/frontend/__init__.py)
- [starling/samplers/__init__.py](starling/samplers/__init__.py)
- [starling/samplers/ddim_sampler.py](starling/samplers/ddim_sampler.py)
- [starling/samplers/ddpm_sampler.py](starling/samplers/ddpm_sampler.py)
- [starling/samplers/plms_sampler.py](starling/samplers/plms_sampler.py)
- [starling/structure/__init__.py](starling/structure/__init__.py)
- [starling/tests/test_sequence_encoder_backend.py](starling/tests/test_sequence_encoder_backend.py)
- [starling/tests/test_sequence_encoder_backend_integration.py](starling/tests/test_sequence_encoder_backend_integration.py)

</details>



Sampling algorithms in STARLING are responsible for the reverse diffusion process, transforming Gaussian noise into structured protein distance map latents. The codebase provides three distinct strategies, ranging from the standard stochastic Markovian process to accelerated non-Markovian and multi-step solvers.

The samplers operate on the latent space defined by the VAE and are conditioned on sequence embeddings and ionic strength values. Each sampler implements a `sample` method that manages the denoising loop, integrates with the `Constraint` framework, and produces distance maps ready for 3D reconstruction.

### Algorithm Selection Logic

The choice of sampler typically involves a trade-off between sampling speed and the fidelity of the generated ensemble.

| Sampler | Strategy | Steps (Typical) | Best For |
| :--- | :--- | :--- | :--- |
| **DDPM** | Stochastic Markovian | 1000 | Maximum diversity, standard reference. |
| **DDIM** | Deterministic/Stochastic | 10 - 100 | Fast generation, deterministic mapping. |
| **PLMS** | Multi-step Adams-Bashforth | 10 - 50 | High-quality accelerated sampling. |

---

## Sampler Architecture and Data Flow

All samplers interact with the `DiffusionModel` to obtain noise predictions and use the `StarlingTokenizer` to process input sequences. They share a common pattern for conditioning the model using `sequence2labels` [starling/samplers/ddpm_sampler.py:86-88]().

### Sampling Data Flow Diagram
This diagram illustrates how the `Sampler` classes bridge the gap between user inputs (Sequence/Ionic Strength) and the internal `DiffusionModel` (ViT backbone).

```mermaid
graph TD
    subgraph "Natural Language & Params"
        Seq["Protein Sequence (str)"]
        IS["Ionic Strength (float)"]
    end

    subgraph "Sampler Entity Space"
        Tokenizer["StarlingTokenizer"]
        S["Sampler Class (DDPM/DDIM/PLMS)"]
        GenLab["generate_labels()"]
    end

    subgraph "Diffusion Core"
        ViT["DiffusionModel (ViT)"]
        Scale["latent_space_scaling_factor"]
    end

    Seq --> Tokenizer
    Tokenizer --> GenLab
    IS --> GenLab
    GenLab -->|"Conditioning Labels"| S
    S -->|"p_sample / denoise"| ViT
    ViT -->|"Noise Prediction"| S
    Scale -.->|"Rescale Latents"| S
```
**Sources:** [starling/samplers/ddpm_sampler.py:44-90](), [starling/samplers/ddim_sampler.py:102-124](), [starling/samplers/plms_sampler.py:124-147]().

---

## DDPM Sampler

The `DDPMSampler` implements the standard Denoising Diffusion Probabilistic Model process. It follows a Markovian chain where each step $t-1$ depends only on step $t$. It is the most computationally expensive but serves as the baseline for ensemble quality.

*   **Key Method:** `p_sample_loop` [starling/samplers/ddpm_sampler.py:150-178]() orchestrates the full sequence of $T$ steps.
*   **Step Logic:** `p_sample` [starling/samplers/ddpm_sampler.py:92-149]() calculates the predicted mean and injects posterior variance noise.
*   **Use Case:** Use this when computational time is not a constraint and you require the most faithful reproduction of the training distribution.

For details, see [DDPM Sampler](#3.1).

**Sources:** [starling/samplers/ddpm_sampler.py:44-230]().

---

## Accelerated Samplers: DDIM and PLMS

STARLING provides two accelerated samplers that significantly reduce the number of required function evaluations (NFE) by discretizing the diffusion process into fewer steps.

### DDIM (Denoising Diffusion Implicit Models)
The `DDIMSampler` allows for non-Markovian sampling, enabling deterministic generation when `ddim_eta` is set to 0.0 [starling/samplers/ddim_sampler.py:27-52](). It supports different discretization schedules such as `uniform` or `quad` [starling/samplers/ddim_sampler.py:73-81]().

### PLMS (Pseudo Linear Multi-Step)
The `PLMSSampler` uses an Adams-Bashforth multi-step method to improve the accuracy of the denoising trajectory at higher speeds [starling/samplers/plms_sampler.py:62-72](). It includes features like `dynamic_thresholding_fn` [starling/samplers/plms_sampler.py:26-46]() to prevent latent collapse during accelerated sampling.

### Discretization and Constraint Interaction
This diagram shows how accelerated samplers map the full training timesteps to a reduced inference schedule while maintaining hooks for physical constraints.

```mermaid
graph LR
    subgraph "Discretization Logic"
        Steps["n_steps (e.g. 50)"]
        FullT["num_timesteps (1000)"]
        Sched["ddim_time_steps"]
    end

    subgraph "Constraint Lifecycle"
        CInit["Constraint.initialize()"]
        CLog["ConstraintLogger"]
        CApply["apply_constraints()"]
    end

    FullT --> Steps
    Steps --> Sched
    Sched -->|"Reverse Loop"| CApply
    CInit --> CApply
    CApply --> CLog
```
**Sources:** [starling/samplers/ddim_sampler.py:72-101](), [starling/samplers/ddim_sampler.py:194-207](), [starling/samplers/plms_sampler.py:89-110]().

For details, see [Accelerated Samplers: DDIM and PLMS](#3.2).

---

## Constraint Integration

All samplers support guided diffusion through the `constraint` parameter in their `sample` methods. When a constraint is provided:
1.  It is initialized with the VAE encoder and latent scaling factors [starling/samplers/ddpm_sampler.py:201-208]().
2.  A `ConstraintLogger` is instantiated to track potential energy and gradients across timesteps [starling/samplers/ddpm_sampler.py:197-200]().
3.  Gradients from physical priors (e.g., `RgConstraint`, `DistanceConstraint`) are applied to the latents during the denoising loop to steer the generation toward specific conformational properties.

**Sources:** [starling/samplers/ddpm_sampler.py:10-15](), [starling/samplers/ddim_sampler.py:11-16](), [starling/samplers/plms_sampler.py:11-16]().

---